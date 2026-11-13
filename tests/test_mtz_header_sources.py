"""The MTZ header is the space-group and wavelength source; Data_1*.cif and
summary.html are optional cross-checks (for example the 2021 iCiU_Ag dataset,
which has only the MTZ files and autoPROC report PDFs)."""

import json
import tempfile
import unittest
from pathlib import Path

from nasolve.automr import AutoMRInputError, prepare_automr
from nasolve.automr_input import discover_dataset
from nasolve.autosol import AutoSolPreparationError, resolve_wavelength
from nasolve.campaigns import plan_campaign
from nasolve.reflection_data import ReflectionDataError, mtz_wavelength

from .helpers import make_dataset, make_mtz_dump, model_text, mtz_dump_columns_text, w_model_text
from .test_autosol_metal import ENV, make_autosol, make_run
from .test_autosol import make_mtz_dump as autosol_mtz_dump
from nasolve.autosol_metal import execute_metal_search


VALID = {"1AP", "DT", "DC", "DG", "DA"}


def _mtz_only_dataset(root: Path) -> Path:
    dataset = make_dataset(root, include_model=False)
    (dataset / "data_1.cif").unlink()
    (dataset / "summary.html").unlink()
    (dataset / "report_staraniso.pdf").write_bytes(b"%PDF-1.4 autoPROC report\n")
    return dataset


def _w_catalogue(root: Path) -> Path:
    catalogue = root / "frames" / "5W6W"
    catalogue.mkdir(parents=True)
    (catalogue / "C_G.pdb").write_text(w_model_text())
    return root / "frames"


class OptionalAutoPROCFilesTests(unittest.TestCase):
    def test_dataset_with_only_mtz_is_discovered(self):
        # Regression: discovery required Data_1*.cif and summary.html.
        with tempfile.TemporaryDirectory() as directory:
            dataset = _mtz_only_dataset(Path(directory) / "iCiU_Ag")
            files = discover_dataset(dataset)
            self.assertEqual(files.reflections.name, "staraniso-alldata.mtz")
            self.assertIsNone(files.metadata)
            self.assertIsNone(files.summary)

    def test_standard_preflight_uses_the_mtz_header_alone(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = _mtz_only_dataset(root / "iCiU_Ag")
            result = prepare_automr(
                dataset, frame_override="W", pair_override="D:T",
                frames_dir=_w_catalogue(root), valid_ligand_codes=VALID,
                mtz_dump_executable=make_mtz_dump(root / "tools"),
            )
            report = json.loads(result.report_path.read_text())
            self.assertIsNone(report["inputs"]["metadata"])
            self.assertIsNone(report["inputs"]["summary"])
            evidence = report["symmetry"]["evidence"]
            self.assertEqual(evidence["compared_sources"], ["MTZ header", "MTZ symmetry matrices"])
            self.assertEqual(evidence["cif_symbols"], [])
            self.assertIsNone(evidence["summary_symbol"])
            log = (result.run_directory / "automr.log").read_text()
            self.assertIn("Data_1*.cif: absent", log)
            self.assertIn("summary.html: absent", log)

    def test_present_cif_or_summary_must_still_agree_with_the_mtz(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frames = _w_catalogue(root)
            dump = make_mtz_dump(root / "tools")
            for name, text in (
                ("data_1.cif", "data_x\n_symmetry.space_group_name_H-M 'P 1'\n"),
                ("summary.html", "Spacegroup name          P1\n"),
            ):
                with self.subTest(name=name):
                    dataset = _mtz_only_dataset(root / name)
                    (dataset / name).write_text(text)
                    with self.assertRaisesRegex(AutoMRInputError, "sources disagree"):
                        prepare_automr(
                            dataset, frame_override="W", pair_override="D:T",
                            frames_dir=frames, valid_ligand_codes=VALID,
                            mtz_dump_executable=dump,
                        )
            dataset = make_dataset(root / "all_present", include_model=False)
            result = prepare_automr(
                dataset, frame_override="W", pair_override="D:T",
                frames_dir=frames, valid_ligand_codes=VALID, mtz_dump_executable=dump,
            )
            evidence = json.loads(result.report_path.read_text())["symmetry"]["evidence"]
            self.assertEqual(
                evidence["compared_sources"],
                ["MTZ header", "MTZ symmetry matrices", "data_1.cif", "summary.html"],
            )

    def test_campaign_plans_a_dataset_with_only_mtz(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            campaign = root / "campaign"
            campaign.mkdir()
            dataset = _mtz_only_dataset(campaign / "iCiU_Ag")
            (dataset / "nasolve.txt").write_text("[automr]\nframe = W\npair = C:G\n")
            frames = root / "MR_frames"
            (frames / "5W6W").mkdir(parents=True)
            (frames / "5W6W" / "C_G.pdb").write_text(model_text())
            (frames / "5W6W" / "seq_base.txt").write_text("ACGT\n")
            entry = plan_campaign(campaign, frames_directory=frames)["datasets"][0]
            self.assertEqual(entry["status"], "DISCOVERED", entry["diagnostic"])
            self.assertNotIn("metadata", entry["inputs"])
            self.assertNotIn("summary", entry["inputs"])


class MtzWavelengthTests(unittest.TestCase):
    def test_data_dataset_wavelength_not_hkl_base_placeholder(self):
        value, dataset = mtz_wavelength(mtz_dump_columns_text(wavelength=1.00743))
        self.assertEqual((value, dataset), (1.00743, "1"))

    def test_missing_implausible_or_disagreeing_wavelengths_stop(self):
        text = mtz_dump_columns_text(wavelength=1.0)
        with self.assertRaisesRegex(ReflectionDataError, "no wavelength"):
            mtz_wavelength(text.replace("    Wavelength: 1\n", ""))
        with self.assertRaisesRegex(ReflectionDataError, "implausible"):
            mtz_wavelength(mtz_dump_columns_text(wavelength=9.999))
        second = mtz_dump_columns_text(
            (("IMEAN", "J"), ("SIGIMEAN", "Q")), wavelength=0.9795
        ).replace("Crystal 1:", "Crystal 3:").replace("Crystal 2:", "Crystal 4:")
        with self.assertRaisesRegex(ReflectionDataError, "disagree"):
            mtz_wavelength(text + second)
        with self.assertRaisesRegex(ReflectionDataError, "no MTZ dataset"):
            mtz_wavelength("Space group symbol from file: H 3\n")

    def test_summary_is_an_optional_cross_check(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mtz = root / "data.mtz"
            dump = mtz_dump_columns_text(wavelength=1.00743)
            value, source = resolve_wavelength(dump, mtz, None)
            self.assertEqual(value, 1.00743)
            self.assertEqual(source["method"], "MTZ header (phenix.mtz.dump)")
            self.assertIsNone(source["summary_cross_check"])

            summary = root / "summary.html"
            summary.write_text("wavelength [A] = 1.007430\n")
            value, source = resolve_wavelength(dump, mtz, summary)
            self.assertEqual(source["summary_cross_check"]["wavelength"], 1.00743)

            summary.write_text("<html>no wavelength here</html>\n")
            value, source = resolve_wavelength(dump, mtz, summary)
            self.assertEqual(value, 1.00743)
            self.assertIsNone(source["summary_cross_check"]["wavelength"])

            summary.write_text("wavelength [A] = 0.979500\n")
            with self.assertRaisesRegex(AutoSolPreparationError, "disagree"):
                resolve_wavelength(dump, mtz, summary)

    def test_metal_search_runs_without_summary_html(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = make_run(root)
            report_path = run / "report.json"
            report = json.loads(report_path.read_text())
            report["inputs"]["summary"] = None
            report_path.write_text(json.dumps(report))
            (root / "dataset" / "summary.html").unlink()
            result = execute_metal_search(
                run, make_autosol(root, []), autosol_mtz_dump(root, wavelength=1.00743),
                ("AG",), environment=ENV, processor_count=8,
            )
            attempt = json.loads(result.report_path.read_text())
            self.assertEqual(attempt["wavelength"], 1.00743)
            self.assertEqual(attempt["wavelength_source"]["method"], "MTZ header (phenix.mtz.dump)")
            self.assertIn("autosol.lambda=1.00743", " ".join(attempt["command"]))
            self.assertIsNotNone(result.status)


if __name__ == "__main__":
    unittest.main()
