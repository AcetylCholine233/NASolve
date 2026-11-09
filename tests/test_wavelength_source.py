"""The X-ray wavelength comes from Data_1*.cif; summary.html is optional."""

import json
import tempfile
import unittest
from pathlib import Path

from nasolve.automr_input import AutoMRIntent, discover_dataset, resolve_automr_input
from nasolve.autosol import (
    AutoSolPreparationError,
    execute_autosol,
    read_cif_wavelength,
    resolve_wavelength,
)
from nasolve.autosol_metal import FOUND, execute_metal_search
from nasolve.campaigns import campaign_status, plan_campaign
from nasolve.symmetry import SymmetryError, assess_standard_symmetry

from . import test_autosol, test_autosol_metal
from .helpers import make_dataset, make_mtz_dump, model_text


ENV = {"PATH": "/usr/bin:/bin"}

# Shape of an autoPROC Data_1_autoPROC_STARANISO_all.cif: several data blocks,
# each with a key-value wavelength entry.
AUTOPROC_CIF = (
    "data_merged\n"
    "_diffrn_radiation.wavelength_id 1\n"
    "_diffrn_radiation_wavelength.id 1\n"
    "_diffrn_radiation_wavelength.wavelength 1.00744\n"
    "_symmetry.space_group_name_H-M 'R 3'\n"
    "data_unmerged\n"
    "_diffrn_radiation_wavelength.id 1\n"
    "_diffrn_radiation_wavelength.wavelength 1.00744\n"
)


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


class CifWavelengthTests(unittest.TestCase):
    def read(self, text):
        with tempfile.TemporaryDirectory() as directory:
            return read_cif_wavelength(write(Path(directory) / "Data_1.cif", text))

    def test_autoproc_key_value_blocks(self):
        value, method = self.read(AUTOPROC_CIF)
        self.assertEqual(value, 1.00744)
        self.assertIn("_diffrn_radiation_wavelength.wavelength", method)

    def test_loop_uses_the_selected_wavelength_id(self):
        value, _ = self.read(
            "data_x\n_diffrn_radiation.wavelength_id 2\nloop_\n"
            "_diffrn_radiation_wavelength.id\n_diffrn_radiation_wavelength.wavelength\n"
            "_diffrn_radiation_wavelength.wt\n1 0.97950 1.0\n2 1.00744 1.0\n#\n"
        )
        self.assertEqual(value, 1.00744)

    def test_several_wavelengths_without_selection_are_refused(self):
        with self.assertRaisesRegex(AutoSolPreparationError, "more than one wavelength"):
            self.read(
                "data_x\nloop_\n_diffrn_radiation_wavelength.id\n"
                "_diffrn_radiation_wavelength.wavelength\n1 0.9795\n2 1.0074\n"
            )

    def test_missing_and_implausible_values_are_refused(self):
        with self.assertRaisesRegex(AutoSolPreparationError, "No _diffrn_radiation"):
            self.read("data_x\n_cell.length_a 10\n")
        with self.assertRaisesRegex(AutoSolPreparationError, "Implausible"):
            self.read("data_x\n_diffrn_radiation_wavelength.wavelength 9.9\n")


class ResolveWavelengthTests(unittest.TestCase):
    def report(self, root: Path, *, cif=None, summary=None):
        inputs = {}
        if cif is not None:
            inputs["metadata"] = str(write(root / "Data_1.cif", cif))
        if summary is not None:
            inputs["summary"] = str(write(root / "summary.html", summary))
        return {"inputs": inputs}

    def test_cif_is_used_without_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            value, source = resolve_wavelength(self.report(root, cif=AUTOPROC_CIF), root)
            self.assertEqual(value, 1.00744)
            self.assertTrue(source["path"].endswith("Data_1.cif"))
            self.assertNotIn("summary_check", source)

    def test_cif_wins_and_summary_must_agree(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = self.report(root, cif=AUTOPROC_CIF, summary="wavelength [A] = 1.00740\n")
            value, source = resolve_wavelength(report, root)
            self.assertEqual(value, 1.00744)
            self.assertIn("agrees", source["summary_check"])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = self.report(root, cif=AUTOPROC_CIF, summary="wavelength [A] = 1.377618\n")
            with self.assertRaisesRegex(AutoSolPreparationError, "disagree"):
                resolve_wavelength(report, root)

    def test_legacy_runs_fall_back_to_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            value, source = resolve_wavelength(
                self.report(root, summary="wavelength [A] = 1.377618\n"), root
            )
            self.assertEqual(value, 1.377618)
            self.assertTrue(source["path"].endswith("summary.html"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = self.report(root, cif="data_x\n", summary="wavelength [A] = 1.377618\n")
            value, source = resolve_wavelength(report, root)
            self.assertEqual(value, 1.377618)
            self.assertIn("not used", source["cif_check"])

    def test_no_source_is_an_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(AutoSolPreparationError, "No wavelength source"):
                resolve_wavelength({"inputs": {}}, root)
            with self.assertRaisesRegex(AutoSolPreparationError, "No _diffrn"):
                resolve_wavelength(self.report(root, cif="data_x\n"), root)


class OptionalSummaryTests(unittest.TestCase):
    def test_dataset_without_summary_is_discovered(self):
        with tempfile.TemporaryDirectory() as directory:
            root = make_dataset(Path(directory) / "PP_Hg")
            (root / "summary.html").unlink()
            files = discover_dataset(root)
            self.assertIsNone(files.summary)
            self.assertEqual(files.metadata.name, "data_1.cif")
            resolved = resolve_automr_input(
                root, AutoMRIntent(), valid_ligand_codes={"DA", "DC", "DG", "DT"}
            )
            self.assertIsNone(resolved.dataset.summary)

    def test_two_summaries_are_still_ambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            root = make_dataset(Path(directory) / "d")
            (root / "SUMMARY.HTML").write_text("other\n")
            if len([p for p in root.iterdir() if p.name.casefold() == "summary.html"]) > 1:
                with self.assertRaises(Exception):
                    discover_dataset(root)

    def test_symmetry_cross_check_without_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mtz = root / "data.mtz"
            mtz.write_bytes(b"mtz")
            cif = write(root / "Data_1.cif", AUTOPROC_CIF)
            assessment = assess_standard_symmetry(
                mtz, cif, None, make_mtz_dump(root, "H 3", 146, matrix_symbol="R 3 :H"),
            )
            self.assertEqual(assessment.evidence.normalized_class, "H3/R3")
            self.assertIsNone(assessment.evidence.summary_symbol)
            bad = write(root / "bad.cif", "data_x\n_symmetry.space_group_name_H-M 'P 1'\n")
            with self.assertRaisesRegex(SymmetryError, "no summary.html"):
                assess_standard_symmetry(mtz, bad, None, make_mtz_dump(root, "H 3", 146))

    def test_campaign_plans_a_dataset_without_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory).resolve()
            campaign = base / "campaign"
            catalogue = base / "MR_frames" / "5W6W"
            catalogue.mkdir(parents=True)
            (catalogue / "C_G.pdb").write_text(model_text())
            (catalogue / "seq_base.txt").write_text("ACGT\n")
            preset = write(base / "preset" / "project.toml",
                           'schema_version = 1\nid = "fixture"\nversion = "1.0"\n'
                           '[automr]\npair = "C:G"\n')
            dataset = make_dataset(campaign / "PP_Hg", include_model=False)
            (dataset / "summary.html").unlink()
            plan = plan_campaign(campaign, preset=preset, frames_directory=base / "MR_frames")
            entry = plan["datasets"][0]
            self.assertEqual(entry["status"], "DISCOVERED", entry.get("diagnostic"))
            self.assertNotIn("summary", entry["inputs"])
            self.assertIn("metadata", entry["inputs"])
            self.assertEqual(campaign_status(campaign)["integrity"], "OK")


class AutoSolWithoutSummaryTests(unittest.TestCase):
    def use_cif_instead_of_summary(self, run: Path, root: Path) -> None:
        report = json.loads((run / "report.json").read_text())
        summary = Path(report["inputs"].pop("summary"))
        summary.unlink()
        report["inputs"]["summary"] = None
        report["inputs"]["metadata"] = str(write(root / "dataset" / "Data_1.cif", AUTOPROC_CIF))
        (run / "report.json").write_text(json.dumps(report))

    def test_metal_search_reads_wavelength_from_cif(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = test_autosol_metal.make_run(root)
            self.use_cif_instead_of_summary(run, root)
            result = execute_metal_search(
                run, test_autosol_metal.make_autosol(root, [("HG", (0.0, 0.0, 0.0))]),
                test_autosol.make_mtz_dump(root), ["HG"], environment=ENV, processor_count=8,
            )
            self.assertEqual(result.status, FOUND)
            args = json.loads((result.attempt_directory / "received_args.json").read_text())
            self.assertIn("autosol.lambda=1.00744", args)
            attempt = json.loads((run / "report.json").read_text())["metal_autosol"]["attempts"][0]
            self.assertTrue(attempt["wavelength_source"]["path"].endswith("Data_1.cif"))

    def test_nucleotide_autosol_reads_wavelength_from_cif(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = test_autosol.hex_cart(0.2, 0.1, 0.3)
            run = test_autosol.make_run(root, model)
            self.use_cif_instead_of_summary(run, root)
            result = execute_autosol(
                run, test_autosol.make_autosol(root, model), test_autosol.make_mtz_dump(root),
                environment=ENV, processor_count=8,
            )
            self.assertEqual(result.status, "AUTOSOL_READY", result.message)
            args = json.loads((result.autosol_directory / "received_args.json").read_text())
            self.assertIn("autosol.lambda=1.00744", args)
            payload = json.loads((run / "report.json").read_text())["autosol"]
            self.assertEqual(payload["wavelength"], 1.00744)


if __name__ == "__main__":
    unittest.main()
