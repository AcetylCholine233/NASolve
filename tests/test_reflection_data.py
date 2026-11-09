import json
import tempfile
import unittest
from pathlib import Path

from nasolve.automr import AutoMRInputError, prepare_automr
from nasolve.automr_input import AutoMRIntent, format_intent, read_intent, resolve_automr_input
from nasolve.phaser import PhaserExecutionError, build_phaser_eff, execute_phaser
from nasolve.reflection_data import (
    MtzColumn,
    ReflectionDataError,
    mean_intensity_label_strings,
    normalize_phaser_data,
    parse_mtz_dump_columns,
    reported_data_kind,
    select_intensity_labels,
)

from .helpers import (
    STARANISO_COLUMNS,
    make_dataset,
    make_mtz_dump,
    make_phaser,
    mtz_dump_columns_text,
)


VALID = {"DA", "DC"}


def _columns(*pairs):
    return tuple(MtzColumn(1, label, kind) for label, kind in pairs)


class ColumnSelectionTests(unittest.TestCase):
    def test_staraniso_mtz_selects_imean_not_amplitudes(self):
        columns = parse_mtz_dump_columns(mtz_dump_columns_text())
        self.assertEqual(len(columns), len(STARANISO_COLUMNS))
        self.assertEqual(columns[5], MtzColumn(1, "IMEAN", "J"))
        self.assertEqual(select_intensity_labels(columns), "IMEAN,SIGIMEAN")

    def test_header_and_summary_lines_are_not_columns(self):
        text = (
            "Number of crystals: 2\n"
            "Number of Miller indices: 1635\n"
            + mtz_dump_columns_text((("H", "H"), ("I", "J"), ("SIGI", "Q")))
        )
        labels = [column.label for column in parse_mtz_dump_columns(text)]
        self.assertEqual(labels, ["H", "I", "SIGI"])

    def test_amplitude_only_mtz_is_rejected_with_its_columns(self):
        columns = _columns(("H", "H"), ("F", "F"), ("SIGF", "Q"), ("FreeR_flag", "I"))
        with self.assertRaisesRegex(ReflectionDataError, r"type J.*F\(F\), SIGF\(Q\)"):
            select_intensity_labels(columns)

    def test_anomalous_only_intensities_are_not_treated_as_merged(self):
        anomalous = _columns(("I(+)", "K"), ("SIGI(+)", "M"), ("I(-)", "K"), ("SIGI(-)", "M"))
        self.assertEqual(mean_intensity_label_strings(anomalous), [])
        mislabelled = _columns(("I(+)", "J"), ("SIGI(+)", "Q"), ("I(-)", "J"), ("SIGI(-)", "Q"))
        self.assertEqual(mean_intensity_label_strings(mislabelled), [])

    def test_iotbx_fqdq_group_does_not_lend_its_sigma_to_an_intensity(self):
        columns = _columns(
            ("F", "F"), ("SIGF", "Q"), ("DANO", "D"), ("SIGDANO", "Q"),
            ("I", "J"), ("SIGI", "Q"),
        )
        self.assertEqual(mean_intensity_label_strings(columns), ["I,SIGI"])

    def test_columns_do_not_pair_across_datasets(self):
        columns = (MtzColumn(1, "I", "J"), MtzColumn(2, "SIGI", "Q"))
        self.assertEqual(mean_intensity_label_strings(columns), [])

    def test_setting_aliases(self):
        for value in ("intensities", "I", "isigi", " I/SIGI "):
            self.assertEqual(normalize_phaser_data(value), "intensities")
        self.assertEqual(normalize_phaser_data(None), "auto")
        with self.assertRaises(ReflectionDataError):
            normalize_phaser_data("F")

    def test_phaser_labin_echo_classification(self):
        self.assertEqual(reported_data_kind("I = I SIGI = SIGI"), "intensities")
        self.assertEqual(reported_data_kind("F = F SIGF = SIGF"), "amplitudes")
        self.assertIsNone(reported_data_kind(None))


class PhaserDataPipelineTests(unittest.TestCase):
    def _preflight(self, root: Path, **kwargs):
        dataset = make_dataset(root / "dataset")
        dump = make_mtz_dump(root, columns=kwargs.pop("columns", STARANISO_COLUMNS))
        return prepare_automr(
            dataset, valid_ligand_codes=VALID, mtz_dump_executable=dump, **kwargs
        )

    def test_intensities_reach_phaser_eff_as_labin(self):
        # Regression: phaser.eff had no labin, so Phenix always chose F,SIGF.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root, phaser_data="intensities")
            report = json.loads(preflight.report_path.read_text())
            self.assertEqual(report["phaser_data"]["requested"], "intensities")
            self.assertEqual(report["phaser_data"]["labin"], "IMEAN,SIGIMEAN")
            snapshot = (preflight.run_directory / "nasolve.input.txt").read_text()
            self.assertIn("phaser_data = intensities", snapshot)
            self.assertIn("Phaser data: intensities (IMEAN,SIGIMEAN)",
                          (preflight.run_directory / "automr.log").read_text())

            result = execute_phaser(
                preflight.report_path,
                make_phaser(root, labin_echo="I = I SIGI = SIGI"),
                phenix_version="2.2.1",
            )
            eff = result.parameter_path.read_text()
            self.assertIn('  labin = "IMEAN,SIGIMEAN"', eff)
            self.assertLess(eff.index("hklin ="), eff.index("labin ="))
            self.assertEqual(result.status, "MR_SUCCESS")
            executed = json.loads(result.report_path.read_text())
            data = executed["execution"]["phaser"]["reflection_data"]
            self.assertEqual(data["labin"], "IMEAN,SIGIMEAN")
            self.assertEqual(data["phaser_reported_kind"], "intensities")
            self.assertIn("Reflection data: intensities (IMEAN,SIGIMEAN)",
                          (result.run_directory / "automr.log").read_text())

    def test_default_keeps_phenix_choice_and_writes_no_labin(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root)
            report = json.loads(preflight.report_path.read_text())
            self.assertEqual(report["phaser_data"]["requested"], "auto")
            self.assertIsNone(report["phaser_data"]["labin"])
            snapshot = (preflight.run_directory / "nasolve.input.txt").read_text()
            self.assertNotIn("phaser_data", snapshot)
            result = execute_phaser(
                preflight.report_path, make_phaser(root, labin_echo="F = F SIGF = SIGF")
            )
            self.assertNotIn("labin", result.parameter_path.read_text())
            self.assertEqual(result.status, "MR_SUCCESS")

    def test_nasolve_txt_setting_is_read_and_cli_overrides_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset")
            config = dataset / "nasolve.txt"
            config.write_text("[automr]\nmode = nonstandard\nmodel = search.pdb\nphaser_data = I\n")
            intent = read_intent(config)
            self.assertEqual(intent.phaser_data, "intensities")
            resolved = resolve_automr_input(dataset, intent, valid_ligand_codes=VALID)
            self.assertEqual(resolved.phaser_data, "intensities")
            self.assertIn("phaser_data = intensities", format_intent(resolved))
            overridden = resolve_automr_input(
                dataset, intent, valid_ligand_codes=VALID, phaser_data_override="auto"
            )
            self.assertEqual(overridden.phaser_data, "auto")
            self.assertNotIn("phaser_data", format_intent(overridden))
            self.assertEqual(
                resolve_automr_input(dataset, AutoMRIntent(), valid_ligand_codes=VALID).phaser_data,
                "auto",
            )
            config.write_text("[automr]\nphaser_data = F\n")
            with self.assertRaisesRegex(AutoMRInputError, "phaser_data"):
                read_intent(config)

    def test_preflight_stops_when_mtz_has_no_merged_intensities(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(AutoMRInputError, "merged intensity column"):
                self._preflight(
                    root, phaser_data="intensities",
                    columns=(("H", "H"), ("K", "H"), ("L", "H"), ("F", "F"), ("SIGF", "Q")),
                )
            self.assertFalse((root / "dataset" / "AutoMR").exists())

    def test_intensities_need_mtz_dump(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = make_dataset(Path(directory) / "dataset")
            with self.assertRaisesRegex(AutoMRInputError, "phenix.mtz.dump"):
                prepare_automr(dataset, valid_ligand_codes=VALID, phaser_data="I")

    def test_phaser_reporting_amplitudes_after_intensity_request_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root, phaser_data="intensities")
            result = execute_phaser(
                preflight.report_path, make_phaser(root, labin_echo="F = F SIGF = SIGF")
            )
            self.assertEqual(result.status, "MR_FAILED")
            self.assertIn("reported LABIN F = F SIGF = SIGF", result.message)

    def test_reports_without_phaser_data_still_run_and_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root, phaser_data="intensities")
            report = json.loads(preflight.report_path.read_text())
            phaser_directory = preflight.run_directory / "Phaser"
            legacy = dict(report)
            legacy.pop("phaser_data")
            self.assertNotIn("labin", build_phaser_eff(legacy, phaser_directory))
            report["phaser_data"]["labin"] = 'IMEAN" }\nphaser { mode = EP_AUTO'
            with self.assertRaisesRegex(PhaserExecutionError, "malformed"):
                build_phaser_eff(report, phaser_directory)


if __name__ == "__main__":
    unittest.main()
