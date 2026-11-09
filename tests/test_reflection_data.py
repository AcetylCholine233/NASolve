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
    merged_arrays,
    normalize_phaser_data,
    parse_mtz_dump_columns,
    reported_data_kind,
    select_labels,
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
    def test_staraniso_mtz_gives_f_sigf_and_imean_sigimean(self):
        columns = parse_mtz_dump_columns(mtz_dump_columns_text())
        self.assertEqual(len(columns), len(STARANISO_COLUMNS))
        self.assertEqual(columns[5].label, "IMEAN")
        self.assertEqual(columns[5].type, "J")
        self.assertEqual(select_labels(columns, "F"), "F,SIGF")
        self.assertEqual(select_labels(columns, "I"), "IMEAN,SIGIMEAN")

    def test_header_and_summary_lines_are_not_columns(self):
        text = (
            "Number of crystals: 2\n"
            "Number of Miller indices: 1635\n"
            + mtz_dump_columns_text((("H", "H"), ("I", "J"), ("SIGI", "Q")))
        )
        labels = [column.label for column in parse_mtz_dump_columns(text)]
        self.assertEqual(labels, ["H", "I", "SIGI"])

    def test_missing_requested_kind_is_rejected_with_the_columns(self):
        amplitudes_only = _columns(("H", "H"), ("F", "F"), ("SIGF", "Q"), ("FreeR_flag", "I"))
        with self.assertRaisesRegex(ReflectionDataError, r"type J.*F\(F\), SIGF\(Q\)"):
            select_labels(amplitudes_only, "I")
        intensities_only = _columns(("H", "H"), ("IMEAN", "J"), ("SIGIMEAN", "Q"))
        with self.assertRaisesRegex(ReflectionDataError, r"phaser_data = F needs .*type F"):
            select_labels(intensities_only, "F")

    def test_anomalous_columns_are_not_treated_as_merged(self):
        anomalous = _columns(("I(+)", "K"), ("SIGI(+)", "M"), ("I(-)", "K"), ("SIGI(-)", "M"),
                             ("F(+)", "G"), ("SIGF(+)", "L"), ("F(-)", "G"), ("SIGF(-)", "L"))
        self.assertEqual(merged_arrays(anomalous), [])
        mislabelled = _columns(("I(+)", "J"), ("SIGI(+)", "Q"), ("I(-)", "J"), ("SIGI(-)", "Q"))
        self.assertEqual(merged_arrays(mislabelled), [])

    def test_iotbx_fqdq_group_is_one_amplitude_array_with_isym(self):
        # ctruncate order: iotbx reads F SIGF DANO SIGDANO (+ISYM) as one array.
        columns = (
            MtzColumn(1, "F", "F"), MtzColumn(1, "SIGF", "Q"),
            MtzColumn(1, "DANO", "D"), MtzColumn(1, "SIGDANO", "Q"),
            MtzColumn(1, "ISYM", "Y", 1635, 0.0, 2.0),
            MtzColumn(1, "IMEAN", "J"), MtzColumn(1, "SIGIMEAN", "Q"),
        )
        self.assertEqual(
            merged_arrays(columns),
            [("F", "F,SIGF,DANO,SIGDANO,ISYM"), ("I", "IMEAN,SIGIMEAN")],
        )
        self.assertEqual(select_labels(columns, "F"), "F,SIGF,DANO,SIGDANO,ISYM")

    def test_amplitude_preference_follows_phenix_phaser(self):
        columns = _columns(("FC", "F"), ("SIGFC", "Q"), ("F_early", "F"), ("SIGF_early", "Q"),
                           ("FOBS", "F"), ("SIGFOBS", "Q"))
        self.assertEqual(select_labels(columns, "F"), "FOBS,SIGFOBS")
        self.assertEqual(select_labels(columns[:4], "F"), "F_early,SIGF_early")

    def test_columns_do_not_pair_across_datasets(self):
        columns = (MtzColumn(1, "I", "J"), MtzColumn(2, "SIGI", "Q"))
        self.assertEqual(merged_arrays(columns), [])

    def test_setting_values(self):
        for value in ("I", "i", "intensities", " I/SIGI "):
            self.assertEqual(normalize_phaser_data(value), "I")
        for value in ("F", "f", "amplitudes", "F,SIGF"):
            self.assertEqual(normalize_phaser_data(value), "F")
        self.assertEqual(normalize_phaser_data(None), "F")
        with self.assertRaisesRegex(ReflectionDataError, "F .*or I"):
            normalize_phaser_data("auto")

    def test_phaser_labin_echo_classification(self):
        self.assertEqual(reported_data_kind("I = I SIGI = SIGI"), "I")
        self.assertEqual(reported_data_kind("F = F SIGF = SIGF"), "F")
        self.assertIsNone(reported_data_kind(None))


class PhaserDataPipelineTests(unittest.TestCase):
    def _preflight(self, root: Path, **kwargs):
        dataset = make_dataset(root / "dataset")
        dump = make_mtz_dump(root / "tools", columns=kwargs.pop("columns", STARANISO_COLUMNS))
        return prepare_automr(
            dataset, valid_ligand_codes=VALID, mtz_dump_executable=dump, **kwargs
        )

    def test_intensities_reach_phaser_eff_as_labin(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root, phaser_data="I")
            report = json.loads(preflight.report_path.read_text())
            self.assertEqual(report["phaser_data"]["requested"], "I")
            self.assertEqual(report["phaser_data"]["labin"], "IMEAN,SIGIMEAN")
            snapshot = (preflight.run_directory / "nasolve.input.txt").read_text()
            self.assertIn("phaser_data = I\n", snapshot)
            self.assertIn("Phaser data: I - intensities (IMEAN,SIGIMEAN)",
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
            self.assertEqual((data["requested"], data["labin"]), ("I", "IMEAN,SIGIMEAN"))
            self.assertEqual(data["phaser_reported_kind"], "I")
            self.assertIn("Reflection data: I - intensities (IMEAN,SIGIMEAN)",
                          (result.run_directory / "automr.log").read_text())

    def test_default_is_an_explicit_f_choice(self):
        # Regression: the default left labin unset ("auto"), so Phaser chose.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root)
            report = json.loads(preflight.report_path.read_text())
            self.assertEqual(report["phaser_data"]["requested"], "F")
            self.assertEqual(report["phaser_data"]["labin"], "F,SIGF")
            snapshot = (preflight.run_directory / "nasolve.input.txt").read_text()
            self.assertIn("phaser_data = F\n", snapshot)
            result = execute_phaser(
                preflight.report_path, make_phaser(root, labin_echo="F = F SIGF = SIGF")
            )
            self.assertIn('  labin = "F,SIGF"', result.parameter_path.read_text())
            self.assertEqual(result.status, "MR_SUCCESS")

    def test_nasolve_txt_setting_is_read_and_cli_overrides_it(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dataset = make_dataset(root / "dataset")
            config = dataset / "nasolve.txt"
            config.write_text("[automr]\nmode = nonstandard\nmodel = search.pdb\nphaser_data = I\n")
            intent = read_intent(config)
            self.assertEqual(intent.phaser_data, "I")
            resolved = resolve_automr_input(dataset, intent, valid_ligand_codes=VALID)
            self.assertEqual(resolved.phaser_data, "I")
            self.assertIn("phaser_data = I\n", format_intent(resolved))
            overridden = resolve_automr_input(
                dataset, intent, valid_ligand_codes=VALID, phaser_data_override="F"
            )
            self.assertEqual(overridden.phaser_data, "F")
            self.assertIn("phaser_data = F\n", format_intent(overridden))
            self.assertEqual(
                resolve_automr_input(dataset, AutoMRIntent(), valid_ligand_codes=VALID).phaser_data,
                "F",
            )
            config.write_text("[automr]\nphaser_data = auto\n")
            with self.assertRaisesRegex(AutoMRInputError, "phaser_data"):
                read_intent(config)

    def test_preflight_stops_when_mtz_lacks_the_requested_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(AutoMRInputError, "merged intensity column"):
                self._preflight(
                    root, phaser_data="I",
                    columns=(("H", "H"), ("K", "H"), ("L", "H"), ("F", "F"), ("SIGF", "Q")),
                )
            with self.assertRaisesRegex(AutoMRInputError, "merged amplitude column"):
                self._preflight(
                    root, columns=(("H", "H"), ("K", "H"), ("L", "H"), ("I", "J"), ("SIGI", "Q")),
                )
            self.assertFalse((root / "dataset" / "AutoMR").exists())

    def test_preflight_needs_mtz_dump(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = make_dataset(Path(directory) / "dataset")
            with self.assertRaisesRegex(AutoMRInputError, "phenix.mtz.dump"):
                prepare_automr(dataset, valid_ligand_codes=VALID)
            self.assertFalse((dataset / "AutoMR").exists())

    def test_phaser_reporting_the_other_kind_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root / "i", phaser_data="I")
            result = execute_phaser(
                preflight.report_path, make_phaser(root / "i", labin_echo="F = F SIGF = SIGF")
            )
            self.assertEqual(result.status, "MR_FAILED")
            self.assertIn("Intensities (IMEAN,SIGIMEAN) were requested", result.message)
            self.assertIn("reported LABIN F = F SIGF = SIGF", result.message)

            preflight = self._preflight(root / "f", phaser_data="F")
            result = execute_phaser(
                preflight.report_path, make_phaser(root / "f", labin_echo="I = I SIGI = SIGI")
            )
            self.assertEqual(result.status, "MR_FAILED")
            self.assertIn("Amplitudes (F,SIGF) were requested", result.message)

    def test_older_reports_still_build_and_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            preflight = self._preflight(root, phaser_data="I")
            report = json.loads(preflight.report_path.read_text())
            phaser_directory = preflight.run_directory / "Phaser"
            legacy = dict(report)
            legacy.pop("phaser_data")
            self.assertNotIn("labin", build_phaser_eff(legacy, phaser_directory))
            schema_1_auto = {**report, "phaser_data": {
                "schema_version": 1, "requested": "auto", "labin": None}}
            self.assertNotIn("labin", build_phaser_eff(schema_1_auto, phaser_directory))
            schema_1_intensities = {**report, "phaser_data": {
                "schema_version": 1, "requested": "intensities", "labin": "IMEAN,SIGIMEAN"}}
            self.assertIn('labin = "IMEAN,SIGIMEAN"',
                          build_phaser_eff(schema_1_intensities, phaser_directory))
            report["phaser_data"]["labin"] = 'IMEAN" }\nphaser { mode = EP_AUTO'
            with self.assertRaisesRegex(PhaserExecutionError, "malformed"):
                build_phaser_eff(report, phaser_directory)
            report["phaser_data"].update(requested="auto", labin=None)
            with self.assertRaisesRegex(PhaserExecutionError, "malformed"):
                build_phaser_eff(report, phaser_directory)


if __name__ == "__main__":
    unittest.main()
