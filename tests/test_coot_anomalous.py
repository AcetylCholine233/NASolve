"""Exercise the Coot startup script with the documented Python map API.

The fake follows the Coot bundled with Phenix 2.2.1 (``coot.py`` SWIG
wrappers): ``valid_labels`` takes ``bool use_weights_flag`` and rejects an
int, and phenix.refine writes the anomalous map as ANOM/PANOM (phenix.maps
writes ANOM/PHANOM). Before 2026-11-23 the script passed 0 and looked only for
PHANOM, so the anomalous map never opened (U-Y_Ag run_004, every round).
The anomalous map opens as an ordinary map, not a difference map (user,
2026-11-23).
"""

import io
import json
import shutil
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from nasolve.coot_view import CootViewError, launch_coot_view, resolve_view_profile
from nasolve.model_assessment import file_sha256
from nasolve.run_context import artifact_reference

from .test_coot_view import add_manual_checkpoint, add_review_checkpoint, make_view_run


class FakeCoot:
    def __init__(self, mtz, labels=("ANOM", "PANOM")):
        self.labels = {str(mtz): set(labels)}
        self.maps = [None, (str(mtz), "2FOFCWT", "PH2FOFCWT", 0),
                     (str(mtz), "FOFCWT", "PHFOFCWT", 1)]
        self.created = []
        self.contours = {}
        self.names = {}
        self.refinement = 1
        self.scroll = 1

    def valid_labels(self, mtz, amplitude, phase, weight, use_weights):
        if not isinstance(use_weights, bool):   # SWIG: bool use_weights_flag
            raise TypeError("in method 'valid_labels', argument 5 of type 'bool'")
        assert (amplitude, phase, weight, use_weights) in (
            ("ANOM", "PANOM", "", False), ("ANOM", "PHANOM", "", False))
        return int({amplitude, phase} <= self.labels.get(mtz, set()))

    def graphics_n_molecules(self):
        return len(self.maps)

    def mtz_hklin_for_map(self, imol):
        return self.maps[imol][0] if self.maps[imol] else None

    def mtz_fp_for_map(self, imol):
        return self.maps[imol][1]

    def mtz_phi_for_map(self, imol):
        return self.maps[imol][2]

    def mtz_use_weight_for_map(self, imol):
        return 0

    def make_and_draw_map(self, mtz, amplitude, phase, weight, use_weights, difference):
        self.created.append((mtz, amplitude, phase, weight, use_weights, difference))
        self.maps.append((mtz, amplitude, phase, difference))
        self.refinement = self.scroll = len(self.maps) - 1
        return self.refinement

    def set_map_is_difference_map(self, imol, difference):
        self.maps[imol] = (*self.maps[imol][:3], difference)

    def set_contour_level_in_sigma(self, imol, level):
        self.contours[imol] = level

    def set_molecule_name(self, imol, name):
        self.names[imol] = name

    def imol_refinement_map(self):
        return self.refinement

    def set_imol_refinement_map(self, imol):
        self.refinement = imol

    def scroll_wheel_map(self):
        return self.scroll

    def set_scroll_wheel_map(self, imol):
        self.scroll = imol


class CootAnomalousTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.run = make_view_run(self.root)
        self.mtz = self.run / "AutoRefine" / "round_001" / "refined_001_map_coeffs.mtz"

    def launch(self, **kwargs):
        result = launch_coot_view(
            self.run, self.root / "coot",
            launcher=lambda *args, **kwargs: SimpleNamespace(pid=1234), **kwargs,
        )
        command = result.command
        self.assertEqual(command[command.index("--auto") + 1], str(result.map_path))
        self.assertIn("--script", command)
        script = Path(command[command.index("--script") + 1])
        self.assertEqual(script.parent, result.working_directory)
        return script, result

    def execute_script(self, script, coot):
        output = io.StringIO()
        with patch.dict("sys.modules", {"coot": coot}), redirect_stdout(output):
            exec(compile(script.read_text(), str(script), "exec"), {})
        return output.getvalue()

    def add_reflections(self, identifier="refine-001"):
        registry_path = self.run / "AutoRefine" / "checkpoints.json"
        registry = json.loads(registry_path.read_text())
        checkpoint = next(c for c in registry["checkpoints"] if c["id"] == identifier)
        path = Path(checkpoint["model"]["path"]).with_suffix(".mtz")
        path.write_bytes(b"refinement reflections including anomalous coefficients")
        reference = artifact_reference(path, self.run)
        reference["sha256"] = file_sha256(path)
        checkpoint["outputs"]["refinement_reflections"] = reference
        registry_path.write_text(json.dumps(registry))
        return path

    def test_show_adds_anomalous_map_without_changing_normal_maps_or_controls(self):
        registry = self.run / "AutoRefine" / "checkpoints.json"
        before_registry = registry.read_bytes()
        before_mtz = self.mtz.read_bytes()
        script, result = self.launch()
        coot = FakeCoot(self.mtz)
        normal_maps = coot.maps[:]
        output = self.execute_script(script, coot)
        self.assertEqual(coot.created, [(str(self.mtz), "ANOM", "PANOM", "", 0, 0)])
        self.assertEqual(coot.maps[3][3], 0)              # not a difference map
        self.assertEqual(coot.maps[:3], normal_maps)
        self.assertEqual(coot.contours, {3: 3.0})
        self.assertEqual(coot.names[3], "Anomalous (ANOM/PANOM)")
        self.assertEqual((coot.refinement, coot.scroll), (1, 1))
        self.assertIn('"status": "loaded"', output)
        self.assertIn('"labels": ["ANOM", "PANOM"]', output)
        launch = json.loads((result.working_directory / "launch.json").read_text())
        self.assertEqual(launch["anomalous_map"]["labels"], [["ANOM", "PANOM"], ["ANOM", "PHANOM"]])
        self.assertEqual(launch["anomalous_map"]["status"], "check-in-coot")
        self.assertEqual(registry.read_bytes(), before_registry)
        self.assertEqual(self.mtz.read_bytes(), before_mtz)

    def test_absent_or_incomplete_pair_does_not_construct_a_map(self):
        script, _ = self.launch()
        for labels in ((), ("ANOM",), ("PANOM",), ("PHANOM",), ("HLA", "HLB", "HLC", "HLD")):
            with self.subTest(labels=labels):
                coot = FakeCoot(self.mtz, labels)
                before = coot.maps[:]
                output = self.execute_script(script, coot)
                self.assertEqual(coot.created, [])
                self.assertEqual(coot.maps, before)
                self.assertEqual((coot.refinement, coot.scroll), (1, 1))
                self.assertIn('"status": "absent"', output)

    def test_phenix_maps_labels_still_open(self):
        # phenix.maps (and Phenix's own Coot GUI) use ANOM/PHANOM.
        script, _ = self.launch()
        coot = FakeCoot(self.mtz, ("ANOM", "PHANOM"))
        output = self.execute_script(script, coot)
        self.assertEqual(coot.created, [(str(self.mtz), "ANOM", "PHANOM", "", 0, 0)])
        self.assertIn("ANOM/PHANOM", coot.names[3])
        self.assertIn('"labels": ["ANOM", "PHANOM"]', output)
        # A file with both: phenix.refine's pair first.
        coot = FakeCoot(self.mtz, ("ANOM", "PANOM", "PHANOM"))
        self.execute_script(script, coot)
        self.assertEqual([entry[2] for entry in coot.created], ["PANOM"])

    def test_existing_anomalous_overlay_is_reused_including_qualified_labels(self):
        script, _ = self.launch()
        coot = FakeCoot(self.mtz)
        coot.maps.append((str(self.mtz), "/crystal/dataset/ANOM",
                          "/crystal/dataset/PANOM", 1))
        for _ in range(2):
            output = self.execute_script(script, coot)
            self.assertIn('"status": "reused"', output)
        self.assertEqual(coot.created, [])
        self.assertEqual(coot.maps[3][3], 0)   # an earlier difference-map copy is switched
        self.assertEqual(coot.contours, {3: 3.0})
        self.assertEqual((coot.refinement, coot.scroll), (1, 1))

    def test_identical_labels_from_other_file_do_not_satisfy_selected_map(self):
        script, _ = self.launch()
        coot = FakeCoot(self.mtz)
        coot.maps.append((str(self.root / self.mtz.name), "ANOM", "PANOM", 1))
        self.execute_script(script, coot)
        self.assertEqual(len(coot.created), 1)
        self.assertEqual(coot.created[0][0], str(self.mtz))

    def test_map_index_zero_is_valid(self):
        script, _ = self.launch()
        coot = FakeCoot(self.mtz)
        coot.maps = []
        coot.refinement = coot.scroll = -1
        output = self.execute_script(script, coot)
        self.assertIn('"status": "loaded"', output)
        self.assertEqual(coot.contours, {0: 3.0})

    def test_optional_map_errors_are_logged_without_closing_ordinary_maps(self):
        script, _ = self.launch()
        for failure in ("reader", "construction", "display", "missing-api"):
            with self.subTest(failure=failure):
                coot = FakeCoot(self.mtz)
                before = coot.maps[:]
                if failure == "reader":
                    coot.valid_labels = Mock(side_effect=ValueError("bad MTZ"))
                elif failure == "construction":
                    coot.make_and_draw_map = lambda *args: -1
                elif failure == "display":
                    coot.set_contour_level_in_sigma = Mock(
                        side_effect=RuntimeError("display failed"))
                else:
                    coot.valid_labels = None
                output = self.execute_script(script, coot)
                self.assertIn('"status": "error"', output)
                self.assertEqual(coot.maps[:3], before)
                self.assertEqual((coot.refinement, coot.scroll), (1, 1))

    def test_separate_map_file_falls_back_to_same_checkpoint_reflections(self):
        reflections = self.add_reflections()
        script, result = self.launch()
        coot = FakeCoot(self.mtz, ())
        coot.labels[str(reflections)] = {"ANOM", "PANOM"}
        self.execute_script(script, coot)
        self.assertEqual(coot.created[0][0], str(reflections))
        self.assertEqual(result.map_path, self.mtz)
        # If both files contain the pair, the selected map file has priority.
        coot = FakeCoot(self.mtz)
        coot.labels[str(reflections)] = {"ANOM", "PANOM"}
        self.execute_script(script, coot)
        self.assertEqual([entry[0] for entry in coot.created], [str(self.mtz)])

    def test_explicit_checkpoint_does_not_borrow_anomalous_maps_from_other_rounds(self):
        reflections = self.add_reflections()
        _, review_map = add_review_checkpoint(self.run)
        script, _ = self.launch(checkpoint="refine-002")
        coot = FakeCoot(review_map, ())
        coot.labels[str(self.mtz)] = coot.labels[str(reflections)] = {"ANOM", "PANOM"}
        output = self.execute_script(script, coot)
        self.assertEqual(coot.created, [])
        self.assertIn('"status": "absent"', output)
        self.assertEqual(json.loads(
            (self.run / "AutoRefine" / "checkpoints.json").read_text()
        )["current"], "refine-001")

    def test_manual_checkpoint_inherits_anomalous_source_after_relocation(self):
        reflections = self.add_reflections()
        add_manual_checkpoint(self.run)
        add_review_checkpoint(self.run)
        old_run = self.run
        relocated = self.root / "relocated" / "dataset"
        shutil.copytree(old_run.parent.parent, relocated)
        shutil.rmtree(old_run.parent.parent)
        self.run = relocated / "AutoMR" / old_run.name
        expected = self.run / reflections.relative_to(old_run)
        script, result = self.launch()
        coot = FakeCoot(result.map_path, ())
        coot.labels[str(expected)] = {"ANOM", "PANOM"}
        self.execute_script(script, coot)
        self.assertEqual(coot.created[0][0], str(expected))
        self.assertEqual(result.checkpoint_id, "manual-001")
        self.assertIn("not recalculated for manual model", result.map_source)
        self.assertNotIn(str(old_run), script.read_text())

    def test_declared_fallback_reflections_require_integrity(self):
        reflections = self.add_reflections()
        original = reflections.read_bytes()
        for mutation in ("corrupt", "missing"):
            with self.subTest(mutation=mutation):
                reflections.write_bytes(original)
                if mutation == "corrupt":
                    reflections.write_bytes(b"changed")
                else:
                    reflections.unlink()
                with self.assertRaisesRegex(CootViewError, "missing or failed checksum"):
                    resolve_view_profile(self.run)

    def test_schema_two_rejects_unchecksummed_fallback(self):
        self.add_reflections()
        path = self.run / "AutoRefine" / "checkpoints.json"
        registry = json.loads(path.read_text())
        registry["schema_version"] = 2
        checkpoint = registry["checkpoints"][0]
        model = self.mtz.with_name("refined_001.pdb")
        checkpoint["model"] = artifact_reference(model, self.run)
        checkpoint["model"]["sha256"] = file_sha256(model)
        map_ref = artifact_reference(self.mtz, self.run)
        map_ref["sha256"] = file_sha256(self.mtz)
        checkpoint["outputs"]["map_coefficients"] = map_ref
        checkpoint["outputs"]["refinement_reflections"].pop("sha256")
        path.write_text(json.dumps(registry))
        with self.assertRaisesRegex(CootViewError, "refinement reflections checksum"):
            resolve_view_profile(self.run)

    def test_startup_script_refuses_symlinks(self):
        script, _ = self.launch()
        script.unlink()
        outside = self.root / "other.py"
        outside.write_text("unchanged\n")
        script.symlink_to(outside)
        with self.assertRaisesRegex(CootViewError, "symlink"):
            self.launch()
        self.assertEqual(outside.read_text(), "unchanged\n")

    def test_stage_views_only_check_their_selected_map(self):
        for stage in ("automr", "postmr", "autosol"):
            with self.subTest(stage=stage):
                script, result = self.launch(stage=stage)
                coot = FakeCoot(result.map_path)
                self.execute_script(script, coot)
                self.assertEqual([entry[0] for entry in coot.created], [str(result.map_path)])
