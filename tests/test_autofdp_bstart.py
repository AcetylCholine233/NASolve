"""AutoFDP step 1: copy an AutoRefine result and set each metal's B from the
nucleotides around it (NARestraints metal-pair search). Pure Python; no
Phenix program runs in this step."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from nasolve.autofdp import AutoFDPError, prepare_start_model, surrounding_b
from nasolve.checkpoints import append_checkpoint, initialize_registry, list_checkpoints
from nasolve.cli import main
from nasolve.model_assessment import file_sha256
from nasolve.run_context import artifact_reference

from .helpers import pdb_record
from .test_autorefine import make_refine_run


def _residue(serial, name, chain, number, b, atoms):
    return [
        pdb_record("ATOM", serial + index, atom, name, chain, number, element=atom[0],
                   x=x, y=y, z=z, b_factor=b)
        for index, (atom, (x, y, z)) in enumerate(atoms)
    ]


def bridged_model_text(metal_b=379.99, occupancy=0.9) -> str:
    """Ag at the origin bridging A:12 DC (N3 2.3 A, B 60) and B:4 DU (O4 2.2 A,
    B 40). B:3 DA (B 10) also touches the Ag with N6 at 3.0 A but is the chain
    neighbour of B:4, so the pair is A:12-B:4. A:20 DG (B 5) is far away."""
    lines = []
    lines += _residue(1, "DC", "A", 12, 60.0, [
        ("P", (6.0, 2.0, 0.0)), ("C1'", (5.0, 1.0, 0.0)), ("O2", (3.0, 1.5, 0.0)),
        ("N3", (2.3, 0.0, 0.0)),
    ])
    lines += _residue(5, "DC", "A", 20, 5.0, [
        ("P", (20.0, 20.0, 20.0)), ("C1'", (21.0, 20.0, 20.0)),
    ])
    lines += _residue(7, "DA", "B", 3, 10.0, [
        ("P", (0.0, 6.0, 0.0)), ("C1'", (0.0, 5.0, 0.0)), ("N6", (0.0, 3.0, 0.0)),
    ])
    lines += _residue(10, "DU", "B", 4, 40.0, [
        ("P", (-6.0, -2.0, 0.0)), ("C1'", (-5.0, -1.0, 0.0)), ("O2", (-4.5, 0.5, 0.0)),
        ("O4", (-2.2, 0.0, 0.0)),
    ])
    lines.append(pdb_record("HETATM", 14, "AG", "AG", "B", 101, occupancy=occupancy,
                            element="AG", b_factor=metal_b))
    return "".join(lines) + "END\n"


def make_refined_run(root: Path, model_text: str = None) -> Path:
    """A PostMR run plus one AutoRefine checkpoint, refine-001, whose result
    is ``model_text``. The current checkpoint stays postmr (a REVIEW round)."""
    run = make_refine_run(root)
    run, registry = initialize_registry(run)
    round_dir = run / "AutoRefine" / "round_001"
    round_dir.mkdir(parents=True)
    model = round_dir / "refined_001.pdb"
    model.write_text(model_text or bridged_model_text())
    root_item = registry["checkpoints"][0]
    append_checkpoint(run, registry, {
        "id": "refine-001",
        "parent": "postmr",
        "kind": "refinement",
        "status": "REVIEW",
        "usable": True,
        "recipe": "AutoRefine/default",
        "label": None,
        "model": {**artifact_reference(model, run), "sha256": file_sha256(model),
                  "size": model.stat().st_size},
        "observations": root_item.get("observations"),
        "phases": root_item.get("phases"),
        "restraints": root_item.get("restraints", []),
        "metrics": {},
    })
    return run


def _b(line: str) -> float:
    return float(line[60:66])


class BStartTests(unittest.TestCase):
    def test_metal_b_set_to_mean_of_bridged_pair_in_a_new_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            run = make_refined_run(Path(directory))
            source = run / "AutoRefine" / "round_001" / "refined_001.pdb"
            registry_before = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            result = prepare_start_model(run, from_checkpoint="refine-001")

            (site,) = result.sites
            self.assertEqual((site.site, site.atom_name, site.element), ("B:101", "AG", "AG"))
            self.assertEqual(site.b_before, 379.99)
            self.assertEqual(site.b_after, 50.0)   # 4 atoms at 60, 4 at 40; B:3 not used
            self.assertEqual(site.occupancy, 0.9)  # occupancy is not touched
            self.assertIn("bridged base pair", site.neighbour_b.method)
            self.assertEqual(
                [(n.residue, n.residue_name, n.mean_b) for n in site.neighbour_b.neighbours],
                [("B:4", "DU", 40.0), ("A:12", "DC", 60.0)],
            )

            # The source copy is byte-identical; only the metal's B column changed.
            self.assertEqual(result.output_directory, run / "AutoFDP" / "fdp_001")
            self.assertEqual(result.source_copy.read_bytes(), source.read_bytes())
            before = source.read_text().splitlines()
            after = result.start_model.read_text().splitlines()
            changed = [i for i, (a, b) in enumerate(zip(before, after)) if a != b]
            self.assertEqual(len(before), len(after))
            self.assertEqual(len(changed), 1)
            line = after[changed[0]]
            self.assertTrue(line.startswith("HETATM") and " AG " in line)
            self.assertEqual(_b(line), 50.0)
            self.assertEqual(line[:60], before[changed[0]][:60])
            self.assertEqual(line[66:], before[changed[0]][66:])

            # Registered as a model-edit checkpoint; current unchanged.
            records, current, _ = list_checkpoints(run)
            self.assertEqual(current, "postmr")
            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            item = next(c for c in registry["checkpoints"] if c["id"] == "autofdp-001")
            self.assertEqual(result.checkpoint_id, "autofdp-001")
            self.assertEqual(
                (item["kind"], item["parent"], item["recipe"], item["status"], item["usable"]),
                ("manual", "refine-001", "AutoFDP/b-start", "REVIEW", True),
            )
            self.assertEqual(item["model"]["relative_path"], "AutoFDP/fdp_001/start_model.pdb")
            self.assertEqual(item["model"]["sha256"], file_sha256(result.start_model))
            self.assertEqual(item["observations"], registry_before["checkpoints"][0]["observations"])
            self.assertEqual(item["autofdp"]["edits"][0]["b_after"], 50.0)

            report = json.loads(result.report_path.read_text())
            self.assertEqual(report["workflow"], "autofdp")
            self.assertEqual(report["status"], "AUTOFDP_PREPARED")
            self.assertEqual(report["source_checkpoint"], "refine-001")
            self.assertEqual(report["start_checkpoint"], "autofdp-001")
            self.assertEqual(report["source_model"]["sha256"], file_sha256(source))
            self.assertEqual(report["sites"][0]["neighbour_mean_b"], 50.0)
            self.assertIn("B 379.99 -> 50.00", (result.output_directory / "autofdp.log").read_text())

            # A second run gets its own numbered folder and checkpoint.
            again = prepare_start_model(run, from_checkpoint="refine-001")
            self.assertEqual(again.output_directory.name, "fdp_002")
            self.assertEqual(again.checkpoint_id, "autofdp-002")

    def test_start_must_be_an_autorefine_result(self):
        # Regression for the 7CU_Ag case: REVIEW rounds are not selected, so
        # the current checkpoint is still postmr.
        with tempfile.TemporaryDirectory() as directory:
            run = make_refined_run(Path(directory))
            with self.assertRaisesRegex(AutoFDPError, r"postmr checkpoint.*--from \(refine-001\)"):
                prepare_start_model(run)
            self.assertFalse((run / "AutoFDP").exists())

    def test_plan_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            run = make_refined_run(Path(directory))
            registry = (run / "AutoRefine" / "checkpoints.json").read_bytes()
            result = prepare_start_model(run, from_checkpoint="refine-001", plan_only=True)
            self.assertEqual(result.sites[0].b_after, 50.0)
            self.assertIsNone(result.start_model)
            self.assertFalse((run / "AutoFDP").exists())
            self.assertEqual((run / "AutoRefine" / "checkpoints.json").read_bytes(), registry)

    def test_site_selection_and_other_b_modes(self):
        with tempfile.TemporaryDirectory() as directory:
            run = make_refined_run(Path(directory))
            with self.assertRaisesRegex(AutoFDPError, "No metal site B:102.*metal sites: B:101"):
                prepare_start_model(run, from_checkpoint="refine-001", sites=["B:102"], plan_only=True)
            kept = prepare_start_model(run, from_checkpoint="refine-001", b_start="model", plan_only=True)
            self.assertEqual(kept.sites[0].b_after, 379.99)
            self.assertIsNone(kept.sites[0].neighbour_b)
            fixed = prepare_start_model(run, from_checkpoint="refine-001", b_start="123.4")
            self.assertEqual(_b(next(line for line in fixed.start_model.read_text().splitlines()
                                     if line.startswith("HETATM"))), 123.4)
            for bad in ("hot", "0", "1000"):
                with self.subTest(b_start=bad), self.assertRaises(AutoFDPError):
                    prepare_start_model(run, from_checkpoint="refine-001", b_start=bad, plan_only=True)

    def test_model_without_a_metal_stops(self):
        with tempfile.TemporaryDirectory() as directory:
            text = "".join(l for l in bridged_model_text().splitlines(keepends=True)
                           if not l.startswith("HETATM"))
            run = make_refined_run(Path(directory), text)
            with self.assertRaisesRegex(AutoFDPError, "No metal HETATM"):
                prepare_start_model(run, from_checkpoint="refine-001")
            self.assertFalse((run / "AutoFDP").exists())


class NeighbourFallbackTests(unittest.TestCase):
    def _model(self, root: Path, metal_xyz) -> tuple[Path, int]:
        text = "".join(_residue(1, "DT", "A", 5, 30.0, [
            ("C1'", (5.0, 0.0, 0.0)), ("O4", (2.5, 0.0, 0.0)), ("P", (7.0, 0.0, 0.0)),
        ]) + _residue(4, "DA", "A", 9, 70.0, [
            ("C1'", (0.0, 0.0, 5.5)), ("P", (0.0, 0.0, 9.0)),
        ]))
        x, y, z = metal_xyz
        text += pdb_record("HETATM", 9, "HG", "HG", "C", 1, element="HG", x=x, y=y, z=z, b_factor=200.0)
        path = root / "model.pdb"
        path.write_text(text + "END\n")
        return path, text.count("\n")

    def test_single_base_contact_then_nearby_nucleotides_then_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path, line = self._model(root, (0.0, 0.0, 0.0))
            single = surrounding_b(path, line)   # O4 of A:5 at 2.5 A; A:9 has no base atom
            self.assertIn("base atoms within 3.6 A", single.method)
            self.assertEqual([n.residue for n in single.neighbours], ["A:5"])
            self.assertEqual(single.b_factor, 30.0)

            path, line = self._model(root, (0.0, 0.0, 3.0))   # O4 3.9 A: no base contact
            nearby = surrounding_b(path, line)
            self.assertIn("nucleotides within 6 A", nearby.method)
            self.assertEqual([n.residue for n in nearby.neighbours], ["A:9", "A:5"])
            self.assertEqual(nearby.b_factor, 46.0)   # atom-weighted: 3 x 30, 2 x 70

            path, line = self._model(root, (40.0, 40.0, 40.0))
            with self.assertRaisesRegex(AutoFDPError, "No nucleotide within 6 A"):
                surrounding_b(path, line)


class AutoFDPCommandTests(unittest.TestCase):
    def test_command_runs_step_one_and_prints_the_change(self):
        with tempfile.TemporaryDirectory() as directory:
            run = make_refined_run(Path(directory))
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main(["autofdp", str(run), "--from", "refine-001"])
            self.assertEqual(code, 0, err.getvalue())
            text = out.getvalue()
            self.assertIn("B:101 AG (AG), occupancy 0.90: B 379.99 -> 50.00", text)
            self.assertIn("B:4 DU  mean B 40.00", text)
            self.assertIn("Checkpoint: autofdp-001 (current checkpoint unchanged)", text)
            self.assertIn("--checkpoint autofdp-001", text)
            self.assertIn("not implemented yet", text)

    def test_command_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            run = make_refined_run(Path(directory))
            err = io.StringIO()
            with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["autofdp", str(run)]), 1)
                self.assertEqual(main(["autofdp", str(run), "--from", "refine-001",
                                       "--element", "B:101=HG"]), 2)
            self.assertIn("not an AutoRefine result", err.getvalue())
            self.assertIn("--element is not implemented yet", err.getvalue())
            self.assertFalse((run / "AutoFDP").exists())


if __name__ == "__main__":
    unittest.main()
