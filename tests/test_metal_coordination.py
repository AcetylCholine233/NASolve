"""Metal coordination restraints: ``nasolve autorefine --coordination``.

User's specification (2026-11-23): one bond per donor with the length from a
lab table (Ag-O 2.5, Ag-N 2.3, Hg-O 2.1, Hg-N 2.1 A), sigma 0.2 A; for linear
coordination (two donors) the angle donor-donor-metal 0 degrees, sigma 3 (not
a physical angle); no preset angle for more donors. Written as phenix.refine
``geometry_restraints.edits`` (the form checked against Phenix 2.2.1's own
definition on 2026-11-23). phenix.refine is simulated.

``--geometry`` (user, 2026-11-24): a preset coordination geometry; for now
``sqp``, square planar: four donors listed in order around the square, and
for each donor-metal-donor angle between neighbours (the last with the first)
90 degrees, sigma 6; and (user, same day) a planarity restraint over the
metal and its four donors (sigma 0.2 A, NASolve's choice: the user gave none).
Without --geometry nothing changes.
"""

import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nasolve.autorefine import AutoRefineError, execute_autorefine
from nasolve.checkpoints import add_checkpoint
from nasolve.metal_coordination import (
    CoordinationError,
    coordination_phil,
    coordination_record,
    describe_site,
    parse_coordination,
    parse_geometry,
    plan_coordination,
)

from .helpers import pdb_record
from .test_autorefine import PHENIX_21, make_mtz_dump, make_refine, make_refine_run


def coordination_model(path: Path, *, metal: str = "AG", altloc: bool = False, n1_z: float = 0.0) -> Path:
    """A metal between N3 of C:1 and O4 of C:2, about on a line."""
    lines = [
        pdb_record("ATOM", 1, "N3", "DT", "C", 1, element="N", x=-2.2),
        pdb_record("ATOM", 2, "C4", "DT", "C", 1, element="C", x=-3.4),
        pdb_record("ATOM", 3, "O4", "DT", "C", 2, element="O", x=2.6, y=0.4),
        pdb_record("ATOM", 4, "S4", "DT", "C", 2, element="S", x=0.0, y=3.0),
        pdb_record("HETATM", 5, metal, metal, "Z", 1, element=metal),
        pdb_record("ATOM", 6, "O6", "DG", "C", 3, element="O", x=0.0, y=-2.5),
        pdb_record("ATOM", 7, "N7", "DG", "C", 3, element="N", x=0.0, y=-9.0),
        pdb_record("ATOM", 8, "N1", "DC", "C", 4, element="N", x=0.0, y=2.3, z=n1_z),
    ]
    if altloc:
        lines.append(lines[0][:16] + "B" + lines[0][17:])
        lines[0] = lines[0][:16] + "A" + lines[0][17:]
    path.write_text("".join(lines) + "END\n")
    return path


class ParseTests(unittest.TestCase):
    def test_specifications(self):
        spec = parse_coordination("B:101=B:4/N2T, A:12/O4")
        self.assertEqual((spec.metal.chain, spec.metal.resid, spec.metal.name), ("B", "101", None))
        self.assertEqual([donor.label for donor in spec.donors], ["B:4/N2T", "A:12/O4"])
        self.assertEqual(parse_coordination("h:1/HG=A:5/N3").metal.label, "h:1/HG")
        for text, message in (
            ("B:101", "give METAL=DONOR"),
            ("B:101=", "no donor atoms"),
            ("B:101=B:4", "needs its atom name"),
            ("B:101=B:4/N2T,B:4/N2T", "given twice"),
            ("B101=B:4/N2T", "not CHAIN:RESID"),
        ):
            with self.subTest(text=text), self.assertRaisesRegex(CoordinationError, message):
                parse_coordination(text)


class GeometryParseTests(unittest.TestCase):
    def test_presets(self):
        bare = parse_geometry("sqp")
        self.assertEqual((bare.metal, bare.preset.code, bare.preset.name), (None, "sqp", "square planar"))
        self.assertEqual((bare.preset.donors, bare.preset.angle_ideal, bare.preset.angle_sigma,
                          bare.preset.plane_sigma), (4, 90.0, 6.0, 0.2))
        named = parse_geometry(" B:101=SQP ")
        self.assertEqual((named.metal.label, named.preset.code), ("B:101", "sqp"))
        for text, message in (
            ("oct", r"unknown geometry 'oct' \(presets: sqp = square planar\)"),
            ("B:101=", "unknown geometry ''"),
            ("B101=sqp", "not CHAIN:RESID"),
            ("B:101=B:4=sqp", "give PRESET or METAL=PRESET"),
        ):
            with self.subTest(text=text), self.assertRaisesRegex(CoordinationError, message):
                parse_geometry(text)


class PlanTests(unittest.TestCase):
    def test_linear_site_gets_two_bonds_and_the_zero_degree_angle(self):
        with tempfile.TemporaryDirectory() as directory:
            model = coordination_model(Path(directory) / "model.pdb")
            site, = plan_coordination(model, [parse_coordination("Z:1=C:1/N3,C:2/O4")])
            self.assertEqual(site.geometry, "linear")
            self.assertEqual([(bond.donor.name, bond.ideal, bond.sigma, bond.current) for bond in site.bonds],
                             [("N3", 2.3, 0.2, 2.2), ("O4", 2.5, 0.2, 2.631)])
            self.assertEqual([atom.name for atom in site.angle.atoms], ["O4", "N3", "AG"])
            self.assertEqual((site.angle.ideal, site.angle.sigma), (0.0, 3.0))
            self.assertAlmostEqual(site.angle.current, 4.76, places=2)
            self.assertEqual(coordination_phil([site]), EXPECTED_PHIL)

    def test_one_or_several_donors_and_mercury(self):
        with tempfile.TemporaryDirectory() as directory:
            model = coordination_model(Path(directory) / "model.pdb", metal="HG")
            single, = plan_coordination(model, [parse_coordination("Z:1=C:2/O4")])
            self.assertEqual((single.geometry, single.angle), ("one donor", None))
            self.assertEqual(single.bonds[0].ideal, 2.1)
            several, = plan_coordination(model, [parse_coordination("Z:1/HG=C:1/N3,C:2/O4,C:3/O6")])
            self.assertEqual((several.geometry, several.angle),
                             ("3 donors (no angle restraint)", None))
            self.assertEqual([bond.ideal for bond in several.bonds], [2.1, 2.1, 2.1])
            self.assertNotIn("angle {", coordination_phil([several]))

    def test_what_is_refused_or_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = coordination_model(root / "model.pdb")
            for text, message in (
                ("Z:1=C:2/S4", r"no Ag-S length in the table \(Ag-O 2.5 A, Ag-N 2.3 A, Hg-O 2.1 A, Hg-N 2.1 A\)"),
                ("Z:2=C:1/N3", "no residue Z:2 in the model"),
                ("C:1=Z:1/AG", "the residue has several atoms .*name one"),
                ("Z:1=C:1/N9", "no atom N9 in DT C:1"),
                ("Z:1=Z:1/AG", "the metal is also a donor"),
            ):
                with self.subTest(text=text), self.assertRaisesRegex(CoordinationError, message):
                    plan_coordination(model, [parse_coordination(text)])
            with self.assertRaisesRegex(CoordinationError, "AG Z:1 AG is given twice"):
                plan_coordination(model, [parse_coordination("Z:1=C:1/N3"), parse_coordination("Z:1/AG=C:2/O4")])
            far, = plan_coordination(model, [parse_coordination("Z:1=C:3/N7")])
            self.assertEqual(far.warnings, ("DG C:3 N7 is 9.00 A from AG Z:1 AG (ideal 2.3 A)",))
            with self.assertRaisesRegex(CoordinationError, "alternate conformations"):
                plan_coordination(coordination_model(root / "alt.pdb", altloc=True),
                                  [parse_coordination("Z:1=C:1/N3")])


class SquarePlanarTests(unittest.TestCase):
    # Around Ag Z:1 in the plane z = 0: C:1 N3 (-x), C:3 O6 (-y), C:2 O4 (+x), C:4 N1 (+y).
    SQUARE = "Z:1=C:1/N3,C:3/O6,C:2/O4,C:4/N1"

    def test_four_bonds_and_four_angles_between_neighbours(self):
        with tempfile.TemporaryDirectory() as directory:
            model = coordination_model(Path(directory) / "model.pdb")
            site, = plan_coordination(model, [parse_coordination(self.SQUARE)], [parse_geometry("sqp")])
            self.assertEqual(site.geometry, "square planar (sqp)")
            self.assertEqual([(bond.donor.name, bond.ideal) for bond in site.bonds],
                             [("N3", 2.3), ("O6", 2.5), ("O4", 2.5), ("N1", 2.3)])
            self.assertIsNone(site.angle)                 # no linear angle
            self.assertEqual([tuple(atom.name for atom in angle.atoms) for angle in site.angles],
                             [("N3", "AG", "O6"), ("O6", "AG", "O4"), ("O4", "AG", "N1"), ("N1", "AG", "N3")])
            self.assertEqual({(angle.ideal, angle.sigma) for angle in site.angles}, {(90.0, 6.0)})
            self.assertEqual([round(angle.current, 1) for angle in site.angles], [90.0, 98.8, 81.2, 90.0])
            self.assertEqual(site.warnings, ())
            self.assertEqual([atom.name for atom in site.plane.atoms], ["AG", "N3", "O6", "O4", "N1"])
            self.assertEqual((site.plane.sigma, site.plane.current_rms), (0.2, 0.0))
            self.assertEqual(coordination_phil([site]), EXPECTED_SQP_PHIL)
            record = coordination_record(site)
            self.assertEqual((record["geometry"], record["preset"], record["angle"]),
                             ("square planar (sqp)", "sqp", None))
            self.assertEqual([item["atoms"] for item in record["angles"]][0],
                             ["DT C:1 N3", "AG Z:1 AG", "DG C:3 O6"])
            self.assertEqual(record["plane"], {
                "atoms": ["AG Z:1 AG", "DT C:1 N3", "DG C:3 O6", "DT C:2 O4", "DC C:4 N1"],
                "sigma": 0.2, "current_rms": 0.0})
            self.assertEqual(describe_site(site), (
                "AG Z:1 AG: C:1 N3 2.3 A (now 2.20), C:3 O6 2.5 A (now 2.50), C:2 O4 2.5 A (now 2.63), "
                "C:4 N1 2.3 A (now 2.30); square planar 90 +/- 6 deg: C:1 N3-AG-C:3 O6 (now 90.0), "
                "C:3 O6-AG-C:2 O4 (now 98.8), C:2 O4-AG-C:4 N1 (now 81.2), C:4 N1-AG-C:1 N3 (now 90.0); "
                "AG and its 4 donors in one plane, sigma 0.2 A (now 0.00 A rms)"))
            # Named metal: the same.
            named, = plan_coordination(model, [parse_coordination(self.SQUARE)], [parse_geometry("Z:1/AG=sqp")])
            self.assertEqual(coordination_phil([named]), EXPECTED_SQP_PHIL)

    def test_current_distance_from_the_plane(self):
        with tempfile.TemporaryDirectory() as directory:
            # N1 1 A above the plane of the others: rms 0.262 A from the best plane (numpy).
            model = coordination_model(Path(directory) / "model.pdb", n1_z=1.0)
            site, = plan_coordination(model, [parse_coordination(self.SQUARE)], [parse_geometry("sqp")])
            self.assertAlmostEqual(site.plane.current_rms, 0.262, places=3)
            # Without --geometry, or for a linear site: no plane.
            linear, = plan_coordination(model, [parse_coordination("Z:1=C:1/N3,C:2/O4")])
            self.assertIsNone(linear.plane)
            self.assertNotIn("planarity", coordination_phil([linear]))

    def test_donors_out_of_order_are_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            model = coordination_model(Path(directory) / "model.pdb")
            # N3 and O4 are across the square, not neighbours.
            site, = plan_coordination(model, [parse_coordination("Z:1=C:1/N3,C:2/O4,C:3/O6,C:4/N1")],
                                      [parse_geometry("sqp")])
            self.assertEqual(site.warnings, (
                "C:1 N3 and C:2 O4 are 171.2 deg apart at AG Z:1 AG: square planar takes the donors "
                "in order around the square, each next to the one before",
                "C:3 O6 and C:4 N1 are 180.0 deg apart at AG Z:1 AG: square planar takes the donors "
                "in order around the square, each next to the one before",
            ))

    def test_what_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            model = coordination_model(Path(directory) / "model.pdb")
            square = parse_coordination(self.SQUARE)
            for specs, geometries, message in (
                ([parse_coordination("Z:1=C:1/N3,C:2/O4")], ["sqp"],
                 r"--geometry sqp \(square planar\) needs 4 donors, listed in order around the "
                 r"square; AG Z:1 AG has 2"),
                ([square, parse_coordination("C:3/N7=C:1/N3")], ["sqp"],
                 "--geometry 'sqp': name the metal when there are several --coordination "
                 r"\(e.g. --geometry Z:1=sqp\)"),
                ([square], ["C:3/O6=sqp"], "--geometry 'C:3/O6=sqp': DG C:3 O6 has no --coordination"),
                ([square], ["sqp", "Z:1=sqp"], "--geometry: AG Z:1 AG is given two geometries"),
                ([], ["sqp"], "--geometry needs --coordination"),
            ):
                with self.subTest(geometries=geometries), self.assertRaisesRegex(CoordinationError, message):
                    plan_coordination(model, specs, [parse_geometry(text) for text in geometries])


EXPECTED_SQP_PHIL = """\
# NASolve metal coordination restraints (nasolve autorefine --coordination)
# Bonds: lab table, sigma 0.2 A. Linear (two donors): angle donor-donor-metal 0 deg, sigma 3 (not a physical angle; keeps the metal on the donor-donor line).
# Square planar (--geometry sqp): angle donor-metal-donor between donors listed next to each other (and the last with the first), 90 deg, sigma 6; the metal and its donors in one plane, sigma 0.2 A.
refinement {
  geometry_restraints.edits {
    # AG Z:1 AG: Z:1=C:1/N3,C:3/O6,C:2/O4,C:4/N1 (square planar (sqp))
    bond {
      action = add
      atom_selection_1 = chain Z and resid 1 and name AG
      atom_selection_2 = chain C and resid 1 and name N3
      distance_ideal = 2.3
      sigma = 0.2
    }
    bond {
      action = add
      atom_selection_1 = chain Z and resid 1 and name AG
      atom_selection_2 = chain C and resid 3 and name O6
      distance_ideal = 2.5
      sigma = 0.2
    }
    bond {
      action = add
      atom_selection_1 = chain Z and resid 1 and name AG
      atom_selection_2 = chain C and resid 2 and name O4
      distance_ideal = 2.5
      sigma = 0.2
    }
    bond {
      action = add
      atom_selection_1 = chain Z and resid 1 and name AG
      atom_selection_2 = chain C and resid 4 and name N1
      distance_ideal = 2.3
      sigma = 0.2
    }
    angle {
      action = add
      atom_selection_1 = chain C and resid 1 and name N3
      atom_selection_2 = chain Z and resid 1 and name AG
      atom_selection_3 = chain C and resid 3 and name O6
      angle_ideal = 90
      sigma = 6
    }
    angle {
      action = add
      atom_selection_1 = chain C and resid 3 and name O6
      atom_selection_2 = chain Z and resid 1 and name AG
      atom_selection_3 = chain C and resid 2 and name O4
      angle_ideal = 90
      sigma = 6
    }
    angle {
      action = add
      atom_selection_1 = chain C and resid 2 and name O4
      atom_selection_2 = chain Z and resid 1 and name AG
      atom_selection_3 = chain C and resid 4 and name N1
      angle_ideal = 90
      sigma = 6
    }
    angle {
      action = add
      atom_selection_1 = chain C and resid 4 and name N1
      atom_selection_2 = chain Z and resid 1 and name AG
      atom_selection_3 = chain C and resid 1 and name N3
      angle_ideal = 90
      sigma = 6
    }
    planarity {
      action = add
      atom_selection = (chain Z and resid 1 and name AG) or (chain C and resid 1 and name N3) or (chain C and resid 3 and name O6) or (chain C and resid 2 and name O4) or (chain C and resid 4 and name N1)
      sigma = 0.2
    }
  }
}
"""


EXPECTED_PHIL = """\
# NASolve metal coordination restraints (nasolve autorefine --coordination)
# Bonds: lab table, sigma 0.2 A. Linear (two donors): angle donor-donor-metal 0 deg, sigma 3 (not a physical angle; keeps the metal on the donor-donor line).
refinement {
  geometry_restraints.edits {
    # AG Z:1 AG: Z:1=C:1/N3,C:2/O4 (linear)
    bond {
      action = add
      atom_selection_1 = chain Z and resid 1 and name AG
      atom_selection_2 = chain C and resid 1 and name N3
      distance_ideal = 2.3
      sigma = 0.2
    }
    bond {
      action = add
      atom_selection_1 = chain Z and resid 1 and name AG
      atom_selection_2 = chain C and resid 2 and name O4
      distance_ideal = 2.5
      sigma = 0.2
    }
    angle {
      action = add
      atom_selection_1 = chain C and resid 2 and name O4
      atom_selection_2 = chain C and resid 1 and name N3
      atom_selection_3 = chain Z and resid 1 and name AG
      angle_ideal = 0
      sigma = 3
    }
  }
}
"""


class AutoRefineCoordinationTests(unittest.TestCase):
    def run_with_metal(self, root: Path) -> Path:
        run = make_refine_run(root)
        model = run / "PostMR" / "Model" / "readyset_model.pdb"
        kept = model.read_text().replace("END\n", "")
        coordination_model(root / "metal.pdb")
        model.write_text(kept + "".join(line + "\n" for line in (root / "metal.pdb").read_text().splitlines()
                                        if line.startswith(("ATOM", "HETATM"))) + "END\n")
        return run

    def test_coordination_is_written_passed_recorded_and_inherited(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = self.run_with_metal(root)

            def refine(**kwargs):
                return execute_autorefine(
                    run, make_refine(root, final_work=0.244, final_free=0.267), make_mtz_dump(root),
                    phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
                    macro_cycles=1, auto_select_success=False, **kwargs,
                )

            first = refine(coordination=("Z:1=C:1/N3,C:2/O4",))
            restraint = (first.report_path.parent / "metal_coordination.eff").resolve()
            self.assertEqual(restraint.read_text(), EXPECTED_PHIL)
            payload = json.loads(first.report_path.read_text())
            self.assertIn(str(restraint), payload["command"])
            recorded = payload["inputs"]["metal_coordination"]
            self.assertEqual(recorded["file"], str(restraint))
            self.assertEqual(recorded["inherited"], [])
            site = recorded["sites"][0]
            self.assertEqual((site["metal"], site["geometry"]), ("AG Z:1 AG", "linear"))
            self.assertEqual([bond["ideal"] for bond in site["bonds"]], [2.3, 2.5])
            self.assertEqual((site["angle"]["ideal"], site["angle"]["sigma"]), (0.0, 3.0))
            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            checkpoint = next(item for item in registry["checkpoints"] if item["id"] == first.checkpoint_id)
            self.assertIn(f"AutoRefine/{first.report_path.parent.name}/metal_coordination.eff",
                          [ref.get("relative_path") for ref in checkpoint["restraints"]])
            self.assertEqual([item["metal"] for item in checkpoint["metal_coordination"]], ["AG Z:1 AG"])

            # A later round from that checkpoint keeps the restraints without the option.
            second = refine(from_checkpoint=first.checkpoint_id)
            second_payload = json.loads(second.report_path.read_text())
            self.assertIn(str(restraint), second_payload["command"])
            self.assertEqual([item["metal"] for item in second_payload["inputs"]["metal_coordination"]["inherited"]],
                             ["AG Z:1 AG"])

            # Restraining the same metal again would add the bonds twice
            # (Phenix stops on that): refused before a round exists.
            rounds = sorted((run / "AutoRefine").glob("round_*"))
            for kwargs, message in (
                ({"from_checkpoint": first.checkpoint_id, "coordination": ("Z:1=C:2/O4",)},
                 f"AG Z:1 AG already has coordination restraints from checkpoint {first.checkpoint_id}"),
                ({"coordination": ("Z:1=C:2/S4",)}, "no Ag-S length in the table"),
                ({"coordination": ("Z:1",)}, "give METAL=DONOR"),
            ):
                with self.subTest(kwargs=kwargs), self.assertRaisesRegex(AutoRefineError, message):
                    refine(**kwargs)
            self.assertEqual(sorted((run / "AutoRefine").glob("round_*")), rounds)

    def test_restraints_are_planned_on_the_given_starting_model(self):
        # --model: the current distances are those of the given model.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = self.run_with_metal(root)
            model = run / "PostMR" / "Model" / "readyset_model.pdb"
            edited = root / "edited.pdb"
            edited.write_text("".join(
                (line[:30] + "%8.3f" % 0.5 + line[38:] if line.startswith("HETATM") and line[12:16].strip() == "AG"
                 else line) + "\n" for line in model.read_text().splitlines()))
            result = execute_autorefine(
                run, make_refine(root), make_mtz_dump(root), phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"}, macro_cycles=1, auto_select_success=False,
                coordination=("Z:1=C:1/N3,C:2/O4",), start_model=edited,
            )
            site = json.loads(result.report_path.read_text())["inputs"]["metal_coordination"]["sites"][0]
            self.assertEqual([bond["current"] for bond in site["bonds"]], [2.7, 2.138])

    def test_an_imported_model_keeps_the_coordination_record(self):
        # Before 2026-11-24 a checkpoints-add import kept the coordination file
        # but not its record, so a second --coordination for the same metal was
        # not refused and phenix.refine got both (Y-T run_002 refine-005).
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = self.run_with_metal(root)
            first = execute_autorefine(
                run, make_refine(root), make_mtz_dump(root), phenix_version=PHENIX_21,
                environment={"PATH": "/usr/bin:/bin"}, macro_cycles=1, auto_select_success=False,
                coordination=("Z:1=C:1/N3,C:2/O4",),
            )
            imported = add_checkpoint(run, name="hand edited", parent=first.checkpoint_id,
                                      model=first.report_path.parent / "refined_001.pdb")
            registry = json.loads((run / "AutoRefine" / "checkpoints.json").read_text())
            by_id = {item["id"]: item for item in registry["checkpoints"]}
            self.assertEqual(by_id[imported.checkpoint_id]["metal_coordination"],
                             by_id[first.checkpoint_id]["metal_coordination"])
            with self.assertRaisesRegex(
                    AutoRefineError, f"AG Z:1 AG already has coordination restraints from checkpoint "
                                     f"{imported.checkpoint_id}"):
                execute_autorefine(
                    run, make_refine(root), make_mtz_dump(root), phenix_version=PHENIX_21,
                    environment={"PATH": "/usr/bin:/bin"}, macro_cycles=1, auto_select_success=False,
                    from_checkpoint="hand edited", coordination=("Z:1=C:2/O4",),
                )

    def test_square_planar_round(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = self.run_with_metal(root)

            def refine(**kwargs):
                return execute_autorefine(
                    run, make_refine(root, final_work=0.244, final_free=0.267), make_mtz_dump(root),
                    phenix_version=PHENIX_21, environment={"PATH": "/usr/bin:/bin"},
                    macro_cycles=1, auto_select_success=False, **kwargs,
                )

            rounds = sorted((run / "AutoRefine").glob("round_*"))
            for kwargs, message in (
                ({"coordination_geometry": ("sqp",)}, "--geometry needs --coordination"),
                ({"coordination": ("Z:1=C:1/N3,C:2/O4",), "coordination_geometry": ("sqp",)},
                 "needs 4 donors"),
            ):
                with self.subTest(kwargs=kwargs), self.assertRaisesRegex(AutoRefineError, message):
                    refine(**kwargs)
            self.assertEqual(sorted((run / "AutoRefine").glob("round_*")), rounds)

            result = refine(coordination=(SquarePlanarTests.SQUARE,), coordination_geometry=("sqp",))
            restraint = result.report_path.parent / "metal_coordination.eff"
            self.assertEqual(restraint.read_text(), EXPECTED_SQP_PHIL)
            site = json.loads(result.report_path.read_text())["inputs"]["metal_coordination"]["sites"][0]
            self.assertEqual((site["geometry"], site["preset"], len(site["angles"]), len(site["plane"]["atoms"])),
                             ("square planar (sqp)", "sqp", 4, 5))

    def test_command_line_forwards_and_prints(self):
        from nasolve.cli import main
        from nasolve.config import AppConfig
        from nasolve.coot_view import CootViewError

        run = Path("/dataset/AutoMR/run_004")
        result = SimpleNamespace(
            status="AUTOREFINE_READY", message="ready", checkpoint_id="refine-002",
            parent_checkpoint="refine-001", report_path=Path("/nonexistent/report.json"),
            statistics={}, selected_as_current=True, model_path=None, model_cif=None,
            reflection_cif=None, map_coefficients=None, log_path=Path("/nonexistent/log"),
            run_directory=run, exit_code=0,
        )
        installation = SimpleNamespace(
            version=PHENIX_21, environment={},
            executables={"phenix.refine": Path("/phenix/phenix.refine"),
                         "phenix.mtz.dump": Path("/phenix/phenix.mtz.dump")},
        )
        payload = {"inputs": {"observation_labels": [], "start_model": {
            "source": "/data/edited.pdb", "copy": "/r/start_model.pdb",
            "instead_of": "/r/round_001/refined_001.pdb"}, "metal_coordination": {
            "file": "/r/metal_coordination.eff", "inherited": [], "sites": [{
                "metal": "AG B:101 AG", "geometry": "linear",
                "summary": "AG B:101 AG: B:4 N2T 2.3 A (now 2.17), A:12 O4 2.5 A (now 2.84); "
                           "linear: O4-N2T-AG 0 +/- 3 deg (now 25.7)",
                "warnings": ["a warning"],
            }]}}, "refinement": {}}
        output = io.StringIO()
        with (patch("nasolve.cli.load_config", return_value=AppConfig()),
              patch("nasolve.cli._workspace_run", return_value=run),
              patch("nasolve.cli.discover_phenix", return_value=installation),
              patch("nasolve.cli.remember_phenix"), patch("nasolve.cli.save_config"),
              patch("nasolve.cli.execute_autorefine", return_value=result) as call,
              patch("pathlib.Path.read_text", return_value=json.dumps(payload)),
              patch("nasolve.cli.resolve_view_profile", side_effect=CootViewError("no view")),
              redirect_stdout(output)):
            code = main(["autorefine", str(run), "--coordination", "B:101=B:4/N2T,A:12/O4",
                         "--coordination", "B:102=B:3/O6", "--geometry", "B:102=sqp",
                         "--model", "/data/edited.pdb"])
        self.assertEqual(code, 0)
        self.assertEqual(call.call_args.kwargs["start_model"], Path("/data/edited.pdb"))
        self.assertIn("Starting model: /data/edited.pdb (instead of the model of refine-001)",
                      output.getvalue())
        self.assertEqual(call.call_args.kwargs["coordination"],
                         ["B:101=B:4/N2T,A:12/O4", "B:102=B:3/O6"])
        self.assertEqual(call.call_args.kwargs["coordination_geometry"], ["B:102=sqp"])
        self.assertIn("Coordination: AG B:101 AG: B:4 N2T 2.3 A (now 2.17), A:12 O4 2.5 A (now 2.84); "
                      "linear: O4-N2T-AG 0 +/- 3 deg (now 25.7)", output.getvalue())
        self.assertIn("Warning: a warning", output.getvalue())


if __name__ == "__main__":
    unittest.main()
