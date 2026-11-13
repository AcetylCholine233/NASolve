import stat
import json
import tempfile
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def symlinked_temporary_directory():
    """Exercise path aliases on every host, including macOS /var -> /private/var."""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        actual = root / "actual"
        actual.mkdir()
        alias = root / "alias"
        alias.symlink_to(actual, target_is_directory=True)
        yield str(alias)


def pdb_record(
    record: str,
    serial: int,
    atom: str,
    residue: str,
    chain: str,
    residue_number: int,
    occupancy: float = 1.0,
    element: str = "P",
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    b_factor: float = 20.0,
) -> str:
    return (
        f"{record:<6}{serial:5d} {atom:>4s} {residue:>3s} {chain:1s}"
        f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{b_factor:6.2f}          {element:>2s}\n"
    )


def model_text() -> str:
    return "".join([
        pdb_record("ATOM", 1, "P", "DA", "A", 1),
        pdb_record("ATOM", 2, "P", "DC", "A", 2),
        pdb_record("HETATM", 3, "MG", "MG", "A", 50, element="MG"),
        "END\n",
    ])


def w_model_text() -> str:
    """Minimal preflight-only inventory including the W recipe's D:1 site."""
    return model_text().replace("END\n", pdb_record("ATOM", 4, "P", "DC", "D", 1) + "END\n")


def make_dataset(root: Path, include_model: bool = True) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "staraniso-alldata.mtz").write_bytes(b"mtz")
    (root / "data_1.cif").write_text(
        "data_test\n_symmetry.space_group_name_H-M 'H 3'\n"
    )
    (root / "summary.html").write_text(
        "<html>\nPossible space group P1\nSpacegroup name          H3\n</html>\n"
    )
    if include_model:
        (root / "search.pdb").write_text(model_text())
    return root


# Column layout of an autoPROC/STARANISO ``staraniso_alldata-unique.mtz``
# (labels and types read from a real file header), in phenix.mtz.dump order.
STARANISO_COLUMNS = (
    ("H", "H"), ("K", "H"), ("L", "H"), ("F", "F"), ("SIGF", "Q"),
    ("IMEAN", "J"), ("SIGIMEAN", "Q"), ("DANO", "D"), ("SIGDANO", "Q"),
    ("ISYM", "Y"), ("I(+)", "K"), ("SIGI(+)", "M"), ("I(-)", "K"),
    ("SIGI(-)", "M"), ("F(+)", "G"), ("SIGF(+)", "L"), ("F(-)", "G"),
    ("SIGF(-)", "L"), ("SA_flag", "I"), ("FreeR_flag", "I"),
    ("F_early", "F"), ("SIGF_early", "Q"), ("F_late", "F"), ("SIGF_late", "Q"),
)


def mtz_dump_columns_text(columns=STARANISO_COLUMNS, wavelength: float = 1.377618) -> str:
    """Return the crystal/dataset part of phenix.mtz.dump's summary (iotbx format).

    As in autoPROC files, H, K, L sit in the HKL_base dataset (placeholder
    wavelength 9.999) and the data columns in one dataset with the real one.
    """
    def block(crystal, name, wave, selected):
        lines = [
            f"Crystal {crystal}:",
            f"  Name: {name}",
            "  Number of datasets: 1",
            "  Dataset 1:",
            f"    Name: {name}",
            f"    Wavelength: {wave:.6g}",
            f"    Number of columns: {len(selected)}",
            "    label        #valid  %valid    min     max type",
        ]
        for label, kind in selected:
            lines.append(f"    {label:<12} 1635 100.00%   -1.00   99.00 {kind}: legend")
        return lines

    indices = [column for column in columns if column[1] == "H"]
    data = [column for column in columns if column[1] != "H"]
    lines = block(1, "HKL_base", 9.999, indices) + block(2, "1", wavelength, data)
    return "\n".join(lines) + "\n"


def make_mtz_dump(
    root: Path,
    symbol: str = "H 3",
    number: int = 146,
    matrix_symbol: str | None = None,
    columns=STARANISO_COLUMNS,
) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    executable = root / "phenix.mtz.dump"
    matrix = symbol if matrix_symbol is None else matrix_symbol
    column_text = mtz_dump_columns_text(columns) if columns else ""
    executable.write_text(
        "#!/bin/sh\n"
        f"echo 'Space group symbol from file: {symbol}'\n"
        f"echo 'Space group number from file: {number}'\n"
        f"echo 'Space group from matrices: {matrix} (No. {number})'\n"
        "cat <<'NASOLVE_COLUMNS'\n"
        + column_text
        + "NASOLVE_COLUMNS\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_phaser(
    root: Path,
    tfz: float | None = 8.4,
    llg: float = 121.0,
    return_code: int = 0,
    create_outputs: bool = True,
    labin_echo: str | None = None,
) -> Path:
    executable = root / "phenix.phaser"
    output_commands = ""
    if create_outputs:
        output_commands = (
            "printf 'MODEL\\nEND\\n' > PHASER.1.pdb\n"
            "printf 'mtz' > PHASER.1.mtz\n"
        )
    score_command = (
        f"echo 'SOLU SET RFZ=9.1 TFZ={tfz} LLG={llg}'\n"
        if tfz is not None else "echo 'No placed solution'\n"
    )
    labin_command = f"echo 'LABIN {labin_echo}'\n" if labin_echo else ""
    executable.write_text(
        "#!/bin/sh\n"
        + labin_command
        + output_commands
        + score_command
        + f"exit {return_code}\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_ready_set(root: Path) -> Path:
    executable = root / "phenix.ready_set"
    executable.write_text(
        "#!/bin/sh\n"
        "cp \"$1\" prepared_model.updated.pdb\n"
        "printf 'data_ligands\\n' > prepared_model.ligands.cif\n"
        "printf 'ready_set {}\\n' > prepared_model.eff\n"
        "echo 'Build ligand and use user provided restraints'\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_coot(root: Path, version: str = "1.1.10") -> Path:
    executable = root / "coot"
    executable.write_text(
        "#!/bin/sh\n"
        "if [ \"$1\" = \"--version\" ]; then\n"
        f"  echo '{version}'\n"
        "  exit 0\n"
        "fi\n"
        "echo 'NASOLVE_COOT_PYTHON_OK'\n"
        "if [ -n \"$NASOLVE_COOT_OUTPUT\" ]; then\n"
        "  cp \"$NASOLVE_COOT_INPUT\" \"$NASOLVE_COOT_OUTPUT\"\n"
        "fi\n"
    )
    executable.chmod(executable.stat().st_mode | stat.S_IXUSR)
    return executable


def make_postmr_report(
    run: Path,
    model: Path,
    first: str = "DE",
    second: str = "DG",
    requested: str = "E:G",
) -> Path:
    phaser = run / "Phaser"
    phaser.mkdir(parents=True)
    solution = phaser / "mr_solution.pdb"
    solution.write_text(model.read_text())
    report = {
        "workflow": "automr",
        "stage": "phaser",
        "status": "MR_SUCCESS",
        "frame": {"name": "W"},
        "post_mr_plan": {
            "sequences": {},
            "standard_pair": {
                "requested": requested,
                "ligand_codes": [first, second],
            },
            "mutations": {},
        },
        "execution": {"phaser": {"solution_pdb": str(solution.resolve())}},
    }
    path = run / "report.json"
    path.write_text(json.dumps(report))
    return path
