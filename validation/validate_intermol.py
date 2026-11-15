#!/usr/bin/env python3
"""Compare fresh AMBER2LAMMPS and InterMol epon energies with LAMMPS."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
THERMO_FIELDS = (
    "E_bond",
    "E_angle",
    "E_dihed",
    "E_impro",
    "E_pair",
    "E_vdwl",
    "E_coul",
    "E_long",
    "E_tail",
    "PotEng",
)
TOLERANCES = {
    "E_bond": 1.0e-6,
    "E_angle": 5.0e-5,
    "E_dihed": 5.0e-5,
    "E_impro": 1.0e-8,
    "E_pair": 4.0e-2,
    "E_vdwl": 1.0e-6,
    "E_coul": 6.0e-2,
    "E_long": 3.0e-2,
    "E_tail": 1.0e-7,
    "PotEng": 4.0e-2,
}


def run(command, cwd, env=None):
    """Run a command and return combined output, raising on failure."""
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode:
        print(result.stdout)
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {' '.join(map(str, command))}"
        )
    return result.stdout


def extract_box_lines(intermol_data):
    """Extract the common cell used to isolate force-field energy differences."""
    wanted = ("xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz")
    found = {}
    for line in intermol_data.read_text().splitlines():
        stripped = line.strip()
        for suffix in wanted:
            if stripped.endswith(suffix):
                found[suffix] = line
    missing = [suffix for suffix in wanted if suffix not in found]
    if missing:
        raise ValueError(f"InterMol data file is missing box fields: {', '.join(missing)}")
    return found


def align_candidate_box(candidate_data, box_lines):
    """Give AMBER2LAMMPS output the InterMol cell without changing coordinates."""
    output = []
    tilt_written = False
    for line in candidate_data.read_text().splitlines():
        stripped = line.strip()
        matched = False
        for suffix in ("xlo xhi", "ylo yhi", "zlo zhi"):
            if stripped.endswith(suffix):
                output.append(box_lines[suffix])
                matched = True
                if suffix == "zlo zhi":
                    output.append(box_lines["xy xz yz"])
                    tilt_written = True
                break
        if not matched and not stripped.endswith("xy xz yz"):
            output.append(line)
    if not tilt_written:
        raise ValueError("AMBER2LAMMPS data file has no zlo/zhi box line")
    candidate_data.write_text("\n".join(output) + "\n")


def parse_thermo(output):
    """Extract the first thermo row containing the energy fields."""
    lines = output.splitlines()
    for index, line in enumerate(lines[:-1]):
        headers = line.split()
        if not all(field in headers for field in THERMO_FIELDS):
            continue
        values = lines[index + 1].split()
        if len(values) != len(headers):
            continue
        parsed = dict(zip(headers, map(float, values)))
        return {field: parsed[field] for field in THERMO_FIELDS}
    raise ValueError("Could not find a complete LAMMPS energy row")


def compare(reference, candidate):
    failures = []
    print("Energy comparison (kcal/mol):")
    print(f"{'Component':<12} {'InterMol':>14} {'AMBER2LAMMPS':>16} {'|Difference|':>14} {'Tolerance':>12}")
    for field in THERMO_FIELDS:
        difference = abs(reference[field] - candidate[field])
        tolerance = TOLERANCES[field]
        print(
            f"{field:<12} {reference[field]:>14.8f} {candidate[field]:>16.8f} "
            f"{difference:>14.8f} {tolerance:>12.8f}"
        )
        if difference > tolerance:
            failures.append(f"{field}: difference {difference} exceeds {tolerance}")
    if failures:
        raise AssertionError("InterMol validation failed:\n" + "\n".join(failures))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--intermol-dir", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    args = parser.parse_args()

    convert_script = args.intermol_dir.resolve() / "intermol" / "convert.py"
    lammps = args.lammps.resolve()
    if not convert_script.is_file():
        parser.error(f"InterMol convert.py not found: {convert_script}")
    if not lammps.is_file():
        parser.error(f"LAMMPS executable not found: {lammps}")

    with tempfile.TemporaryDirectory(prefix="amber2lammps_intermol_") as tmp_name:
        tmp = Path(tmp_name)
        intermol_out = tmp / "intermol"
        candidate_out = tmp / "amber2lammps"
        intermol_out.mkdir()
        candidate_out.mkdir()
        intermol_prmtop = intermol_out / "epon.prmtop"
        intermol_crd = intermol_out / "epon.crd"
        shutil.copy2(ROOT / "epon.prmtop", intermol_prmtop)
        shutil.copy2(ROOT / "epon.crd", intermol_crd)

        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            filter(None, (str(args.intermol_dir.resolve()), env.get("PYTHONPATH", "")))
        )
        run(
            [
                sys.executable,
                str(convert_script),
                "--amb_in",
                str(intermol_prmtop),
                str(intermol_crd),
                "--lammps",
                "--odir",
                str(intermol_out),
                "--oname",
                "epon_intermol",
                "--force",
            ],
            cwd=intermol_out,
            env=env,
        )

        sys.path.insert(0, str(ROOT))
        from amber_to_lammps import amber2lammps

        candidate_data = candidate_out / "epon_validation.data"
        candidate_param = candidate_out / "epon_validation.parm"
        amber2lammps(
            data_file=str(candidate_data),
            param_file=str(candidate_param),
            topologies=[str(ROOT / "epon.prmtop")],
            molecule_counts=[1],
            pdb_file=str(ROOT / "epon.pdb"),
            charges_target=[0],
        )

        intermol_data = intermol_out / "epon_intermol.lmp"
        intermol_input = intermol_out / "epon_intermol.input"
        align_candidate_box(candidate_data, extract_box_lines(intermol_data))
        candidate_input = candidate_out / "epon_validation.input"
        candidate_input.write_text((ROOT / "epon_validation.input").read_text())

        intermol_energy = parse_thermo(
            run([str(lammps), "-in", intermol_input.name, "-log", "none"], cwd=intermol_out)
        )
        candidate_energy = parse_thermo(
            run([str(lammps), "-in", candidate_input.name, "-log", "none"], cwd=candidate_out)
        )
        compare(intermol_energy, candidate_energy)

    print("InterMol/LAMMPS validation passed.")


if __name__ == "__main__":
    main()
