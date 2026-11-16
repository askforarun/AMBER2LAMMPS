#!/usr/bin/env python3
"""Generate and execute a minimal AMBER 12-6-4 C4 LAMMPS fixture."""

import argparse
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import parmed as pmd


ROOT = Path(__file__).resolve().parents[2]


def run(command, cwd):
    result = subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed with status {result.returncode}: {' '.join(map(str, command))}\n"
            f"{result.stdout}"
        )
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lammps", type=Path, required=True)
    args = parser.parse_args()
    lammps = args.lammps.resolve()
    if not lammps.is_file():
        parser.error(f"LAMMPS executable not found: {lammps}")

    sys.path.insert(0, str(ROOT))
    from amber_to_lammps import amber2lammps

    with tempfile.TemporaryDirectory(prefix="amber2lammps-c4-") as tmp_name:
        work = Path(tmp_name)
        topology = pmd.load_file(str(ROOT / "ethanol.prmtop"))
        coefficient_count = len(topology.parm_data["LENNARD_JONES_ACOEF"])
        # A small nonzero like-like C4 coefficient is sufficient to exercise the
        # converter's complete 12-6-4 path and LAMMPS's LEPTON pair sub-style.
        topology.add_flag(
            "LENNARD_JONES_CCOEF",
            "5E16.8",
            data=[5.0] + [0.0] * (coefficient_count - 1),
        )
        prmtop = work / "c4.prmtop"
        topology.write_parm(str(prmtop))

        data_file = work / "c4.data"
        parameter_file = work / "c4.parm"
        amber2lammps(
            str(data_file),
            str(parameter_file),
            [str(prmtop)],
            [1],
            str(ROOT / "ethanol.pdb"),
            [0],
        )
        parameters = parameter_file.read_text()
        required = (
            "pair_style hybrid/overlay",
            " lepton ",
            'lepton "-5/r^4"',
        )
        for marker in required:
            if marker not in parameters:
                raise AssertionError(f"Generated C4 parameters omit {marker!r}")

        input_file = work / "run0.in"
        input_file.write_text(
            "units real\n"
            "atom_style full\n"
            "boundary p p p\n"
            "read_data c4.data\n"
            "include c4.parm\n"
            "kspace_style pppm 1.0e-6\n"
            "thermo 1\n"
            "thermo_style custom step pe evdwl ecoul elong\n"
            "run 0\n"
        )
        output = run([lammps, "-in", input_file.name, "-log", "none"], work)
        if "LEPTON" not in run([lammps, "-h"], work):
            raise AssertionError("LAMMPS ran the fixture but did not report LEPTON support")
        matches = re.findall(
            r"^\s*0\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+"
            r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$",
            output,
            flags=re.MULTILINE,
        )
        if not matches:
            raise AssertionError(f"Could not find run-0 energy row in LAMMPS output:\n{output}")
        energies = tuple(float(value) for value in matches[-1])
        if not all(math.isfinite(value) for value in energies):
            raise AssertionError(f"Non-finite C4 run-0 energies: {energies}")
        print("12-6-4 C4 LEPTON run 0 passed")
        print(f"PotEng={energies[0]:.9g} E_vdwl={energies[1]:.9g} "
              f"E_coul={energies[2]:.9g} E_long={energies[3]:.9g}")


if __name__ == "__main__":
    main()
