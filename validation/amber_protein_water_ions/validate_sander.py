#!/usr/bin/env python3
"""Validate ff14SB + TIP3P + NaCl conversion against AmberTools sander."""

import argparse
import math
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import parmed as pmd
from scipy.io import netcdf_file


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = HERE / "source"
AMBER = HERE / "sander"
LAMMPS = HERE / "amber2lammps"
TOLERANCES = {
    "bond": 5.0e-4,
    "angle": 5.0e-4,
    "dihedral": 5.0e-4,
    "vdw": 2.0e-2,
    "coulomb": 1.0e-2,
    "potential": 3.0e-2,
}
FORCE_MAX_TOLERANCE = 1.0e-2
FORCE_RMSE_TOLERANCE = 1.0e-3


def run(command, cwd=HERE):
    result = subprocess.run(
        [str(item) for item in command], cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    if result.returncode:
        print(result.stdout)
        raise RuntimeError(f"Command failed: {' '.join(map(str, command))}")
    return result.stdout


def amber_energies(path):
    text = path.read_text()
    match = re.search(
        r"NSTEP =\s*0.*?EPtot\s*=\s*([-+0-9.E]+).*?"
        r"BOND\s*=\s*([-+0-9.E]+)\s+ANGLE\s*=\s*([-+0-9.E]+)\s+"
        r"DIHED\s*=\s*([-+0-9.E]+).*?1-4 NB\s*=\s*([-+0-9.E]+)\s+"
        r"1-4 EEL\s*=\s*([-+0-9.E]+)\s+VDWAALS\s*=\s*([-+0-9.E]+).*?"
        r"EELEC\s*=\s*([-+0-9.E]+)", text, re.DOTALL,
    )
    if not match:
        raise ValueError("Could not parse the sander NSTEP 0 energy block")
    potential, bond, angle, dihedral, vdw14, coul14, vdw, coul = map(float, match.groups())
    return {
        "bond": bond, "angle": angle, "dihedral": dihedral,
        "vdw": vdw + vdw14, "coulomb": coul + coul14,
        "potential": potential,
    }


def lammps_energies(output):
    wanted = ("E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "E_long", "PotEng")
    lines = output.splitlines()
    for index, line in enumerate(lines[:-1]):
        headers = line.split()
        if not all(field in headers for field in wanted):
            continue
        parsed = dict(zip(headers, map(float, lines[index + 1].split())))
        return {
            "bond": parsed["E_bond"], "angle": parsed["E_angle"],
            "dihedral": parsed["E_dihed"], "vdw": parsed["E_vdwl"],
            "coulomb": parsed["E_coul"] + parsed["E_long"],
            "potential": parsed["PotEng"],
        }
    raise ValueError("Could not parse the LAMMPS energy row")


def force_statistics():
    with netcdf_file(AMBER / "sander.mdfrc", "r", mmap=False) as force_file:
        amber = force_file.variables["forces"].data.copy()[0]
    lines = (LAMMPS / "lammps_forces.dump").read_text().splitlines()
    header = lines.index("ITEM: ATOMS id fx fy fz")
    rows = np.array([[float(value) for value in line.split()] for line in lines[header + 1:]])
    lammps = rows[np.argsort(rows[:, 0]), 1:]
    if amber.shape != lammps.shape:
        raise ValueError(f"Force shape mismatch: AMBER {amber.shape}, LAMMPS {lammps.shape}")
    difference = amber - lammps
    return float(np.max(np.abs(difference))), float(math.sqrt(np.mean(difference**2)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sander", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    args = parser.parse_args()
    sander = args.sander.absolute()
    lammps = args.lammps.absolute()
    for executable in (sander, lammps):
        if not executable.is_file():
            parser.error(f"Executable not found: {executable}")

    topology = pmd.load_file(str(SOURCE / "system.prmtop"))
    raw_charge = sum(atom.charge for atom in topology.atoms)
    sys.path.insert(0, str(ROOT))
    from amber_to_lammps import amber2lammps
    amber2lammps(
        data_file=str(LAMMPS / "system.data"),
        param_file=str(LAMMPS / "system.parm"),
        topologies=[str(SOURCE / "system.prmtop")], molecule_counts=[1],
        pdb_file=str(SOURCE / "system.pdb"), charges_target=[raw_charge],
    )
    run([sys.executable, "build_reference.py"], cwd=AMBER)
    run([
        sander, "-O", "-i", "sander_force.in", "-o", "sander_force.out",
        "-p", "reference.prmtop", "-c", "reference.rst7",
        "-r", "sander_force.rst7", "-inf", "sander_force.mdinfo",
        "-frc", "sander.mdfrc",
    ], cwd=AMBER)
    output = run([lammps, "-in", "system.input", "-log", "none"], cwd=LAMMPS)

    amber = amber_energies(AMBER / "sander_force.out")
    lammps_values = lammps_energies(output)
    failures = []
    print(f"{'Component':<12} {'sander':>14} {'LAMMPS':>14} {'difference':>14} {'tolerance':>12}")
    for component, reference in amber.items():
        candidate = lammps_values[component]
        difference = abs(reference - candidate)
        tolerance = TOLERANCES[component]
        print(f"{component:<12} {reference:>14.7f} {candidate:>14.7f} {difference:>14.7f} {tolerance:>12.7f}")
        if difference > tolerance:
            failures.append(f"{component} difference {difference:g} > {tolerance:g}")

    force_max, force_rmse = force_statistics()
    print(f"Force max difference: {force_max:.9f} kcal/mol/Angstrom")
    print(f"Force RMSE:           {force_rmse:.9f} kcal/mol/Angstrom")
    if force_max > FORCE_MAX_TOLERANCE:
        failures.append(f"force maximum {force_max:g} > {FORCE_MAX_TOLERANCE:g}")
    if force_rmse > FORCE_RMSE_TOLERANCE:
        failures.append(f"force RMSE {force_rmse:g} > {FORCE_RMSE_TOLERANCE:g}")
    if failures:
        raise AssertionError("protein/water/ion validation failed:\n" + "\n".join(failures))
    run([sys.executable, "prepare_shake_data.py"], cwd=HERE)
    run([lammps, "-in", "system_shake.input", "-log", "none"], cwd=LAMMPS)
    print("ff14SB/TIP3P/NaCl sander validation passed.")


if __name__ == "__main__":
    main()
