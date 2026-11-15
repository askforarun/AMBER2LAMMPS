#!/usr/bin/env python3
"""Validate the mixed GAFF2/GLYCAM LAMMPS output against sander."""

import argparse
import math
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
from scipy.io import netcdf_file


HERE = Path(__file__).resolve().parent
ENERGY_TOLERANCE = 5.0e-4
FORCE_MAX_TOLERANCE = 0.1
FORCE_RMSE_TOLERANCE = 0.02


def run(command):
    result = subprocess.run(
        [str(item) for item in command], cwd=HERE, text=True,
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
        r"EELEC\s*=\s*([-+0-9.E]+)",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError("Could not parse the sander NSTEP 0 energy block")
    potential, bond, angle, dihedral, vdw14, coul14, vdw, coul = map(float, match.groups())
    return {
        "bond": bond,
        "angle": angle,
        "dihedral": dihedral,
        "vdw": vdw + vdw14,
        "coulomb": coul + coul14,
        "potential": potential,
    }


def lammps_energies(output):
    fields = ("E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "E_long", "PotEng")
    lines = output.splitlines()
    for index, line in enumerate(lines[:-1]):
        headers = line.split()
        if not all(field in headers for field in fields):
            continue
        values = lines[index + 1].split()
        parsed = dict(zip(headers, map(float, values)))
        return {
            "bond": parsed["E_bond"],
            "angle": parsed["E_angle"],
            "dihedral": parsed["E_dihed"],
            "vdw": parsed["E_vdwl"],
            "coulomb": parsed["E_coul"] + parsed["E_long"],
            "potential": parsed["PotEng"],
        }
    raise ValueError("Could not parse the LAMMPS energy row")


def force_statistics():
    with netcdf_file(HERE / "sander.mdfrc", "r", mmap=False) as force_file:
        amber = force_file.variables["forces"].data.copy()[0]

    lines = (HERE / "lammps_forces.dump").read_text().splitlines()
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
    for executable in (args.sander, args.lammps):
        if not executable.resolve().is_file():
            parser.error(f"Executable not found: {executable}")

    run([sys.executable, "build_amber_reference.py"])
    run([
        args.sander, "-O", "-i", "sander_force.in", "-o", "sander_force.out",
        "-p", "amber_reference.prmtop", "-c", "amber_reference.rst7",
        "-r", "sander_force.rst7", "-inf", "sander_force.mdinfo", "-frc", "sander.mdfrc",
    ])
    lammps_output = run([args.lammps, "-in", "ethanol_glucose.input", "-log", "none"])

    amber = amber_energies(HERE / "sander_force.out")
    lammps = lammps_energies(lammps_output)
    failures = []
    print(f"{'Component':<12} {'sander':>14} {'LAMMPS':>14} {'difference':>14}")
    for component in amber:
        difference = abs(amber[component] - lammps[component])
        print(f"{component:<12} {amber[component]:>14.7f} {lammps[component]:>14.7f} {difference:>14.7f}")
        if difference > ENERGY_TOLERANCE:
            failures.append(f"{component} difference {difference:g} > {ENERGY_TOLERANCE:g}")

    force_max, force_rmse = force_statistics()
    print(f"Force max difference: {force_max:.9f} kcal/mol/Angstrom")
    print(f"Force RMSE:           {force_rmse:.9f} kcal/mol/Angstrom")
    if force_max > FORCE_MAX_TOLERANCE:
        failures.append(f"force maximum {force_max:g} > {FORCE_MAX_TOLERANCE:g}")
    if force_rmse > FORCE_RMSE_TOLERANCE:
        failures.append(f"force RMSE {force_rmse:g} > {FORCE_RMSE_TOLERANCE:g}")
    if failures:
        raise AssertionError("sander validation failed:\n" + "\n".join(failures))
    print("sander energy and force validation passed.")


if __name__ == "__main__":
    main()
