#!/usr/bin/env python3
"""Validate a non-covalent ff14SB peptide + GAFF2 ligand system."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import parmed as pmd
from scipy.io import netcdf_file

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from validation.glycam_glucose.validate_sander import amber_energies, lammps_energies


WORK = HERE / "work"
EXPECTED = json.loads((HERE / "expected_results.json").read_text())


def run(command, cwd=HERE):
    result = subprocess.run(
        [str(item) for item in command], cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    if result.returncode:
        print(result.stdout)
        raise RuntimeError(f"Command failed: {' '.join(map(str, command))}")
    return result.stdout


def centered(coordinates, target):
    coordinates = np.asarray(coordinates).copy()
    return coordinates - coordinates.mean(axis=0) + np.asarray(target)


def prepare_inputs():
    WORK.mkdir(exist_ok=True)
    peptide = pmd.load_file(str(WORK / "peptide.prmtop"))
    peptide_pdb = pmd.load_file(str(WORK / "peptide-source.pdb"))
    peptide.coordinates = centered(peptide_pdb.coordinates, [10.0, 16.0, 16.0])
    peptide.box = None
    peptide.save(str(WORK / "peptide.prmtop"), overwrite=True)

    ligand = pmd.load_file(str(ROOT / "ethanol.prmtop"))
    ligand_pdb = pmd.load_file(str(ROOT / "ethanol.pdb"))
    ligand.coordinates = centered(ligand_pdb.coordinates, [22.0, 16.0, 16.0])
    ligand.box = None

    combined = peptide + ligand
    combined.box = [32.0, 32.0, 32.0, 90.0, 90.0, 90.0]
    combined.save(str(WORK / "system.pdb"), overwrite=True)
    combined.save(str(WORK / "reference.prmtop"), overwrite=True)
    combined.save(str(WORK / "reference.rst7"), overwrite=True)
    return peptide, ligand


def force_statistics():
    with netcdf_file(WORK / "sander.mdfrc", "r", mmap=False) as force_file:
        amber = force_file.variables["forces"].data.copy()[0]
    lines = (WORK / "lammps_forces.dump").read_text().splitlines()
    header = lines.index("ITEM: ATOMS id fx fy fz")
    rows = np.array([[float(value) for value in row.split()] for row in lines[header + 1:]])
    lammps = rows[np.argsort(rows[:, 0]), 1:]
    difference = amber - lammps
    return float(np.max(np.abs(difference))), float(math.sqrt(np.mean(difference**2)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sander", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    args = parser.parse_args()
    sander, lammps = args.sander.absolute(), args.lammps.absolute()
    for executable in (sander, lammps):
        if not executable.is_file():
            parser.error(f"Executable not found: {executable}")

    tleap = sander.parent / "tleap"
    if not tleap.is_file():
        parser.error(f"tleap not found beside sander: {tleap}")
    WORK.mkdir(exist_ok=True)
    run([tleap, "-f", HERE / "tleap.in"])
    peptide, ligand = prepare_inputs()
    sys.path.insert(0, str(ROOT))
    from amber_to_lammps import amber2lammps
    amber2lammps(
        str(WORK / "system.data"), str(WORK / "system.parm"),
        [str(WORK / "peptide.prmtop"), str(ROOT / "ethanol.prmtop")], [1, 1],
        str(WORK / "system.pdb"),
        [sum(a.charge for a in peptide.atoms), sum(a.charge for a in ligand.atoms)],
    )
    run([
        sander, "-O", "-i", HERE / "sander_force.in", "-o", "sander_force.out",
        "-p", "reference.prmtop", "-c", "reference.rst7",
        "-r", "sander_force.rst7", "-inf", "sander_force.mdinfo",
        "-frc", "sander.mdfrc",
    ], WORK)
    output = run([lammps, "-in", "system.input", "-log", "none"])
    reference = amber_energies(WORK / "sander_force.out")
    candidate = lammps_energies(output)
    failures = []
    print(f"{'Component':<12} {'sander':>14} {'LAMMPS':>14} {'difference':>14}")
    for name, value in reference.items():
        difference = abs(value - candidate[name])
        print(f"{name:<12} {value:>14.7f} {candidate[name]:>14.7f} {difference:>14.7f}")
        if difference > EXPECTED["energy_tolerance_kcal_per_mol"]:
            failures.append(f"{name} difference {difference:g}")
    force_max, force_rmse = force_statistics()
    print(f"Force max difference: {force_max:.9f} kcal/mol/Angstrom")
    print(f"Force RMSE:           {force_rmse:.9f} kcal/mol/Angstrom")
    if force_max > EXPECTED["force_max_tolerance_kcal_per_mol_angstrom"]:
        failures.append(f"force maximum {force_max:g}")
    if force_rmse > EXPECTED["force_rmse_tolerance_kcal_per_mol_angstrom"]:
        failures.append(f"force RMSE {force_rmse:g}")
    if failures:
        raise AssertionError("protein-ligand validation failed:\n" + "\n".join(failures))
    print("ff14SB + GAFF2 protein-ligand validation passed.")


if __name__ == "__main__":
    main()
