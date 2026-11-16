#!/usr/bin/env python3
"""Validate pure GLYCAM06j-1 glucose against AmberTools sander."""

import argparse
import json
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
SHARED_SOURCE = ROOT / "validation" / "amber_family_mix" / "source"
WORK = HERE / "work"
EXPECTED = json.loads((HERE / "expected_results.json").read_text())
ENERGY_TOLERANCE = EXPECTED["energy_tolerance_kcal_per_mol"]
FORCE_MAX_TOLERANCE = EXPECTED["force_max_tolerance_kcal_per_mol_angstrom"]
FORCE_RMSE_TOLERANCE = EXPECTED["force_rmse_tolerance_kcal_per_mol_angstrom"]


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


def prepare_inputs():
    WORK.mkdir(exist_ok=True)
    topology = pmd.load_file(str(SHARED_SOURCE / "glucose.prmtop"))
    coordinates = pmd.load_file(str(SHARED_SOURCE / "glucose.pdb"))
    topology.coordinates = coordinates.coordinates
    topology.box = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]
    topology.save(str(WORK / "reference.prmtop"), overwrite=True)
    topology.save(str(WORK / "reference.rst7"), overwrite=True)
    topology.save(str(WORK / "glucose.pdb"), overwrite=True)
    return topology


def force_statistics():
    with netcdf_file(WORK / "sander.mdfrc", "r", mmap=False) as force_file:
        amber = force_file.variables["forces"].data.copy()[0]
    lines = (WORK / "lammps_forces.dump").read_text().splitlines()
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

    topology = prepare_inputs()
    raw_charge = sum(atom.charge for atom in topology.atoms)
    sys.path.insert(0, str(ROOT))
    from amber_to_lammps import amber2lammps
    amber2lammps(
        data_file=str(WORK / "glucose.data"),
        param_file=str(WORK / "glucose.parm"),
        topologies=[str(SHARED_SOURCE / "glucose.prmtop")], molecule_counts=[1],
        pdb_file=str(WORK / "glucose.pdb"), charges_target=[raw_charge],
    )
    run([
        sander, "-O", "-i", HERE / "sander_force.in", "-o", "sander_force.out",
        "-p", "reference.prmtop", "-c", "reference.rst7",
        "-r", "sander_force.rst7", "-inf", "sander_force.mdinfo",
        "-frc", "sander.mdfrc",
    ], cwd=WORK)
    output = run([lammps, "-in", "glucose.input", "-log", "none"], cwd=HERE)

    amber = amber_energies(WORK / "sander_force.out")
    candidate = lammps_energies(output)
    failures = []
    print(f"{'Component':<12} {'sander':>14} {'LAMMPS':>14} {'difference':>14}")
    for component, reference in amber.items():
        difference = abs(reference - candidate[component])
        print(f"{component:<12} {reference:>14.7f} {candidate[component]:>14.7f} {difference:>14.7f}")
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
        raise AssertionError("pure GLYCAM validation failed:\n" + "\n".join(failures))
    print("Pure GLYCAM sander energy and force validation passed.")


if __name__ == "__main__":
    main()
