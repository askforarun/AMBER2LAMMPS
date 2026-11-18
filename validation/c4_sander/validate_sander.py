#!/usr/bin/env python3
"""Validate a Li/Merz Zn(II) 12-6-4 TIP3P system against sander."""

import argparse
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
import parmed as pmd
from scipy.io import netcdf_file


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = HERE / "source"
WORK = HERE / "work"
EXPECTED = json.loads((HERE / "expected_results.json").read_text())
TOLERANCES = EXPECTED["tolerances"]


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
        "bond": bond,
        "angle": angle,
        "dihedral": dihedral,
        "vdw": vdw + vdw14,
        "coulomb": coul + coul14,
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
            "bond": parsed["E_bond"],
            "angle": parsed["E_angle"],
            "dihedral": parsed["E_dihed"],
            "vdw": parsed["E_vdwl"],
            "coulomb": parsed["E_coul"] + parsed["E_long"],
            "potential": parsed["PotEng"],
        }
    raise ValueError("Could not parse the LAMMPS energy row")


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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sander", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    args = parser.parse_args()
    sander = args.sander.resolve()
    lammps = args.lammps.resolve()
    for executable in (sander, lammps):
        if not executable.is_file():
            parser.error(f"Executable not found: {executable}")
    if "LEPTON" not in run([lammps, "-h"]):
        parser.error(f"LAMMPS executable lacks the LEPTON package: {lammps}")

    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir()

    topology = pmd.load_file(str(SOURCE / "system_1264.prmtop"), xyz=str(SOURCE / "system.pdb"))
    c4_coefficients = topology.parm_data.get("LENNARD_JONES_CCOEF", [])
    nonzero_c4 = [value for value in c4_coefficients if abs(value) > 1.0e-12]
    if not nonzero_c4:
        raise AssertionError("Reference topology contains no nonzero C4 coefficients")
    topology.save(str(WORK / "reference.rst7"), format="rst7")
    raw_charge = sum(atom.charge for atom in topology.atoms)
    sys.path.insert(0, str(ROOT))
    from amber_to_lammps import amber2lammps
    amber2lammps(
        data_file=str(WORK / "system.data"),
        param_file=str(WORK / "system.parm"),
        topologies=[str(SOURCE / "system_1264.prmtop")],
        molecule_counts=[1],
        pdb_file=str(SOURCE / "system.pdb"),
        charges_target=[raw_charge],
    )

    run([
        sander, "-O", "-i", HERE / "sander_force.in", "-o", "sander_force.out",
        "-p", SOURCE / "system_1264.prmtop", "-c", "reference.rst7",
        "-r", "sander_force.rst7", "-inf", "sander_force.mdinfo",
        "-frc", "sander.mdfrc",
    ], cwd=WORK)
    output = run([lammps, "-in", HERE / "system.input", "-log", "none"], cwd=HERE)

    amber = amber_energies(WORK / "sander_force.out")
    candidate = lammps_energies(output)
    failures = []
    print(f"{'Component':<12} {'sander':>14} {'LAMMPS':>14} {'difference':>14} {'tolerance':>12}")
    for component, reference in amber.items():
        difference = abs(reference - candidate[component])
        tolerance = TOLERANCES[f"{component}_energy_kcal_per_mol"]
        print(
            f"{component:<12} {reference:>14.7f} {candidate[component]:>14.7f} "
            f"{difference:>14.7f} {tolerance:>12.7f}"
        )
        if difference > tolerance:
            failures.append(f"{component} difference {difference:g} > {tolerance:g}")
    force_max, force_rmse = force_statistics()
    print(f"Force max difference: {force_max:.9f} kcal/mol/Angstrom")
    print(f"Force RMSE:           {force_rmse:.9f} kcal/mol/Angstrom")
    if force_max > TOLERANCES["force_max_kcal_per_mol_angstrom"]:
        failures.append(
            f"force maximum {force_max:g} > "
            f"{TOLERANCES['force_max_kcal_per_mol_angstrom']:g}"
        )
    if force_rmse > TOLERANCES["force_rmse_kcal_per_mol_angstrom"]:
        failures.append(
            f"force RMSE {force_rmse:g} > "
            f"{TOLERANCES['force_rmse_kcal_per_mol_angstrom']:g}"
        )
    if failures:
        raise AssertionError("12-6-4 sander validation failed:\n" + "\n".join(failures))
    print(f"Validated {len(nonzero_c4)} nonzero C4 coefficient(s).")
    print("Zn(II)/TIP3P 12-6-4 sander energy and force validation passed.")


if __name__ == "__main__":
    main()
