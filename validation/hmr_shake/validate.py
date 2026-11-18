#!/usr/bin/env python3
"""Validate HMR masses, energies, forces, and TIP3P SHAKE execution."""

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
BASE = ROOT / "validation" / "amber_protein_water_ions"
SOURCE = BASE / "source"
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


def section_records(path, name, next_names):
    lines = path.read_text().splitlines()
    start = next(index for index, line in enumerate(lines) if line.strip() == name)
    end = len(lines)
    for index in range(start + 1, len(lines)):
        if lines[index].strip() in next_names:
            end = index
            break
    return [line.split() for line in lines[start + 1:end] if line.strip()]


def lammps_mass_mapping(data_file):
    masses = {
        int(row[0]): float(row[1])
        for row in section_records(data_file, "Masses", {"Atoms"})
    }
    atom_rows = section_records(data_file, "Atoms", {"Bonds"})
    atom_types = [int(row[2]) for row in atom_rows]
    return masses, atom_types


def verify_masses(original, hmr, data_file):
    masses, atom_types = lammps_mass_mapping(data_file)
    if len(atom_types) != len(hmr.atoms):
        raise AssertionError(f"Atom count mismatch: {len(atom_types)} != {len(hmr.atoms)}")
    differences = []
    changed = 0
    for index, (before, after, type_id) in enumerate(zip(original.atoms, hmr.atoms, atom_types), 1):
        emitted = masses[type_id]
        difference = abs(emitted - after.mass)
        if difference > 5.0e-7:
            differences.append(f"atom {index}: AMBER {after.mass:.8f}, LAMMPS {emitted:.8f}")
        if abs(before.mass - after.mass) > 1.0e-8:
            changed += 1
        if after.residue.name == "WAT" and abs(before.mass - after.mass) > 1.0e-8:
            differences.append(f"water atom {index} was repartitioned")
    if differences:
        raise AssertionError("HMR mass verification failed:\n" + "\n".join(differences[:20]))
    before_total = sum(atom.mass for atom in original.atoms)
    after_total = sum(atom.mass for atom in hmr.atoms)
    emitted_total = sum(masses[type_id] for type_id in atom_types)
    if abs(before_total - after_total) > 1.0e-6 or abs(after_total - emitted_total) > 5.0e-4:
        raise AssertionError(
            f"Total mass mismatch: original={before_total}, HMR={after_total}, LAMMPS={emitted_total}"
        )
    if changed == 0:
        raise AssertionError("HMR changed no atom masses")
    return changed, before_total, masses, atom_types


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
    with netcdf_file(WORK / "sander.mdfrc", "r", mmap=False) as force_file:
        amber = force_file.variables["forces"].data.copy()[0]
    lines = (WORK / "lammps_forces.dump").read_text().splitlines()
    header = lines.index("ITEM: ATOMS id fx fy fz")
    rows = np.array([[float(value) for value in line.split()] for line in lines[header + 1:]])
    lammps = rows[np.argsort(rows[:, 0]), 1:]
    difference = amber - lammps
    return float(np.max(np.abs(difference))), float(math.sqrt(np.mean(difference**2)))


def write_lammps_input(path, data_name, parameter_name, dump_forces=False, shake=None):
    text = (
        "units real\natom_style full\nboundary p p p\n"
        "dielectric 1.0000346111225924\n"
        "bond_style harmonic\nangle_style harmonic\ndihedral_style fourier\n"
        "pair_style lj/cut/coul/long 9.0 9.0\n"
        "pair_modify mix arithmetic tail yes table 0\n"
        f"read_data {data_name}\ninclude {parameter_name}\n"
        "kspace_style ewald 1.0e-8\nkspace_modify gewald 0.42001\n"
    )
    if shake is not None:
        oxygen_type, hydrogen_type, oh_bond_type, water_angle_type = shake
        text += (
            f"group tip3p type {oxygen_type} {hydrogen_type}\n"
            f"fix rigid_tip3p tip3p shake 1.0e-6 200 0 "
            f"b {oh_bond_type} a {water_angle_type}\n"
        )
    text += "thermo 1\nthermo_style custom ebond eangle edihed eimp evdwl ecoul elong etail pe\n"
    if dump_forces:
        text += (
            "dump validation_forces all custom 1 lammps_forces.dump id fx fy fz\n"
            "dump_modify validation_forces sort id\n"
        )
    text += "run 0\n"
    path.write_text(text)


def water_type_ids(topology, atom_types, data_file):
    water_atoms = [atom for atom in topology.atoms if atom.residue.name == "WAT"]
    oxygen_types = {atom_types[atom.idx] for atom in water_atoms if atom.element == 8}
    hydrogen_types = {atom_types[atom.idx] for atom in water_atoms if atom.element == 1}
    if len(oxygen_types) != 1 or len(hydrogen_types) != 1:
        raise AssertionError(f"Unexpected water atom types: O={oxygen_types}, H={hydrogen_types}")
    oxygen_type = oxygen_types.pop()
    hydrogen_type = hydrogen_types.pop()

    atom_type_by_id = {index + 1: type_id for index, type_id in enumerate(atom_types)}
    bond_rows = section_records(data_file, "Bonds", {"Angles"})
    oh_types = set()
    hh_types = set()
    for row in bond_rows:
        bond_type, atom1, atom2 = int(row[1]), int(row[2]), int(row[3])
        pair = {atom_type_by_id[atom1], atom_type_by_id[atom2]}
        if pair == {oxygen_type, hydrogen_type}:
            oh_types.add(bond_type)
        elif atom_type_by_id[atom1] == atom_type_by_id[atom2] == hydrogen_type:
            hh_types.add(bond_type)
    if len(oh_types) != 1 or len(hh_types) != 1:
        raise AssertionError(f"Unexpected water bond types: O-H={oh_types}, H-H={hh_types}")
    return oxygen_type, hydrogen_type, oh_types.pop(), hh_types.pop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sander", type=Path, required=True)
    parser.add_argument("--lammps", type=Path, required=True)
    args = parser.parse_args()
    sander, lammps = args.sander.resolve(), args.lammps.resolve()
    for executable in (sander, lammps):
        if not executable.is_file():
            parser.error(f"Executable not found: {executable}")

    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir()
    original = pmd.load_file(str(SOURCE / "system.prmtop"))
    hmr = pmd.load_file(str(SOURCE / "system.prmtop"), xyz=str(SOURCE / "system.pdb"))
    pmd.tools.actions.HMassRepartition(hmr).execute()
    hmr.write_parm(str(WORK / "system_hmr.prmtop"))
    hmr.save(str(WORK / "reference.rst7"), format="rst7")

    raw_charge = sum(atom.charge for atom in hmr.atoms)
    sys.path.insert(0, str(ROOT))
    from amber_to_lammps import amber2lammps
    amber2lammps(
        data_file=str(WORK / "system.data"), param_file=str(WORK / "system.parm"),
        topologies=[str(WORK / "system_hmr.prmtop")], molecule_counts=[1],
        pdb_file=str(SOURCE / "system.pdb"), charges_target=[raw_charge],
    )
    changed, total_mass, _, atom_types = verify_masses(original, hmr, WORK / "system.data")
    reference_metadata = EXPECTED["reference"]
    if changed != reference_metadata["changed_atom_masses"]:
        raise AssertionError(
            f"Expected {reference_metadata['changed_atom_masses']} changed masses, got {changed}"
        )
    if abs(total_mass - reference_metadata["total_mass_dalton"]) > 1.0e-6:
        raise AssertionError(
            f"Expected total mass {reference_metadata['total_mass_dalton']}, got {total_mass}"
        )

    run([
        sander, "-O", "-i", BASE / "sander" / "sander_force.in", "-o", "sander_force.out",
        "-p", "system_hmr.prmtop", "-c", "reference.rst7",
        "-r", "sander_force.rst7", "-inf", "sander_force.mdinfo", "-frc", "sander.mdfrc",
    ], cwd=WORK)
    write_lammps_input(WORK / "system.input", "system.data", "system.parm", dump_forces=True)
    output = run([lammps, "-in", "system.input", "-log", "none"], cwd=WORK)

    amber, candidate = amber_energies(WORK / "sander_force.out"), lammps_energies(output)
    failures = []
    print(f"{'Component':<12} {'sander':>14} {'LAMMPS':>14} {'difference':>14} {'tolerance':>12}")
    for component, reference in amber.items():
        difference = abs(reference - candidate[component])
        tolerance = TOLERANCES[f"{component}_energy_kcal_per_mol"]
        print(f"{component:<12} {reference:>14.7f} {candidate[component]:>14.7f} {difference:>14.7f} {tolerance:>12.7f}")
        if difference > tolerance:
            failures.append(f"{component} difference {difference:g} > {tolerance:g}")
    force_max, force_rmse = force_statistics()
    if force_max > TOLERANCES["force_max_kcal_per_mol_angstrom"]:
        failures.append(f"force maximum {force_max:g} exceeds tolerance")
    if force_rmse > TOLERANCES["force_rmse_kcal_per_mol_angstrom"]:
        failures.append(f"force RMSE {force_rmse:g} exceeds tolerance")

    from validation.amber_protein_water_ions.prepare_shake_data import prepare
    oxygen_type, hydrogen_type, oh_bond_type, hh_bond_type = water_type_ids(
        hmr, atom_types, WORK / "system.data"
    )
    angle_rows = section_records(WORK / "system.data", "Angles", {"Dihedrals"})
    existing_angle_types = [int(row[1]) for row in angle_rows]
    water_angle_type = max(existing_angle_types, default=0) + 1
    water_count = prepare(
        WORK / "system.data", WORK / "system.parm",
        WORK / "system_shake.data", WORK / "system_shake.parm",
        oxygen_type=oxygen_type, hydrogen_type=hydrogen_type,
        oh_bond_type=oh_bond_type, hh_bond_type=hh_bond_type,
        water_angle_type=water_angle_type,
    )
    if water_count != reference_metadata["shake_water_count"]:
        raise AssertionError(
            f"Expected {reference_metadata['shake_water_count']} SHAKE waters, got {water_count}"
        )
    write_lammps_input(
        WORK / "system_shake.input", "system_shake.data", "system_shake.parm",
        shake=(oxygen_type, hydrogen_type, oh_bond_type, water_angle_type),
    )
    run([lammps, "-in", "system_shake.input", "-log", "none"], cwd=WORK)

    print(f"HMR changed {changed} atom masses; total mass preserved at {total_mass:.6f} Da.")
    print(f"Force max difference: {force_max:.9f} kcal/mol/Angstrom")
    print(f"Force RMSE:           {force_rmse:.9f} kcal/mol/Angstrom")
    print(f"LAMMPS SHAKE initialized for {water_count} TIP3P waters.")
    if failures:
        raise AssertionError("HMR/SHAKE validation failed:\n" + "\n".join(failures))
    print("HMR mass, sander energy/force, and LAMMPS SHAKE validation passed.")


if __name__ == "__main__":
    main()
