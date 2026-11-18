#!/usr/bin/env python3
"""Convert AMBER's three-bond TIP3P triangles to LAMMPS SHAKE topology.

AMBER encodes each rigid water with two O-H bonds and one H-H bond.  LAMMPS
SHAKE instead requires two O-H bonds and an H-O-H angle.  This fixture-specific
preparation preserves the geometry and exclusions while changing only that
topological representation.
"""

import argparse
import math
import re
from pathlib import Path


def section_index(lines, name):
    for index, line in enumerate(lines):
        if line.strip() == name:
            return index
    raise ValueError(f"Missing {name} section")


def records(lines, start, end):
    return [line.split() for line in lines[start + 1:end] if line.strip()]


def replace_header_count(lines, label, value):
    pattern = re.compile(rf"^\s*\d+\s+{re.escape(label)}\s*$")
    matches = [index for index, line in enumerate(lines) if pattern.match(line)]
    if len(matches) != 1:
        raise ValueError(f"Expected one '{label}' count, found {len(matches)}")
    lines[matches[0]] = f"{value} {label} \n"


def header_count(lines, label):
    pattern = re.compile(rf"^\s*(\d+)\s+{re.escape(label)}\s*$")
    matches = [pattern.match(line) for line in lines if pattern.match(line)]
    if len(matches) != 1:
        raise ValueError(f"Expected one '{label}' count, found {len(matches)}")
    return int(matches[0].group(1))


def prepare(data_input, parameters_input, data_output, parameters_output,
            oxygen_type=11, hydrogen_type=12, oh_bond_type=8, hh_bond_type=9,
            water_angle_type=12, oh_distance=0.9572, hh_distance=1.5136):
    lines = Path(data_input).read_text().splitlines(keepends=True)
    declared_angle_types = header_count(lines, "angle types")
    atoms_at = section_index(lines, "Atoms")
    bonds_at = section_index(lines, "Bonds")
    angles_at = section_index(lines, "Angles")
    dihedrals_at = section_index(lines, "Dihedrals")

    atom_types = {
        int(fields[0]): int(fields[2])
        for fields in records(lines, atoms_at, bonds_at)
    }
    bonds = [[int(value) for value in fields[:4]] for fields in records(lines, bonds_at, angles_at)]
    angles = [[int(value) for value in fields[:5]] for fields in records(lines, angles_at, dihedrals_at)]

    oh_neighbors = {}
    for _, bond_type, atom1, atom2 in bonds:
        if bond_type != oh_bond_type:
            continue
        types = (atom_types[atom1], atom_types[atom2])
        if types == (oxygen_type, hydrogen_type):
            oxygen, hydrogen = atom1, atom2
        elif types == (hydrogen_type, oxygen_type):
            oxygen, hydrogen = atom2, atom1
        else:
            raise ValueError(f"Bond type {oh_bond_type} is not exclusively water O-H")
        oh_neighbors.setdefault(hydrogen, set()).add(oxygen)

    retained_bonds = []
    water_angles = []
    for _, bond_type, atom1, atom2 in bonds:
        if bond_type != hh_bond_type:
            retained_bonds.append((bond_type, atom1, atom2))
            continue
        if atom_types[atom1] != hydrogen_type or atom_types[atom2] != hydrogen_type:
            raise ValueError(f"Bond type {hh_bond_type} is not exclusively water H-H")
        common = oh_neighbors.get(atom1, set()) & oh_neighbors.get(atom2, set())
        if len(common) != 1:
            raise ValueError(f"H-H bond {atom1}-{atom2} has {len(common)} common water oxygens")
        water_angles.append((water_angle_type, atom1, common.pop(), atom2))

    if not water_angles:
        raise ValueError("No AMBER three-bond water triangles were found")

    bond_lines = [
        f"{index} {bond_type} {atom1} {atom2}\n"
        for index, (bond_type, atom1, atom2) in enumerate(retained_bonds, 1)
    ]
    next_angle = max((angle[0] for angle in angles), default=0) + 1
    angle_lines = [" ".join(map(str, angle)) + "\n" for angle in angles]
    angle_lines.extend(
        f"{index} {angle_type} {atom1} {oxygen} {atom2}\n"
        for index, (angle_type, atom1, oxygen, atom2) in enumerate(water_angles, next_angle)
    )

    output = lines[:bonds_at + 1] + ["\n"] + bond_lines + ["\n"]
    output += lines[angles_at:angles_at + 1] + ["\n"] + angle_lines + ["\n"]
    output += lines[dihedrals_at:]
    replace_header_count(output, "bonds", len(retained_bonds))
    replace_header_count(output, "angles", len(angles) + len(water_angles))
    replace_header_count(
        output, "bond types",
        max((bond_type for bond_type, _, _ in retained_bonds), default=0),
    )
    replace_header_count(output, "angle types", max(declared_angle_types, water_angle_type))
    Path(data_output).write_text("".join(output))

    cosine = (2.0 * oh_distance**2 - hh_distance**2) / (2.0 * oh_distance**2)
    angle_degrees = math.degrees(math.acos(cosine))
    parameters = Path(parameters_input).read_text()
    parameters, removed_coefficients = re.subn(
        rf"^\s*bond_coeff\s+{hh_bond_type}(?:\s|$).*\n?",
        "", parameters, flags=re.MULTILINE,
    )
    if removed_coefficients != 1:
        raise ValueError(
            f"Expected one bond_coeff for removed H-H type {hh_bond_type}, "
            f"found {removed_coefficients}"
        )
    parameters += (
        "\n# SHAKE representation of AMBER three-bond TIP3P; K=0 because fix shake constrains it\n"
        f"angle_coeff {water_angle_type} 0.0000 {angle_degrees:.6f}  # HW-OW-HW\n"
    )
    Path(parameters_output).write_text(parameters)
    return len(water_angles)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-input", default="amber2lammps/system.data")
    parser.add_argument("--parameters-input", default="amber2lammps/system.parm")
    parser.add_argument("--data-output", default="amber2lammps/system_shake.data")
    parser.add_argument("--parameters-output", default="amber2lammps/system_shake.parm")
    args = parser.parse_args()
    count = prepare(args.data_input, args.parameters_input, args.data_output, args.parameters_output)
    print(f"Prepared {count} TIP3P waters for LAMMPS SHAKE")


if __name__ == "__main__":
    main()
