#!/usr/bin/env python3
"""Measure how a per-template residual charge accumulates under replication."""

import argparse
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile

import parmed as pmd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

from amber_to_lammps import amber2lammps


def lammps_atom_charge(path):
    """Sum charges in the Atoms section of a LAMMPS full-style data file."""
    total = 0.0
    in_atoms = False
    seen_atom = False
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if line == "Atoms" or line.startswith("Atoms #"):
            in_atoms = True
            continue
        if not in_atoms or not line or line.startswith("#"):
            continue
        fields = line.split()
        if not fields[0].isdigit():
            if seen_atom:
                break
            continue
        if len(fields) < 4:
            raise ValueError(f"Malformed Atoms row in {path}: {line}")
        total += float(fields[3])
        seen_atom = True
    if not seen_atom:
        raise ValueError(f"No Atoms rows found in {path}")
    return total


def generate_results():
    topology_path = ROOT / "epon.prmtop"
    pdb_path = ROOT / "epon.pdb"
    intermol_path = ROOT / "validation/intermol_epon/intermol/epon_intermol.lmp"
    topology = pmd.load_file(str(topology_path))
    source_charge = float(sum(atom.charge for atom in topology.atoms))

    with tempfile.TemporaryDirectory(prefix="a2l-charge-benchmark-") as directory:
        directory = Path(directory)
        raw_data = directory / "raw.data"
        normalized_data = directory / "normalized.data"
        with contextlib.redirect_stdout(io.StringIO()):
            amber2lammps(
                str(raw_data), str(directory / "raw.parm"), [str(topology_path)],
                [1], str(pdb_path), [source_charge],
            )
            amber2lammps(
                str(normalized_data), str(directory / "normalized.parm"),
                [str(topology_path)], [1], str(pdb_path), [0.0],
            )
        amber_raw_charge = lammps_atom_charge(raw_data)
        amber_normalized_charge = lammps_atom_charge(normalized_data)

    intermol_charge = lammps_atom_charge(intermol_path)
    counts = (1, 10, 100, 1000)
    rows = [
        {
            "copies": count,
            "raw_topology_total_e": source_charge * count,
            "intermol_total_e": intermol_charge * count,
            "amber2lammps_raw_target_total_e": amber_raw_charge * count,
            "amber2lammps_neutral_target_total_e": amber_normalized_charge * count,
        }
        for count in counts
    ]
    return {
        "fixture": "GAFF2 epon",
        "atoms_per_template": len(topology.atoms),
        "target_charge_e": 0.0,
        "source_topology_charge_e": source_charge,
        "intermol_emitted_charge_e": intermol_charge,
        "amber2lammps_raw_target_charge_e": amber_raw_charge,
        "amber2lammps_neutral_target_charge_e": amber_normalized_charge,
        "amber2lammps_uniform_shift_per_atom_e": -source_charge / len(topology.atoms),
        "replication": rows,
        "interpretation": (
            "Charge-preserving conversion retains the per-template residual, so its "
            "system total grows linearly with the number of copies. AMBER2LAMMPS "
            "normalization removes that residual when the requested target is zero."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = generate_results()
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
