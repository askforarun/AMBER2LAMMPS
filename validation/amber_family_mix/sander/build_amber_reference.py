#!/usr/bin/env python3
"""Build one AMBER reference topology for the GAFF2/GLYCAM mixture."""

from pathlib import Path

import parmed as pmd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = HERE.parent / "source"

ethanol = pmd.load_file(str(ROOT / "ethanol.prmtop"))
glucose = pmd.load_file(str(SOURCE / "glucose.prmtop"))
pdb = pmd.load_file(str(SOURCE / "ethanol_glucose.pdb"))
coordinates = pdb.coordinates

system = ethanol + glucose
system.coordinates = coordinates

# Use the same CRYST1 cell that AMBER2LAMMPS reads from the packed PDB.
if pdb.box is None:
    raise ValueError("ethanol_glucose.pdb must contain a CRYST1 record")
system.box = pdb.box

system.save(str(HERE / "amber_reference.prmtop"), overwrite=True)
system.save(str(HERE / "amber_reference.rst7"), overwrite=True)

print(f"Wrote {len(system.atoms)} atoms with box {system.box.tolist()}")
