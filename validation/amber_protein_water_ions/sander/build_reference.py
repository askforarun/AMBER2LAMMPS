#!/usr/bin/env python3
"""Build the sander reference from the exact PDB coordinates used by LAMMPS."""

from pathlib import Path

import parmed as pmd


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "source"

topology = pmd.load_file(str(SOURCE / "system.prmtop"))
pdb = pmd.load_file(str(SOURCE / "system.pdb"))
if pdb.box is None:
    raise ValueError("system.pdb must contain the periodic CRYST1 cell")
if len(topology.atoms) != len(pdb.atoms):
    raise ValueError("PDB and topology atom counts differ")

topology.coordinates = pdb.coordinates
topology.box = pdb.box
topology.save(str(HERE / "reference.prmtop"), overwrite=True)
topology.save(str(HERE / "reference.rst7"), overwrite=True)
print(f"Wrote {len(topology.atoms)} atoms with box {topology.box.tolist()}")
