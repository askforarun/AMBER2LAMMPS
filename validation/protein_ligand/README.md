# ff14SB peptide + GAFF2 ligand validation

This non-covalent fixture combines the existing ACE-ALA-NME ff14SB peptide
with the existing GAFF2 ethanol ligand in a 32 Angstrom periodic cell. The
driver derives a dry peptide topology from the first three residues of the
committed protein/water/ion fixture, so source force-field files are not
duplicated.

```bash
python validate_sander.py \
  --sander "$CONDA_PREFIX/bin/sander" \
  --lammps /path/to/lmp
```

It regenerates separate peptide and ligand templates, the combined PDB, the
AMBER reference and LAMMPS files in ignored `work/`, then compares all reported
energy components and every Cartesian force component. This validates a
non-covalent protein-ligand system; covalent protein-ligand linkages remain a
separate unsupported capability.

The committed run has a total-energy difference of 0.000104 kcal/mol, maximum
force-component difference of 0.000289 kcal/mol/Angstrom, and force RMSE of
0.0000688 kcal/mol/Angstrom.
