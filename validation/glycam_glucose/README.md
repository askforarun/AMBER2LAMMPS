# Pure GLYCAM06j-1 glucose validation

This fixture validates a neutral, standalone alpha-D-glucopyranose molecule
parameterized with GLYCAM06j-1. It reuses the topology and coordinates generated
by the documented `tleap` recipe in `validation/amber_family_mix/source/` and
places them in a 30 Angstrom orthorhombic periodic cell.

The validation regenerates the AMBER2LAMMPS files, runs AmberTools `sander` and
LAMMPS on identical coordinates and cell, and compares bond, angle, dihedral,
van der Waals, electrostatic and total potential energies plus every Cartesian
force component. Generated files are kept in the ignored `work/` directory.

## Run

Requirements: AmberTools 25, ParmEd, NumPy, SciPy, and LAMMPS with the
MOLECULE, EXTRA-MOLECULE, and KSPACE packages.

```bash
python validate_sander.py \
  --sander "$CONDA_PREFIX/bin/sander" \
  --lammps /path/to/lmp
```

The committed tolerances and observed reference values are recorded in
`expected_results.json`. This fixture validates the uniform GLYCAM SCEE/SCNB
profile of 1.0/1.0 independently of the heterogeneous GAFF2 + GLYCAM fixture.
The committed run has a total-energy difference of 0.00018 kcal/mol, maximum
force-component difference of 0.000530 kcal/mol/Angstrom, and force RMSE of
0.000118 kcal/mol/Angstrom.
