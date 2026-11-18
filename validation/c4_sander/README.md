# Zn(II) 12-6-4 validation against `sander`

This fixture validates AMBER2LAMMPS 12-6-4 output against AmberTools `sander`
for a neutral periodic system containing one Zn(II) ion, 339 TIP3P waters, and
two chloride counterions. The Zn(II) Lennard-Jones parameters are the Li/Merz
TIP3P 12-6-4 set distributed with AmberTools. ParmEd's `add12_6_4` action
generates the topology's `LENNARD_JONES_CCOEF` array. The divalent-ion model is
described by Li and Merz, *J. Chem. Theory Comput.* **10**, 289–297 (2014),
[doi:10.1021/ct400751u](https://doi.org/10.1021/ct400751u).

Both engines use the same PDB-rounded coordinates, cell, 9 Angstrom cutoff,
Ewald coefficient, and no dispersion-tail correction. The validation compares
bonded, van der Waals, electrostatic, and total potential energies and every
Cartesian force component. Explicit tolerances and observed reference
differences are stored in `expected_results.json`.

Run the committed fixture with:

```console
python validate_sander.py \
  --sander "$CONDA_PREFIX/bin/sander" \
  --lammps /path/to/lepton-enabled/lmp
```

LAMMPS must include `MOLECULE`, `EXTRA-MOLECULE`, `KSPACE`, and `LEPTON`.

## Rebuilding the source topology

From `source/`, with AmberTools activated:

```console
tleap -f build_tleap.in
parmed -i parmed.in
```

The build recipes use the AmberTools-distributed
`frcmod.ions234lm_1264_tip3p` and `lj_1264_pol.dat` files. `system.pdb` and
`system_1264.prmtop` are committed so CI does not depend on stochastic ion
placement during `tleap` regeneration.
