# 12-6-4 C4 LEPTON smoke validation

This validation creates a temporary AMBER topology containing one nonzero
`LENNARD_JONES_CCOEF`, converts it, and executes a LAMMPS `run 0`. It verifies
that the generated `hybrid/overlay` and Lepton `-C4/r^4` coefficient are
accepted by a LAMMPS build with `PKG_LEPTON=on` and produce finite energies.

The fixture is deliberately synthetic. It tests executable compatibility of
the converter's 12-6-4 output; it is not an independent scientific validation
against an AMBER 12-6-4 ion model.

Run it with:

```console
python validation/c4_lepton/validate.py --lammps /path/to/lmp
```

CI builds the pinned LAMMPS release with LEPTON and runs this check on Linux.
