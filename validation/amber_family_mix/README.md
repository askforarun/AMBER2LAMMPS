# GAFF2 ethanol + GLYCAM glucose mixture

This experimental fixture exercises two AMBER-family 1–4 scaling profiles in
one non-covalent periodic system. Ethanol uses GAFF2 SCEE/SCNB 1.2/2.0;
glucose uses GLYCAM06j-1 1.0/1.0. The generated LAMMPS include uses two
instances of `lj/cut/coul/long` under `pair_style hybrid`, per-instance special
weights, and all 66 I–J assignments for its 11 atom types.

## Rebuild the inputs

Requirements: AmberTools 25 (including GLYCAM06j-1), Packmol, ParmEd, and
LAMMPS with the KSPACE package.

```bash
conda run -n AmberTools25 tleap -f tleap_glucose.in
packmol < packmol.in
python add_cryst1.py

python ../../amber_to_lammps.py ethanol_glucose.data ethanol_glucose.parm \
  ethanol_glucose.pdb \
  -t ../../ethanol.prmtop glucose.prmtop -c 1 1 \
  --charges 0.0000010021786493 0

lmp_serial -in ethanol_glucose.input
```

The nonzero ethanol target preserves the raw charge from the committed
topology; it prevents charge normalization from contaminating the comparison.
`ethanol_glucose.pdb` has a 30 Angstrom orthorhombic `CRYST1` cell shared by
both engines.

## Independent AMBER reference

```bash
python build_amber_reference.py
python validate_sander.py \
  --sander "$CONDA_PREFIX/bin/sander" \
  --lammps /path/to/lmp
```

The script rebuilds the combined AMBER reference, runs `sander` and LAMMPS on
the same coordinates and cell, then compares bonded, van der Waals,
electrostatic and total energies plus every Cartesian force component. The
committed result has a total-energy difference of about 0.00007 kcal/mol, a
maximum force-component difference of 0.0437 kcal/mol/Angstrom, and a force
RMSE of 0.0135 kcal/mol/Angstrom. See `expected_results.json` for the numerical
values and pass/fail tolerances.

The LAMMPS input matches AMBER's reported 30x30x30 PME mesh, interpolation
order 4 and Ewald coefficient 0.42001. It also uses `dielectric
1.0000346111225924`, the ratio between the LAMMPS real-units Coulomb factor
(332.06371) and AMBER's factor (332.0522173). Without that explicit conversion,
the small difference between physical constants obscures the force-field
comparison.

The example is intentionally non-covalent. Cross-family covalent linkages and
mixed 12-6-4 systems remain unsupported and are rejected rather than silently
approximated.
