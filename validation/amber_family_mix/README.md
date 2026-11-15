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
conda run -n AmberTools25 sander -O \
  -i sander.in -o sander.out \
  -p amber_reference.prmtop -c amber_reference.rst7 \
  -r sander_restart.rst7 -inf sander.mdinfo
```

For the committed coordinates, both programs report 2.5296 kcal/mol for bonds,
3.2474 for angles, and 2.9723 for dihedrals. The combined AMBER 1–4 and
ordinary van der Waals energy agrees with LAMMPS at the printed precision.
The current total-energy difference is about 0.049 kcal/mol, isolated to the
real/reciprocal electrostatic split. This fixture therefore proves the hybrid
representation runs and applies both scaling profiles, but it is not yet the
final publication-grade validation. The next step is matching the Ewald/PPPM
settings and adding per-atom force comparison in CI.

The example is intentionally non-covalent. Cross-family covalent linkages and
mixed 12-6-4 systems remain unsupported and are rejected rather than silently
approximated.
