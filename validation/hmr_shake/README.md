# HMR mass and TIP3P SHAKE validation

This validation reuses the public ff14SB/TIP3P/NaCl fixture in
`validation/amber_protein_water_ions/`. At runtime it applies ParmEd's default
hydrogen mass repartitioning to the ACE-ALA-NME peptide while leaving TIP3P
water unchanged, converts that topology, and checks every emitted LAMMPS atom
mass against the HMR AMBER topology.

The check then:

1. verifies that repartitioning changed 18 peptide atom masses while preserving
   the system's total mass;
2. compares bonded, van der Waals, electrostatic, and total energies against
   AmberTools `sander`;
3. compares every Cartesian force component;
4. derives the TIP3P atom and bond type IDs from the generated data rather than
   hard-coding them;
5. converts AMBER's three-bond water triangles to the representation required
   by LAMMPS `fix shake`; and
6. executes a LAMMPS `run 0` with SHAKE initialized for all 380 waters.

Run it with:

```console
python validation/hmr_shake/validate.py \
  --sander "$CONDA_PREFIX/bin/sander" \
  --lammps /path/to/lmp
```

LAMMPS must include `MOLECULE`, `EXTRA-MOLECULE`, `KSPACE`, and `RIGID`.
The fixture validates mass preservation, potential-energy/force equivalence,
and constraint initialization. It does not claim to benchmark long-timescale
integration stability or justify a particular HMR timestep.
