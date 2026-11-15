# ff14SB peptide with TIP3P water and NaCl

This validation contains a capped ACE-ALA-NME peptide parameterized with
ff14SB, 380 explicit TIP3P waters, one Na+ and one Cl-. The 1,164-atom
orthorhombic periodic system exercises protein bonded terms, rigid-water bond
topology, ion parameters, exclusions, 1-4 scaling, LJ mixing and long-range
electrostatics.

```text
source/          tleap recipe and generated AMBER topology/PDB
sander/          exact-coordinate reference builder and force input
amber2lammps/    generated data/parameter files and LAMMPS inputs
validate_sander.py
expected_results.json
```

## Run the validation

Requirements are AmberTools 25, ParmEd, NumPy, SciPy and LAMMPS with MOLECULE,
EXTRA-MOLECULE, KSPACE and RIGID packages. RIGID supplies `fix shake` for the
TIP3P dynamics check.

```bash
cd validation/amber_protein_water_ions
python validate_sander.py \
  --sander "$CONDA_PREFIX/bin/sander" \
  --lammps /path/to/lmp
```

The validator regenerates the AMBER2LAMMPS files from the committed topology
and PDB, rebuilds a `sander` reference using those exact PDB coordinates and
cell, runs both engines, and compares six energy groups plus all 3,492 force
components. `sander` uses a converged 96x96x90 PME grid (order 6); LAMMPS uses
converged Ewald. Both use coefficient 0.42001 and the LAMMPS input compensates
for the small difference between AMBER and LAMMPS Coulomb conversion factors.

Observed differences for the committed fixture are:

| Quantity | Difference |
|---|---:|
| Bond energy | 0.0000352 kcal/mol |
| Angle energy | 0.0000211 kcal/mol |
| Dihedral energy | 0.0000308 kcal/mol |
| van der Waals energy | 0.01145 kcal/mol |
| Coulomb energy | 0.00620 kcal/mol |
| Total potential energy | 0.01820 kcal/mol |
| Maximum force component | 0.005188 kcal/mol/Angstrom |
| Force RMSE | 0.000164 kcal/mol/Angstrom |

The small van der Waals energy difference is consistent with the engines'
analytical dispersion-tail implementations. Numerical tolerances are committed
in `expected_results.json` and enforced by the validator.

## TIP3P constraints for dynamics

AMBER represents TIP3P as a rigid triangle. The topology contains HW-OW and
HW-HW bonds but no water angle. A LAMMPS dynamics input must constrain those
bonds with `fix shake`; otherwise it runs a flexible model that is not AMBER
TIP3P. `amber2lammps/system_shake.input` demonstrates the correct command for
this fixture:

```lammps
fix rigid_tip3p all shake 1.0e-6 200 0 b 8 9
```

Bond type IDs depend on the complete input topology, so users must verify the
generated `bond_coeff` comments rather than copying `8 9` to unrelated systems.

## Regenerate the AMBER source

```bash
cd validation/amber_protein_water_ions/source
conda run -n AmberTools25 tleap -f tleap_protein_water_ions.in
```

tleap's `addionsrand` placement may vary between invocations. The committed PDB
and topology define the reference used by CI.
