# Supported systems and validation status

AMBER2LAMMPS targets classical, fixed-charge members of the AMBER force-field
family. It does not claim to translate arbitrary LAMMPS force fields. This
matrix separates code paths that merely exist from those validated against an
independent AMBER engine.

## Status definitions

- **Validated**: energies and Cartesian forces are compared against AmberTools
  `sander` with explicit numerical tolerances in an automated or reproducible
  validation.
- **Regression-tested**: behavior and generated syntax are covered by automated
  tests, but no complete AMBER energy-and-force reference exists yet.
- **Implemented**: the converter has an explicit representation and safety
  checks, but broader reference validation is still required.
- **Unsupported**: conversion stops with a clear error rather than silently
  substituting a different model.

## Current matrix

| Capability | Status | Evidence or boundary |
|---|---|---|
| GAFF2 + GLYCAM non-covalent mixture | **Validated** | `validation/amber_family_mix/`; required `sander/LAMMPS energy and force validation` CI job |
| GAFF2 small molecules | Regression-tested | Ethanol, epon, aspirin, benzene and mixed-molecule fixtures; GAFF2 terms are also exercised by the validated mixed system |
| GLYCAM carbohydrates | **Validated** for standalone glucose | Pure GLYCAM06j-1 energy/force comparison in `validation/glycam_glucose/`; GLYCAM terms are also exercised by the validated mixed system |
| Mixed SCEE/SCNB profiles | **Validated** for separate GAFF2 + GLYCAM templates | Repeated `pair_style hybrid` instances, per-instance `pair_modify ... special`, complete I-J matrix; a combined mixed-profile `prmtop` fails closed |
| Uniform SCEE/SCNB profile | **Validated** for GLYCAM06j-1 and ff14SB | Active proper-dihedral values are extracted from the topology and emitted as `special_bonds` |
| Lorentz-Berthelot LJ mixing | **Validated** in the mixed fixture | Explicit like-like and cross-family energy/force comparison |
| AMBER off-diagonal A/B overrides | Regression-tested | Detected and emitted as explicit I-J coefficients |
| 12-6-4 C4 terms | **Validated** for Zn(II) in periodic TIP3P water | Li/Merz Zn(II) 12-6-4 energy and all Cartesian forces compared with `sander` in `validation/c4_sander/`; multiple topology files fail closed because cross-topology C4 terms are undefined; requires LAMMPS LEPTON |
| Mixed SCEE/SCNB plus 12-6-4 | **Unsupported** | Fails before output because the combined representation is not yet validated |
| Hydrogen mass repartitioning | **Validated** for an ff14SB peptide in TIP3P/NaCl | Every emitted mass is checked against the HMR topology; total mass, `sander` energies/forces, and LAMMPS TIP3P SHAKE initialization are validated in `validation/hmr_shake/`. Long-timescale integration and timestep selection are not benchmarked |
| Orthorhombic periodic `CRYST1` cells | **Validated** in the mixed fixture | Same 30 Angstrom cell used by `sander` and LAMMPS |
| Triclinic `CRYST1` cells | **Unsupported** | Clear error; tilt factors are not inferred from PDB coordinates |
| ff14SB capped peptide | **Validated** | ACE-ALA-NME energy and all Cartesian forces compared with `sander` |
| Larger folded proteins | Planned | The peptide covers ff14SB terms, but a representative folded-protein fixture remains desirable |
| Protein + GAFF2 ligand | **Validated** for non-covalent ff14SB peptide + ethanol | Energy and all Cartesian forces compared with `sander` in `validation/protein_ligand/` |
| Protein + GLYCAM carbohydrate | Planned | Requires protein fixture and linkage/non-covalent cases |
| TIP3P water with Na+/Cl- | **Validated** | 380-water periodic ff14SB fixture; energy/force validation and SHAKE example |
| Covalent bonds between separate input topologies | **Not represented or automatically detected** | Input topologies are independent molecule templates; PDB proximity alone cannot establish intended covalent connectivity |
| Reactive, polarizable or many-body models | Out of scope | The converter targets classical fixed-charge AMBER-family topologies |

## Validation policy

A capability moves to **Validated** only when its fixture:

1. can be regenerated from documented source inputs;
2. runs on both AMBER and LAMMPS using the same coordinates and cell;
3. compares bonded, van der Waals, electrostatic and total energies;
4. compares every Cartesian force component;
5. uses committed, justified tolerances; and
6. runs in CI when its dependencies can be installed reproducibly.

Passing a LAMMPS syntax check alone is not scientific validation.

Topologies without an active proper dihedral return no SCEE/SCNB profile. This
is expected for water, monatomic ions, and other models with no 1-4 owner; no
`special_bonds` command is needed for interactions that do not exist. A
topology that contains active 1-4 owners with conflicting profiles is rejected.
