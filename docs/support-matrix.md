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
| GLYCAM carbohydrates | Regression-tested | Neutral glucose fixture; GLYCAM terms are exercised by the validated mixed system |
| Mixed SCEE/SCNB profiles | **Validated** for GAFF2 + GLYCAM | Repeated `pair_style hybrid` instances, per-instance `pair_modify ... special`, complete I-J matrix |
| Uniform SCEE/SCNB profile | Regression-tested | Active proper-dihedral values are extracted from the topology and emitted as `special_bonds` |
| Lorentz-Berthelot LJ mixing | **Validated** in the mixed fixture | Explicit like-like and cross-family energy/force comparison |
| AMBER off-diagonal A/B overrides | Regression-tested | Detected and emitted as explicit I-J coefficients |
| 12-6-4 C4 terms | Implemented | `hybrid/overlay` + Lepton output and unit tests; requires LAMMPS LEPTON package |
| Mixed SCEE/SCNB plus 12-6-4 | **Unsupported** | Fails before output because the combined representation is not yet validated |
| Hydrogen mass repartitioning | Regression-tested | Per-atom mass variants and unchanged LJ parameters are tested |
| Orthorhombic periodic `CRYST1` cells | **Validated** in the mixed fixture | Same 30 Angstrom cell used by `sander` and LAMMPS |
| Triclinic `CRYST1` cells | **Unsupported** | Clear error; tilt factors are not inferred from PDB coordinates |
| ff14SB capped peptide | **Validated** | ACE-ALA-NME energy and all Cartesian forces compared with `sander` |
| Larger folded proteins | Planned | The peptide covers ff14SB terms, but a representative folded-protein fixture remains desirable |
| Protein + GAFF2 ligand | Planned | Requires protein fixture and mixed-system validation |
| Protein + GLYCAM carbohydrate | Planned | Requires protein fixture and linkage/non-covalent cases |
| TIP3P water with Na+/Cl- | **Validated** | 380-water periodic ff14SB fixture; energy/force validation and SHAKE example |
| Covalent bonds between separate input topologies | **Unsupported** | Input topologies are independent molecule templates; cross-topology connectivity is not invented |
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
