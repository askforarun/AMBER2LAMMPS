# Changelog

All notable changes to AMBER2LAMMPS are documented here. The project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.0] - 2026-11-18

### Added

- Standalone GLYCAM06j-1 glucose energy-and-force validation against
  AmberTools `sander`.
- A non-covalent ff14SB peptide + GAFF2 ligand energy-and-force validation
  fixture against AmberTools `sander`.
- A Li/Merz Zn(II) 12-6-4 energy-and-force fixture against `sander`, using a
  pinned LEPTON-enabled LAMMPS build.
- Default ParmEd hydrogen-mass-repartitioning checks covering every emitted
  mass, total-mass conservation, and TIP3P SHAKE initialization reported by
  LAMMPS.
- The first internal package-split boundary for force-field model objects and
  topology scaling analysis, preserving the public CLI and Python API.
- AMBER-family conversion with topology-derived SCEE/SCNB 1-4 scaling and
  explicit fail-closed boundaries.
- Exact non-covalent GAFF2/GLYCAM mixed scaling using repeated LAMMPS
  `pair_style hybrid` instances, per-instance `pair_modify ... special`
  settings, and an explicit I-J coefficient matrix.
- A reproducible GAFF2 ethanol + GLYCAM glucose validation bundle, separated
  into source, `sander`, and AMBER2LAMMPS inputs and outputs.
- Direct AmberTools `sander` energy-and-force validation as a required CI job.
- A reproducible ff14SB peptide + TIP3P water + NaCl energy-and-force
  validation, including a LAMMPS SHAKE example for rigid water dynamics.
- A supported-system matrix, scientific design document, JOSS-readiness map,
  contribution guide, governance policy, research-use record, and AI-use
  disclosure draft.

### Changed

- Uniform-system output now selects AMBER's Lorentz-Berthelot mixing rule
  explicitly.
- Nonbonded epsilon and sigma values now derive from the full-precision AMBER
  A/B coefficient arrays instead of four-decimal formatted diagnostic text.
- Atom types with otherwise identical parameters are separated when their
  topology-specific 1-4 scaling profiles differ.
- InterMol validation is retained as an optional secondary cross-check;
  `sander` is now the primary AMBER reference.

### Fixed

- Scientific CI validations now run as independent fail-fast steps, and the
  protein-ligand tleap recipe is tracked in the repository.
- MCP `run 0` validation serializes workspace access, validates provenance
  fields, restricts output-path characters, and rejects executable commands,
  unapproved pair styles, or non-generated nonbonded-control arguments.
- Removed machine-specific absolute paths from generated data and parameter
  comments, and report unsupported 10-12 H-bond diagonals explicitly.
- Removed a false warning for correctly sized combined PDB files.
- Corrected the GLYCAM reducing-end fixture so its cap is covalently bonded.

## [1.1.0] - 2026-11-15

### Added

- Optional MCP server for invoking conversions from compatible AI clients.
- MCP input validation and regression tests.
- Narrated demonstration video and animated README preview.
- GitHub Actions testing across supported Python versions.
- Pinned InterMol/LAMMPS energy validation with runnable validation artifacts.
- Optional, sanitized LAMMPS type labels.
- Orthogonal PDB `CRYST1` cell handling.

### Changed

- Raised the minimum supported Python version to 3.9.
- Debug files now use output-specific names beside the data file.

### Fixed

- Prevented debug cleanup from deleting unrelated user files.
- Preserved per-atom masses for hydrogen mass repartitioned topologies.
- Replaced silent unknown-type fallback with explicit errors.
- Hardened invalid type labels and placeholder or malformed `CRYST1` records.

## [1.0.0] - 2026-11-14

### Added

- Initial packaged PyPI release of `amber-to-lammps`.
- Command-line and Python conversion interfaces.
- Multiple-topology and replicated-molecule conversion.
- AMBER topology, PDB coordinate, charge, bonded-term, and nonbonded-term
  conversion to LAMMPS data and parameter files.
- Initial regression tests and packaging metadata.

[Unreleased]: https://github.com/askforarun/AMBER2LAMMPS/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/askforarun/AMBER2LAMMPS/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/askforarun/AMBER2LAMMPS/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/askforarun/AMBER2LAMMPS/releases/tag/v1.0.0
