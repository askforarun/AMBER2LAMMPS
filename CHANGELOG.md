# Changelog

All notable changes to AMBER2LAMMPS are documented here. The project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Experimental AMBER-family conversion with topology-derived SCEE/SCNB 1-4
  scaling.
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

- Atom types with otherwise identical parameters are separated when their
  topology-specific 1-4 scaling profiles differ.
- InterMol validation is retained as an optional secondary cross-check;
  `sander` is now the primary AMBER reference.

### Fixed

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

[Unreleased]: https://github.com/askforarun/AMBER2LAMMPS/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/askforarun/AMBER2LAMMPS/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/askforarun/AMBER2LAMMPS/releases/tag/v1.0.0
