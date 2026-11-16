# Contributing to AMBER2LAMMPS

Contributions through GitHub issues and pull requests are welcome. Scientific
changes must preserve AMBER semantics or fail clearly when exact conversion is
not possible.

## Development workflow

1. Open an issue describing the problem, affected force field and expected
   AMBER behavior for scientific changes.
2. Create a focused branch from `main`.
3. Add or update regression tests.
4. Run `python -m unittest -v` and `git diff --check`.
5. For changes affecting energies or forces, add a reproducible AMBER
   reference fixture and compare against `sander`.
6. Open a pull request explaining the scientific mapping, validation results,
   limitations and compatibility impact.

## Scientific acceptance criteria

A pull request that adds a force field, combination rule or interaction must
include:

- provenance and generation instructions for redistributable fixtures;
- an explicit mapping from AMBER terms to LAMMPS commands;
- energy-component and Cartesian-force comparisons where applicable;
- justified numerical tolerances;
- a negative test showing unsupported combinations fail clearly; and
- documentation in `docs/support-matrix.md`.

Do not commit private molecular structures, credentials, PyPI tokens, licensed
files that cannot be redistributed, or generated build environments.

## Reporting problems

Use a GitHub issue and include the command, AMBER2LAMMPS version, Python and
LAMMPS versions, relevant error text, and a minimal redistributable input when
possible. Security-sensitive reports should not include secrets in an issue;
contact the maintainer privately through the address in the package metadata.

## AI-assisted contributions

Contributors must disclose material generative-AI assistance in the pull
request and confirm that a human reviewed, tested and takes responsibility for
the submitted code and scientific claims.
