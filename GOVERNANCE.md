# Governance and support

AMBER2LAMMPS is currently a maintainer-led project. Arun Sridhar is the lead
maintainer and has final responsibility for releases, scientific support
claims and repository administration.

## Decisions

- Routine fixes and documentation changes are reviewed through pull requests.
- New scientific capabilities require tests, reference validation and an
  update to the support matrix.
- When contributors disagree about a scientific mapping, the AMBER reference
  implementation and reproducible energy/force evidence take precedence.
- Backward-incompatible CLI or output changes require a documented migration
  path and a major-version decision under semantic versioning.

## Releases

Releases use annotated `vMAJOR.MINOR.PATCH` Git tags, matching Python package
versions and GitHub Releases. User-visible changes are recorded in
`CHANGELOG.md`. Release candidates must pass all required CI checks.

## Support expectations

The project is provided without warranty under the MIT license. Maintainers
aim to acknowledge reproducible bug reports within two weeks, but no response
or resolution time is guaranteed. Supported behavior is limited to the matrix
in `docs/support-matrix.md`; entries marked planned or unsupported are not
promises of immediate implementation.

## Adding maintainers

Sustained contributors may be invited as maintainers based on technical and
scientific contributions, constructive review, and demonstrated commitment to
the fail-closed validation policy. Governance changes are made publicly by
pull request.
