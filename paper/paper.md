---
title: 'AMBER2LAMMPS: A validated, provenance-aware, and agent-ready framework for AMBER-family simulations in LAMMPS'
tags:
  - Python
  - molecular dynamics
  - AMBER
  - LAMMPS
  - interoperability
  - provenance
authors:
  - name: Arun Srikanth Sridhar
    orcid: 0000-0002-1984-3513
    affiliation: 1
    corresponding: true
    email: askforarun@gmail.com
affiliations:
  - index: 1
    name: School of Energy, Geoscience, Infrastructure and Society, Heriot-Watt University, Edinburgh, United Kingdom, EH14 4AS
date: 16 November 2026
bibliography: paper.bib
---

# Summary

Molecular simulations are commonly prepared with one software ecosystem and
executed with another. AMBER topology files describe not only atoms and bonds,
but also masses, partial charges, nonbonded coefficients, exclusions, and
force-field-specific scaling of interactions between atoms separated by three
bonds. Translating these choices incorrectly can produce a LAMMPS input that
runs successfully while representing a different molecular model.

AMBER2LAMMPS is an open-source Python package that converts fully parameterized
AMBER topology and PDB coordinate files into LAMMPS data and parameter files.
Its scope is classical fixed-charge AMBER-family models, including validated
GAFF2 small-molecule, GLYCAM06 carbohydrate, ff14SB peptide, TIP3P water, and
ion fixtures. It provides a command-line interface, Python API, and a bounded
Model Context Protocol (MCP) interface for external software agents. Unsupported
scientific combinations fail before output rather than silently adopting a
fallback model.

# Statement of need

LAMMPS [@thompson2022lammps] is flexible enough to represent many AMBER-family
potentials, but that flexibility places responsibility for an exact mapping on
the user. A recurring difficulty is heterogeneous simulation: GAFF2 conventionally
uses SCEE/SCNB values of 1.2/2.0, whereas GLYCAM uses 1.0/1.0. LAMMPS normally
expresses these weights through one global `special_bonds` command. Additional
complications include hydrogen mass repartitioning, off-diagonal Lennard-Jones
overrides, residual molecular charges, periodic cells, and 12-6-4 C4 terms.

AMBER2LAMMPS targets researchers who already have parameterized AMBER systems
and need a traceable LAMMPS representation. It reads topology-native parameters
rather than assigning atom types or charges. For separate GAFF2 and GLYCAM
templates, it constructs repeated hybrid pair-style instances with per-instance
1-4 weights and a complete cross-type coefficient matrix. A combined topology
containing conflicting active scaling profiles is rejected because ownership of
individual interactions cannot be inferred safely from a global profile.

# State of the field

Moltemplate is a broad, extensible LAMMPS model builder with libraries and
workflows spanning atomistic and coarse-grained models [@jewett2021moltemplate].
Its bundled AMBER-oriented material includes GAFF/GAFF2 templates and import
utilities, while atom typing and charge assignment remain external. No bundled
automatic GAFF2--GLYCAM heterogeneous-scaling workflow was identified; a user
could construct an appropriate LAMMPS hybrid model manually.

InterMol is a multi-engine format converter and single-point energy-comparison
tool [@shirts2017intermol]. Its native conversions cover Desmond, GROMACS, and
LAMMPS; its AMBER route uses ParmEd to create an intermediate GROMACS model.
That route preserves input charges in the tested epon fixture and rejects a
combined topology with heterogeneous SCEE values. This is appropriate behavior
for a general converter, not a defect. AMBER2LAMMPS makes a narrower trade-off:
it directly consumes parameterized AMBER topologies and specializes in explicit,
fail-closed AMBER-family mappings.

| Property | Moltemplate | InterMol | AMBER2LAMMPS |
|---|---|---|---|
| Primary role | General LAMMPS model builder | Multi-engine conversion and energy comparison | Direct AMBER-family to LAMMPS conversion |
| AMBER input | GAFF/GAFF2 libraries and external parameter imports; no `prmtop` reader | ParmEd/GROMACS intermediate route | Reads parameterized `prmtop` files through ParmEd |
| Mixed SCEE/SCNB | Manual hybrid construction required; no bundled example found | Combined mixed-SCEE AMBER topology rejected on the GROMACS route | Automatic for separate GAFF2 and GLYCAM templates; combined mixed-profile topology rejected |
| Reference checks | AMBER equivalence is not claimed by its documentation | Cross-engine single-point energies | Required `sander` energy and Cartesian-force CI for validated fixtures |

# Software design

Conversion is separated into topology interpretation, type and interaction
assembly, and deterministic output. Type identities include physical parameters,
so equal AMBER names cannot merge different masses or Lennard-Jones values.
SCEE/SCNB values are derived from active proper dihedrals. Unknown parameters,
unsupported scaling combinations, unsafe cells, and undefined cross-topology C4
interactions raise explicit errors. Orthorhombic `CRYST1` records are supported;
triclinic cells are rejected rather than inferred from insufficient coordinate
metadata.

Scientific validation compares component energies and every Cartesian force
against AmberTools `sander` [@case2023ambertools]. Required CI fixtures cover
standalone GLYCAM06j-1 glucose [@kirschner2008glycam], non-covalent
GAFF2--GLYCAM, ff14SB peptide--GAFF2 ligand [@maier2015ff14sb], and an
ff14SB/TIP3P/NaCl periodic system. Hydrogen mass repartitioning is
regression-tested but is not claimed as independently `sander`-validated.
Synthetic 12-6-4 output completes a LEPTON-enabled LAMMPS `run 0`; this proves
executable compatibility, not AMBER-reference scientific equivalence.

Each optional CLI provenance record—and every MCP conversion record—contains
input and output hashes, software versions, charge decisions, detected scaling
profiles, box metadata, warnings, and the LAMMPS execution contract. MCP tools
return compact structured summaries while retaining complete molecular artifacts
and logs on disk. This “agent-ready” description concerns bounded execution,
traceability, and context efficiency; it does not imply that a language model
selects force fields or establishes scientific correctness.

# Research impact statement

The project supplies installable releases, documented examples, machine-readable
provenance, fail-closed boundaries, and reproducible AMBER-to-LAMMPS energy and
force fixtures. These materials support independent examination of a conversion
rather than requiring trust in a generated input file. A charge-policy benchmark
also measures one epon template and projects its approximately -0.003 e residual
arithmetically to larger copy counts; it does not claim replicated-box conversions
or characterize faithful charge preservation as an InterMol error.

At the time of this manuscript draft, the repository does not yet document a
genuine research calculation using AMBER2LAMMPS or an independent user report.
Consequently, this draft is not ready for JOSS submission. Those evidence gates
must be satisfied with verifiable use rather than replaced by anticipated impact.

# AI usage disclosure

OpenAI Codex (GPT-5 family; exact hosted revision not exposed) assisted with
implementation, refactoring, tests, CI, validation scripts, documentation, and
manuscript drafting. Anthropic Claude and xAI Grok were used for independent
repository-based reviews; their exact hosted model revisions were not recorded.
The human author directed the scientific scope and design, reviewed proposed
changes, and required automated tests plus AMBER and LAMMPS execution before
accepting scientific claims. The author remains responsible for accuracy,
originality, licensing, and the final manuscript. Exact product/model metadata
and the completeness of this disclosure must be confirmed before submission.

# Acknowledgements

No funding statement has yet been supplied. Before submission, the author must
identify all financial support and state the sponsor's role, or explicitly state
that the work received no specific funding. The author reports no conflict of
interest for this draft; this statement must be reconfirmed at submission.

# References
