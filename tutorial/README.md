# AMBER2LAMMPS Beamer tutorial

This folder contains a presentation-style tutorial for AMBER2LAMMPS 1.2.0.
It covers installation, single- and multiple-template conversion, charge
targets, AMBER 1-4 scaling, LAMMPS input construction, the Python API,
provenance, the bounded MCP workflow, validation evidence, and supported
boundaries.

Every command in the slides runs from the repository root on committed files:
`ethanol.prmtop`/`ethanol.pdb`, `multi_ethanol.pdb`, and the GAFF2 + GLYCAM
mixture in `validation/amber_family_mix/source/`. The LAMMPS input and the
`run 0` energies shown were checked with LAMMPS 22 Jul 2025.

Build the slides with:

```bash
cd tutorial
pdflatex -interaction=nonstopmode -halt-on-error amber2lammps-tutorial.tex
pdflatex -interaction=nonstopmode -halt-on-error amber2lammps-tutorial.tex
```

The compiled PDF is intentionally kept with the source so users can read the
tutorial without a LaTeX installation.
