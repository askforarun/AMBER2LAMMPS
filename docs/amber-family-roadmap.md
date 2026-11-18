# AMBER-family support roadmap

## Scope

AMBER2LAMMPS will target validated classical molecular-mechanics systems built
from AMBER protein force fields, GAFF2 small molecules, GLYCAM carbohydrates,
AMBER water models, and compatible ion parameters. It will not claim generic
support for unrelated LAMMPS force fields or reactive and many-body potentials.

## Safety rule

Every AMBER term that affects energy or forces must either be represented
exactly in the generated LAMMPS model or produce a clear error before output is
written. The converter must not silently substitute a default force-field
behavior.

## Milestones

1. Extract active proper-dihedral SCEE/SCNB values from every topology, ignore
   improper/repeated-term placeholders, emit uniform scaling automatically, and
   reject incompatible profiles.
2. Add redistributable GLYCAM and AMBER protein fixtures with documented
   generation recipes and provenance. (GLYCAM and ff14SB peptide fixtures added.)
3. Represent different GAFF2 and GLYCAM 1–4 scaling exactly in a non-covalent
   mixed system. (Implemented and validated with repeated hybrid instances.)
4. Validate pure GAFF2, pure GLYCAM, pure protein, protein–ligand,
   protein–carbohydrate, GAFF2–GLYCAM, and three-family mixtures.
5. Add explicit solvent, ions, periodic cells, and HMR variants. (TIP3P,
   Na+/Cl-, an orthorhombic cell, HMR mass mapping, and TIP3P SHAKE
   initialization are validated.)
6. Compare AMBER and LAMMPS bonded terms, van der Waals, direct and reciprocal
   electrostatics, total potential energy, and per-atom forces in CI.
7. Treat covalent protein–GLYCAM and GAFF2–GLYCAM linkages as a separate
   capability requiring linkage-specific reference systems.

## Current branch status

The converter derives SCEE/SCNB scaling from active proper dihedrals. Uniform
systems write one `special_bonds` command. A non-covalent mixed GAFF2/GLYCAM
system now writes repeated `pair_style hybrid` instances with per-instance 1–4
weights and a complete explicit I–J matrix. The ethanol/glucose fixture is
validated directly against AmberTools `sander`, including energy components
and all Cartesian forces. With matched PME settings and Coulomb conversion
factors, total energy differs by 0.00010 kcal/mol, maximum force-component
difference is 0.000543 kcal/mol/Angstrom, and force RMSE is 0.0000817
kcal/mol/Angstrom. Standalone GLYCAM06j-1 glucose is independently validated.
A periodic Li/Merz Zn(II) 12-6-4 TIP3P system is also validated against
`sander` energies and forces using LEPTON-enabled LAMMPS. Mixed scaling with
12-6-4 terms remains fail-closed.

## Publication target

The intended JOSS title is **AMBER2LAMMPS: A validated, provenance-aware, and
agent-ready framework for AMBER-family simulations in LAMMPS**. Here,
"agent-ready" means a structured and bounded interface that exposes compact,
machine-readable conversion and validation results while retaining complete
provenance and large molecular artifacts on disk. It does not mean that the
converter contains an LLM or autonomously designs simulations. Submission follows only after the supported-system
matrix has energy-and-force validation, documented research use, contribution
guidelines, a changelog, governance/support expectations, and a tagged archive.
