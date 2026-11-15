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
   generation recipes and provenance.
3. Represent different GAFF2 and GLYCAM 1–4 scaling exactly in a non-covalent
   mixed system.
4. Validate pure GAFF2, pure GLYCAM, pure protein, protein–ligand,
   protein–carbohydrate, GAFF2–GLYCAM, and three-family mixtures.
5. Add explicit solvent, ions, periodic cells, and HMR variants.
6. Compare AMBER and LAMMPS bonded terms, van der Waals, direct and reciprocal
   electrostatics, total potential energy, and per-atom forces in CI.
7. Treat covalent protein–GLYCAM and GAFF2–GLYCAM linkages as a separate
   capability requiring linkage-specific reference systems.

## Current branch status

The converter now derives uniform SCEE/SCNB scaling from active proper
dihedrals and writes the corresponding `special_bonds` command. A topology with
multiple active profiles, or a mixture whose topologies require different
profiles, is detected and rejected. Exact mixed-profile output is the next
implementation milestone; detection alone is not support.

## Publication target

The intended JOSS title is **AMBER2LAMMPS: Validated AMBER-family force-field
conversion for LAMMPS**. Submission follows only after the supported-system
matrix has energy-and-force validation, documented research use, contribution
guidelines, a changelog, governance/support expectations, and a tagged archive.
