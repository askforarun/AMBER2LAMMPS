# Scientific design and extension model

## Interoperability problem

AMBER topology files store a molecular model, not merely coordinates. A useful
LAMMPS conversion must preserve bonded parameters, atom masses and charges,
Lennard-Jones tables, exclusions, topology-specific 1-4 scaling and the
periodic cell. The difficult case is a heterogeneous system whose AMBER-family
components use different conventions—for example GAFF2 SCEE/SCNB 1.2/2.0 and
GLYCAM 1.0/1.0—while LAMMPS normally exposes global `special_bonds` weights.

AMBER2LAMMPS solves this recurring interoperability problem by deriving the
model from each `.prmtop`, selecting an explicit LAMMPS representation, and
failing before output when it cannot preserve a term.

Charge handling is explicit rather than implicit. A required per-template
target equal to the topology sum preserves its written charges for direct
engine comparisons; a different target applies a uniform per-atom correction.
The reproducible
[`validation/charge_replication/`](../validation/charge_replication/README.md)
benchmark converts one template and arithmetically projects the measured charge
to larger copy counts. It documents the linear consequence of copying an
unchanged residual and the policy trade-off between preservation and a requested
target; it is not a replicated-box conversion benchmark.

## Scientific mapping

| AMBER concept | LAMMPS representation |
|---|---|
| Atom mass, charge and LJ type | `Masses`, `Atoms`, and `pair_coeff` |
| Harmonic bond and angle | `bond_style harmonic`, `angle_style harmonic` |
| Multi-term proper/improper torsion | `dihedral_style fourier` |
| Uniform SCEE/SCNB | `special_bonds` with reciprocal scale factors |
| Multiple topology-level SCEE/SCNB profiles | Repeated instances under `pair_style hybrid` plus per-instance `pair_modify ... special` |
| Standard LJ cross interactions | Lorentz-Berthelot arithmetic sigma/geometric epsilon |
| Off-diagonal AMBER A/B pair | Explicit I-J coefficient derived from the raw table |
| 12-6-4 C4 term | `hybrid/overlay` with a Lepton `-C4/r^4` term |
| Orthorhombic PDB cell | LAMMPS `0..L` bounds from `CRYST1` |

Mixed-profile output assigns every I-J coefficient explicitly. A pair of atom
types from the same topology profile uses that profile's hybrid instance.
Cross-family atom types use the standard cross coefficient and do not require
1-4 scaling because separate topology templates have no covalent 1-4 path.
Consequently, mixed-profile conversion accepts separate topology templates;
a combined topology containing multiple active profiles fails closed. A
topology with no active proper-dihedral 1-4 owner—such as water or a monatomic
ion—has no scaling profile because there is no corresponding 1-4 interaction.

## Fail-closed rules

The converter raises an error for unsupported or ambiguous physics, including:

- unresolved atom parameters;
- conflicting raw pair tables;
- multiple active scaling profiles inside one topology;
- mixed scaling combined with 12-6-4 terms;
- non-orthogonal PDB cells; and
- output/debug path collisions.

It never falls back to atom type 1 or silently replaces an unsupported
nonbonded model.

## Extension points

The current implementation uses explicit intermediate contexts rather than
embedding system-specific conditions in file-writing loops:

- `SpecialBondScaling` represents topology-level SCEE/SCNB behavior;
- `TopologyContext` owns global type mappings, raw pair tables and scaling
  profiles;
- `ConnectivityContext` owns deduplicated bonded types and connectivity;
- `BoxBounds` and `MoleculeSpan` isolate coordinate/cell and replication logic;
- `topology_14_scaling()` and `find_pair_specific_nonbonded_terms()` analyze
  source physics;
- `build_explicit_pair_coeff_lines()` and
  `build_mixed_scaling_pair_coeff_lines()` select nonbonded representations;
- section writers consume the contexts without reinterpreting the topology.

New AMBER-family capabilities should add a topology analysis, an explicit
intermediate representation, a fail-closed compatibility rule, and an
AMBER/LAMMPS reference fixture. They should not add molecule-name checks or
special cases for a particular PDB.

The package split has begun behind the unchanged public CLI and Python API:
`amber_to_lammps_core.model` owns scaling model objects and
`amber_to_lammps_core.topology` owns topology-level scaling analysis. Later
steps will extract nonbonded, connectivity, writer, and CLI boundaries only
after equivalent characterization coverage; scientific behavior must not be
changed merely to complete the reorganization.
