
#!/usr/bin/env python3
"""
AMBER to LAMMPS conversion script

This script converts multiple AMBER topology files with molecule counts and a PDB coordinate file 
to LAMMPS data format. It requires parmed, numpy, and the AMBER force field files.

Usage:
    python amber_to_lammps.py <data_file> <param_file> <pdb_file> -t <top1.prmtop> [top2.prmtop ...] -c <count1> [count2 ...] --charges <q1> [q2 ...]

Arguments:
    data_file       Output LAMMPS data file name
    param_file      Output LAMMPS parameter file name
    pdb_file        PDB coordinate file (typically from packmol) containing all molecules
    -t / --topologies   One or more AMBER topology files (.prmtop)
    -c / --counts       Number of molecules for each topology file (same order/length as topologies)
    --charges           Target net charge per topology (same order/length as topologies)
"""

import argparse
import os
import sys
from dataclasses import dataclass
import parmed as pmd
import numpy as np


class UnsupportedNonbondedError(ValueError):
    """Raised when pair-specific parameters cannot be combined safely."""


class UnsupportedScalingError(ValueError):
    """Raised when AMBER 1-4 scaling cannot be represented exactly."""


@dataclass(frozen=True)
class SpecialBondScaling:
    """AMBER divisors and their equivalent LAMMPS 1-4 weights."""

    scee: float
    scnb: float

    @property
    def coulomb_weight(self):
        return 1.0 / self.scee

    @property
    def lj_weight(self):
        return 1.0 / self.scnb


def topology_14_scaling(parm, topology_name="topology"):
    """Return the unique active proper-dihedral SCEE/SCNB profile.

    AMBER stores zero scaling placeholders for some impropers and repeated
    torsion terms. Only proper dihedrals that own their end-to-end interaction
    (``ignore_end`` is false) determine the LAMMPS 1-4 weights.
    """
    profiles = set()
    for dihedral in parm.dihedrals:
        if dihedral.improper or dihedral.ignore_end:
            continue
        terms = dihedral.type if isinstance(dihedral.type, (list, tuple)) else [dihedral.type]
        for term in terms:
            if term is None:
                continue
            scee = float(term.scee)
            scnb = float(term.scnb)
            if scee <= 0.0 or scnb <= 0.0:
                raise UnsupportedScalingError(
                    f"{topology_name} has non-positive SCEE/SCNB values on an active "
                    f"proper dihedral: SCEE={scee:g}, SCNB={scnb:g}."
                )
            profiles.add((round(scee, 10), round(scnb, 10)))

    if not profiles:
        return None
    if len(profiles) != 1:
        details = ", ".join(f"SCEE={scee:g}/SCNB={scnb:g}" for scee, scnb in sorted(profiles))
        raise UnsupportedScalingError(
            f"{topology_name} uses multiple active 1-4 scaling profiles ({details}); "
            "one global LAMMPS special_bonds setting cannot represent them exactly."
        )
    scee, scnb = next(iter(profiles))
    return SpecialBondScaling(scee=scee, scnb=scnb)


def resolve_system_14_scaling(profiles, topologies):
    """Return a common profile or fail closed on mixed AMBER-family scaling."""
    used = {}
    for topology, profile in zip(topologies, profiles):
        if profile is not None:
            used.setdefault(profile, []).append(topology)
    if len(used) <= 1:
        return next(iter(used), None)

    details = "; ".join(
        f"SCEE={profile.scee:g}/SCNB={profile.scnb:g}: {', '.join(names)}"
        for profile, names in used.items()
    )
    raise UnsupportedScalingError(
        "The input topologies use different AMBER 1-4 scaling profiles. "
        f"Current LAMMPS output has one global special_bonds setting ({details}). "
        "Mixed GAFF2/GLYCAM-style scaling is detected but not yet representable exactly."
    )


def _amber_lj_pair(parm, type_i, type_j):
    """Return the raw (A, B) coefficients for two zero-based AMBER LJ types."""
    ntypes = parm.ptr("ntypes")
    table_index = parm.parm_data["NONBONDED_PARM_INDEX"][ntypes * type_i + type_j] - 1
    if table_index < 0:
        return None
    return (
        parm.parm_data["LENNARD_JONES_ACOEF"][table_index],
        parm.parm_data["LENNARD_JONES_BCOEF"][table_index],
    )


def _close_coefficient(actual, expected, rel_tol=1e-6, abs_tol=1e-12):
    return abs(actual - expected) <= max(abs_tol, rel_tol * max(abs(actual), abs(expected)))


def find_pair_specific_nonbonded_terms(parm):
    """Describe terms that require explicit I-J pair coefficients."""
    problems = []
    ccoef = parm.parm_data.get("LENNARD_JONES_CCOEF", [])
    nonzero_c4 = sum(not _close_coefficient(value, 0.0) for value in ccoef)
    if nonzero_c4:
        problems.append(f"{nonzero_c4} nonzero 12-6-4 C4 coefficient(s)")

    ntypes = parm.ptr("ntypes")
    diagonal = [_amber_lj_pair(parm, i, i) for i in range(ntypes)]
    modified_pairs = []
    for i in range(ntypes):
        for j in range(i + 1, ntypes):
            actual = _amber_lj_pair(parm, i, j)
            if actual is None or diagonal[i] is None or diagonal[j] is None:
                continue

            ai, bi = diagonal[i]
            aj, bj = diagonal[j]
            if ai <= 0.0 or bi <= 0.0 or aj <= 0.0 or bj <= 0.0:
                expected = (0.0, 0.0)
            else:
                rmin_i = (2.0 * ai / bi) ** (1.0 / 6.0)
                rmin_j = (2.0 * aj / bj) ** (1.0 / 6.0)
                epsilon_i = bi * bi / (4.0 * ai)
                epsilon_j = bj * bj / (4.0 * aj)
                rmin_ij = 0.5 * (rmin_i + rmin_j)
                epsilon_ij = (epsilon_i * epsilon_j) ** 0.5
                expected = (epsilon_ij * rmin_ij**12, 2.0 * epsilon_ij * rmin_ij**6)

            if not (
                _close_coefficient(actual[0], expected[0])
                and _close_coefficient(actual[1], expected[1])
            ):
                modified_pairs.append(f"{i + 1}-{j + 1}")

    if modified_pairs:
        preview = ", ".join(modified_pairs[:8])
        if len(modified_pairs) > 8:
            preview += f", and {len(modified_pairs) - 8} more"
        problems.append(
            f"{len(modified_pairs)} off-diagonal LJ pair override(s) "
            f"(AMBER type pairs {preview})"
        )
    return problems


def _parse_pdb_serial(serial_field, fallback=None):
    """Parse PDB atom serials, accepting Packmol base-36 overflow values."""
    token = serial_field.strip()
    if not token:
        if fallback is None:
            raise ValueError("Empty PDB atom serial field")
        return fallback

    try:
        return int(token)
    except ValueError:
        try:
            return int(token.upper(), 36)
        except ValueError:
            if fallback is None:
                raise
            return fallback

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Convert AMBER files to LAMMPS data format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Single molecule type
    python amber_to_lammps.py data.lammps parm.lammps combined.pdb -t epon.prmtop -c 1 --charges 0
    
    # Multiple molecule types
    python amber_to_lammps.py system.data system.parm combined.pdb -t mol1.prmtop mol2.prmtop -c 5 3 --charges 0 0 --verbose
    
    # With custom buffer
    python amber_to_lammps.py output.data output.parm all_molecules.pdb -t topo1.prmtop topo2.prmtop -c 10 20 --charges 0 0 -b 5.0
        """
    )
    
    parser.add_argument('data_file', help='Output LAMMPS data file name')
    parser.add_argument('param_file', help='Output LAMMPS parameter file name')
    parser.add_argument('pdb_file', help='PDB coordinate file (typically from packmol) containing all molecules')
    parser.add_argument('-t', '--topologies', nargs='+', required=True,
                        help='AMBER topology files (.prmtop) - specify one or more')
    parser.add_argument('-c', '--counts', type=int, nargs='+', required=True,
                        help='Number of molecules for each topology file (same order as --topologies)')
    parser.add_argument('-b', '--buffer', type=float, default=3.8,
                        help='Buffer size around molecule (default: 3.8)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Write debug files <data stem>.pairs.txt, .bonds.txt, .angles.txt, .dihedrals.txt beside the data file')
    parser.add_argument('--charges', type=float, nargs='+', required=True,
                        help='Target net charge per topology (same length/order as --topologies). Use 0 0 ... for neutral molecules.')
    parser.add_argument('--c4-cutoff', type=float, default=9.0,
                        help='Cutoff in Angstrom for generated 12-6-4 hybrid/overlay pair styles (default: 9.0)')
    parser.add_argument('--c4-coulomb', choices=('long', 'cut'), default='long',
                        help='Coulomb treatment for generated 12-6-4 pair style (default: long)')
    parser.add_argument('--typelabels', action='store_true',
                        help='Write out LAMMPS typelabels. Separate bond, angle, etc types will be created for separate types even if the parameters are the same.')
    
    return parser.parse_args()

def validate_files(topologies, molecule_counts, pdb_file):
    """Validate input files exist and are readable"""
    if len(topologies) != len(molecule_counts):
        print(f"Error: Number of topology files ({len(topologies)}) must match number of molecule counts ({len(molecule_counts)})")
        sys.exit(1)
    
    # Validate topology files
    for i, topology in enumerate(topologies):
        if not os.path.exists(topology):
            print(f"Error: Topology file '{topology}' (type {i+1}) not found")
            sys.exit(1)
        if not os.access(topology, os.R_OK):
            print(f"Error: Topology file '{topology}' (type {i+1}) is not readable")
            sys.exit(1)
        if molecule_counts[i] <= 0:
            print(f"Error: Molecule count for topology {i+1} must be positive (got {molecule_counts[i]})")
            sys.exit(1)
    
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found")
        sys.exit(1)
    if not os.access(pdb_file, os.R_OK):
        print(f"Error: PDB file '{pdb_file}' is not readable")
        sys.exit(1)

DEBUG_FILES = ('pairs.txt', 'bonds.txt', 'angles.txt', 'dihedrals.txt')


def debug_file_paths(data_file, param_file=None):
    """Return {debug name: path} for --keep-temp output, prefixed with the data file's stem.

    Writing ``<stem>.pairs.txt`` etc. beside the data file (rather than bare
    ``pairs.txt`` in the working directory) keeps --keep-temp from clobbering
    unrelated files that happen to share those common names.
    """
    stem = os.path.splitext(data_file)[0]
    paths = {name: f"{stem}.{name}" for name in DEBUG_FILES}
    outputs = {os.path.realpath(p) for p in (data_file, param_file) if p}
    for path in paths.values():
        if os.path.realpath(path) in outputs:
            raise ValueError(
                f"Debug file '{path}' would overwrite an output file; "
                "choose a different data or parameter file name."
            )
    return paths

# Charge corrections larger than this (in e per molecule) are reported even
# without --verbose. Charge-fitting round-off (e.g. AM1-BCC) is typically
# ~1e-3 e per molecule, so 0.01 e flags a genuinely requested change in net
# charge or a mistyped target without warning on routine neutralization.
CHARGE_WARNING_THRESHOLD = 0.01

def parse_pdb_coordinates(pdb_file, verbose=False):
    """Parse coordinates and atom information from PDB file"""
    if verbose:
        print(f"Parsing coordinates from PDB file: {pdb_file}")
    
    atoms = []
    x, y, z = [], [], []
    box_bounds = None
    
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('CRYST1'):
                try:
                    xhi = float(line[6:15].strip())
                    yhi = float(line[15:24].strip())
                    zhi = float(line[24:33].strip())
                    alpha = float(line[33:40].strip())
                    beta = float(line[40:47].strip())
                    gamma = float(line[47:54].strip())
                    if min(xhi, yhi, zhi) <= 0:
                        raise ValueError("cell lengths must be positive")
                    if not all(abs(angle - 90.0) <= 1.0e-3 for angle in (alpha, beta, gamma)):
                        raise ValueError(
                            "Non-orthogonal CRYST1 cells are not supported; "
                            f"got alpha={alpha:g}, beta={beta:g}, gamma={gamma:g}"
                        )
                    if all(abs(length - 1.0) <= 1.0e-3 for length in (xhi, yhi, zhi)):
                        print(
                            "Warning: ignoring placeholder 1 x 1 x 1 Angstrom "
                            "CRYST1 cell; using coordinate bounds and --buffer instead."
                        )
                        continue
                    box_bounds = BoxBounds(0,xhi,0,yhi,0,zhi)
                except ValueError as exc:
                    if "Non-orthogonal CRYST1" in str(exc):
                        raise
                    print(f"Warning: Could not parse PDB CRYST1 record: {line.strip()} ({exc})")
                    continue
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atom_num = _parse_pdb_serial(line[6:11], fallback=len(atoms) + 1)
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    x_coord = float(line[30:38].strip())
                    y_coord = float(line[38:46].strip())
                    z_coord = float(line[46:54].strip())

                    atoms.append({
                        'num': atom_num,
                        'name': atom_name,
                        'res_name': res_name,
                        'res_num': res_num,
                        'x': x_coord,
                        'y': y_coord,
                        'z': z_coord
                    })

                    x.append(x_coord)
                    y.append(y_coord)
                    z.append(z_coord)

                except (ValueError, IndexError):
                    if verbose:
                        print(f"Warning: Could not parse PDB line: {line.strip()}")
                    continue
    
    if verbose:
        print(f"Parsed {len(atoms)} atoms from PDB file")
    
    if len(atoms) == 0:
        raise ValueError(f"No atoms found in PDB file '{pdb_file}'. File may be empty or malformed.")

    if box_bounds is not None:
        cell_lengths = (box_bounds.xhi, box_bounds.yhi, box_bounds.zhi)
        coordinate_spans = (max(x) - min(x), max(y) - min(y), max(z) - min(z))
        axes = ("x", "y", "z")
        too_small = [
            f"{axis}: atom span {span:.6g} > cell length {length:.6g}"
            for axis, span, length in zip(axes, coordinate_spans, cell_lengths)
            if span > length + 1.0e-6
        ]
        if too_small:
            raise ValueError(
                "PDB CRYST1 cell is smaller than the coordinate span ("
                + "; ".join(too_small)
                + ")"
            )
    
    return atoms, x, y, z, box_bounds

@dataclass
class TopologyContext:
    all_parms: list
    atom_type_mapping: dict
    mass_list: list
    nonbonded_params: dict
    total_atoms_per_topology: list
    atom_types_per_topology: list
    type_remaps: list
    type_origins: dict
    raw_pair_params: dict
    has_pair_specific_terms: bool
    has_c4: bool
    scaling_profiles: list
    system_14_scaling: object


@dataclass
class BoxBounds:
    xlo: float
    xhi: float
    ylo: float
    yhi: float
    zlo: float
    zhi: float


@dataclass
class MoleculeSpan:
    mol_id: int
    topo_idx: int
    start: int
    end: int
    atoms_per_molecule: int


@dataclass
class ConnectivityContext:
    bond_count: int
    angle_count: int
    dihedral_count: int
    bond_type_count: int
    angle_type_count: int
    dihedral_type_count: int
    bond_lines: list
    angle_lines: list
    dihedral_lines: list
    bond_coeff_lines: list
    angle_coeff_lines: list
    dihedral_coeff_lines: list
    bond_type_label_lines: list
    angle_type_label_lines: list
    dihedral_type_label_lines: list
    bond_debug: list
    angle_debug: list
    dihedral_debug: list


def detect_multi_mode(topologies, molecule_counts, pdb_file, verbose):
    """Detect single vs multi-topology/copy workflows."""
    multi_mode = False

    if len(topologies) > 1:
        if verbose:
            print(f"✓ Auto-detected multi-topology system: {len(topologies)} topology files")
        multi_mode = True
    elif molecule_counts[0] > 1:
        if verbose:
            print(f"✓ Auto-detected multi-copy system: {molecule_counts[0]} copies of single topology")
        multi_mode = True
    elif len(topologies) == 1 and molecule_counts[0] == 1:
        try:
            with open(pdb_file, "r") as f:
                atom_count = sum(1 for line in f if line.startswith(("ATOM", "HETATM")))

            temp_parm = pmd.load_file(topologies[0])
            single_molecule_atoms = len(temp_parm.atoms)

            if atom_count > single_molecule_atoms:
                raise ValueError(
                    f"PDB file '{pdb_file}' contains {atom_count} atoms, but single topology expects {single_molecule_atoms} atoms per molecule. "
                    f"This appears to be a combined PDB file with multiple molecules. "
                    f"Please update your molecule count to {atom_count // single_molecule_atoms} to match the PDB content."
                )
        except Exception as e:
            if verbose:
                print(f"Warning: Could not analyze PDB file for combined detection: {e}")
                print("         Proceeding with single molecule mode")

    if not multi_mode and verbose:
        print("✓ Single molecule mode: 1 topology, 1 copy")

    return multi_mode


def initialize_output_files(data_file, param_file, topologies, pdb_file):
    with open(data_file, "w") as f:
        f.write("LAMMPS data file from AMBER conversion\n")
        f.write(f"#Source: {', '.join(topologies)}, {pdb_file}\n\n")

    with open(param_file, "w") as f:
        f.write("# Force field parameters generated from AMBER topologies\n\n")


def load_topology_context(topologies, verbose, amber_parm_cls, print_details):
    all_parms = []
    atom_type_mapping = {}
    mass_list = []
    nonbonded_params = {}
    total_atoms_per_topology = []
    atom_types_per_topology = []
    type_remaps = []
    type_origins = {}
    raw_pair_params = {}
    has_pair_specific_terms = False
    has_c4 = False
    param_tol = 1e-6
    variants_by_amber_type = {}
    scaling_profiles = []

    for i, topology in enumerate(topologies):
        if verbose:
            print(f"Loading topology {i+1}: {topology}")

        parm = amber_parm_cls(topology)
        scaling_profiles.append(topology_14_scaling(parm, topology))
        pair_specific = find_pair_specific_nonbonded_terms(parm)
        has_c4 = has_c4 or any("C4" in item for item in pair_specific)
        has_pair_specific_terms = has_pair_specific_terms or bool(pair_specific)
        all_parms.append(parm)
        total_atoms_per_topology.append(len(parm.atoms))
        atom_types_per_topology.append(set())
        type_remap = {}  # zero-based atom index -> canonical LAMMPS type name
        type_remaps.append(type_remap)

        if verbose:
            print(f"  Found {len(parm.atoms)} atoms, {len(parm.bonds)} bonds, {len(parm.angles)} angles, {len(parm.dihedrals)} dihedrals")

        lj_details = print_details(parm, "@1-{}".format(len(parm.atoms)))

        for line in str(lj_details).split("\n"):
            line = line.strip()
            if not line or not line[0].isdigit():
                continue

            parts = line.split()
            if len(parts) < 10:
                continue

            try:
                atom_index = int(parts[0]) - 1
                atom_type = parts[4]
                atom_mass = float(parts[8])
                lj_radius_amber = float(parts[6])
                lj_depth_amber = float(parts[7])
            except (ValueError, IndexError):
                continue
            lj_sigma = lj_radius_amber * (1 / (2 ** (1 / 6))) * 2

            # One LAMMPS type per distinct (AMBER type, LJ, mass) combination.
            # Atoms sharing an AMBER type name can differ in mass (hydrogen mass
            # repartitioning) or LJ parameters (conflicts between topologies).
            canonical_name = None
            for candidate in variants_by_amber_type.get(atom_type, []):
                existing = nonbonded_params[candidate]
                existing_mass = mass_list[atom_type_mapping[candidate] - 1]
                if (
                    abs(existing["lj_epsilon"] - lj_depth_amber) <= param_tol
                    and abs(existing["lj_sigma"] - lj_sigma) <= param_tol
                    and abs(existing_mass - atom_mass) <= param_tol
                ):
                    canonical_name = candidate
                    break

            if canonical_name is None:
                if atom_type not in atom_type_mapping:
                    canonical_name = atom_type
                else:
                    base_name = f"{atom_type}_top{i+1}"
                    canonical_name = base_name
                    suffix = 2
                    while canonical_name in atom_type_mapping:
                        canonical_name = f"{base_name}_{suffix}"
                        suffix += 1
                    print(
                        f"Atom type '{atom_type}' has differing LJ or mass parameters "
                        f"(e.g. hydrogen mass repartitioning or a conflict between topologies); "
                        f"using '{canonical_name}' for the affected atoms in topology {i+1}"
                    )
                variants_by_amber_type.setdefault(atom_type, []).append(canonical_name)
                atom_type_mapping[canonical_name] = len(atom_type_mapping) + 1
                mass_list.append(atom_mass)
                nonbonded_params[canonical_name] = {
                    "lj_epsilon": lj_depth_amber,
                    "lj_sigma": lj_sigma,
                }

            type_remap[atom_index] = canonical_name
            atom_types_per_topology[-1].add(canonical_name)
            type_origins.setdefault(canonical_name, set()).add(i + 1)

        missing = [idx + 1 for idx in range(len(parm.atoms)) if idx not in type_remap]
        if missing:
            preview = ", ".join(str(idx) for idx in missing[:8])
            if len(missing) > 8:
                preview += f", and {len(missing) - 8} more"
            raise ValueError(
                f"Could not read atom type, LJ, or mass parameters for {len(missing)} atom(s) "
                f"in topology {i+1} ({topology}): atom(s) {preview}"
            )

        nb_index_by_type = {}
        for atom in parm.atoms:
            nb_index_by_type.setdefault(type_remap[atom.idx], atom.nb_idx - 1)

        topology_types = sorted(atom_types_per_topology[-1])
        ccoef = parm.parm_data.get("LENNARD_JONES_CCOEF")
        for left_pos, left_type in enumerate(topology_types):
            for right_type in topology_types[left_pos:]:
                left_nb = nb_index_by_type[left_type]
                right_nb = nb_index_by_type[right_type]
                ab = _amber_lj_pair(parm, left_nb, right_nb)
                if ab is None:
                    continue
                ntypes = parm.ptr("ntypes")
                table_index = parm.parm_data["NONBONDED_PARM_INDEX"][ntypes * left_nb + right_nb] - 1
                c4 = ccoef[table_index] if ccoef is not None else 0.0
                key = tuple(sorted((left_type, right_type)))
                left_lj = nonbonded_params[left_type]
                right_lj = nonbonded_params[right_type]
                mixed_epsilon = (left_lj["lj_epsilon"] * right_lj["lj_epsilon"]) ** 0.5
                mixed_sigma = 0.5 * (left_lj["lj_sigma"] + right_lj["lj_sigma"])
                expected_a = 4.0 * mixed_epsilon * mixed_sigma**12
                expected_b = 4.0 * mixed_epsilon * mixed_sigma**6
                is_override = left_type != right_type and not (
                    _close_coefficient(ab[0], expected_a)
                    and _close_coefficient(ab[1], expected_b)
                )
                value = (float(ab[0]), float(ab[1]), float(c4), is_override)
                previous = raw_pair_params.get(key)
                if previous is not None and (
                    previous[3] != value[3]
                    or not all(_close_coefficient(old, new) for old, new in zip(previous[:3], value[:3]))
                ):
                    raise UnsupportedNonbondedError(
                        f"Pair parameters for atom types {key[0]}-{key[1]} conflict between topologies. "
                        "Rename the atom types or convert these topologies separately."
                    )
                raw_pair_params[key] = value

    system_14_scaling = resolve_system_14_scaling(scaling_profiles, topologies)

    if verbose:
        print(f"Found {len(atom_type_mapping)} unique atom types")
        if system_14_scaling is not None:
            print(
                "AMBER 1-4 scaling: "
                f"SCEE={system_14_scaling.scee:g}, SCNB={system_14_scaling.scnb:g} "
                f"(LAMMPS Coulomb={system_14_scaling.coulomb_weight:.10g}, "
                f"LJ={system_14_scaling.lj_weight:.10g})"
            )

    return TopologyContext(
        all_parms=all_parms,
        atom_type_mapping=atom_type_mapping,
        mass_list=mass_list,
        nonbonded_params=nonbonded_params,
        total_atoms_per_topology=total_atoms_per_topology,
        atom_types_per_topology=atom_types_per_topology,
        type_remaps=type_remaps,
        type_origins=type_origins,
        raw_pair_params=raw_pair_params,
        has_pair_specific_terms=has_pair_specific_terms,
        has_c4=has_c4,
        scaling_profiles=scaling_profiles,
        system_14_scaling=system_14_scaling,
    )


def validate_pdb_atom_count(pdb_atoms, molecule_counts, total_atoms_per_topology, multi_mode, topologies, pdb_file, verbose):
    contrib = [f"{count}*{atoms_per_topo}={count * atoms_per_topo}" for count, atoms_per_topo in zip(molecule_counts, total_atoms_per_topology)]
    breakdown = " + ".join(contrib)
    expected_total_atoms = sum(count * atoms_per_topo for count, atoms_per_topo in zip(molecule_counts, total_atoms_per_topology))

    if len(pdb_atoms) != expected_total_atoms:
        raise ValueError(
            f"Atom count mismatch: PDB has {len(pdb_atoms)} atoms but expected {expected_total_atoms} "
            f"({breakdown}). Check molecule counts, topology order, and PDB atom ordering."
        )
    if verbose:
        print(f"Atom count check passed: PDB={len(pdb_atoms)} matches expected {expected_total_atoms} ({breakdown})")

    if multi_mode:
        single_molecule_atoms = total_atoms_per_topology[0] if len(topologies) == 1 else sum(total_atoms_per_topology)
        if len(pdb_atoms) <= single_molecule_atoms:
            print(f"Warning: PDB file '{pdb_file}' contains {len(pdb_atoms)} atoms, but multi-mode expects")
            print(f"         a combined PDB file with {expected_total_atoms} atoms.")
            print("         Ensure you're using a combined PDB file from PackMol or similar tool.")
            print(f"         Expected: {expected_total_atoms} atoms ({breakdown})")
        else:
            print(f"✓ PDB file '{pdb_file}' contains {len(pdb_atoms)} atoms (appears to be combined)")

    return expected_total_atoms


def compute_box_bounds(x_coords, y_coords, z_coords, buffer):
    return BoxBounds(
        xlo=np.min(x_coords) - buffer,
        xhi=np.max(x_coords) + buffer,
        ylo=np.min(y_coords) - buffer,
        yhi=np.max(y_coords) + buffer,
        zlo=np.min(z_coords) - buffer,
        zhi=np.max(z_coords) + buffer,
    )


def write_data_header(data_file, expected_total_atoms, atom_type_mapping, conn_ctx, box_bounds):
    with open(data_file, "a") as f:
        f.write(f"{expected_total_atoms} atoms \n")
        f.write(f"{len(atom_type_mapping)} atom types \n")
        f.write(f"{conn_ctx.bond_count} bonds \n")
        f.write(f"{conn_ctx.bond_type_count} bond types \n")
        f.write(f"{conn_ctx.angle_count} angles \n")
        f.write(f"{conn_ctx.angle_type_count} angle types \n")
        f.write(f"{conn_ctx.dihedral_count} dihedrals \n")
        f.write(f"{conn_ctx.dihedral_type_count} dihedral types \n\n")
        f.write(f"{box_bounds.xlo} {box_bounds.xhi} xlo xhi \n")
        f.write(f"{box_bounds.ylo} {box_bounds.yhi} ylo yhi \n")
        f.write(f"{box_bounds.zlo} {box_bounds.zhi} zlo zhi \n")


def build_charges(all_parms, molecule_counts, charges_target, topologies, verbose):
    if len(charges_target) != len(topologies):
        raise ValueError(f"Number of charges provided ({len(charges_target)}) must match number of topologies ({len(topologies)})")

    charges = []
    charge_tol = 1e-6

    for topo_idx, (parm, count, target_charge_per_mol) in enumerate(zip(all_parms, molecule_counts, charges_target)):
        topo_charges = np.array([atom.charge for atom in parm.atoms], dtype=float)
        actual_charge_per_mol = float(np.sum(topo_charges))
        diff = target_charge_per_mol - actual_charge_per_mol

        if abs(diff) > charge_tol:
            shift = diff / len(topo_charges)
            topo_charges = topo_charges + shift
            if abs(diff) > CHARGE_WARNING_THRESHOLD:
                print(
                    f"Warning: topology {topo_idx+1} ({topologies[topo_idx]}) has net charge "
                    f"{actual_charge_per_mol:.6f} e but --charges requests {target_charge_per_mol:.6f} e; "
                    f"shifting every atom by {shift:+.6f} e ({diff:+.6f} e per molecule). "
                    f"Check --charges if this change was not intended."
                )
            elif verbose:
                print(
                    f"Charge adjust topo {topo_idx+1}: actual {actual_charge_per_mol:.6f} -> "
                    f"target {target_charge_per_mol:.6f} (shift {shift:.6f}/atom)"
                )
        elif verbose:
            print(f"Charge check topo {topo_idx+1}: {actual_charge_per_mol:.6f} matches target {target_charge_per_mol:.6f}")

        for _ in range(count):
            charges.extend(topo_charges.tolist())

    net_charge = float(np.sum(charges))
    target_total_charge = float(np.sum([c * q for c, q in zip(molecule_counts, charges_target)]))
    if abs(net_charge - target_total_charge) > 1e-4:
        raise ValueError(f"Total charge mismatch after adjustment: got {net_charge:.6f}, expected {target_total_charge:.6f}")
    if verbose:
        print(f"Total charge check passed: {net_charge:.6f} matches expected {target_total_charge:.6f}")

    return charges


def build_molecule_spans(all_parms, molecule_counts):
    molecule_spans = []
    running_offset = 0
    for topo_idx, (parm, count) in enumerate(zip(all_parms, molecule_counts)):
        atoms_per_molecule = len(parm.atoms)
        for _ in range(count):
            start = running_offset
            end = start + atoms_per_molecule - 1
            molecule_spans.append(
                MoleculeSpan(
                    mol_id=len(molecule_spans) + 1,
                    topo_idx=topo_idx,
                    start=start,
                    end=end,
                    atoms_per_molecule=atoms_per_molecule,
                )
            )
            running_offset += atoms_per_molecule
    return molecule_spans


def build_pair_coeff_map(atom_type_mapping, nonbonded_params):
    pair_coeff_map = {}
    for atom_type in atom_type_mapping.keys():
        if atom_type in nonbonded_params:
            lj_epsilon = nonbonded_params[atom_type]["lj_epsilon"]
            lj_sigma = nonbonded_params[atom_type]["lj_sigma"]
            type_id = atom_type_mapping[atom_type]
            pair_coeff_map[atom_type] = f"pair_coeff {type_id} {type_id} {lj_epsilon} {lj_sigma} # {atom_type}"
    return pair_coeff_map


def write_masses_section(data_file, top_ctx, verbose):
    if verbose:
        print("Writing masses section...")

    with open(data_file, "a") as f:
        f.write("\nMasses \n\n")
        for i in range(len(top_ctx.atom_type_mapping)):
            f.write(f"{i+1} {top_ctx.mass_list[i]} \n")

def _sanitize_type_label(label):
    """Make an AMBER type name valid as a LAMMPS type label."""
    sanitized = "".join("_" if char.isspace() else char for char in str(label))
    if not sanitized:
        sanitized = "type"
    if sanitized[0].isdigit() or sanitized[0] in "*#":
        sanitized = f"type_{sanitized}"
    return sanitized


def build_atom_type_labels(atom_type_mapping):
    """Map canonical atom types to stable, valid, unique LAMMPS labels."""
    labels = {}
    used = set()
    for atom_type in atom_type_mapping:
        labels[atom_type] = _unique_type_label(_sanitize_type_label(atom_type), used)
    return labels


def write_atom_type_labels_section(data_file, atom_type_mapping, atom_type_labels, verbose):
    if verbose:
        print("Writing atom labels section...")

    with open(data_file, "a") as flammps:
        flammps.write("\nAtom Type Labels\n\n")

        for atom_type in atom_type_mapping.keys():
            type_id = atom_type_mapping[atom_type]
            flammps.write(f"{type_id} {atom_type_labels[atom_type]}\n")

def write_connectivity_labels_section(data_file, conn_ctx, verbose):
    if verbose:
        print("Writing connectivity type labels section...")

    with open(data_file, "a") as flammps:
        flammps.write("\nBond Type Labels\n\n")
        flammps.write("\n".join(conn_ctx.bond_type_label_lines) + "\n")
        flammps.write("\nAngle Type Labels\n\n")
        flammps.write("\n".join(conn_ctx.angle_type_label_lines) + "\n")
        flammps.write("\nDihedral Type Labels\n\n")
        flammps.write("\n".join(conn_ctx.dihedral_type_label_lines) + "\n")

def _ab_to_epsilon_sigma(acoef, bcoef):
    if acoef <= 0.0 or bcoef <= 0.0:
        return 0.0, 0.0
    return bcoef * bcoef / (4.0 * acoef), (acoef / bcoef) ** (1.0 / 6.0)


def build_explicit_pair_coeff_lines(top_ctx, c4_cutoff, c4_coulomb):
    """Build a complete I-J matrix, preserving AMBER A, B, and C4 coefficients."""
    if c4_cutoff <= 0.0:
        raise ValueError(f"--c4-cutoff must be positive (got {c4_cutoff})")

    atom_types = sorted(top_ctx.atom_type_mapping, key=top_ctx.atom_type_mapping.get)
    lj_style = f"lj/cut/coul/{c4_coulomb}"
    lines = []
    if top_ctx.has_c4:
        lines.extend([
            "# 12-6-4 topology: this include configures the required LEPTON overlay.",
            "# LAMMPS must be built with the LEPTON package.",
            f"pair_style hybrid/overlay {lj_style} {c4_cutoff:.10g} {c4_cutoff:.10g} lepton {c4_cutoff:.10g}",
            "pair_modify tail no",
            "# AMBER's C4 CustomNonbondedForce excludes 1-2, 1-3, and 1-4 pairs.",
            "pair_modify pair lepton special lj 0.0 0.0 0.0",
        ])

    for left_pos, left_type in enumerate(atom_types):
        for right_type in atom_types[left_pos:]:
            left_id = top_ctx.atom_type_mapping[left_type]
            right_id = top_ctx.atom_type_mapping[right_type]
            raw = top_ctx.raw_pair_params.get(tuple(sorted((left_type, right_type))))
            left = top_ctx.nonbonded_params[left_type]
            right = top_ctx.nonbonded_params[right_type]
            epsilon = (left["lj_epsilon"] * right["lj_epsilon"]) ** 0.5
            sigma = 0.5 * (left["lj_sigma"] + right["lj_sigma"])
            c4 = raw[2] if raw is not None else 0.0
            if raw is not None and raw[3]:
                epsilon, sigma = _ab_to_epsilon_sigma(raw[0], raw[1])

            style_field = f" {lj_style}" if top_ctx.has_c4 else ""
            lines.append(
                f"pair_coeff {left_id} {right_id}{style_field} "
                f"{epsilon:.16g} {sigma:.16g} # {left_type}-{right_type}"
            )
            if top_ctx.has_c4 and not _close_coefficient(c4, 0.0):
                lines.append(
                    f'pair_coeff {left_id} {right_id} lepton "-{c4:.16g}/r^4" '
                    f"# C4 {left_type}-{right_type}"
                )
    return lines


def write_atoms_section(data_file, pdb_atoms, molecule_spans, all_parms, type_remaps, atom_type_mapping, charges, verbose):
    if verbose:
        print("Writing atoms section...")

    with open(data_file, "a") as flammps:
        flammps.write("Atoms\n\n")
        span_idx = 0
        for i, pdb_atom in enumerate(pdb_atoms):
            while span_idx < len(molecule_spans) and i > molecule_spans[span_idx].end:
                span_idx += 1

            if span_idx >= len(molecule_spans):
                raise ValueError(f"Atom index {i} exceeds computed molecule spans; check PDB ordering.")

            span = molecule_spans[span_idx]
            atom_idx_in_mol = i - span.start
            canonical_atom_type = type_remaps[span.topo_idx].get(atom_idx_in_mol)
            atom_type_id = atom_type_mapping.get(canonical_atom_type)
            if atom_type_id is None:
                amber_type = all_parms[span.topo_idx].atoms[atom_idx_in_mol].type
                raise ValueError(
                    f"No LAMMPS atom type for atom {i+1} (topology {span.topo_idx+1}, "
                    f"atom {atom_idx_in_mol+1}, AMBER type '{amber_type}')"
                )

            flammps.write(
                f"{i+1} {span.mol_id} {atom_type_id} {charges[i]:.10f} "
                f"{pdb_atom['x']:.4f} {pdb_atom['y']:.4f} {pdb_atom['z']:.4f} 0 0 0\n"
            )


def _unique_type_label(base, used_labels):
    """Return a stable unique LAMMPS label without merging parameter conflicts."""
    label = base
    suffix = 2
    while label in used_labels:
        label = f"{base}_{suffix}"
        suffix += 1
    used_labels.add(label)
    return label


def process_connectivity(all_parms, molecule_counts, topologies, use_typelabels,
                         type_remaps=None, atom_type_labels=None):
    bond_type_registry = {}
    angle_type_registry = {}
    dihedral_type_registry = {}
    bond_coeff_lines = []
    angle_coeff_lines = []
    dihedral_coeff_lines = []
    bond_type_label_lines = []
    angle_type_label_lines = []
    dihedral_type_label_lines = []
    used_bond_labels = set()
    used_angle_labels = set()
    used_dihedral_labels = set()
    bond_type_ids_per_topo = []
    angle_type_ids_per_topo = []
    dih_entries_per_topo = []

    def dihedral_signature(terms):
        result = []
        for term in terms:
            phi_k = float(term.phi_k)
            per = int(round(float(term.per)))
            phase = float(term.phase)
            if abs(phase) <= 2 * np.pi + 0.1:
                phase = np.degrees(phase)
            result.append((round(phi_k, 4), per, round(phase, 4)))
        return tuple(sorted(result))

    def labelled_atom_type(topo_idx, atom):
        if type_remaps is None:
            return _sanitize_type_label(atom.type)
        canonical = type_remaps[topo_idx][atom.idx]
        return atom_type_labels[canonical]

    for topo_idx, parm in enumerate(all_parms):
        bond_type_ids = []
        for bond in parm.bonds:
            if bond.type is None:
                raise ValueError(f"Bond parameters missing for atoms {bond.atom1.idx}-{bond.atom2.idx} in topology {topo_idx+1}")
            k = float(bond.type.k)
            req = float(bond.type.req)
            display_key = (bond.atom1.type, bond.atom2.type)
            type_key = tuple(labelled_atom_type(topo_idx, atom) for atom in (bond.atom1, bond.atom2)) if use_typelabels else display_key
            if type_key[0] > type_key[1]:
                type_key = tuple(reversed(type_key))
            param_sig = (round(k, 6), round(req, 6))
            sig = (type_key, param_sig) if use_typelabels else param_sig
            if sig not in bond_type_registry:
                type_id = len(bond_type_registry) + 1
                bond_type_registry[sig] = type_id
                comment_key = type_key if use_typelabels else display_key
                bond_coeff_lines.append(
                    f"bond_coeff {type_id} {k:.4f} {req:.4f}  # {'-'.join(comment_key)}"
                )
                label = _unique_type_label("-".join(type_key), used_bond_labels)
                bond_type_label_lines.append(f"{type_id} {label}")
            bond_type_ids.append(bond_type_registry[sig])
        bond_type_ids_per_topo.append(bond_type_ids)

        angle_type_ids = []
        for angle in parm.angles:
            if angle.type is None:
                raise ValueError(
                    f"Angle parameters missing for atoms {angle.atom1.idx}-{angle.atom2.idx}-{angle.atom3.idx} in topology {topo_idx+1}"
                )
            k = float(angle.type.k)
            theteq = float(angle.type.theteq)
            display_key = (angle.atom1.type, angle.atom2.type, angle.atom3.type)
            type_key = tuple(labelled_atom_type(topo_idx, atom) for atom in (angle.atom1, angle.atom2, angle.atom3)) if use_typelabels else display_key
            if type_key[0] > type_key[2]:
                type_key = tuple(reversed(type_key))
            param_sig = (round(k, 6), round(theteq, 6))
            sig = (type_key, param_sig) if use_typelabels else param_sig
            if sig not in angle_type_registry:
                type_id = len(angle_type_registry) + 1
                angle_type_registry[sig] = type_id
                comment_key = type_key if use_typelabels else display_key
                angle_coeff_lines.append(
                    f"angle_coeff {type_id} {k:.4f} {theteq:.4f}  # {'-'.join(comment_key)}"
                )
                label = _unique_type_label("-".join(type_key), used_angle_labels)
                angle_type_label_lines.append(f"{type_id} {label}")
            angle_type_ids.append(angle_type_registry[sig])
        angle_type_ids_per_topo.append(angle_type_ids)

        dih_by_idx = {}
        dih_first = {}
        for dih in parm.dihedrals:
            if dih.type is None:
                raise ValueError(
                    f"Dihedral parameters missing for atoms {dih.atom1.idx}-{dih.atom2.idx}-{dih.atom3.idx}-{dih.atom4.idx} in topology {topo_idx+1}"
                )
            idx_key = (dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx)
            if isinstance(dih.type, (list, tuple)):
                terms = [term for term in dih.type if term is not None]
            else:
                terms = [dih.type]
            if idx_key not in dih_by_idx:
                dih_by_idx[idx_key] = []
                dih_first[idx_key] = dih
            dih_by_idx[idx_key].extend(terms)

        dih_entries = []
        for idx_key, terms in dih_by_idx.items():
            dih = dih_first[idx_key]
            sig_terms = dihedral_signature(terms)
            display_key = (dih.atom1.type, dih.atom2.type, dih.atom3.type, dih.atom4.type)
            type_key = tuple(
                labelled_atom_type(topo_idx, atom)
                for atom in (dih.atom1, dih.atom2, dih.atom3, dih.atom4)
            ) if use_typelabels else display_key
            if type_key[0] > type_key[3] or (type_key[0] == type_key[3] and type_key[1] > type_key[2]):
                type_key = tuple(reversed(type_key))
            sig = (type_key, sig_terms) if use_typelabels else sig_terms
            if sig not in dihedral_type_registry:
                type_id = len(dihedral_type_registry) + 1
                dihedral_type_registry[sig] = type_id
                coeff_str = " ".join(f"{phi_k:.4f} {per} {phase:.4f}" for phi_k, per, phase in sig_terms)
                comment_key = type_key if use_typelabels else display_key
                dihedral_coeff_lines.append(
                    f"dihedral_coeff {type_id} {len(sig_terms)} {coeff_str}  # {'-'.join(comment_key)}"
                )
                label = _unique_type_label("-".join(type_key), used_dihedral_labels)
                dihedral_type_label_lines.append(f"{type_id} {label}")
            dih_entries.append((dihedral_type_registry[sig], dih))
        dih_entries_per_topo.append(dih_entries)

    bond_lines = []
    angle_lines = []
    dihedral_lines = []
    bond_debug = []
    angle_debug = []
    dihedral_debug = []
    bond_count = 0
    angle_count = 0
    dihedral_count = 0
    atom_offset = 0

    for topo_idx, (parm, count) in enumerate(zip(all_parms, molecule_counts)):
        atoms_per_molecule = len(parm.atoms)
        bond_type_ids = bond_type_ids_per_topo[topo_idx]
        angle_type_ids = angle_type_ids_per_topo[topo_idx]
        dih_entries = dih_entries_per_topo[topo_idx]

        for mol_idx in range(count):
            base_offset = atom_offset + mol_idx * atoms_per_molecule

            for bond_idx, bond in enumerate(parm.bonds):
                bond_count += 1
                type_id = bond_type_ids[bond_idx]
                atom1 = bond.atom1.idx + 1 + base_offset
                atom2 = bond.atom2.idx + 1 + base_offset
                bond_lines.append(f"{bond_count} {type_id} {atom1} {atom2}")
                if mol_idx == 0:
                    bond_debug.append(
                        f"{bond_count}\t{type_id}\t{atom1}\t{atom2}\t{bond.type.k:.6f}\t{bond.type.req:.6f}\t{topo_idx+1}\t{mol_idx+1}"
                    )

            for angle_idx, angle in enumerate(parm.angles):
                angle_count += 1
                type_id = angle_type_ids[angle_idx]
                atom1 = angle.atom1.idx + 1 + base_offset
                atom2 = angle.atom2.idx + 1 + base_offset
                atom3 = angle.atom3.idx + 1 + base_offset
                angle_lines.append(f"{angle_count} {type_id} {atom1} {atom2} {atom3}")
                if mol_idx == 0:
                    angle_debug.append(
                        f"{angle_count}\t{type_id}\t{atom1}\t{atom2}\t{atom3}\t{angle.type.k:.6f}\t{angle.type.theteq:.6f}\t{topo_idx+1}\t{mol_idx+1}"
                    )

            for type_id, dih in dih_entries:
                dihedral_count += 1
                atom1 = dih.atom1.idx + 1 + base_offset
                atom2 = dih.atom2.idx + 1 + base_offset
                atom3 = dih.atom3.idx + 1 + base_offset
                atom4 = dih.atom4.idx + 1 + base_offset

                dihedral_lines.append(f"{dihedral_count} {type_id} {atom1} {atom2} {atom3} {atom4}")
                if mol_idx == 0:
                    dihedral_debug.append(
                        f"{dihedral_count}\t{type_id}\t{atom1}\t{atom2}\t{atom3}\t{atom4}\t{topo_idx+1}\t{mol_idx+1}"
                    )

        atom_offset += count * atoms_per_molecule

    return ConnectivityContext(
        bond_count=bond_count,
        angle_count=angle_count,
        dihedral_count=dihedral_count,
        bond_type_count=len(bond_type_registry),
        angle_type_count=len(angle_type_registry),
        dihedral_type_count=len(dihedral_type_registry),
        bond_lines=bond_lines,
        angle_lines=angle_lines,
        dihedral_lines=dihedral_lines,
        bond_coeff_lines=bond_coeff_lines,
        angle_coeff_lines=angle_coeff_lines,
        dihedral_coeff_lines=dihedral_coeff_lines,
        bond_type_label_lines=bond_type_label_lines,
        angle_type_label_lines=angle_type_label_lines,
        dihedral_type_label_lines=dihedral_type_label_lines,
        bond_debug=bond_debug,
        angle_debug=angle_debug,
        dihedral_debug=dihedral_debug,
    )


def write_connectivity_sections(data_file, conn_ctx, verbose):
    if verbose:
        print("Writing bonds section...")
    with open(data_file, "a") as flammps:
        flammps.write("\nBonds \n\n")
        flammps.write("\n".join(conn_ctx.bond_lines) + "\n")

    if verbose:
        print("Writing angles section...")
    with open(data_file, "a") as flammps:
        flammps.write("\nAngles \n\n")
        flammps.write("\n".join(conn_ctx.angle_lines) + "\n")

    if verbose:
        print("Writing dihedrals section...")
    with open(data_file, "a") as flammps:
        flammps.write("\nDihedrals \n\n")
        flammps.write("\n".join(conn_ctx.dihedral_lines) + "\n")



def write_debug_files(keep_temp, topologies, atom_types_per_topology, atom_type_mapping, nonbonded_params, type_origins, conn_ctx,
                      data_file, param_file=None):
    if not keep_temp:
        return

    paths = debug_file_paths(data_file, param_file)

    with open(paths["pairs.txt"], "w") as ftemp:
        ftemp.write(f"# {os.path.basename(paths['pairs.txt'])} generated by amber_to_lammps.py\n")
        ftemp.write("# columns: type_id\tname\tepsilon\tsigma\ttopology_sources\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            for atom_type in sorted(atom_types_per_topology[topo_idx]):
                type_id = atom_type_mapping[atom_type]
                params = nonbonded_params[atom_type]
                origins = ",".join(str(t) for t in sorted(type_origins.get(atom_type, {topo_idx+1})))
                ftemp.write(f"{type_id}\t{atom_type}\t{params['lj_epsilon']:.6f}\t{params['lj_sigma']:.6f}\t{origins}\n")
            ftemp.write("\n")

    with open(paths["bonds.txt"], "w") as ftemp:
        ftemp.write(f"# {os.path.basename(paths['bonds.txt'])} generated by amber_to_lammps.py\n")
        ftemp.write("# columns: instance_id\ttype_id\tatom1\tatom2\tk\treq\ttopology_idx\tmolecule_idx\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            topo_lines = [line for line in conn_ctx.bond_debug if line.split("\t")[-2] == str(topo_idx + 1)]
            if topo_lines:
                ftemp.write("\n".join(topo_lines) + "\n")
            else:
                ftemp.write("# (none)\n")
            ftemp.write("\n")

    with open(paths["angles.txt"], "w") as ftemp:
        ftemp.write(f"# {os.path.basename(paths['angles.txt'])} generated by amber_to_lammps.py\n")
        ftemp.write("# columns: instance_id\ttype_id\ta1\ta2\ta3\tk\ttheta\ttopology_idx\tmolecule_idx\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            topo_lines = [line for line in conn_ctx.angle_debug if line.split("\t")[-2] == str(topo_idx + 1)]
            if topo_lines:
                ftemp.write("\n".join(topo_lines) + "\n")
            else:
                ftemp.write("# (none)\n")
            ftemp.write("\n")

    with open(paths["dihedrals.txt"], "w") as ftemp:
        ftemp.write(f"# {os.path.basename(paths['dihedrals.txt'])} generated by amber_to_lammps.py\n")
        ftemp.write("# columns: instance_id\ttype_id\ta1\ta2\ta3\ta4\ttopology_idx\tmolecule_idx\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            topo_lines = [line for line in conn_ctx.dihedral_debug if line.split("\t")[-2] == str(topo_idx + 1)]
            if topo_lines:
                ftemp.write("\n".join(topo_lines) + "\n")
            else:
                ftemp.write("# (none)\n")
            ftemp.write("\n")

    print(f"Wrote debug files (--keep-temp): {', '.join(paths.values())}")


def write_parameter_file(param_file, topologies, atom_types_per_topology, pair_coeff_map, conn_ctx,
                         top_ctx, c4_cutoff, c4_coulomb):
    printed_pair_types = set()
    with open(param_file, "a") as flammpsparm:
        if top_ctx.system_14_scaling is not None:
            scaling = top_ctx.system_14_scaling
            flammpsparm.write("# AMBER 1-4 nonbonded scaling read from the topology\n")
            flammpsparm.write(
                "special_bonds "
                f"lj 0.0 0.0 {scaling.lj_weight:.10g} "
                f"coul 0.0 0.0 {scaling.coulomb_weight:.10g}\n"
            )
        if top_ctx.has_pair_specific_terms:
            flammpsparm.write("\n# Explicit AMBER nonbonded I-J coefficient matrix\n")
            flammpsparm.write("\n".join(build_explicit_pair_coeff_lines(top_ctx, c4_cutoff, c4_coulomb)) + "\n")
        else:
            for topo_idx, topo_name in enumerate(topologies):
                flammpsparm.write(f"\n# Topology {topo_idx+1}: {topo_name}\n")
                flammpsparm.write("# Nonbonded\n")
                for atom_type in sorted(atom_types_per_topology[topo_idx]):
                    line = pair_coeff_map.get(atom_type)
                    if line is None:
                        continue
                    if atom_type in printed_pair_types:
                        flammpsparm.write(f"# pair_coeff for {atom_type} already defined above\n")
                    else:
                        flammpsparm.write(line + "\n")
                        printed_pair_types.add(atom_type)

        flammpsparm.write("\n# Bond Coefficients (deduplicated across all topologies and copies)\n")
        flammpsparm.write("\n".join(conn_ctx.bond_coeff_lines) + "\n")

        flammpsparm.write("\n# Angle Coefficients (deduplicated across all topologies and copies)\n")
        flammpsparm.write("\n".join(conn_ctx.angle_coeff_lines) + "\n")

        flammpsparm.write("\n# Dihedral Coefficients (deduplicated across all topologies and copies)\n")
        flammpsparm.write("\n".join(conn_ctx.dihedral_coeff_lines) + "\n")


def amber2lammps(data_file, param_file, topologies, molecule_counts, pdb_file, charges_target,
                 buffer=3.8, verbose=False, keep_temp=False, c4_cutoff=9.0, c4_coulomb="long", use_typelabels=False):
    amber_parm_cls = pmd.amber.AmberParm
    print_details = pmd.tools.actions.printDetails

    if keep_temp:
        debug_file_paths(data_file, param_file)  # fail before writing anything on a name collision

    multi_mode = detect_multi_mode(topologies, molecule_counts, pdb_file, verbose)

    if verbose:
        print("Converting multiple AMBER topologies to LAMMPS format...")
        print(f"Output files: {data_file}, {param_file}")
        for i, (topo, count) in enumerate(zip(topologies, molecule_counts)):
            print(f"  Topology {i+1}: {topo} ({count} molecules)")

    pdb_atoms, x_coords, y_coords, z_coords, box_bounds = parse_pdb_coordinates(pdb_file, verbose)
    top_ctx = load_topology_context(topologies, verbose, amber_parm_cls, print_details)
    initialize_output_files(data_file, param_file, topologies, pdb_file)

    expected_total_atoms = validate_pdb_atom_count(
        pdb_atoms=pdb_atoms,
        molecule_counts=molecule_counts,
        total_atoms_per_topology=top_ctx.total_atoms_per_topology,
        multi_mode=multi_mode,
        topologies=topologies,
        pdb_file=pdb_file,
        verbose=verbose,
    )

    if not box_bounds:
        box_bounds = compute_box_bounds(x_coords, y_coords, z_coords, buffer)
    if verbose:
        print(
            f"Box dimensions: X[{box_bounds.xlo:.3f}, {box_bounds.xhi:.3f}], "
            f"Y[{box_bounds.ylo:.3f}, {box_bounds.yhi:.3f}], "
            f"Z[{box_bounds.zlo:.3f}, {box_bounds.zhi:.3f}]"
        )

    atom_type_labels = build_atom_type_labels(top_ctx.atom_type_mapping) if use_typelabels else None
    conn_ctx = process_connectivity(
        top_ctx.all_parms,
        molecule_counts,
        topologies,
        use_typelabels,
        type_remaps=top_ctx.type_remaps,
        atom_type_labels=atom_type_labels,
    )

    write_data_header(
        data_file=data_file,
        expected_total_atoms=expected_total_atoms,
        atom_type_mapping=top_ctx.atom_type_mapping,
        conn_ctx=conn_ctx,
        box_bounds=box_bounds,
    )

    if use_typelabels:
        write_atom_type_labels_section(
            data_file=data_file,
            atom_type_mapping=top_ctx.atom_type_mapping,
            atom_type_labels=atom_type_labels,
            verbose=verbose,
        )

        write_connectivity_labels_section(
            data_file=data_file,
            conn_ctx=conn_ctx,
            verbose=verbose,
        )

    write_masses_section(
        data_file=data_file,
        top_ctx=top_ctx,
        verbose=verbose,
    )

    charges = build_charges(top_ctx.all_parms, molecule_counts, charges_target, topologies, verbose)
    molecule_spans = build_molecule_spans(top_ctx.all_parms, molecule_counts)
    pair_coeff_map = build_pair_coeff_map(top_ctx.atom_type_mapping, top_ctx.nonbonded_params)

    write_atoms_section(
        data_file=data_file,
        pdb_atoms=pdb_atoms,
        molecule_spans=molecule_spans,
        all_parms=top_ctx.all_parms,
        type_remaps=top_ctx.type_remaps,
        atom_type_mapping=top_ctx.atom_type_mapping,
        charges=charges,
        verbose=verbose,
    )

    write_connectivity_sections(data_file, conn_ctx, verbose)

    write_debug_files(
        keep_temp=keep_temp,
        topologies=topologies,
        atom_types_per_topology=top_ctx.atom_types_per_topology,
        atom_type_mapping=top_ctx.atom_type_mapping,
        nonbonded_params=top_ctx.nonbonded_params,
        type_origins=top_ctx.type_origins,
        conn_ctx=conn_ctx,
        data_file=data_file,
        param_file=param_file,
    )

    write_parameter_file(
        param_file=param_file,
        topologies=topologies,
        atom_types_per_topology=top_ctx.atom_types_per_topology,
        pair_coeff_map=pair_coeff_map,
        conn_ctx=conn_ctx,
        top_ctx=top_ctx,
        c4_cutoff=c4_cutoff,
        c4_coulomb=c4_coulomb,
    )

    if verbose:
        print("Conversion complete!")
        print("Generated files:")
        print(f"  - {data_file} (LAMMPS data file)")
        print(f"  - {param_file} (LAMMPS parameters)")
        print("Summary:")
        print(f"  - {expected_total_atoms} atoms")
        print(f"  - {conn_ctx.bond_count} bonds ({conn_ctx.bond_type_count} unique types)")
        print(f"  - {conn_ctx.angle_count} angles ({conn_ctx.angle_type_count} unique types)")
        print(f"  - {conn_ctx.dihedral_count} dihedrals ({conn_ctx.dihedral_type_count} unique types)")

def main():
    """Main function"""
    args = parse_arguments()
    
    # Validate input files and get multi-mode status
    validate_files(args.topologies, args.counts, args.pdb_file)
    
    # Run conversion
    try:
        amber2lammps(
            data_file=args.data_file,
            param_file=args.param_file,
            topologies=args.topologies,
            molecule_counts=args.counts,
            pdb_file=args.pdb_file,
            charges_target=args.charges,
            buffer=args.buffer,
            verbose=args.verbose,
            keep_temp=args.keep_temp,
            c4_cutoff=args.c4_cutoff,
            c4_coulomb=args.c4_coulomb,
            use_typelabels=args.typelabels,
        )
        print(f"✓ Conversion completed successfully!")
        print(f"Output files: {args.data_file}, {args.param_file}")
        return True
        
    except Exception as e:
        print(f"✗ Error during conversion: {e}")
        return False


def entry_point():
    """Run the installed console command with conventional process exit codes."""
    return 0 if main() else 1


if __name__ == "__main__":
    sys.exit(entry_point())

    
    
