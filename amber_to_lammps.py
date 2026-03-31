
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
                        help='Keep temporary files (pairs.txt, bonds.txt, angles.txt, dihedrals.txt) after conversion')
    parser.add_argument('--charges', type=float, nargs='+', required=True,
                        help='Target net charge per topology (same length/order as --topologies). Use 0 0 ... for neutral molecules.')
    
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

def cleanup_temp_files(verbose=False, keep_temp=False):
    """Remove temporary files if they exist"""
    if keep_temp:
        if verbose:
            print("Keeping temporary files as requested:")
            temp_files = ['bonds.txt', 'angles.txt', 'dihedrals.txt', 'pairs.txt']
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    print(f"  - {temp_file}")
        return
        
    temp_files = ['bonds.txt', 'angles.txt', 'dihedrals.txt', 'pairs.txt']
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
                if verbose:
                    print(f"Removed existing temporary file: {temp_file}")
            except OSError as e:
                print(f"Warning: Could not remove {temp_file}: {e}")

def parse_pdb_coordinates(pdb_file, verbose=False):
    """Parse coordinates and atom information from PDB file"""
    if verbose:
        print(f"Parsing coordinates from PDB file: {pdb_file}")
    
    atoms = []
    x, y, z = [], [], []
    
    with open(pdb_file) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
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
    
    return atoms, x, y, z

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
    param_tol = 1e-6

    for i, topology in enumerate(topologies):
        if verbose:
            print(f"Loading topology {i+1}: {topology}")

        parm = amber_parm_cls(topology)
        all_parms.append(parm)
        total_atoms_per_topology.append(len(parm.atoms))
        atom_types_per_topology.append(set())
        type_remap = {}
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
                atom_type = parts[4]
                atom_mass = float(parts[8])
                lj_radius_amber = float(parts[6])
                lj_depth_amber = float(parts[7])
                lj_sigma = lj_radius_amber * (1 / (2 ** (1 / 6))) * 2

                canonical_name = atom_type

                if canonical_name in nonbonded_params:
                    existing = nonbonded_params[canonical_name]
                    existing_mass = mass_list[atom_type_mapping[canonical_name] - 1]
                    if (
                        abs(existing["lj_epsilon"] - lj_depth_amber) > param_tol
                        or abs(existing["lj_sigma"] - lj_sigma) > param_tol
                        or abs(existing_mass - atom_mass) > param_tol
                    ):
                        base_name = f"{atom_type}_top{i+1}"
                        canonical_name = base_name
                        suffix = 2
                        while canonical_name in nonbonded_params:
                            canonical_name = f"{base_name}_{suffix}"
                            suffix += 1
                        print(f"Atom type conflict for '{atom_type}' between topologies; renaming to '{canonical_name}' for topology {i+1}")

                type_remap[atom_type] = canonical_name
                atom_types_per_topology[-1].add(canonical_name)
                type_origins.setdefault(canonical_name, set()).add(i + 1)

                if canonical_name not in atom_type_mapping:
                    atom_type_mapping[canonical_name] = len(atom_type_mapping) + 1
                    mass_list.append(atom_mass)

                if canonical_name not in nonbonded_params:
                    nonbonded_params[canonical_name] = {
                        "lj_epsilon": lj_depth_amber,
                        "lj_sigma": lj_sigma,
                    }
            except (ValueError, IndexError):
                continue

    if verbose:
        print(f"Found {len(atom_type_mapping)} unique atom types")

    return TopologyContext(
        all_parms=all_parms,
        atom_type_mapping=atom_type_mapping,
        mass_list=mass_list,
        nonbonded_params=nonbonded_params,
        total_atoms_per_topology=total_atoms_per_topology,
        atom_types_per_topology=atom_types_per_topology,
        type_remaps=type_remaps,
        type_origins=type_origins,
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


def write_data_header(data_file, expected_total_atoms, top_ctx, conn_ctx, box_bounds):
    with open(data_file, "a") as f:
        f.write(f"{expected_total_atoms} atoms \n")
        f.write(f"{len(top_ctx.atom_type_mapping)} atom types \n")
        f.write(f"{conn_ctx.bond_count} bonds \n")
        f.write(f"{conn_ctx.bond_type_count} bond types \n")
        f.write(f"{conn_ctx.angle_count} angles \n")
        f.write(f"{conn_ctx.angle_type_count} angle types \n")
        f.write(f"{conn_ctx.dihedral_count} dihedrals \n")
        f.write(f"{conn_ctx.dihedral_type_count} dihedral types \n\n")
        f.write(f"{box_bounds.xlo} {box_bounds.xhi} xlo xhi \n")
        f.write(f"{box_bounds.ylo} {box_bounds.yhi} ylo yhi \n")
        f.write(f"{box_bounds.zlo} {box_bounds.zhi} zlo zhi \n\n")
        f.write("Masses \n\n")
        for i in range(len(top_ctx.atom_type_mapping)):
            f.write(f"{i+1} {top_ctx.mass_list[i]} \n")


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
            if verbose:
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
            parm = all_parms[span.topo_idx]
            atom_type_str = parm.atoms[atom_idx_in_mol].type
            canonical_atom_type = type_remaps[span.topo_idx].get(atom_type_str, atom_type_str)
            atom_type_id = atom_type_mapping.get(canonical_atom_type, 1)

            flammps.write(
                f"{i+1} {span.mol_id} {atom_type_id} {charges[i]:.10f} "
                f"{pdb_atom['x']:.4f} {pdb_atom['y']:.4f} {pdb_atom['z']:.4f} 0 0 0\n"
            )


def process_connectivity(all_parms, molecule_counts, topologies):
    bond_type_registry = {}
    angle_type_registry = {}
    dihedral_type_registry = {}
    bond_coeff_lines = []
    angle_coeff_lines = []
    dihedral_coeff_lines = []
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

    for topo_idx, parm in enumerate(all_parms):
        bond_type_ids = []
        for bond in parm.bonds:
            if bond.type is None:
                raise ValueError(f"Bond parameters missing for atoms {bond.atom1.idx}-{bond.atom2.idx} in topology {topo_idx+1}")
            k = float(bond.type.k)
            req = float(bond.type.req)
            sig = (round(k, 6), round(req, 6))
            if sig not in bond_type_registry:
                type_id = len(bond_type_registry) + 1
                bond_type_registry[sig] = type_id
                bond_coeff_lines.append(
                    f"bond_coeff {type_id} {k:.4f} {req:.4f}  # {bond.atom1.type}-{bond.atom2.type}"
                )
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
            sig = (round(k, 6), round(theteq, 6))
            if sig not in angle_type_registry:
                type_id = len(angle_type_registry) + 1
                angle_type_registry[sig] = type_id
                angle_coeff_lines.append(
                    f"angle_coeff {type_id} {k:.4f} {theteq:.4f}  # {angle.atom1.type}-{angle.atom2.type}-{angle.atom3.type}"
                )
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
            sig = dihedral_signature(terms)
            if sig not in dihedral_type_registry:
                type_id = len(dihedral_type_registry) + 1
                dihedral_type_registry[sig] = type_id
                coeff_str = " ".join(f"{phi_k:.4f} {per} {phase:.4f}" for phi_k, per, phase in sig)
                dihedral_coeff_lines.append(
                    f"dihedral_coeff {type_id} {len(sig)} {coeff_str}  # {dih.atom1.type}-{dih.atom2.type}-{dih.atom3.type}-{dih.atom4.type}"
                )
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



def write_debug_files(keep_temp, topologies, atom_types_per_topology, atom_type_mapping, nonbonded_params, type_origins, conn_ctx):
    if not keep_temp:
        return

    with open("pairs.txt", "w") as ftemp:
        ftemp.write("# pairs.txt generated by amber_to_lammps.py\n")
        ftemp.write("# columns: type_id\tname\tepsilon\tsigma\ttopology_sources\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            for atom_type in sorted(atom_types_per_topology[topo_idx]):
                type_id = atom_type_mapping[atom_type]
                params = nonbonded_params[atom_type]
                origins = ",".join(str(t) for t in sorted(type_origins.get(atom_type, {topo_idx+1})))
                ftemp.write(f"{type_id}\t{atom_type}\t{params['lj_epsilon']:.6f}\t{params['lj_sigma']:.6f}\t{origins}\n")
            ftemp.write("\n")

    with open("bonds.txt", "w") as ftemp:
        ftemp.write("# bonds.txt generated by amber_to_lammps.py\n")
        ftemp.write("# columns: instance_id\ttype_id\tatom1\tatom2\tk\treq\ttopology_idx\tmolecule_idx\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            topo_lines = [line for line in conn_ctx.bond_debug if line.split("\t")[-2] == str(topo_idx + 1)]
            if topo_lines:
                ftemp.write("\n".join(topo_lines) + "\n")
            else:
                ftemp.write("# (none)\n")
            ftemp.write("\n")

    with open("angles.txt", "w") as ftemp:
        ftemp.write("# angles.txt generated by amber_to_lammps.py\n")
        ftemp.write("# columns: instance_id\ttype_id\ta1\ta2\ta3\tk\ttheta\ttopology_idx\tmolecule_idx\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            topo_lines = [line for line in conn_ctx.angle_debug if line.split("\t")[-2] == str(topo_idx + 1)]
            if topo_lines:
                ftemp.write("\n".join(topo_lines) + "\n")
            else:
                ftemp.write("# (none)\n")
            ftemp.write("\n")

    with open("dihedrals.txt", "w") as ftemp:
        ftemp.write("# dihedrals.txt generated by amber_to_lammps.py\n")
        ftemp.write("# columns: instance_id\ttype_id\ta1\ta2\ta3\ta4\ttopology_idx\tmolecule_idx\n\n")
        for topo_idx, topo_name in enumerate(topologies):
            ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
            topo_lines = [line for line in conn_ctx.dihedral_debug if line.split("\t")[-2] == str(topo_idx + 1)]
            if topo_lines:
                ftemp.write("\n".join(topo_lines) + "\n")
            else:
                ftemp.write("# (none)\n")
            ftemp.write("\n")



def write_parameter_file(param_file, topologies, atom_types_per_topology, pair_coeff_map, conn_ctx):
    printed_pair_types = set()
    with open(param_file, "a") as flammpsparm:
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


def amber2lammps(data_file, param_file, topologies, molecule_counts, pdb_file, charges_target, buffer=3.8, verbose=False, keep_temp=False):
    amber_parm_cls = pmd.amber.AmberParm
    print_details = pmd.tools.actions.printDetails

    multi_mode = detect_multi_mode(topologies, molecule_counts, pdb_file, verbose)
    cleanup_temp_files(verbose, keep_temp)

    if verbose:
        print("Converting multiple AMBER topologies to LAMMPS format...")
        print(f"Output files: {data_file}, {param_file}")
        for i, (topo, count) in enumerate(zip(topologies, molecule_counts)):
            print(f"  Topology {i+1}: {topo} ({count} molecules)")

    initialize_output_files(data_file, param_file, topologies, pdb_file)
    pdb_atoms, x_coords, y_coords, z_coords = parse_pdb_coordinates(pdb_file, verbose)
    top_ctx = load_topology_context(topologies, verbose, amber_parm_cls, print_details)

    expected_total_atoms = validate_pdb_atom_count(
        pdb_atoms=pdb_atoms,
        molecule_counts=molecule_counts,
        total_atoms_per_topology=top_ctx.total_atoms_per_topology,
        multi_mode=multi_mode,
        topologies=topologies,
        pdb_file=pdb_file,
        verbose=verbose,
    )

    box_bounds = compute_box_bounds(x_coords, y_coords, z_coords, buffer)
    if verbose:
        print(
            f"Box dimensions: X[{box_bounds.xlo:.3f}, {box_bounds.xhi:.3f}], "
            f"Y[{box_bounds.ylo:.3f}, {box_bounds.yhi:.3f}], "
            f"Z[{box_bounds.zlo:.3f}, {box_bounds.zhi:.3f}]"
        )

    conn_ctx = process_connectivity(top_ctx.all_parms, molecule_counts, topologies)

    write_data_header(
        data_file=data_file,
        expected_total_atoms=expected_total_atoms,
        top_ctx=top_ctx,
        conn_ctx=conn_ctx,
        box_bounds=box_bounds,
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
    )

    write_parameter_file(
        param_file=param_file,
        topologies=topologies,
        atom_types_per_topology=top_ctx.atom_types_per_topology,
        pair_coeff_map=pair_coeff_map,
        conn_ctx=conn_ctx,
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

    cleanup_temp_files(verbose, keep_temp)

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
            keep_temp=args.keep_temp
        )
        print(f"✓ Conversion completed successfully!")
        print(f"Output files: {args.data_file}, {args.param_file}")
        return True
        
    except Exception as e:
        print(f"✗ Error during conversion: {e}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)

    
    
