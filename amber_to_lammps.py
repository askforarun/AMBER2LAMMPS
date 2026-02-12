
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
import parmed as pmd
import numpy as np

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
                    atom_num = int(line[6:11].strip())
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
                    
                except (ValueError, IndexError) as e:
                    if verbose:
                        print(f"Warning: Could not parse PDB line: {line.strip()}")
                    continue
    
    if verbose:
        print(f"Parsed {len(atoms)} atoms from PDB file")
    
    if len(atoms) == 0:
        raise ValueError(f"No atoms found in PDB file '{pdb_file}'. File may be empty or malformed.")
    
    return atoms, x, y, z

def amber2lammps(data_file, param_file, topologies, molecule_counts, pdb_file, charges_target, buffer=3.8, verbose=False, keep_temp=False):
    AmberParm=pmd.amber.AmberParm
    printDetails = pmd.tools.actions.printDetails

    # Auto-detect multi-molecule scenarios
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
        # Check if PDB appears to be a combined file even for single molecule request
        try:
            with open(pdb_file, 'r') as f:
                atom_count = 0
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        atom_count += 1
            
            # Load topology to get expected single molecule atom count
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
        print(f"✓ Single molecule mode: 1 topology, 1 copy")

    # Clean up temporary files if they exist
    cleanup_temp_files(verbose, keep_temp)

    # Setup output files
    if verbose:
        print(f"Converting multiple AMBER topologies to LAMMPS format...")
        print(f"Output files: {data_file}, {param_file}")
        for i, (topo, count) in enumerate(zip(topologies, molecule_counts)):
            print(f"  Topology {i+1}: {topo} ({count} molecules)")
    
    # Initialize output files
    with open(data_file, "w") as f:
        f.write(f"LAMMPS data file from AMBER conversion\n")
        f.write(f"#Source: {', '.join(topologies)}, {pdb_file}\n\n")

    with open(param_file, "w") as f:
        f.write(f"# Force field parameters generated from AMBER topologies\n\n")

    # Parse PDB coordinates
    pdb_atoms, x_coords, y_coords, z_coords = parse_pdb_coordinates(pdb_file, verbose)
    
    # Load all topologies and collect atom type information
    all_parms = []
    atom_type_mapping = {}
    mass_list = []
    nonbonded_params = {}
    total_atoms_per_topology = []
    atom_types_per_topology = []
    type_remaps = []  # per-topology map: original atom type -> canonical unique name
    param_tol = 1e-6
    type_origins = {}  # canonical atom type -> set of topology indices that defined it
    
    for i, topology in enumerate(topologies):
        if verbose:
            print(f"Loading topology {i+1}: {topology}")
        
        parm = AmberParm(topology)
        all_parms.append(parm)
        total_atoms_per_topology.append(len(parm.atoms))
        atom_types_per_topology.append(set())
        type_remap = {}
        type_remaps.append(type_remap)
        
        if verbose:
            print(f"  Found {len(parm.atoms)} atoms, {len(parm.bonds)} bonds, {len(parm.angles)} angles, {len(parm.dihedrals)} dihedrals")
        
        # Extract LJ coefficients and atom information using printDetails
        lj_details = printDetails(parm, "@1-{}".format(len(parm.atoms)))
        
        # Parse atom types and masses
        for line in str(lj_details).split('\n'):
            line = line.strip()
            if not line or not line[0].isdigit():
                continue
                
            parts = line.split()
            if len(parts) >= 10:
                try:
                    atom_type = parts[4]
                    atom_mass = float(parts[8])
                    lj_radius_amber = float(parts[6])
                    lj_depth_amber = float(parts[7])
                    lj_sigma = lj_radius_amber * (1/(2**(1/6))) * 2

                    canonical_name = atom_type

                    if canonical_name in nonbonded_params:
                        existing = nonbonded_params[canonical_name]
                        existing_mass = mass_list[atom_type_mapping[canonical_name]-1]
                        if (abs(existing['lj_epsilon'] - lj_depth_amber) > param_tol or
                            abs(existing['lj_sigma'] - lj_sigma) > param_tol or
                            abs(existing_mass - atom_mass) > param_tol):
                            # Same atom type label but different parameters: namespace per topology
                            base_name = f"{atom_type}_top{i+1}"
                            canonical_name = base_name
                            suffix = 2
                            while canonical_name in nonbonded_params:
                                canonical_name = f"{base_name}_{suffix}"
                                suffix += 1
                            print(f"Atom type conflict for '{atom_type}' between topologies; renaming to '{canonical_name}' for topology {i+1}")

                    type_remap[atom_type] = canonical_name
                    atom_types_per_topology[-1].add(canonical_name)

                    # Add atom type if not already present
                    type_origins.setdefault(canonical_name, set()).add(i + 1)

                    if canonical_name not in atom_type_mapping:
                        atom_type_mapping[canonical_name] = len(atom_type_mapping) + 1
                        mass_list.append(atom_mass)
                    
                    # Store LJ parameters if not already stored
                    if canonical_name not in nonbonded_params:
                        nonbonded_params[canonical_name] = {
                            'lj_epsilon': lj_depth_amber,
                            'lj_sigma': lj_sigma
                        }
                        
                except (ValueError, IndexError) as e:
                    continue
    
    if verbose:
        print(f"Found {len(atom_type_mapping)} unique atom types")
    
    # Calculate expected total atoms and sanity-check vs PDB
    contrib = [f"{count}*{atoms_per_topo}={count * atoms_per_topo}"
               for count, atoms_per_topo in zip(molecule_counts, total_atoms_per_topology)]
    expected_total_atoms = sum(count * atoms_per_topo for count, atoms_per_topo in zip(molecule_counts, total_atoms_per_topology))
    
    if len(pdb_atoms) != expected_total_atoms:
        breakdown = " + ".join(contrib)
        raise ValueError(
            f"Atom count mismatch: PDB has {len(pdb_atoms)} atoms but expected {expected_total_atoms} "
            f"({breakdown}). Check molecule counts, topology order, and PDB atom ordering."
        )
    elif verbose:
        breakdown = " + ".join(contrib)
        print(f"Atom count check passed: PDB={len(pdb_atoms)} matches expected {expected_total_atoms} ({breakdown})")
    
    # Additional validation for multi-mode: ensure PDB appears to be combined
    if multi_mode:
        single_molecule_atoms = total_atoms_per_topology[0] if len(topologies) == 1 else sum(total_atoms_per_topology)
        if len(pdb_atoms) <= single_molecule_atoms:
            print(f"Warning: PDB file '{pdb_file}' contains {len(pdb_atoms)} atoms, but multi-mode expects")
            print(f"         a combined PDB file with {expected_total_atoms} atoms.")
            print(f"         Ensure you're using a combined PDB file from PackMol or similar tool.")
            print(f"         Expected: {expected_total_atoms} atoms ({breakdown})")
        else:
            print(f"✓ PDB file '{pdb_file}' contains {len(pdb_atoms)} atoms (appears to be combined)")
    
    # Calculate box dimensions with buffer
    xlo = np.min(x_coords) - buffer
    xhi = np.max(x_coords) + buffer
    ylo = np.min(y_coords) - buffer
    yhi = np.max(y_coords) + buffer
    zlo = np.min(z_coords) - buffer
    zhi = np.max(z_coords) + buffer
    
    if verbose:
        print(f"Box dimensions: X[{xlo:.3f}, {xhi:.3f}], Y[{ylo:.3f}, {yhi:.3f}], Z[{zlo:.3f}, {zhi:.3f}]")
    
    # Calculate total bonds, angles, dihedrals
    total_bonds = sum(count * len(parm.bonds) for count, parm in zip(molecule_counts, all_parms))
    total_angles = sum(count * len(parm.angles) for count, parm in zip(molecule_counts, all_parms))
    total_dihedrals = sum(count * len(parm.dihedrals) for count, parm in zip(molecule_counts, all_parms))
    
    # Write header information to data file
    with open(data_file, "a") as f:
        f.write(f"{expected_total_atoms} atoms \n")
        f.write(f"{len(atom_type_mapping)} atom types \n")
        f.write(f"{total_bonds} bonds \n")
        f.write(f"{total_bonds} bond types \n")
        f.write(f"{total_angles} angles \n")
        f.write(f"{total_angles} angle types \n")
        f.write(f"{total_dihedrals} dihedrals \n")
        f.write(f"{total_dihedrals} dihedral types \n\n")
        f.write(f"{xlo} {xhi} xlo xhi \n")
        f.write(f"{ylo} {yhi} ylo yhi \n")
        f.write(f"{zlo} {zhi} zlo zhi \n\n")
        f.write("Masses \n\n")
        for i in range(len(atom_type_mapping)):
            f.write(f"{i+1} {mass_list[i]} \n")
    
    # Extract charges from all topologies and adjust per-topology to user targets
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
                print(f"Charge adjust topo {topo_idx+1}: actual {actual_charge_per_mol:.6f} -> target {target_charge_per_mol:.6f} (shift {shift:.6f}/atom)")
        elif verbose:
            print(f"Charge check topo {topo_idx+1}: {actual_charge_per_mol:.6f} matches target {target_charge_per_mol:.6f}")
        
        # Repeat adjusted charges for each molecule of this type
        for _ in range(count):
            charges.extend(topo_charges.tolist())
    
    # Final sanity on total charge
    net_charge = float(np.sum(charges))
    target_total_charge = float(np.sum([c * q for c, q in zip(molecule_counts, charges_target)]))
    if abs(net_charge - target_total_charge) > 1e-4:
        raise ValueError(f"Total charge mismatch after adjustment: got {net_charge:.6f}, expected {target_total_charge:.6f}")
    elif verbose:
        print(f"Total charge check passed: {net_charge:.6f} matches expected {target_total_charge:.6f}")
    
    # Precompute molecule spans to map atoms to molecule IDs and topologies
    molecule_spans = []
    running_offset = 0
    for topo_idx, (parm, count) in enumerate(zip(all_parms, molecule_counts)):
        atoms_per_molecule = len(parm.atoms)
        for _ in range(count):
            start = running_offset
            end = start + atoms_per_molecule - 1
            molecule_spans.append({
                'id': len(molecule_spans) + 1,
                'topo_idx': topo_idx,
                'start': start,
                'end': end,
                'atoms_per_molecule': atoms_per_molecule
            })
            running_offset += atoms_per_molecule
    
    # Generate pair coefficients
    pair_coeff_map = {}
    for atom_type in atom_type_mapping.keys():
        if atom_type in nonbonded_params:
            lj_epsilon = nonbonded_params[atom_type]['lj_epsilon']
            lj_sigma = nonbonded_params[atom_type]['lj_sigma']
            type_id = atom_type_mapping[atom_type]
            pair_coeff_map[atom_type] = f"pair_coeff {type_id} {type_id} {lj_epsilon} {lj_sigma} # {atom_type}"
    
    # Write atoms section
    if verbose:
        print("Writing atoms section...")
    
    with open(data_file, "a") as flammps:
        flammps.write("Atoms\n\n")
        span_idx = 0
        for i, pdb_atom in enumerate(pdb_atoms):
            # Advance span index until the current atom falls inside the molecule span
            while span_idx < len(molecule_spans) and i > molecule_spans[span_idx]['end']:
                span_idx += 1
            
            if span_idx >= len(molecule_spans):
                raise ValueError(f"Atom index {i} exceeds computed molecule spans; check PDB ordering.")
            
            span = molecule_spans[span_idx]
            atom_idx_in_mol = i - span['start']
            parm = all_parms[span['topo_idx']]
            atom_type_str = parm.atoms[atom_idx_in_mol].type
            canonical_atom_type = type_remaps[span['topo_idx']].get(atom_type_str, atom_type_str)
            atom_type_id = atom_type_mapping.get(canonical_atom_type, 1)
            molecule_id = span['id']
            
            flammps.write(
                f"{i+1} {molecule_id} {atom_type_id} {charges[i]:.10f} "
                f"{pdb_atom['x']:.4f} {pdb_atom['y']:.4f} {pdb_atom['z']:.4f} 0 0 0\n"
            )
    
    # Process bonds, angles, dihedrals for all topologies
    bond_count = 0
    bond_lines = []
    bond_coeffs_by_topo = [[] for _ in topologies]
    bond_debug = []  # type_id atom1 atom2 k req topo_idx mol_idx
    
    angle_count = 0
    angle_lines = []
    angle_coeffs_by_topo = [[] for _ in topologies]
    angle_debug = []  # type_id a1 a2 a3 k theta topo_idx mol_idx
    
    dihedral_count = 0
    dihedral_lines = []
    dihedral_coeffs_by_topo = [[] for _ in topologies]
    dihedral_debug = []  # type_id a1 a2 a3 a4 n_terms (k n phase)* topo_idx mol_idx
    
    atom_offset = 0
    
    for topo_idx, (parm, count) in enumerate(zip(all_parms, molecule_counts)):
        atoms_per_molecule = len(parm.atoms)
        
        # Process each replica of this topology
        for mol_idx in range(count):
            base_offset = atom_offset + mol_idx * atoms_per_molecule
            
            # Bonds
            for bond in parm.bonds:
                if bond.type is None:
                    raise ValueError(f"Bond parameters missing for atoms {bond.atom1.idx}-{bond.atom2.idx} in topology {topo_idx+1}")
                bond_count += 1
                atom1 = bond.atom1.idx + 1 + base_offset
                atom2 = bond.atom2.idx + 1 + base_offset
                k = float(bond.type.k)
                r0 = float(bond.type.req)
                bond_lines.append(f"{bond_count} {bond_count} {atom1} {atom2}")
                bond_coeffs_by_topo[topo_idx].append(f"bond_coeff {bond_count} {k:.4f} {r0:.4f}")
                bond_debug.append(f"{bond_count}\t{bond_count}\t{atom1}\t{atom2}\t{k:.6f}\t{r0:.6f}\t{topo_idx+1}\t{mol_idx+1}")
            
            # Angles
            for angle in parm.angles:
                if angle.type is None:
                    raise ValueError(f"Angle parameters missing for atoms {angle.atom1.idx}-{angle.atom2.idx}-{angle.atom3.idx} in topology {topo_idx+1}")
                angle_count += 1
                atom1 = angle.atom1.idx + 1 + base_offset
                atom2 = angle.atom2.idx + 1 + base_offset
                atom3 = angle.atom3.idx + 1 + base_offset
                k = float(angle.type.k)
                theta0 = float(angle.type.theteq)
                angle_lines.append(f"{angle_count} {angle_count} {atom1} {atom2} {atom3}")
                angle_coeffs_by_topo[topo_idx].append(f"angle_coeff {angle_count} {k:.4f} {theta0:.4f}")
                angle_debug.append(f"{angle_count}\t{angle_count}\t{atom1}\t{atom2}\t{atom3}\t{k:.6f}\t{theta0:.6f}\t{topo_idx+1}\t{mol_idx+1}")
            
            # Dihedrals
            for dih in parm.dihedrals:
                # Skip impropers that may be stored separately
                if dih.type is None:
                    raise ValueError(f"Dihedral parameters missing for atoms {dih.atom1.idx}-{dih.atom2.idx}-{dih.atom3.idx}-{dih.atom4.idx} in topology {topo_idx+1}")
                
                # Collect all multi-term dihedral components (AMBER stores Fourier series as a list)
                if isinstance(dih.type, (list, tuple)):
                    terms = [t for t in dih.type if t is not None]
                else:
                    terms = [dih.type]
                
                if not terms:
                    raise ValueError(f"Dihedral parameters missing for atoms {dih.atom1.idx}-{dih.atom2.idx}-{dih.atom3.idx}-{dih.atom4.idx} in topology {topo_idx+1}")
                
                dihedral_count += 1
                atom1 = dih.atom1.idx + 1 + base_offset
                atom2 = dih.atom2.idx + 1 + base_offset
                atom3 = dih.atom3.idx + 1 + base_offset
                atom4 = dih.atom4.idx + 1 + base_offset
                
                coeff_parts = []
                for term in terms:
                    phi_k = float(term.phi_k)
                    per = int(round(float(term.per)))
                    phase = float(term.phase)
                    # Heuristic: convert to degrees if value looks like radians
                    if abs(phase) <= 2 * np.pi + 0.1:
                        phase = np.degrees(phase)
                    coeff_parts.append(f"{phi_k:.4f} {per} {phase:.4f}")
                
                n_terms = len(coeff_parts)
                dihedral_lines.append(f"{dihedral_count} {dihedral_count} {atom1} {atom2} {atom3} {atom4}")
                dihedral_coeffs_by_topo[topo_idx].append(f"dihedral_coeff {dihedral_count} {n_terms} " + " ".join(coeff_parts))
                dihedral_debug.append(
                    f"{dihedral_count}\t{dihedral_count}\t{atom1}\t{atom2}\t{atom3}\t{atom4}\t{n_terms}\t"
                    + "\t".join(coeff_parts)
                    + f"\t{topo_idx+1}\t{mol_idx+1}"
                )
        
        atom_offset += count * atoms_per_molecule
    
    # Write bonds section
    if verbose:
        print("Writing bonds section...")
    
    with open(data_file, "a") as flammps:
        flammps.write("\nBonds \n\n")
        flammps.write("\n".join(bond_lines) + "\n")
    
    # Write angles section
    if verbose:
        print("Writing angles section...")
    
    with open(data_file, "a") as flammps:
        flammps.write("\nAngles \n\n")
        flammps.write("\n".join(angle_lines) + "\n")
    
    # Write dihedrals section
    if verbose:
        print("Writing dihedrals section...")
    
    with open(data_file, "a") as flammps:
        flammps.write("\nDihedrals \n\n")
        flammps.write("\n".join(dihedral_lines) + "\n")

    # Optionally write debug temp files grouped per topology
    if keep_temp:
        # pairs.txt: per-topology atom type params (no cross terms)
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

        # bonds.txt
        with open("bonds.txt", "w") as ftemp:
            ftemp.write("# bonds.txt generated by amber_to_lammps.py\n")
            ftemp.write("# columns: bond_id\tbond_type\tatom1\tatom2\tk\treq\ttopology_idx\tmolecule_idx\n\n")
            for topo_idx, topo_name in enumerate(topologies):
                ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
                topo_lines = [line for line in bond_debug if line.split("\t")[-2] == str(topo_idx + 1)]
                if topo_lines:
                    ftemp.write("\n".join(topo_lines) + "\n")
                else:
                    ftemp.write("# (none)\n")
                ftemp.write("\n")

        # angles.txt
        with open("angles.txt", "w") as ftemp:
            ftemp.write("# angles.txt generated by amber_to_lammps.py\n")
            ftemp.write("# columns: angle_id\tangle_type\ta1\ta2\ta3\tk\ttheta\ttopology_idx\tmolecule_idx\n\n")
            for topo_idx, topo_name in enumerate(topologies):
                ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
                topo_lines = [line for line in angle_debug if line.split("\t")[-2] == str(topo_idx + 1)]
                if topo_lines:
                    ftemp.write("\n".join(topo_lines) + "\n")
                else:
                    ftemp.write("# (none)\n")
                ftemp.write("\n")

        # dihedrals.txt
        with open("dihedrals.txt", "w") as ftemp:
            ftemp.write("# dihedrals.txt generated by amber_to_lammps.py\n")
            ftemp.write("# columns: dih_id\tdih_type\ta1\ta2\ta3\ta4\tn_terms\t(k n phase)*\ttopology_idx\tmolecule_idx\n\n")
            for topo_idx, topo_name in enumerate(topologies):
                ftemp.write(f"# Topology {topo_idx+1}: {topo_name}\n")
                topo_lines = [line for line in dihedral_debug if line.split("\t")[-2] == str(topo_idx + 1)]
                if topo_lines:
                    ftemp.write("\n".join(topo_lines) + "\n")
                else:
                    ftemp.write("# (none)\n")
                ftemp.write("\n")

    # Write grouped parameter file sections per topology
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
            
            flammpsparm.write("# Bonds\n")
            if bond_coeffs_by_topo[topo_idx]:
                flammpsparm.write("\n".join(bond_coeffs_by_topo[topo_idx]) + "\n")
            else:
                flammpsparm.write("# (none)\n")
            
            flammpsparm.write("# Angles\n")
            if angle_coeffs_by_topo[topo_idx]:
                flammpsparm.write("\n".join(angle_coeffs_by_topo[topo_idx]) + "\n")
            else:
                flammpsparm.write("# (none)\n")
            
            flammpsparm.write("# Dihedrals\n")
            if dihedral_coeffs_by_topo[topo_idx]:
                flammpsparm.write("\n".join(dihedral_coeffs_by_topo[topo_idx]) + "\n")
            else:
                flammpsparm.write("# (none)\n")
    
    # Clean up temporary files
    if not keep_temp:
        temp_files = ["bonds.txt", "angles.txt", "dihedrals.txt", "pairs.txt"]
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
    
    if verbose:
        print(f"Conversion complete!")
        print(f"Generated files:")
        print(f"  - {data_file} (LAMMPS data file)")
        print(f"  - {param_file} (LAMMPS parameters)")
        print(f"Summary:")
        print(f"  - {expected_total_atoms} atoms")
        print(f"  - {bond_count} bonds") 
        print(f"  - {angle_count} angles")
        print(f"  - {dihedral_count} dihedrals")
    
    # Clean up temporary files at the end
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

    
    
