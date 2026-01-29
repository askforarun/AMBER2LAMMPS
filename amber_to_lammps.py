
#!/usr/bin/env python3
"""
AMBER to LAMMPS conversion script

This script converts AMBER topology and coordinate files to LAMMPS data format.
It requires parmed, numpy, and the AMBER force field files.

Usage:
    python amber_to_lammps.py <data_file> <param_file> <topology> <crd>

Arguments:
    data_file    Output LAMMPS data file name
    param_file   Output LAMMPS parameter file name  
    topology     AMBER topology file (.prmtop)
    crd          AMBER coordinate file (.crd)
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
    python amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.crd
    python amber_to_lammps.py system.data system.parm system.prmtop system.crd
        """
    )
    
    parser.add_argument('data_file', help='Output LAMMPS data file name')
    parser.add_argument('param_file', help='Output LAMMPS parameter file name')
    parser.add_argument('topology', help='AMBER topology file (.prmtop)')
    parser.add_argument('crd', help='AMBER coordinate file (.crd)')
    parser.add_argument('-b', '--buffer', type=float, default=3.8,
                        help='Buffer size around molecule (default: 3.8)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep temporary files (bonds.txt, angles.txt, dihedrals.txt) after conversion')
    parser.add_argument('--charge', type=int, required=True,
                        help='Target net charge (integer). If 0, charges will be normalized to neutrality. If non-zero, system will maintain this charge')
    
    return parser.parse_args()

def validate_files(topology, crd):
    """Validate input files exist and are readable"""
    files = [topology, crd]
    file_types = ["topology", "coordinate"]
    
    for file, file_type in zip(files, file_types):
        if not os.path.exists(file):
            print(f"Error: {file_type} file '{file}' not found")
            sys.exit(1)
        if not os.access(file, os.R_OK):
            print(f"Error: {file_type} file '{file}' is not readable")
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

def amber2lammps(data_file, param_file, topology, crd, charge, buffer=3.8, verbose=False, keep_temp=False):
    AmberParm=pmd.amber.AmberParm
    printBonds=pmd.tools.actions.printBonds
    printAngles=pmd.tools.actions.printAngles
    printDihedrals=pmd.tools.actions.printDihedrals

    # Clean up temporary files if they exist
    cleanup_temp_files(verbose, keep_temp)

    # Setup output files
    if verbose:
        print(f"Converting AMBER files to LAMMPS format...")
        print(f"Output files: {data_file}, {param_file}")
    
    # Initialize output files
    with open(data_file, "w") as f:
        f.write(f"LAMMPS data file from AMBER conversion\n")
        f.write(f"#Source: {topology}, {crd}\n\n")

    with open(param_file, "w") as f:
        f.write(f"# Force field parameters generated from AMBER topology\n\n")

    
    # Load AMBER topology
    if verbose:
        print(f"Loading topology from {topology}...")
    
    parm = AmberParm(topology)
    natoms = len(parm.atoms)
    
    if verbose:
        print(f"Found {natoms} atoms")
        print(f"Found {len(parm.bonds)} bonds")
        print(f"Found {len(parm.angles)} angles") 
        print(f"Found {len(parm.dihedrals)} dihedrals")

    # Extract bonds, angles, dihedrals to temporary files
    Bonds = printBonds(parm, "@1-{}".format(natoms))
    with open("bonds.txt", "w") as f:
        print(Bonds, file=f)

    Angles = printAngles(parm, "@1-{}".format(natoms))
    with open("angles.txt", "w") as f:
        print(Angles, file=f)

    Dihedrals = printDihedrals(parm, "@1-{}".format(natoms))
    with open("dihedrals.txt", "w") as f:
        print(Dihedrals, file=f)

    # Extract LJ coefficients and atom information using printDetails
    if verbose:
        print("Extracting LJ coefficients and atom information from topology...")
    
    printDetails = pmd.tools.actions.printDetails
    lj_details = printDetails(parm, "@1-{}".format(natoms))

    # Parse atom types and masses from topology using printDetails
    if verbose:
        print(f"Parsing atom types and masses from topology...")
    
    atype = {}
    mass = []
    
    # Use printDetails to get atom information including masses
    for line in str(lj_details).split('\n'):
        line = line.strip()
        if not line or not line[0].isdigit():
            continue
            
        parts = line.split()
        if len(parts) >= 10:
            try:
                atom_num = int(parts[0])
                atom_name = parts[3]
                atom_type = parts[4]
                atom_mass = float(parts[8])  # Mass is in column 8
                
                # Add atom type if not already present
                if atom_type not in atype:
                    atype[atom_type] = len(atype) + 1
                    mass.append(atom_mass)
                    
            except (ValueError, IndexError) as e:
                continue
    
    if verbose:
        print(f"Found {len(atype)} atom types")

    # Parse coordinates from CRD file
    if verbose:
        print(f"Parsing coordinates from {crd}...")
    
    x, y, z = [], [], []
    charges = []
    with open(crd) as f:
        lines = f.readlines()
        # Skip header lines (first 2 lines)
        for line in lines[2:]:
            if line.strip():
                # CRD format has coordinates in groups of 3 per line
                coords = line.split()
                for i in range(0, len(coords), 3):
                    if i + 2 < len(coords):
                        try:
                            x.append(float(coords[i]))
                            y.append(float(coords[i+1]))
                            z.append(float(coords[i+2]))
                        except (ValueError, IndexError):
                            continue
    
    # Extract charges from topology
    if verbose:
        print("Extracting charges from topology...")
    
    charges = [atom.charge for atom in parm.atoms]

    # Validate coordinate parsing
    if verbose:
        print(f"Parsed {len(x)} coordinates from CRD file")
    
    if len(x) == 0:
        raise ValueError(f"No coordinates found in CRD file '{crd}'. File may be empty or malformed.")
    
    if len(x) != natoms:
        raise ValueError(f"Coordinate count mismatch: Found {len(x)} coordinates but expected {natoms} atoms from topology. CRD file '{crd}' may be incomplete or corrupted.")

    # Calculate box dimensions with buffer
    xlo = np.min(x) - buffer
    xhi = np.max(x) + buffer
    ylo = np.min(y) - buffer
    yhi = np.max(y) + buffer
    zlo = np.min(z) - buffer
    zhi = np.max(z) + buffer
    
    if verbose:
        print(f"Box dimensions: X[{xlo:.3f}, {xhi:.3f}], Y[{ylo:.3f}, {yhi:.3f}], Z[{zlo:.3f}, {zhi:.3f}]")

    # Write header information to data file
    with open(data_file, "a") as f:
        f.write(f"{natoms} atoms \n")
        f.write(f"{len(atype)} atom types \n")
        f.write(f"{len(parm.bonds)} bonds \n")
        f.write(f"{len(parm.bonds)} bond types \n")
        f.write(f"{len(parm.angles)} angles \n")
        f.write(f"{len(parm.angles)} angle types \n")
        f.write(f"{len(parm.dihedrals)} dihedrals \n")
        f.write(f"{len(parm.dihedrals)} dihedral types \n\n")
        f.write(f"{xlo} {xhi} xlo xhi \n")
        f.write(f"{ylo} {yhi} ylo yhi \n")
        f.write(f"{zlo} {zhi} zlo zhi \n\n")
        f.write("Masses \n\n")
        for i in range(len(atype)):
            f.write(f"{i+1} {mass[i]} \n")
                
    # Normalize charges to achieve target charge
    net_charge = float(np.sum(charges))
    if verbose:
        print(f"Current net charge: {net_charge:.6f}")
        print(f"Target net charge: {charge}")
    
    # Calculate required shift: target - current
    charge_difference = charge - net_charge
    charge_tol = 1e-6  # Same tolerance as neutrality check
    
    if abs(charge_difference) > charge_tol:
        shift = charge_difference / len(charges)
        # Apply the shift to each atom
        charges = [c + shift for c in charges]
        
        if verbose:
            print(f"Applied uniform shift of {shift:.6f} per atom to achieve target charge")
            print(f"New net charge: {np.sum(charges):.6f}")
    elif verbose:
        print(f"Net charge already matches target ({target_charge}) within tolerance ({charge_tol}); no shift applied")

    # Parse LJ coefficients from printDetails output
    nonbonded_params = {}
    lines = str(lj_details).split('\n')
    
    # Skip header lines and find the actual data
    for line in lines:
        line = line.strip()
        if not line or not line[0].isdigit():
            continue
            
        parts = line.split()
        if len(parts) >= 10:  # LJ details line format
            try:
                atom_num = int(parts[0])
                atom_name = parts[3]
                atom_type = parts[4]
                
                # Extract LJ parameters from specific columns
                # Column 6: LJ Radius (AMBER units), Column 7: LJ Depth (AMBER units)
                lj_radius_amber = float(parts[6])  # AMBER LJ radius
                lj_depth_amber = float(parts[7])   # AMBER LJ depth (epsilon)
                
                # Convert AMBER LJ radius to LAMMPS sigma
                # Conversion factor: (1/2^(1/6))*2 = 1.781798
                lj_sigma = lj_radius_amber * (1/(2**(1/6))) * 2
                lj_epsilon = lj_depth_amber  # AMBER epsilon is already in LAMMPS units
                
                # Store the LJ parameters if not already stored
                if atom_type not in nonbonded_params:
                    nonbonded_params[atom_type] = {
                        'lj_epsilon': lj_epsilon,
                        'lj_sigma': lj_sigma
                    }
                    
            except (ValueError, IndexError) as e:
                if verbose:
                    print(f"Warning: Could not parse line: {line}")
                continue

    # Generate pair coefficients from topology
    pair_coeffs = []
    for atom_type in atype.keys():
        if atom_type in nonbonded_params:
            lj_epsilon = nonbonded_params[atom_type]['lj_epsilon']
            lj_sigma = nonbonded_params[atom_type]['lj_sigma']
            type_id = atype[atom_type]
            pair_coeffs.append(f"pair_coeff {type_id} {type_id} {lj_epsilon} {lj_sigma} # {atom_type}")
    
    # Write all pair coefficients at once
    with open(param_file, "a") as flammps:
        flammps.write("\n".join(pair_coeffs) + "\n")

    # Generate pairs.txt with nonbonded entries
    with open("pairs.txt", "w") as f:
        f.write("# Nonbonded pairs: atom_number, atom_name, atom_type, epsilon, sigma\n")
        f.write(f"{'Atom':>6} {'Name':>6} {'Type':>6} {'Epsilon':>10} {'Sigma':>10}\n")
        f.write(f"{'-'*6} {'-'*6} {'-'*6} {'-'*10} {'-'*10}\n")
        atom_number = 1
        # Get atom types and names from topology
        for atom in parm.atoms:
            atom_name = atom.name
            atom_type_str = atom.type
            if atom_type_str in nonbonded_params:
                epsilon = nonbonded_params[atom_type_str]['lj_epsilon']
                sigma = nonbonded_params[atom_type_str]['lj_sigma']
                f.write(f"{atom_number:6d} {atom_name:>6s} {atom_type_str:>6s} {epsilon:10.6f} {sigma:10.6f}\n")
            else:
                f.write(f"{atom_number:6d} {atom_name:>6s} {atom_type_str:>6s} {0:10.6f} {0:10.6f} # No LJ parameters found\n")
            atom_number += 1

    # Write atoms section
    if verbose:
        print("Writing atoms section...")
    
    with open(data_file, "a") as flammps:
        flammps.write("Atoms\n\n")
        for i, atom in enumerate(parm.atoms):
            atom_type_str = atom.type
            if atom_type_str in atype:
                type_id = atype[atom_type_str]
            else:
                type_id = 1  # Default type if not found
            
            flammps.write(f"{i+1} {type_id} {type_id} {charges[i]:.10f} {x[i]:.4f} {y[i]:.4f} {z[i]:.4f} 0 0 0\n")

    # Write bonds section
    if verbose:
        print("Writing bonds section...")
    
    bond_count = 0
    bond_lines = []
    bond_coeffs = []
    
    with open("bonds.txt") as f:
        for line in f:
            if "Atom" not in line and line.split():
                bond_count += 1
                bond_lines.append(f"{bond_count} {bond_count} {int(line.split()[0])} {int(line.split()[4])}")
                bond_coeffs.append(f"bond_coeff {bond_count} {line.split()[9]} {line.split()[8]}")
    
    # Write all bonds at once
    with open(data_file, "a") as flammps:
        flammps.write("\nBonds \n\n")
        flammps.write("\n".join(bond_lines) + "\n")
    
    # Write all bond coefficients at once
    with open(param_file, "a") as flammpsparm:
        flammpsparm.write("\n".join(bond_coeffs) + "\n")

    # Write angles section
    if verbose:
        print("Writing angles section...")
    
    angle_count = 0
    angle_lines = []
    angle_coeffs = []
    
    with open("angles.txt") as f:
        for line in f:
            if "Atom" not in line and line.split():
                angle_count += 1
                angle_lines.append(f"{angle_count} {angle_count} {int(line.split()[0])} {int(line.split()[4])} {int(line.split()[8])}")
                angle_coeffs.append(f"angle_coeff {angle_count} {line.split()[12]} {line.split()[13]}")
    
    # Write all angles at once
    with open(data_file, "a") as flammps:
        flammps.write("\nAngles \n\n")
        flammps.write("\n".join(angle_lines) + "\n")
    
    # Write all angle coefficients at once
    with open(param_file, "a") as flammpsparm:
        flammpsparm.write("\n".join(angle_coeffs) + "\n")

    # Write dihedrals section
    if verbose:
        print("Writing dihedrals section...")
    
    dihedral_count = 0
    dihedral_lines = []
    dihedral_coeffs = []
    
    with open("dihedrals.txt") as f:
        for line in f:
            if "Atom" not in line and line.split():
                dihedral_count += 1
                if str(line.split()[0]) == 'M' or str(line.split()[0]) == 'I':
                    dihedral_lines.append(f"{dihedral_count} {dihedral_count} {int(line.split()[1])} {int(line.split()[5])} {int(line.split()[9])} {int(line.split()[13])}")
                    dihedral_coeffs.append(f"dihedral_coeff {dihedral_count} 1 {line.split()[17]} {int(float(line.split()[18]))} {line.split()[19]}")
                else:
                    dihedral_lines.append(f"{dihedral_count} {dihedral_count} {int(line.split()[0])} {int(line.split()[4])} {int(line.split()[8])} {int(line.split()[12])}")
                    dihedral_coeffs.append(f"dihedral_coeff {dihedral_count} 1 {line.split()[16]} {int(float(line.split()[17]))} {line.split()[18]}")
    
    # Write all dihedrals at once
    with open(data_file, "a") as flammps:
        flammps.write("\nDihedrals \n\n")
        flammps.write("\n".join(dihedral_lines) + "\n")
    
    # Write all dihedral coefficients at once
    with open(param_file, "a") as flammpsparm:
        flammpsparm.write("\n".join(dihedral_coeffs) + "\n")

    # Clean up temporary files
    if verbose:
        if keep_temp:
            print("Temporary files preserved:")
            temp_files = ["bonds.txt", "angles.txt", "dihedrals.txt", "pairs.txt"]
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    print(f"  - {temp_file}")
        else:
            print("Cleaning up temporary files...")
    
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
        print(f"  - {natoms} atoms")
        print(f"  - {bond_count} bonds") 
        print(f"  - {angle_count} angles")
        print(f"  - {dihedral_count} dihedrals")
    
    # Clean up temporary files at the end
    cleanup_temp_files(verbose, keep_temp)

def main():
    """Main function"""
    args = parse_arguments()
    
    # Validate input files
    validate_files(args.topology, args.crd)
    
    # Run conversion
    try:
        amber2lammps(
            data_file=args.data_file,
            param_file=args.param_file,
            topology=args.topology,
            crd=args.crd, 
            charge=args.charge,
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

    
    






