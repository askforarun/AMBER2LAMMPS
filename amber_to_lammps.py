
#!/usr/bin/env python3
"""
AMBER to LAMMPS conversion script

This script converts AMBER topology and coordinate files to LAMMPS data format.
It requires parmed, numpy, and the AMBER force field files.

Usage:
    python amber_to_lammps.py <data_file> <param_file> <topology> <mol2> <frcmod>

Arguments:
    data_file    Output LAMMPS data file name
    param_file   Output LAMMPS parameter file name  
    topology     AMBER topology file (.prmtop)
    mol2         MOL2 coordinate file
    frcmod       Force field parameter file (.frcmod)
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
    python amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod
    python amber_to_lammps.py system.data system.parm system.prmtop system.mol2 system.frcmod
        """
    )
    
    parser.add_argument('data_file', help='Output LAMMPS data file name')
    parser.add_argument('param_file', help='Output LAMMPS parameter file name')
    parser.add_argument('topology', help='AMBER topology file (.prmtop)')
    parser.add_argument('mol2', help='MOL2 coordinate file')
    parser.add_argument('frcmod', help='Force field parameter file (.frcmod)')
    parser.add_argument('-b', '--buffer', type=float, default=3.8,
                        help='Buffer size around molecule (default: 3.8)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    
    return parser.parse_args()

def validate_files(topology, mol2, frcmod):
    """Validate input files exist and are readable"""
    for filepath, filetype in [(topology, 'topology'), (mol2, 'MOL2'), (frcmod, 'frcmod')]:
        if not os.path.exists(filepath):
            print(f"Error: {filetype} file '{filepath}' not found")
            sys.exit(1)
        if not os.access(filepath, os.R_OK):
            print(f"Error: Cannot read {filetype} file '{filepath}'")
            sys.exit(1)

def cleanup_temp_files(verbose=False):
    """Remove temporary files if they exist"""
    temp_files = ['bonds.txt', 'angles.txt', 'dihedrals.txt']
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
                if verbose:
                    print(f"Removed existing temporary file: {temp_file}")
            except OSError as e:
                print(f"Warning: Could not remove {temp_file}: {e}")

def amber2lammps(data_file, param_file, topology, mol2, frcmod, buffer=3.8, verbose=False):
    AmberParm=pmd.amber.AmberParm
    printBonds=pmd.tools.actions.printBonds
    printAngles=pmd.tools.actions.printAngles
    printDihedrals=pmd.tools.actions.printDihedrals

    # Clean up temporary files if they exist
    cleanup_temp_files(verbose)

    # Setup output files
    if verbose:
        print(f"Converting AMBER files to LAMMPS format...")
        print(f"Output files: {data_file}, {param_file}")
    
    # Initialize output files
    with open(data_file, "w") as f:
        f.write(f"LAMMPS data file from AMBER conversion\n")
        f.write(f"#Source: {topology}, {mol2}, {frcmod}\n\n")

    with open(param_file, "w") as f:
        f.write(f"# Force field parameters generated from {frcmod}\n\n")

    
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

    # Parse atom types and masses from frcmod
    if verbose:
        print(f"Parsing atom types from {frcmod}...")
    
    atype = {}
    mass = []
    with open(frcmod) as f:
        for line in f:
            if "MASS" in line:
                for i, line_1 in enumerate(f):
                    if "BOND" in line_1:
                        break
                    if line_1.split():
                        atype[line_1.split()[0]] = i + 1
                        mass.append(line_1.split()[1])
                break
    
    if verbose:
        print(f"Found {len(atype)} atom types")

    # Parse coordinates and charges from mol2
    if verbose:
        print(f"Parsing coordinates from {mol2}...")
    
    x, y, z = [], [], []
    charges = []
    with open(mol2) as f: 
        for line in f:
            if "ATOM" in line:
                for i, line_1 in enumerate(f):
                    if "BOND" in line_1:
                        break
                    if line_1.split():
                        charges.append(float(line_1.split()[8]))
                        x.append(float(line_1.split()[2]))
                        y.append(float(line_1.split()[3]))
                        z.append(float(line_1.split()[4]))
                break

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

        f.write("\nAtoms\n\n")
                
    # Normalize charges to ensure neutrality
    if verbose:
        print(f"Normalizing charges (total charge: {np.sum(charges):.6f})...")
    
    if np.sum(charges) >= 0:
        charges = charges - (abs(np.sum(charges))/len(charges))
    elif np.sum(charges) < 0:
        charges = charges + (abs(np.sum(charges))/len(charges))
    
    if verbose:
        print(f"Normalized total charge: {np.sum(charges):.6f}")

    # Parse nonbonded parameters from frcmod
    if verbose:
        print("Parsing nonbonded parameters...")
    
    with open(frcmod) as f:
        for line in f:
            if "NONBON" in line:
                for i, line_1 in enumerate(f):
                    if line_1.split():
                        with open(param_file, "a") as flammps:
                            flammps.write(f"pair_coeff {i+1} {i+1} {line_1.split()[2]} {0.890898718140339*2*float(line_1.split()[1])} # {line_1.split()[0]}\n")

    # Write atoms section
    if verbose:
        print("Writing atoms section...")
    
    with open(mol2) as f:  
        for line in f:
            if "ATOM" in line:
                with open(data_file, "a") as flammps:
                    for i, line_1 in enumerate(f):
                        if "BOND" in line_1:
                            break
                        if line_1.split():
                            atom_type = atype.get(str(line_1.split()[5]))
                            if atom_type is None:
                                print(f"Warning: Atom type {line_1.split()[5]} not found in frcmod")
                                atom_type = 1
                            flammps.write(f"{i+1} 1 {atom_type} {charges[i]} {line_1.split()[2]} {line_1.split()[3]} {line_1.split()[4]} 0 0 0\n")
                break

    # Write bonds section
    if verbose:
        print("Writing bonds section...")
    
    bond_count = 0
    with open(data_file, "a") as flammps:
        flammps.write("\nBonds \n\n")
        with open("bonds.txt") as f:
            for line in f:
                if "Atom" not in line and line.split():
                    bond_count += 1
                    flammps.write(f"{bond_count} {bond_count} {int(line.split()[0])} {int(line.split()[4])}\n")
                    with open(param_file, "a") as flammpsparm:
                        flammpsparm.write(f"bond_coeff {bond_count} {line.split()[9]} {line.split()[8]}\n")

    # Write angles section
    if verbose:
        print("Writing angles section...")
    
    angle_count = 0
    with open(data_file, "a") as flammps:
        flammps.write("\nAngles \n\n")
        with open("angles.txt") as f:
            for line in f:
                if "Atom" not in line and line.split():
                    angle_count += 1
                    flammps.write(f"{angle_count} {angle_count} {int(line.split()[0])} {int(line.split()[4])} {int(line.split()[8])}\n")
                    with open(param_file, "a") as flammpsparm:
                        flammpsparm.write(f"angle_coeff {angle_count} {line.split()[12]} {line.split()[13]}\n")

    # Write dihedrals section
    if verbose:
        print("Writing dihedrals section...")
    
    dihedral_count = 0
    with open(data_file, "a") as flammps:
        flammps.write("\nDihedrals \n\n")
        with open("dihedrals.txt") as f:
            for line in f:
                if "Atom" not in line and line.split():
                    dihedral_count += 1
                    if str(line.split()[0]) == 'M' or str(line.split()[0]) == 'I':
                        flammps.write(f"{dihedral_count} {dihedral_count} {int(line.split()[1])} {int(line.split()[5])} {int(line.split()[9])} {int(line.split()[13])}\n")
                        with open(param_file, "a") as flammpsparm:
                            flammpsparm.write(f"dihedral_coeff {dihedral_count} 1 {line.split()[17]} {int(float(line.split()[18]))} {line.split()[19]}\n")
                    else:
                        flammps.write(f"{dihedral_count} {dihedral_count} {int(line.split()[0])} {int(line.split()[4])} {int(line.split()[8])} {int(line.split()[12])}\n")
                        with open(param_file, "a") as flammpsparm:
                            flammpsparm.write(f"dihedral_coeff {dihedral_count} 1 {line.split()[16]} {int(float(line.split()[17]))} {line.split()[18]}\n")

    # Clean up temporary files
    if verbose:
        print("Cleaning up temporary files...")
    
    temp_files = ["bonds.txt", "angles.txt", "dihedrals.txt"]
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
    cleanup_temp_files(verbose)

def main():
    """Main function"""
    args = parse_arguments()
    
    # Validate input files
    validate_files(args.topology, args.mol2, args.frcmod)
    
    # Run conversion
    try:
        amber2lammps(
            data_file=args.data_file,
            param_file=args.param_file,
            topology=args.topology,
            mol2=args.mol2, 
            frcmod=args.frcmod,
            buffer=args.buffer,
            verbose=args.verbose
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

    
    







