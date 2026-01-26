# AMBER to LAMMPS Tutorial

A comprehensive Python utility to convert AMBER topology, forcefield, and coordinate files to LAMMPS data format. This tool provides both CLI and Python API interfaces with validation-oriented error handling.

## Table of Contents

- [Introduction](#introduction)
- [Project and Download](#project-and-download)
- [Requirements](#requirements)
- [Install Dependencies](#install-dependencies)
- [Command Reference](#command-reference)
- [Conversion Process](#conversion-process)
- [Workflow Examples](#workflow-examples)
- [Troubleshooting](#troubleshooting)
- [Validation](#validation)
- [Getting Help](#getting-help)
- [Citation](#citation)
- [License](#license)

## Introduction

This tool converts molecular simulations from AMBER format to LAMMPS format, enabling you to:
- Convert AMBER topology files (`.prmtop`) to LAMMPS data files
- Transform AMBER force field parameters (`.frcmod`) to LAMMPS parameter files
- Process coordinate files (`.mol2`) with atomic charges
- Generate ready-to-run LAMMPS simulations

**Typical workflow**: PDB structure → AMBER files → LAMMPS files → LAMMPS simulation

## Project and Download

### GitHub Clone

```bash
git clone https://github.com/your-repo/AMBER2LAMMPS.git
cd AMBER2LAMMPS
```

### External Tool Information

This utility integrates with several external tools:
- **AmberTools**: For generating AMBER topology and parameter files
- **ParmEd**: For reading AMBER topology structures
- **LAMMPS**: For running molecular dynamics simulations
- **Open Babel** (optional): For converting between molecular formats

## Requirements

### Platform Compatibility

AMBER2LAMMPS has been validated and tested on:

- **Linux** (Ubuntu, CentOS, Red Hat, Debian) - Fully tested and validated
- **macOS** (Intel and Apple Silicon) - Fully tested
- **Windows** (Windows 10/11 with WSL2 and native Python) - Tested with WSL2 and native Python
- **WSL2 is recommended for Windows users for best compatibility.**

### Structure Files
- **Input**: PDB file of your molecule (or SMILES string convertible to PDB)
- **Output**: AMBER files (`.prmtop`, `.mol2`, `.frcmod`) → LAMMPS files (`.data`, `.lammps`)

### AmberTools Utilities
- `tleap`: For generating AMBER topology files
- `antechamber`: For generating force field parameters
- `parmchk2`: For validating force field parameters

### Python Packages
- `parmed`: For reading AMBER topology structures
- `numpy`: For numerical operations and array handling

### LAMMPS Packages
- LAMMPS executable (`lmp` or `lammps`)
- Molecular package (for standard molecular dynamics)

### Optional Open Babel
- For SMILES to PDB conversion (alternative to other structure sources)

## Install Dependencies

### AmberTools

```bash
# Download AmberTools (free for academic use)
wget https://ambermd.org/cgi-bin/amber.py?action=download_ambertools

# Follow installation instructions
# Typical installation:
tar -xzf ambertools*.tar.gz
cd ambertools*/
./configure --prefix=/path/to/install/amber
make install

# Add to environment
export PATH=/path/to/install/amber/bin:$PATH
source /path/to/install/amber/amber.sh
```

### Python Packages

```bash
# Install required packages
pip install parmed numpy

# Or using conda
conda install -c conda-forge parmed numpy
```

### Open Babel (Optional)

```bash
# Ubuntu/Debian
sudo apt-get install openbabel

# macOS with Homebrew
brew install open-babel

# Conda
conda install -c conda-forge openbabel
```

## Command Reference

### CLI Arguments

| Argument | Required? | Description |
| --- | --- | --- |
| `data_file` | yes | Output LAMMPS data filename |
| `param_file` | yes | Output LAMMPS parameter filename |
| `topology` | yes | AMBER topology file (`.prmtop`) |
| `mol2` | yes | MOL2 coordinate file (with charges) |
| `frcmod` | yes | AMBER force field parameter file (`.frcmod`) |
| `-b, --buffer` | optional | Vacuum padding (Å) for the simulation box. Default: `3.8`. |
| `--verbose` | optional | Print step-by-step progress, counts, and box size. Default: `False`. |
| `-h, --help` | optional | Show help message. |

### CLI Usage

```bash
python3 amber_to_lammps.py <data_file> <param_file> <topology> <mol2> <frcmod> [options]
```

**Outputs**: a LAMMPS data file (`<data_file>`, e.g., `data.lammps`) and a separate parameter file
(`<param_file>`, e.g., `parm.lammps`).

**CLI Examples**:
```bash
# Basic conversion
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod

# With verbose output and custom buffer
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod --verbose -b 5.0

# Using absolute paths
python3 amber_to_lammps.py /home/user/lammps/output/data.lammps /home/user/lammps/output/param.lammps /home/user/amber/topology.prmtop /home/user/amber/coords.mol2 /home/user/amber/params.frcmod -b 4.5
```

### Python API Usage

```python
from amber_to_lammps import amber2lammps, validate_files

# Define files
data_file = 'data.lammps'
param_file = 'parm.lammps'
topology = 'epon.prmtop'
mol2 = 'epon.mol2'
frcmod = 'epon.frcmod'

# Optional validation
validate_files(topology, mol2, frcmod)

# Run conversion
amber2lammps(
    data_file=data_file,
    param_file=param_file,
    topology=topology,
    mol2=mol2,
    frcmod=frcmod,
    buffer=3.8,
    verbose=True,
)
```

### LAMMPS Input Script

Create a LAMMPS input file to use the generated files:

```lammps
# LAMMPS input script
units real
atom_style full

# Include force field parameters
include parm.lammps

# Read molecular structure
read_data data.lammps

# Define simulation box
boundary p p p

# Settings
pair_style lj/cut/coul/long 10.0
bond_style harmonic
angle_style harmonic
dihedral_style harmonic

# Run minimization
minimize 1.0e-4 1.0e-6 1000 10000

# Output
thermo 100
thermo_style custom step temp pe ke etotal press density
```

Run LAMMPS:
```bash
lmp < input.lammps
```

## Conversion Process

### Input Validation Behavior
- **CLI**: Automatically validates that input files (`topology`, `mol2`, `frcmod`) exist and are readable before conversion
- **Python API**: Validation is optional via `validate_files()` function

### AMBER Topology Loading
1. **ParmEd integration**: Uses `pmd.amber.AmberParm` to load `.prmtop` files
2. **Structure extraction**: Reads atoms, bonds, angles, and dihedrals from topology
3. **Force field parsing**: Extracts bonded and nonbonded parameters

### Atom Typing and Mass Mapping
1. **MASS entries**: Parses `MASS` section from `.frcmod` file
2. **Type mapping**: Maps AMBER atom types to sequential LAMMPS type IDs
3. **Fallback handling**: Missing atom types default to type 1 with warning

### Coordinate and Charge Processing
1. **MOL2 parsing**: Extracts coordinates and charges from `ATOM` section
2. **Charge validation**: Checks total system charge
3. **Normalization**: Applies uniform charge shift if `|net charge| > 1e-6`

### Box Creation Algorithm
1. **Bounding box**: Calculates min/max coordinates from molecular structure
2. **Buffer addition**: Adds specified padding (default: 3.8 Å) around molecule
3. **Output format**: Writes LAMMPS box dimensions to data file

### Charge Normalization
- **Threshold**: `1e-6` elementary charge
- **Method**: Uniform shift across all atoms to achieve neutrality
- **Verification**: Reports final net charge after normalization

### Nonbonded Parameter Mapping
1. **NONBON parsing**: Extracts LJ parameters from `.frcmod` file
2. **Conversion**: Transforms AMBER epsilon/sigma to LAMMPS format
3. **Pair coefficients**: Writes `pair_coeff` entries to parameter file

### Topology Term Export
1. **Bonds**: Converts AMBER bond definitions to LAMMPS format
2. **Angles**: Maps AMBER angle parameters to LAMMPS coefficients
3. **Dihedrals**: Handles proper and improper dihedrals with appropriate coefficients

### File Cleanup
- **Temporary files**: Removes intermediate files (`bonds.txt`, `angles.txt`, `dihedrals.txt`)
- **Error handling**: Cleans up even if conversion fails
- **Verbose output**: Reports cleanup status when enabled

### Verbose Diagnostics
When `--verbose` is enabled, reports:
- Atom, bond, angle, and dihedral counts
- Box dimensions with buffer
- Charge normalization details
- File processing status
- Temporary file cleanup

## Workflow Examples

### Prepare Input Files (Ethanol Example)

#### Step 1: Get Structure (Optional SMILES to PDB)

```bash
# Using Open Babel
obabel -:"CCO" -opdb -O ethanol.pdb --gen3d

# Or use existing PDB file
cp your_ethanol.pdb ethanol.pdb
```

#### Step 2: Generate AMBER Topology and Parameters

```bash
# Create tleap input file
cat > tleap.in << EOF
source leaprc.gaff2
ethanol = loadmol2 ethanol.mol2
check ethanol
saveamberparm ethanol ethanol.prmtop ethanol.inpcrd
saveoff ethanol ethanol.lib
quit
EOF

# Generate MOL2 with charges
antechamber -i ethanol.pdb -fi pdb -o ethanol.mol2 -fo mol2 -c bcc -s 2

# Generate force field parameters
parmchk2 -i ethanol.mol2 -f mol2 -o ethanol.frcmod

# Run tleap to create topology
tleap -f tleap.in
```

### Basic Conversion Workflow

#### CLI Usage with LAMMPS Execution

```bash
# Convert AMBER to LAMMPS
python3 amber_to_lammps.py ethanol.data ethanol.parm ethanol.prmtop ethanol.mol2 ethanol.frcmod --verbose

# Create LAMMPS input script
cat > ethanol.lmp << EOF
units real
atom_style full
boundary p p p

include ethanol.parm
read_data ethanol.data

pair_style lj/cut/coul/long 10.0
bond_style harmonic
angle_style harmonic
dihedral_style harmonic

minimize 1.0e-4 1.0e-6 1000 10000

thermo 100
thermo_style custom step temp pe ke etotal press density
EOF

# Run LAMMPS
lmp < ethanol.lmp
```

#### Python API Usage with LAMMPS Execution

```python
from amber_to_lammps import amber2lammps, validate_files
import subprocess
import os

# Define files
data_file = 'ethanol.data'
param_file = 'ethanol.parm'
topology = 'ethanol.prmtop'
mol2 = 'ethanol.mol2'
frcmod = 'ethanol.frcmod'

# Validate inputs
validate_files(topology, mol2, frcmod)

# Convert to LAMMPS
amber2lammps(
    data_file=data_file,
    param_file=param_file,
    topology=topology,
    mol2=mol2,
    frcmod=frcmod,
    buffer=3.8,
    verbose=True,
)

# Create LAMMPS input script
lammps_input = """
units real
atom_style full
boundary p p p

include ethanol.parm
read_data ethanol.data

pair_style lj/cut/coul/long 10.0
bond_style harmonic
angle_style harmonic
dihedral_style harmonic

minimize 1.0e-4 1.0e-6 1000 10000

thermo 100
thermo_style custom step temp pe ke etotal press density
"""

with open('ethanol.lmp', 'w') as f:
    f.write(lammps_input)

# Run LAMMPS
subprocess.run(['lmp', '<', 'ethanol.lmp'], shell=True)
```

### Batch Processing Multiple Molecules

```python
from amber_to_lammps import amber2lammps, validate_files

molecules = ['ethanol', 'benzene', 'aspirin']

for mol in molecules:
    print(f"Processing {mol}...")
    
    # Validate files
    validate_files(f'{mol}.prmtop', f'{mol}.mol2', f'{mol}.frcmod')
    
    # Convert
    amber2lammps(
        data_file=f'{mol}.data',
        param_file=f'{mol}.parm',
        topology=f'{mol}.prmtop',
        mol2=f'{mol}.mol2',
        frcmod=f'{mol}.frcmod',
        verbose=True
    )
    
    print(f"✓ {mol} conversion complete")
```

## Troubleshooting

### Input File Issues

**Problem**: `Error: topology file not found`
- **Solution**: Check file paths and ensure files exist in specified locations
- **CLI**: Automatically validates file existence before conversion

**Problem**: `Error: Cannot read topology file`
- **Solution**: Check file permissions and ensure files are not corrupted
- **Test**: `file your.prmtop` should show readable file type

**Problem**: Missing atom types in frcmod
- **Solution**: Ensure all atom types in MOL2 file are defined in frcmod MASS section
- **Fallback**: Script defaults missing types to type 1 with warning

### Parameter Problems

**Problem**: `Warning: Atom type X not found in frcmod`
- **Solution**: Add missing atom type to frcmod MASS section
- **Format**: `MASS\n X 12.01  # Carbon`

**Problem**: Charge normalization issues
- **Solution**: Check MOL2 file for correct charge values
- **Verification**: Net charge should be reasonable (< 10e)

### LAMMPS Integration

**Problem**: LAMMPS cannot read data file
- **Solution**: Ensure LAMMPS version supports molecular atom style
- **Check**: `atom_style full` in LAMMPS input

**Problem**: Parameter file not found
- **Solution**: Use correct path in `include` statement
- **Example**: `include ../params/ethanol.parm`

### Performance Optimization

**Large molecules** (>10,000 atoms):
- Increase buffer size: `-b 5.0`
- Use verbose mode to monitor progress: `--verbose`

**Memory issues**:
- Close other applications
- Use machine with more RAM
- Consider splitting large systems

## Validation

The AMBER2LAMMPS conversion has been validated using **InterMol** to ensure energy accuracy.

### InterMol Validation

**Installation**: https://github.com/shirtsgroup/InterMol

**Usage Example**:
```bash
python convert.py --amb_in epon.prmtop epon.crd --lammps
```

**Generated Files**:
- `epon_converted.input` (LAMMPS input file)
- `epon_converted.lmp` (LAMMPS data file)

## Getting Help

### GitHub Issues
- Report bugs: https://github.com/your-repo/AMBER2LAMMPS/issues
- Feature requests: https://github.com/your-repo/AMBER2LAMMPS/issues/new

### Community Support
- Discussions: https://github.com/your-repo/AMBER2LAMMPS/discussions
- Email: your-email@example.com

## Citation

If you use this software in your research, please cite it as:

**DOI**: [10.5281/zenodo.18114886](https://doi.org/10.5281/zenodo.18114886)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

**Open Source**: This is an open-source project. The source code is freely available for use, modification, and distribution under the MIT License.
