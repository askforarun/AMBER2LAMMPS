# AMBER2LAMMPS - Enhanced Script

A Python utility to convert AMBER molecular dynamics files to LAMMPS data format with enhanced command-line interface and error handling.

## Features

- Convert AMBER topology (.prmtop) and MOL2 coordinates to LAMMPS data format
- Extract force field parameters from AMBER frcmod files
- Command-line interface with comprehensive options
- Verbose output for debugging and monitoring
- Automatic file validation and error handling
- Configurable output file naming
- Automatic cleanup of temporary files

## Installation

### Prerequisites

- Python 3.6 or higher
- conda (recommended) or pip

### Required Dependencies

```bash
# Using conda (recommended)
conda install numpy
pip install parmed

# Or using pip only
pip install numpy parmed
```

### AMBERTools (Optional)

If you have AMBERTools installed, you can activate the environment:

```bash
conda activate Ambertools23  # or your AMBERTools version
```

## Usage

### 1. Command Line Interface

#### Basic Usage

```bash
python3 amber_to_lammps.py <data_file> <param_file> <topology> <mol2> <frcmod>
```

#### Arguments

**Positional Arguments:**
- `data_file`: Output LAMMPS data file name
- `param_file`: Output LAMMPS parameter file name  
- `topology`: AMBER topology file (.prmtop)
- `mol2`: MOL2 coordinate file  
- `frcmod`: Force field parameter file (.frcmod)

**Optional Arguments:**
- `-b, --buffer`: Buffer size around molecule in Å (default: 3.8)
- `--verbose`: Enable verbose output
- `-h, --help`: Show help message

#### Examples

```bash
# Standard conversion
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod

# Custom output filenames
python3 amber_to_lammps.py system.data system.parm system.prmtop system.mol2 system.frcmod

# Verbose mode with custom buffer
python3 amber_to_lammps.py my_data.lammps my_params.lammps epon.prmtop epon.mol2 epon.frcmod --verbose -b 5.0

# Get help
python3 amber_to_lammps.py --help
```

#### Output Files

The script generates two files with your specified names:
- `{data_file}` - LAMMPS data file
- `{param_file}` - LAMMPS parameters file

### 2. Python Function Usage

You can also import and use the conversion function directly in your Python scripts.

#### Basic Function Import

```python
from amber_to_lammps import amber2lammps

# Basic conversion
amber2lammps(
    data_file='data.lammps',
    param_file='parm.lammps',
    topology='epon.prmtop',
    mol2='epon.mol2', 
    frcmod='epon.frcmod'
)

# With custom options
amber2lammps(
    data_file='my_system.lammps',
    param_file='my_system_parm.lammps',
    topology='system.prmtop',
    mol2='system.mol2',
    frcmod='epon.frcmod',
    buffer=5.0,
    verbose=True
)
```

#### Advanced Usage with Error Handling

```python
from amber_to_lammps import amber2lammps, validate_files
import sys

def convert_system(data_file, param_file, topology, mol2, frcmod):
    """Convert AMBER system with error handling"""
    
    # Validate input files
    validate_files(topology, mol2, frcmod)
    
    try:
        amber2lammps(
            data_file=data_file,
            param_file=param_file,
            topology=topology,
            mol2=mol2,
            frcmod=frcmod,
            verbose=True
        )
        print(f"✓ Conversion completed: {data_file}")
        return True
    except Exception as e:
        print(f"✗ Conversion failed: {e}")
        return False

# Usage
if __name__ == "__main__":
    success = convert_system('epon.lammps', 'epon_parm.lammps', 
                           'epon.prmtop', 'epon.mol2', 'epon.frcmod')
    sys.exit(0 if success else 1)
```

#### Simple Usage Example

```python
#!/usr/bin/env python3
from amber_to_lammps import amber2lammps, validate_files
import sys

def convert_single_system():
    """Convert a single AMBER system"""
    
    # Input files
    data_file = 'system.lammps'
    param_file = 'system_parm.lammps'
    topology = 'system.prmtop'
    mol2 = 'system.mol2'
    frcmod = 'epon.frcmod'
    
    # Validate input files
    validate_files(topology, mol2, frcmod)
    
    try:
        amber2lammps(
            data_file=data_file,
            param_file=param_file,
            topology=topology,
            mol2=mol2,
            frcmod=frcmod,
            verbose=True
        )
        print(f"✓ Conversion completed: {data_file}")
        return True
    except Exception as e:
        print(f"✗ Conversion failed: {e}")
        return False

if __name__ == "__main__":
    success = convert_single_system()
    sys.exit(0 if success else 1)
```

## Function Parameters

### `amber2lammps()` Function

```python
amber2lammps(data_file, param_file, topology, mol2, frcmod, buffer=3.8, verbose=False)
```

**Parameters:**
- `data_file` (str): Output LAMMPS data file name
- `param_file` (str): Output LAMMPS parameter file name
- `topology` (str): Path to AMBER topology file (.prmtop)
- `mol2` (str): Path to MOL2 coordinate file
- `frcmod` (str): Path to force field parameter file (.frcmod)
- `buffer` (float, optional): Buffer size around molecule in Å (default: 3.8)
- `verbose` (bool, optional): Enable verbose output (default: False)

### `validate_files()` Function

```python
validate_files(topology, mol2, frcmod)
```

**Parameters:**
- `topology` (str): Path to AMBER topology file
- `mol2` (str): Path to MOL2 file
- `frcmod` (str): Path to frcmod file

**Behavior:** Exits with error code 1 if any file is missing or unreadable.

## Example Workflow

### Complete Example with Test Data

```bash
# 1. Install dependencies
pip install numpy parmed

# 2. Run conversion (specifying output filenames)
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod

# 3. Check output
ls -la data.lammps parm.lammps
# Output: data.lammps  parm.lammps
```

### Using in LAMMPS

After conversion, you can use the generated files in your LAMMPS input:

```lammps
# Include the parameters
include "parm.lammps"

# Read the data file
read_data "data.lammps"

# Your LAMMPS simulation commands here
```

## Troubleshooting

### Common Issues

1. **ModuleNotFoundError: No module named 'parmed'**
   ```bash
   pip install parmed
   ```

2. **File not found errors**
   - Check that all input files exist and are readable
   - Use absolute paths if files are in different directories

3. **Atom type not found warnings**
   - Ensure your frcmod file contains all atom types from your MOL2 file
   - Check that atom type names match exactly

### Verbose Mode

Use `--verbose` to see detailed progress:
```bash
python3 amber_to_lammps.py -t epon.prmtop -m epon.mol2 -f epon.frcmod --verbose
```

This will show:
- File loading progress
- Number of atoms, bonds, angles, dihedrals found
- Box dimensions
- Charge normalization
- Section-by-section writing progress

## License

This enhanced script maintains the same license as the original AMBER2LAMMPS project.

## Contributing

When making modifications:
1. Create a new branch: `git checkout -b feature-name`
2. Make your changes
3. Test thoroughly with different input files
4. Commit and push your changes

## Original vs Enhanced

| Feature | Original | Enhanced |
|---------|----------|----------|
| Command-line interface | No | Yes (argparse) |
| Error handling | Basic | Comprehensive |
| File validation | No | Yes |
| Verbose output | No | Yes |
| Configurable output | Fixed names | Custom prefixes |
| Temporary file cleanup | No | Automatic |
| Documentation | Minimal | Comprehensive |
