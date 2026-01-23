# AMBER2LAMMPS

A Python utility to convert AMBER molecular dynamics files to LAMMPS data format with enhanced command-line interface and error handling. Also mentioned on the [LAMMPS pre/post processing tools page](https://www.lammps.org/prepost.html).

## Citation

If you use this software in your research, please cite it as:

**DOI:** [10.5281/zenodo.18114886](https://doi.org/10.5281/zenodo.18114886)


### Key Features 

- Convert AMBER topology (.prmtop) and MOL2 coordinates to LAMMPS data format
- Extract force field parameters from AMBER frcmod files
- **Charge normalization** - Automatically normalizes atomic charges to ensure zero net charge (improvement over InterMol)
- **Separate data and parameter files** - Generates distinct LAMMPS data file and parameter file for better organization and flexibility (improvement over InterMol's mixed output)
- Command-line interface with comprehensive options
- Verbose output for debugging and monitoring
- Automatic file validation and error handling
- Input/output file management and error checking
- **Configurable output file naming** - Custom names for data and parameter files (improvement over fixed InterMol naming)
- Automatic cleanup of temporary files

## Command Reference

### CLI

```bash
python3 amber_to_lammps.py <data_file> <param_file> <topology> <mol2> <frcmod> [options]
```

| Argument | Required? | Description |
| --- | --- | --- |
| `data_file` | yes | Output LAMMPS data filename |
| `param_file` | yes | Output LAMMPS parameter filename |
| `topology` | yes | AMBER topology file (`.prmtop`) - filename or path to input topology |
| `mol2` | yes | MOL2 coordinate file (with charges) - filename or path to input coordinates |
| `frcmod` | yes | AMBER force field parameter file (`.frcmod`) - filename or path to input parameters |
| `-b, --buffer` | optional | Vacuum padding (Å) used to set the simulation box. Larger values add more empty space around the molecule. Default: `3.8`. |
| `--verbose` | optional | Print step-by-step progress, counts, and derived values (box size, charge shifts). |
| `-h, --help` | optional | Show help message. |

#### Examples

```bash
# Standard conversion
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod

# Custom buffer and verbose logging (adds 5 Å padding around molecule)
python3 amber_to_lammps.py my_data.lammps my_params.lammps epon.prmtop epon.mol2 epon.frcmod --verbose -b 5.0

# Custom output names
python3 amber_to_lammps.py system.data system.parm system.prmtop system.mol2 system.frcmod
```

### Python API

```python
from amber_to_lammps import amber2lammps, validate_files

amber2lammps(
    data_file='data.lammps',
    param_file='parm.lammps',
    topology='epon.prmtop',
    mol2='epon.mol2',
    frcmod='epon.frcmod',
    buffer=3.8,       # optional
    verbose=False,    # optional
)
```

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

### LAMMPS Installation

LAMMPS is required to run the generated input files. Install using one of these methods:
#### From source
Download from https://lammps.org/ and follow the build instructions. **Include these packages while compiling:** `MOLECULE KSPACE EXTRA-MOLECULE`

## Conversion Logic and Option Effects

Below is the sequence implemented in `amber_to_lammps.py` and how options influence it:

1. **Input validation**: `validate_files` checks that `topology`, `mol2`, and `frcmod` exist and are readable. Fail-fast prevents partial outputs and exits with code 1 if any file is missing or unreadable.
2. **AMBER topology load**: ParmEd reads atoms, bonds, angles, and dihedrals from `.prmtop`. Counts are echoed in `--verbose` mode.
3. **Atom typing and masses**: The script parses `MASS` entries in `frcmod` to map atom types to sequential IDs and masses. Missing types warn and fall back to type 1.
4. **Coordinates and charges**: Coordinates and per-atom charges are read from the MOL2 `ATOM` block.
5. **Box creation (`--buffer`)**: The simulation box is a bounding box around the coordinates expanded equally by the buffer. Increase `--buffer` if atoms are near the box edge; decrease for tighter, smaller boxes.
6. **Charge normalization**: Total charge is shifted uniformly across atoms so the system is neutral. This avoids LAMMPS warnings from slight charge drift.
7. **Nonbonded parameters**: `NONBON` terms in `frcmod` become `pair_coeff` entries in the parameter file.
8. **Topology terms**: Bonds/angles/dihedrals are exported via ParmEd to temporary files and then written to the LAMMPS data/parameter files with matching coefficients.
9. **Cleanup**: Temporary helper files (`bonds.txt`, `angles.txt`, `dihedrals.txt`) are removed; verbose mode reports all generated counts.
10. **Verbose diagnostics**: With `--verbose`, the script echoes counts, box extents (after buffer), the applied charge shift to achieve neutrality, and warnings if atom types/parameters are missing so you can intervene early.

## Guided Tutorial (end-to-end)

Use the included EPON example to see the full workflow and rationale.

### 1) Input Files  
   - `epon.prmtop`: topology/connectivity from AMBER  
   - `epon.mol2`: coordinates **with per-atom charges** from Antechamber  
   - `epon.frcmod`: GAFF/GAFF2 parameters not in the base force field  
   Together they contain all structure + charge + parameters needed for LAMMPS.

### 2) Run the conversion  
   ```bash
   python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod --verbose -b 4.5
   ```
   - `--verbose` shows the parsed counts, box extents, and charge shift so you can spot issues early.
   - `-b 4.5` pads the box by 4.5 Å in each direction; increase if you see atoms too close to periodic images.

### 3) Inspect the outputs  
   - `data.lammps`: contains atom coordinates, types, charges, and topology sections. Confirm atom/bond counts match `--verbose` output.  
   - `parm.lammps`: contains `pair_coeff`, `bond_coeff`, `angle_coeff`, and `dihedral_coeff`. Spot-check a few entries against your frcmod if parameters seem off.

### 4) Run LAMMPS  
   Use the provided `example_lammps_input.lmp`:
   ```bash
   lmp < example_lammps_input.lmp
   ```
   The input `include`s `parm.lammps` first, then `read_data` loads `data.lammps`. Keeping parameters separate makes it easier to swap parameter sets without regenerating coordinates.



## AMBERTools Setup

If you have AMBERTools installed, you can activate the environment (or install from https://ambermd.org/GetAmber.php#ambertools):

```bash
conda activate Ambertools23  # or your AMBERTools version
```

## Preparing AMBER Inputs (optional refresher)

### 1. AMBER Workflow (Generating Input Files)

If you need to generate AMBER input files from PDB, use this workflow:

```python
import subprocess

# Step 1: Generate MOL2 file with charges
cmd1 = "antechamber -j 4 -at gaff2 -dr no -fi pdb -fo mol2 -i epon.pdb -o epon.mol2 -c bcc"
subprocess.run(cmd1, shell=True)

# Step 2: Generate force field parameters (ensure -Y option is activated)
cmd2 = "parmchk2 -i epon.mol2 -o epon.frcmod -f mol2 -a Y"
subprocess.run(cmd2, shell=True)

# Step 3: Create tleap input file
with open("tleap.in", "w") as f:
    f.write("source leaprc.gaff2\n")
    f.write("SUS = loadmol2 epon.mol2\n")
    f.write("check SUS\n")
    f.write("loadamberparams epon.frcmod\n")
    f.write("saveamberparm SUS epon.prmtop epon.crd\n")
    f.write("quit")

# Step 4: Run tleap to generate AMBER files
cmd3 = "tleap -f tleap.in"
subprocess.run(cmd3, shell=True)

# Check log file for any errors
file_path = './leap.log'
```

**Prerequisites for AMBER Workflow:**
- AMBERTools must be installed and in PATH
- Input PDB file: `epon.pdb`
- Output files: `epon.prmtop`, `epon.crd`, `epon.mol2`, `epon.frcmod`

### Python Usage Patterns

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
    success = convert_system('epon.lammps', 'epon_parm.lammps', 'epon.prmtop', 'epon.mol2', 'epon.frcmod')
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

### Common Issues and Solutions When running lammps after input generation 

- **Bad charges or non-neutral system**: A MOL2 generated without charges (or with mixed charge methods) can lead to unexpected shifts when neutrality is enforced. Regenerate MOL2 with the same charge model (e.g., `antechamber -c bcc`) and confirm charges sum to the intended value.  
- **Box too small**: If atoms are interacting with periodic images, increase `--buffer`. For larger or elongated molecules, visualize the bounding box or print verbose logs to confirm min/max extents.  
- **Parameter gaps**: If LAMMPS errors on missing `pair_coeff`/`bond_coeff` entries, the required terms were not present in `frcmod`. Re-run `parmchk2` on the current MOL2 or add missing terms manually.  
- **Wrong file versions**: Changing the molecule but reusing old `prmtop`/`mol2`/`frcmod` leads to inconsistent topology vs coordinates. Regenerate the full AMBER set whenever chemistry changes.  
- **Energy mismatches vs reference**: Compare per-term energies (bond, angle, dihedral, nonbonded) against AMBER or InterMol. Large deviations usually trace back to parameter gaps or atom-type mismatches rather than box size.


## Validation

The AMBER2LAMMPS conversion has been validated using **InterMol** to ensure energy accuracy.

**GitHub Repository:** https://github.com/shirtsgroup/InterMol

**Usage Example:**
```bash
python convert.py --amb_in epon.prmtop epon.crd --lammps
```

**Generated Files:**
- `epon_converted.input` (LAMMPS input file)
- `epon_converted.lmp` (LAMMPS data file)

**Example LAMMPS Output:**
```text
E_bond        E_angle        E_dihed        E_impro         E_pair         E_vdwl         E_coul         E_long         E_tail         PotEng
2.3161274      6.0940384      12.475809      0             -8.8739005      10.824738      97.869973     -117.56861     -0.0044166818   12.012074
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
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod --verbose
```

This will show:
- File loading progress
- Number of atoms, bonds, angles, dihedrals found
- Box dimensions derived from coordinates and buffer
- Charge shift applied to neutralize total charge
- Any warnings about missing atom types or parameters

## Features Comparison

| Feature | AMBER2LAMMPS | InterMol |
|---------|-------------|----------|
| Configurable output | Custom names | Fixed names |
| Charge normalization | Yes | No |
| Separate data/parameter files | Yes | No |

## Contributing

When making modifications:
1. Create a new branch: `git checkout -b feature-name`
2. Make your changes
3. Test thoroughly with different input files
4. Commit and push your changes

## License

**Original Project:** The AMBER2LAMMPS project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


