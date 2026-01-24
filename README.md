# Overview

A Python utility to convert AMBER topology, forcefield, and coordinate files to LAMMPS data format. Features enhanced command-line interface, Python API, and comprehensive error handling.

## Open Source

This is an open-source project. The source code is freely available for use, modification, and distribution under the MIT License.


## Citation

If you use this software in your research, please cite it as:

**DOI:** [10.5281/zenodo.18114886](https://doi.org/10.5281/zenodo.18114886)

## What You Need

- **Structure**: A PDB file of your molecule (or a SMILES string you can convert to PDB; see SMILES to PDB workflow below).
- **AMBER prep tools**: `.prmtop`, `.mol2`, `.frcmod` files. AmberTools (`antechamber`, `parmchk2`, `tleap`) is needed - See instructions below on how to install and generate these files.
- **Python packages**: `parmed` and `numpy` (installation instructions below).
- **LAMMPS**: Installed and available on your `PATH` (`which lmp` or `lmp -help` to confirm) - installation instructions below.

### SMILES to PDB Workflow

You can use `obabel` to generate a PDB file from SMILES.

#### Installation

```bash
# macOS
brew install open-babel

# conda
conda install -c openbabel openbabel

# Ubuntu/Debian
sudo apt-get install openbabel
```

#### Generate PDB from SMILES

```bash
# Basic conversion
obabel -:CCO -opdb -O ethanol.pdb --gen3d

# With explicit hydrogens (recommended)
obabel -:CCO -h -opdb -O ethanol.pdb --gen3d

# Other examples
obabel -:c1ccccc1 -opdb -O benzene.pdb --gen3d  # Benzene
obabel -:"CC(=O)OC1=CC=CC=C1C(=O)O" -opdb -O aspirin.pdb --gen3d  # Aspirin
```


### AmberTools Installation

Install AmberTools from https://ambermd.org/GetAmber.php#ambertools and activate the environment:

```bash
conda activate Ambertools23  # or your AMBERTools version

```

### Generate AMBER input files (`prmtop`, `mol2`, `frcmod` files)

Use the following python code to generate AMBER input files:

```python
import subprocess

# Set your PDB filename here
pdb_file = "epon.pdb"  # Replace with your PDB filename
base_name = pdb_file.replace(".pdb", "")  # Extract base name from PDB filename

# Generate MOL2 file with charges from PDB file
cmd1 = f"antechamber -j 4 -at gaff2 -dr yes -fi pdb -fo mol2 -i {pdb_file} -o {base_name}.mol2 -c bcc"
subprocess.run(cmd1, shell=True)
# For details on the antechamber command, see AMBER Manual https://ambermd.org/Manuals.php

# Generate force field parameters (ensure -Y option is activated)
cmd2 = f"parmchk2 -i {base_name}.mol2 -o {base_name}.frcmod -f mol2 -a Y"
subprocess.run(cmd2, shell=True)

# Create tleap input file
with open("tleap.in", "w") as f:
    f.write("source leaprc.gaff2\n")
    f.write(f"SUS = loadmol2 {base_name}.mol2\n")
    f.write("check SUS\n")
    f.write(f"loadamberparams {base_name}.frcmod\n")
    f.write(f"saveamberparm SUS {base_name}.prmtop {base_name}.crd\n")
    f.write("quit")

# Run tleap to generate AMBER files
cmd3 = "tleap -f tleap.in"
subprocess.run(cmd3, shell=True)

# Check log file for any errors
file_path = './leap.log'

# Output files generated: epon.prmtop, epon.crd, epon.mol2, epon.frcmod
```

### Python Package Installation

### Prerequisites

- Python 3.6 or higher
- conda (recommended) or pip

### Required Dependencies

**Using conda (recommended)**

```bash
conda install -c conda-forge parmed numpy
```

**Using pip**

```bash
pip install parmed numpy
```

### LAMMPS Installation

#### Install from source

Download from https://lammps.org/ and follow the build instructions. **Include these packages while compiling:** `MOLECULE KSPACE EXTRA-MOLECULE` (otherwise the generated input files will not run).


## Command Reference

### CLI

```bash
python3 amber_to_lammps.py <data_file> <param_file> <topology> <mol2> <frcmod> [options]
```

Outputs: a LAMMPS data file (`<data_file>`, e.g., `data.lammps`) and a separate parameter file (`<param_file>`, e.g., `parm.lammps`).

| Argument | Required? | Description |
| --- | --- | --- |
| `data_file` | yes | Output LAMMPS data filename |
| `param_file` | yes | Output LAMMPS parameter filename |
| `topology` | yes | AMBER topology file (`.prmtop`) - filename or path to input topology |
| `mol2` | yes | MOL2 coordinate file (with charges) - filename or path to input coordinates |
| `frcmod` | yes | AMBER force field parameter file (`.frcmod`) - filename or path to input parameters |
| `-b, --buffer` | optional | Vacuum padding (Å) used to set the simulation box. Larger values add more empty space around the molecule. Default: `3.8`. |
| `--verbose` | optional | Print step-by-step progress, counts, and box size. |
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




## Conversion Logic and Options

Below is the sequence implemented in `amber_to_lammps.py` and how various options influence the output:

1. **Input validation**: `validate_files` checks that `topology`, `mol2`, and `frcmod` exist and are readable. Fail-fast prevents partial outputs and exits with code 1 if any file is missing or unreadable.
2. **AMBER topology load**: ParmEd reads atoms, bonds, angles, and dihedrals from `.prmtop`. Counts are echoed in `--verbose` mode.
3. **Atom typing and masses**: The script parses `MASS` entries in `frcmod` to map atom types to sequential IDs and masses. Missing types warn and fall back to type 1.
4. **Coordinates and charges**: Coordinates and per-atom charges are read from the MOL2 `ATOM` block.
5. **Box creation (`--buffer`)**: The simulation box is a bounding box around the coordinates expanded equally by the buffer. Increase `--buffer` if atoms are near the box edge; decrease for tighter, smaller boxes.
6. **Charge normalization**: Total charge is shifted uniformly across atoms so the system is neutral, but only if |net charge| > 1e-6. This introduces a small (typically negligible) offset to whatever charge method you used (e.g., RESP, AM1-BCC); if the net charge is within tolerance, no shift is applied.
7. **Nonbonded parameters**: `NONBON` terms in `frcmod` become `pair_coeff` entries in the parameter file.
8. **Topology terms**: Bonds/angles/dihedrals are exported via ParmEd to temporary files and then written to the LAMMPS data/parameter files with matching coefficients.
9. **Cleanup**: Temporary helper files (`bonds.txt`, `angles.txt`, `dihedrals.txt`) are removed; verbose mode reports all generated counts.
10. **Verbose diagnostics**: With `--verbose`, the script echoes counts, box extents (after buffer), and the normalized total charge; missing atom types are reported if encountered.

## Tutorial (End-to-End Workflow)



After generating the PDB file, you can follow steps 1-3 of the tutorial below.

---

### Usage With CLI

#### 1. AMBER Workflow (Generating Input Files)


#### 2. Run the Conversion
```bash
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod --verbose -b 4.5
```

#### 3. Run LAMMPS
Use the provided `example_lammps_input.lmp`:

```bash
lmp < example_lammps_input.lmp
```

**Note**: The LAMMPS input file `example_lammps_input.lmp` `include`s `parm.lammps` and the `read_data` loads `data.lammps`. Keeping parameters separate makes it easier to swap parameter sets without regenerating coordinates.

### Usage with Python API

```python
from amber_to_lammps import amber2lammps, validate_files

# Convert AMBER system
data_file = 'data.lammps'
param_file = 'parm.lammps'
topology = 'epon.prmtop'
mol2 = 'epon.mol2'
frcmod = 'epon.frcmod'

# Validate input files
validate_files(topology, mol2, frcmod)

# Run conversion
amber2lammps(
    data_file=data_file,
    param_file=param_file,
    topology=topology,
    mol2=mol2,
    frcmod=frcmod,
    verbose=True
)
```

## Validation

The AMBER2LAMMPS conversion has been validated using **InterMol** to ensure energy accuracy.

You can install **InterMol** from https://github.com/shirtsgroup/InterMol

**Usage Example (InterMol):**
```bash
python convert.py --amb_in epon.prmtop epon.crd --lammps
```

**Generated Files:**
- `epon_converted.input` (LAMMPS input file)
- `epon_converted.lmp` (LAMMPS data file)

```bash
lmp < epon_converted.input
```

**LAMMPS Output:**
```text
E_bond        E_angle        E_dihed        E_impro         E_pair         E_vdwl         E_coul         E_long         E_tail         PotEng
2.3161274      6.0940384      12.475809      0             -8.8739005      10.824738      97.869973     -117.56861     -0.0044166818   12.012074
```
You will get a WARNING that the system is not charge neutral like the following

WARNING: System is not charge neutral, net charge = -0.002996

This is because charges generated from antechamber AM1-BCC method carry slight excess charge and the maginitude of the charge is greater than 1e-6

**Output From amber2lammps Conversion (Step 3 of Tutorial)**

```
lmp < example_lammps_input.lmp
```
**LAMMPS Output:** 

```text
E_bond        E_angle        E_dihed        E_impro         E_pair         E_vdwl         E_coul         E_long         E_tail         PotEng
2.3161274      6.0940126      12.475827      0             -9.1202197      10.511316      108.01561     -127.64715     -0.31768789     11.765747
```
No WARNING messages. Charge normalization is applied to make the system charge neutral.

## Contributing

When making modifications:
1. Create a new branch: `git checkout -b feature-name`
2. Make your changes
3. Test thoroughly with different input files
4. Commit and push your changes
