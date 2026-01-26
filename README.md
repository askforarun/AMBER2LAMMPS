# Overview

A Python utility to convert AMBER topology, forcefield, and coordinate files to LAMMPS data format.
Provides an enhanced CLI, a Python API, and validation-oriented error handling.

## What This Tool Does

This tool helps you run molecular dynamics simulations in **LAMMPS** when you have your molecular system set up in **AMBER** format.

**Typical workflow:**
1. Start with a molecular structure (PDB file or SMILES string)
2. Use AMBER tools to create AMBER files (`.prmtop`, `.mol2`, `.frcmod`)
3. **Convert AMBER files to LAMMPS format** using this tool
4. Run your simulation in LAMMPS

**What gets converted:**
- **`.prmtop`** (topology file): Contains bonds, angles, atom types → LAMMPS data file
- **`.mol2`** (coordinates file): Contains atomic positions and charges → LAMMPS coordinates
- **`.frcmod`** (force field file): Contains interaction parameters → LAMMPS parameters

**What you get as output:**
- **LAMMPS data file** (e.g., `data.lammps`): Contains atomic coordinates, box dimensions, and molecular topology
- **LAMMPS parameter file** (e.g., `parm.lammps`): Contains force field parameters for bonds, angles, and nonbonded interactions


## Table of Contents
- [Open Source](#open-source)
- [Citation](#citation)
- [Platform Compatibility](#platform-compatibility)
- [What You Need](#what-you-need)
- [Installation](#installation)
  - [AmberTools](#ambertools)
  - [Python Packages](#python-packages)
  - [LAMMPS](#lammps)
- [SMILES to PDB Workflow](#smiles-to-pdb-workflow)
- [Command Reference](#command-reference)
- [How It Works](#how-it-works)
- [Tutorial](#tutorial)
  - [CLI with LAMMPS execution](#cli-with-lammps-execution)
  - [Python API with LAMMPS execution](#python-api-with-lammps-execution)
- [Validation with InterMol](#validation-with-intermol)
- [Contributing](#contributing)
- [Acknowledgements](#acknowledgements)

## Open Source

This is an open-source project. The source code is freely available for use, modification, and
distribution under the MIT License.

## Citation

If you use this software in your research, please cite it as:

**DOI:** [10.5281/zenodo.18114886](https://doi.org/10.5281/zenodo.18114886)

## Platform Compatibility

AMBER2LAMMPS has been validated and tested on:

- **Linux** (Ubuntu, CentOS, Red Hat, Debian) - Fully tested and validated
- **macOS** (Intel and Apple Silicon) - Fully tested
- **Windows** (Windows 10/11 with WSL2 and native Python) - Tested with WSL2 and native Python
- **WSL2 is recommended for Windows users for best compatibility.**

## What You Need

- **Structure input**: A PDB file of your molecule (or a SMILES string you can convert to PDB; see the workflow below).
- **AMBER prep tools**: AmberTools (`antechamber`, `parmchk2`, `tleap`) to generate AMBER files:
  - **`.prmtop`**: Topology file defining molecular structure (bonds, angles, atom types)
  - **`.mol2`**: Coordinate file with atomic positions and partial charges
  - **`.frcmod`**: Force field parameters defining how atoms interact
- **Python**: Python 3.8+ with `parmed` (for reading AMBER files) and `numpy` (for calculations).
- **LAMMPS**: A build that includes `MOLECULE`, `KSPACE`, and `EXTRA-MOLECULE` packages on your `PATH` (`lmp -h` to confirm).

## Installation

### AmberTools

Install AmberTools (e.g., AmberTools23) from https://ambermd.org/GetAmber.php#ambertools and
activate the environment:

```bash
conda activate AmberTools23  # or your AmberTools version
```

### Python Packages

```bash
# Recommended
conda install -c conda-forge parmed numpy

# Alternative
pip install parmed numpy
```

### LAMMPS

Download LAMMPS (https://lammps.org/) and build a recent stable version with packages `MOLECULE`,
`KSPACE`, and `EXTRA-MOLECULE` enabled. Verify your install with `lmp -h` or `which lmp`.

## SMILES to PDB Workflow

Use `obabel` to generate a PDB file from SMILES.

### Installation

```bash
# macOS
brew install open-babel

# conda
conda install -c openbabel openbabel

# Ubuntu/Debian
sudo apt-get install openbabel
```

### Generate PDB from SMILES

```bash
# Basic conversion
obabel -:CCO -opdb -O ethanol.pdb --gen3d

# With explicit hydrogens (recommended)
obabel -:CCO -h -opdb -O ethanol.pdb --gen3d

# Other examples
obabel -:c1ccccc1 -opdb -O benzene.pdb --gen3d  # Benzene
obabel -:"CC(=O)OC1=CC=CC=C1C(=O)O" -opdb -O aspirin.pdb --gen3d  # Aspirin
```

## Command Reference

Outputs: a LAMMPS data file (`<data_file>`, e.g., `data.lammps`) and a separate parameter file
(`<param_file>`, e.g., `parm.lammps`).

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


## How It Works

1. **Input validation**: `validate_files` checks that `topology`, `mol2`, and `frcmod` exist and are readable. The CLI performs this validation automatically, while it's optional for the Python API.
2. **Load topology and parameters**: ParmEd reads atoms, bonds, angles, and dihedrals from `.prmtop`; `frcmod` provides masses and nonbonded terms.
3. **Atom typing and masses**: `MASS` entries in `frcmod` map atom types to sequential IDs and masses; missing types warn and fall back to type 1.
4. **Coordinates and box**: Coordinates and charges come from the MOL2 `ATOM` block. The box is a bounding box expanded by the `buffer`.
5. **Charge normalization**: If `|net charge| > 1e-6`, charges are shifted uniformly to make the system neutral.
6. **Nonbonded coefficients**: `NONBON` entries in `frcmod` become `pair_coeff` terms in `parm.lammps`.
7. **Bonded coefficients**: Bond, angle, and dihedral coefficients come from the ParmEd-exported AMBER terms (via `frcmod`/`.prmtop`) and are written into the data/parameter files.
8. **Export and cleanup**: Topology terms are written to LAMMPS data/parameter files; temporary helper files are removed.

## Tutorial

The tutorial assumes you are running from an AMBER2LAMMPS checkout or otherwise have access to
`amber_to_lammps.py`.

### CLI with LAMMPS execution

#### 1. Generate AMBER Input Files (`prmtop`, `mol2`, `frcmod`)

Activate the AmberTools environment before running:

```python
import subprocess

pdb_file = "epon.pdb"  # Replace with your PDB filename
base_name = pdb_file.replace(".pdb", "")

# Generate MOL2 file with charges from PDB file
cmd1 = f"antechamber -j 4 -at gaff2 -dr yes -fi pdb -fo mol2 -i {pdb_file} -o {base_name}.mol2 -c bcc"
subprocess.run(cmd1, shell=True)

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

# Output files generated: epon.prmtop, epon.crd, epon.mol2, epon.frcmod
```

#### 2. Run the Conversion

```bash
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.mol2 epon.frcmod --verbose -b 4.5
```

#### 3. Run LAMMPS

Use the provided `example_lammps_input.lmp`:

```bash
lmp < example_lammps_input.lmp
```

**Note**: The LAMMPS input file `example_lammps_input.lmp` `include`s `parm.lammps` (force field
parameters) and `read_data` loads `data.lammps` (coordinates). Keeping parameters separate makes it
easier to swap parameter sets without regenerating coordinates.

#### Additional CLI examples

```bash
# Custom buffer and verbose logging (adds 5 Å padding)
python3 amber_to_lammps.py my_data.lammps my_params.lammps \
   ethanol.prmtop ethanol.mol2 ethanol.frcmod --verbose -b 5.0

# Custom output names
python3 amber_to_lammps.py system.data system.parm system.prmtop \
   system.mol2 system.frcmod

# Minimal output without verbose logging
python3 amber_to_lammps.py small.data small.parm \
   ethanol.prmtop ethanol.mol2 ethanol.frcmod -b 3.0

# Using absolute paths with custom buffer
python3 amber_to_lammps.py /home/user/lammps/output/data.lammps /home/user/lammps/output/param.lammps /home/user/amber/topology.prmtop /home/user/amber/coords.mol2 /home/user/amber/params.frcmod -b 4.5
```

### Python API with LAMMPS execution

```python
from amber_to_lammps import amber2lammps, validate_files

data_file = 'data.lammps'
param_file = 'parm.lammps'
topology = 'epon.prmtop'
mol2 = 'epon.mol2'
frcmod = 'epon.frcmod'

# Optional validation
validate_files(topology, mol2, frcmod)

# Convert
amber2lammps(
    data_file=data_file,
    param_file=param_file,
    topology=topology,
    mol2=mol2,
    frcmod=frcmod,
    buffer=5.0,      # 5 Å padding
    verbose=True     # Show progress
)

print("Completed conversion")
```

#### Additional API examples

```python
# Example 1: Basic conversion without validation
from amber_to_lammps import amber2lammps

amber2lammps(
    data_file='system.data',
    param_file='system.parm',
    topology='system.prmtop',
    mol2='system.mol2',
    frcmod='system.frcmod'
)

# Example 2: Custom buffer and verbose output
amber2lammps(
    data_file='molecule.data',
    param_file='molecule.parm',
    topology='molecule.prmtop',
    mol2='molecule.mol2',
    frcmod='molecule.frcmod',
    buffer=5.0,
    verbose=True
)

# Example 3: Batch processing multiple molecules
from amber_to_lammps import amber2lammps, validate_files

molecules = ['ethanol', 'benzene', 'aspirin']

for mol in molecules:
    validate_files(f'{mol}.prmtop', f'{mol}.mol2', f'{mol}.frcmod')
    amber2lammps(
        data_file=f'{mol}.data',
        param_file=f'{mol}.parm',
        topology=f'{mol}.prmtop',
        mol2=f'{mol}.mol2',
        frcmod=f'{mol}.frcmod',
        verbose=True
    )
```


## Validation with InterMol

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

This is because charges generated from antechamber AM1-BCC method carry slight excess charge and the
magnitude of the charge is greater than 1e-6

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

## Troubleshooting

- **Missing MASS/atom types**: Check `ethanol.frcmod` warnings; add types or correct `frcmod` before
converting.
- **Net charge not zero**: Charge normalization shifts charges; verify source charges and rerun if
unintended.
- **Atoms too close or outside box**: Increase `--buffer` and rerun; inspect verbose box extents.
- **LAMMPS run errors about packages**: Rebuild LAMMPS with `MOLECULE`, `KSPACE`, `EXTRA-MOLECULE`.

## Contributing

When making modifications:
1. Create a new branch: `git checkout -b feature-name`
2. Make your changes
3. Test thoroughly with different input files
4. Commit and push your changes

## Acknowledgements

Thanks to the following people for their suggestions in improving the tutorial:

- Dr. Axel Kohlmeyer
- Dr. Germain Clavier
