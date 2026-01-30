# Overview

Python utility that converts AMBER `.prmtop`/`.crd` files to LAMMPS data/parameter formats through CLI or Python API, featuring validation-first error handling.

## What This Tool Does

Run molecular dynamics simulations in **LAMMPS** using systems prepared in **AMBER** format.

**Typical workflow:**
1. Start with a molecular structure (PDB file or SMILES string).
2. Generate AMBER files (`.prmtop`, `.crd`) with AmberTools.
3. Convert the AMBER files to LAMMPS format with this tool.
4. Run your simulation in LAMMPS.

**Input → Output Mapping**
- **`.prmtop`**: Molecular topology (bonds, angles, atom types, force-field parameters) → LAMMPS data and parameter files.
- **`.crd`**: atomic positions → LAMMPS coordinates.
- **Outputs**: `data.lammps` (coordinates, box, topology) and `parm.lammps` (bonded and nonbonded parameters).


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
    - [Additional CLI examples](#additional-cli-examples)
  - [Python API with LAMMPS execution](#python-api-with-lammps-execution)
    - [Additional API examples](#additional-api-examples)
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

Tested and validated on:
- **Linux** (Ubuntu, CentOS, Red Hat, Debian)
- **macOS** (Intel and Apple Silicon)
- **Windows 10/11** (native Python or WSL2; WSL2 recommended)

## What You Need

- **Structure input**: A PDB file (or a SMILES string that you convert to PDB; see workflow below).
- **AMBER prep**: AmberTools (`antechamber`, `tleap`) to generate `.prmtop`.
- **Python**: Python 3.8+ with `parmed` and `numpy`.
- **LAMMPS**: Build with `MOLECULE`, `KSPACE`, and `EXTRA-MOLECULE` packages on your `PATH` (`lmp -h` to confirm).

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
`KSPACE`, and `EXTRA-MOLECULE` enabled. Verify your installation with `lmp` or `which lmp`.

## SMILES to PDB Workflow

Use Open Babel (`obabel`) to generate a 3D PDB from SMILES.

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

| Argument | Required | Description |
|-----------|----------|-------------|
| `data_file` | yes | Output LAMMPS data file name |
| `param_file` | yes | Output LAMMPS parameter file name |
| `topology` | yes | AMBER topology file (`.prmtop`) |
| `crd` | yes | AMBER coordinate file (`.crd`) |
| `-b, --buffer` | optional | Vacuum padding (Å) for the simulation box. Default: `3.8`. |
| `--verbose` | optional | Print step-by-step progress, counts, and box size. Default: `False`. |
| `--keep-temp` | optional | Keep temporary files (bonds.txt, angles.txt, dihedrals.txt, pairs.txt) after conversion. Default: `False`. |
| `--charge` | yes | Target net charge (integer). Applies a uniform offset to every atom to reach this charge (1e-6 tolerance). Does not add counterions—add them upstream for charged systems. |
| `-h, --help` | optional | Show help message. |

## How It Works

1. **Input validation**: `validate_files` checks that `topology` and `crd` exist and are readable. The CLI performs this validation automatically, while the Python API requires manual validation if desired.
2. **Load topology and parameters**: ParmEd reads atoms, bonds, angles, dihedrals, masses, and LJ parameters directly from `.prmtop`.
3. **Atom types and masses**: ParmEd extracts atom types and masses from `.prmtop` using printDetails.
4. **Coordinates and box**: Coordinates come from the CRD file and charges from topology. The box is a bounding box expanded by the `buffer`.
5. **Charge normalization**: Charges are shifted uniformly to achieve the target net charge specified by `--charge`. A tolerance of 1e-6 is used for charge matching.
6. **Nonbonded coefficients**: LJ parameters extracted from topology become `pair_coeff` terms in `parm.lammps`.
7. **Bonded coefficients**: Bond, angle, and dihedral coefficients come from the ParmEd-exported AMBER terms in `.prmtop` and are written into the data/parameter files. Each bond/angle/dihedral gets a unique type ID equal to its index, so the number of types matches the number of terms.
8. **Export and cleanup**: Topology terms are written to LAMMPS data/parameter files; temporary helper files (pairs.txt, bonds.txt, angles.txt, dihedrals.txt) are removed unless `--keep-temp` is specified. These files can aid debugging, feed other scripts, or serve as supporting information in publications.

### Understanding Charge Schemes and Normalization

**AMBER Charge Schemes:**
AMBER uses various charge methods including:
- **RESP** (Restrained Electrostatic Potential): Derived from quantum mechanical electrostatic potential
- **AM1-BCC**: Semi-empirical charges with bond charge corrections
- **CM5** or **CM1A**: Charge models based on atomic charges


**Charge Normalization in AMBER2LAMMPS:**
- Most systems should typically be neutral for PME convergence; AM1-BCC charge calculations often leave small residuals (±0.003).
- For neutral molecules, run with `--charge 0` so AMBER2LAMMPS applies a uniform offset that removes the residual charge; this prevents the error from scaling up when the system is replicated in LAMMPS.
- For intentionally charged species (e.g., protonated or deprotonated), add counterions in tleap/packmol before conversion. AMBER2LAMMPS never adds ions; it only shifts existing charges to your requested total.
- How the flag works:
  - `--charge 0`: adds a uniform offset so total charge is 0 within 1e-6.
  - `--charge +1` (or any integer): adds a uniform offset so the total matches that integer within 1e-6.
  - The same constant is added to every atom, so relative charge differences are preserved.

**Example:** If your system has net charge +0.003 and you specify `--charge 0`, each atom's charge will be reduced by `(0.003 ÷ number_of_atoms)` to reach neutrality.


## Tutorial

The tutorial assumes you are running from an AMBER2LAMMPS checkout or otherwise have access to
`amber_to_lammps.py`.


### CLI with LAMMPS execution

#### 1. Generate AMBER Input Files (`prmtop`, `crd`)

Activate the AmberTools environment before running:

```python
import subprocess

pdb_file = "epon.pdb"  # Replace with your PDB filename
base_name = pdb_file.replace(".pdb", "")

# Generate MOL2 file with charges from PDB file
cmd1 = f"antechamber -j 4 -at gaff2 -dr yes -fi pdb -fo mol2 -i {pdb_file} -o {base_name}.mol2 -c bcc"
subprocess.run(cmd1, shell=True)
# -c flag specifies the charge method. bcc is used for AM1-BCC charges
# -nc flag specifies net charge (default: 0.0). Use for charged molecules (e.g., -nc 1 for +1 charge)

# Create tleap input file
with open("tleap.in", "w") as f:
    f.write("source leaprc.gaff2\n") # Load the Gaff2 force field
    f.write(f"SUS = loadmol2 {base_name}.mol2\n") # Load the molecule
    f.write("check SUS\n") # Check the molecule
    f.write(f"saveamberparm SUS {base_name}.prmtop {base_name}.crd\n") # Save the AMBER files
    f.write("quit")

# Run tleap to generate AMBER files
cmd3 = "tleap -f tleap.in"
subprocess.run(cmd3, shell=True)

# Output files generated: epon.prmtop, epon.crd, epon.mol2
```

#### 2. Run the Conversion

```bash
python3 amber_to_lammps.py data.lammps parm.lammps epon.prmtop epon.crd --charge 0 --verbose -b 4.5
```

#### 3. Run LAMMPS

Use the provided `example_lammps_input.lmp`:

```bash
lmp < example_lammps_input.lmp
```

**Note**: The LAMMPS input file `example_lammps_input.lmp` includes `parm.lammps` (force-field
parameters) and uses `read_data` to load `data.lammps` (coordinates). Keeping parameters separate
makes it easy to swap parameter sets without regenerating coordinates.

#### Additional CLI examples

```bash
# Custom buffer and verbose logging (adds 5 Å padding)
python3 amber_to_lammps.py my_data.lammps my_params.lammps \
   ethanol.prmtop ethanol.crd --charge 0 --verbose -b 5.0

# Custom output names
python3 amber_to_lammps.py system.data system.parm system.prmtop \
   system.crd --charge 0

# Minimal output without verbose logging
python3 amber_to_lammps.py data.lammps parm.lammps ethanol.prmtop ethanol.crd --charge 0 -b 3.0

# Using absolute paths with custom buffer
python3 amber_to_lammps.py /home/user/lammps/output/data.lammps /home/user/lammps/output/param.lammps /home/user/amber/topology.prmtop /home/user/amber/coords.crd --charge 0 -b 4.5

# Keep temporary files for debugging
python3 amber_to_lammps.py debug_data.lammps debug_parm.lammps molecule.prmtop molecule.crd --charge 0 --keep-temp --verbose
```
`example_lammps_input.lmp` includes `parm.lammps` and reads `data.lammps`; adjust those paths if you rename outputs.

### Python API with LAMMPS execution

```python
from amber_to_lammps import amber2lammps, validate_files
import subprocess

data_file = 'data.lammps'
param_file = 'parm.lammps'
topology = 'epon.prmtop'
crd = 'epon.crd'

# Optional validation
validate_files(topology, crd)

# Convert
amber2lammps(
    data_file=data_file,
    param_file=param_file,
    topology=topology,
    crd=crd,
    charge=0,
    buffer=3.8,
    verbose=True,
    keep_temp=False
)

# Run LAMMPS
subprocess.run("lmp < example_lammps_input.lmp", shell=True)
print("Completed conversion")
```


#### Additional API examples

```python
# Example 1: Custom buffer and keep temporary files
from amber_to_lammps import amber2lammps

amber2lammps(
    data_file='epon.data',
    param_file='epon.parm',
    topology='epon.prmtop',
    crd='epon.crd',
    charge=0,  # Target neutral charge
    buffer=5.0,
    verbose=True,
    keep_temp=True  # Keep temporary files for inspection
)

# Example 2: Batch processing multiple molecules
from amber_to_lammps import amber2lammps, validate_files

molecules = ['ethanol', 'benzene', 'aspirin']

for mol in molecules:
    validate_files(f'{mol}.prmtop', f'{mol}.crd')
    amber2lammps(
        data_file=f'{mol}.data',
        param_file=f'{mol}.parm',
        topology=f'{mol}.prmtop',
        crd=f'{mol}.crd',
        charge=0,
        verbose=True,
        keep_temp=False
    )
```


## Validation with InterMol

Conversion results have been cross-checked with **InterMol**.

Install InterMol: https://github.com/shirtsgroup/InterMol

**Usage (InterMol):**
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
Results from InterMol conversion exhibit slight non-neutrality (e.g., net charge = -0.002996) because AM1-BCC charges are not renormalized.

**Output from AMBER2LAMMPS (Tutorial Step 3)**

```bash
lmp < example_lammps_input.lmp
```

```text
E_bond        E_angle        E_dihed        E_impro         E_pair         E_vdwl         E_coul         E_long         E_tail         PotEng
2.3161274      6.0940126      12.475827      0             -9.1202197      10.511316      108.01561     -127.64715     -0.31768789     11.765747
```

## Contributing

When making changes:
1. Create a branch: `git checkout -b feature-name`
2. Implement the change
3. Test with varied inputs
4. Commit and push

## Acknowledgements

Thanks to Dr. Axel Kohlmeyer and Dr. Germain Clavier for tutorial feedback, and to Dr. Andrew Jewitt
(author of moltemplate) for discussions via the mailing list and email.
