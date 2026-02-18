# Overview

Python utility that converts AMBER `.prmtop`/PDB files to LAMMPS data/parameter formats through CLI or Python API, featuring validation-first error handling. **Can handle multiple AMBER topology files for mixed molecular systems (e.g., drug + solvent mixtures).**

## Key Features

- **🔧 Multi-Topology Conversion**: Convert mixed systems using multiple AMBER `.prmtop` files in a single run
- **🧪 Mixed & Replicated Systems**: Handle mixtures (e.g., drug + solvent) and multiple copies of one molecule via `-c`
- **✅ Validation-First**: Clear errors for missing files, mismatched list lengths, and PDB atom-count/order issues
- **🔄 CLI + Python API**: Same conversion engine from the command line or `amber2lammps(...)` in Python
- **📦 Debug Artifacts (`--keep-temp`)**: Optionally write `pairs.txt`, `bonds.txt`, `angles.txt`, `dihedrals.txt` grouped by topology

## What This Tool Does

Run molecular dynamics simulations in **LAMMPS** using systems prepared in **AMBER** format. Supports both single-molecule and mixed molecular systems with multiple topology files.

**Typical workflow:**
1. Start with a molecular structure (PDB file or SMILES string).
2. Generate AMBER topology files (`.prmtop`) with AmberTools (and optionally `.crd` for other tools/validation).
3. Convert (`.prmtop` + PDB coordinates) to LAMMPS format with this tool.
4. Run your simulation in LAMMPS.

**Input → Output Mapping**
- **`.prmtop`**: Molecular topology (bonds, angles, atom types, force-field parameters) → LAMMPS data and parameter files. Supports multiple topologies for mixed systems.
- **PDB**: Atomic coordinates (single molecule or combined from PackMol) → LAMMPS coordinates with proper atom indexing.
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
  - [Multiple Copies and Mixed Molecular System Workflow](#multiple-copies-and-mixed-molecular-system-workflow)
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

- **Structure input**: A PDB file of the molecular structure (or a SMILES string that you convert to PDB; see workflow below). For mixed systems, use a combined PDB file from PACKMOL containing all molecules.
- **AMBER prep**: AmberTools (`antechamber`, `tleap`) to generate `.prmtop`. For mixed systems, generate separate `.prmtop` files for each molecule type.
- **Python**: Python 3.8+ with `parmed` and `numpy`.
- **LAMMPS**: Build with `MOLECULE`, `KSPACE`, and `EXTRA-MOLECULE` packages on your `PATH` (`lmp -h` to confirm).

**Note**: If you are using a PDB file generated from SMILES or another source, it is recommended that you pass in the PDB file with the `antechamber -dr yes` option. PDB format note: The converter uses fixed-column parsing of ATOM/HETATM records. PackMol-style PDBs work well; non-standard PDBs (e.g., unusual formatting, multi-model files) may fail or require cleanup.

## Installation

### PACKMOL 

Required for running LAMMPS simulation multiple copies of the molecule or mixed molecular systems. Install from:
https://m3g.github.io/packmol/

PACKMOL creates combined PDB files with multiple molecules positioned in a simulation box.

Verify installation by running:
```bash
which packmol
```

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
| `pdb_file` | yes | Combined PDB containing all molecules (e.g., from PackMol) |
| `-t, --topologies` | yes | One or more AMBER topology files (`.prmtop`) |
| `-c, --counts` | yes | Molecule counts for each topology (same order/length as `--topologies`) |
| `--charges` | yes | Target net charge per topology (list matching `--topologies`; use `0 0 ...` for neutrals) |
| `-b, --buffer` | optional | Vacuum padding (Å) for the simulation box. Default: `3.8`. |
| `--verbose` | optional | Print step-by-step progress, counts, and box size. Default: `False`. |
| `--keep-temp` | optional | Keep temporary files (`pairs.txt`, `bonds.txt`, `angles.txt`, `dihedrals.txt`) grouped by topology. Default: `False`. |
| `-h, --help` | optional | Show help message. |

## How It Works

1. **Input validation**: `validate_files` checks that `topologies`, `counts`, and `pdb_file` exist and lengths match. CLI does this automatically; Python API requires calling it yourself.
2. **Load topology and parameters**: ParmEd reads atoms, bonds, angles, dihedrals, masses, and LJ parameters directly from each `.prmtop` (supports multiple topologies).
3. **Atom types and masses**: Atom types are namespaced on conflict across topologies; masses and LJ params are captured per canonical type.
4. **Coordinates and box**: Coordinates come from the combined PDB; the box is min/max of those coordinates expanded by `buffer`.
5. **Charge normalization**: For each topology, charges are uniformly shifted to hit the user-provided `--charges` target (tolerance 1e-6). Totals are checked across all molecules.
6. **Nonbonded coefficients**: Like–like `pair_coeff` lines are emitted per atom type; cross terms are left to LAMMPS mixing rules.
7. **Bonded coefficients**: Bonds, angles, dihedrals are written from ParmEd data; each instance gets a unique type ID. Multi-term torsions are preserved for `dihedral_style fourier`.
8. **Export and cleanup**: Data/parameter files are written; debug files (`pairs.txt`, `bonds.txt`, `angles.txt`, `dihedrals.txt`) are kept if `--keep-temp` is set.

### Understanding Charge Schemes and Normalization

**AMBER Charge Schemes:**
AMBER uses various charge methods including:
- **RESP** (Restrained Electrostatic Potential): Derived from quantum mechanical electrostatic potential
- **AM1-BCC**: Semi-empirical charges with bond charge corrections
- **CM5** or **CM1A**: Charge models based on atomic charges


**Charge Normalization in AMBER2LAMMPS:**
- Most systems should typically be neutral for PME convergence; AM1-BCC charge calculations often leave small residuals (±0.003).
- For neutral molecules, run with `--charges 0` so AMBER2LAMMPS applies a uniform offset that removes the residual charge; this prevents the error from scaling up when the system is replicated in LAMMPS or when you have multiple copies of the same molecule.
- For intentionally charged species (e.g., protonated or deprotonated), add counterions in tleap/packmol before conversion. AMBER2LAMMPS never adds ions; it only shifts existing charges to your requested total.
- How the flag works:
  - `--charges 0`: adds a uniform offset so total charge is 0 within 1e-6.
  - `--charges +1` (or any integer): adds a uniform offset so the total matches that integer within 1e-6.
  - The same constant is added to every atom, so relative charge differences are preserved.

**Example:** If your system has net charge +0.003 and you specify `--charges 0`, each atom's charge will be reduced by `(0.003 ÷ number_of_atoms)` to reach neutrality.


## Tutorial

The tutorial assumes you are running from an AMBER2LAMMPS checkout or otherwise have access to
`amber_to_lammps.py`.


### CLI with LAMMPS execution

#### 1. Generate AMBER Topology (`.prmtop`)

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
    f.write(f"MOLECULE = loadmol2 {base_name}.mol2\n") # Load the molecule
    f.write("check MOLECULE\n") # Check the molecule
    f.write(f"saveamberparm MOLECULE {base_name}.prmtop {base_name}.crd\n") # Save the AMBER files
    f.write("quit")

# Run tleap to generate AMBER files
cmd3 = "tleap -f tleap.in"
subprocess.run(cmd3, shell=True)

# Output files generated: epon.prmtop, epon.crd, epon.mol2
# Note: amber_to_lammps.py uses epon.prmtop + epon.pdb; the .crd is optional (kept for validation/other tools).
```

#### 2. Run the Conversion

```bash
python3 amber_to_lammps.py data.lammps parm.lammps epon.pdb \
  -t epon.prmtop -c 1 --charges 0 --verbose -b 4.5
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
python3 amber_to_lammps.py my_data.lammps my_params.lammps ethanol.pdb \
  -t ethanol.prmtop -c 1 --charges 0 --verbose -b 5.0

# Custom output names
python3 amber_to_lammps.py system.data system.parm system.pdb \
  -t system.prmtop -c 1 --charges 0

# Minimal output without verbose logging
python3 amber_to_lammps.py data.lammps parm.lammps ethanol.pdb \
  -t ethanol.prmtop -c 1 --charges 0 -b 3.0

# Using absolute paths with custom buffer
python3 amber_to_lammps.py /home/user/lammps/output/data.lammps /home/user/lammps/output/param.lammps /home/user/amber/packmol/combined.pdb \
  -t /home/user/amber/topology.prmtop -c 1 --charges 0 -b 4.5

# Single topology with multiple copies (e.g., 10 ethanol molecules)
python3 amber_to_lammps.py multi_ethanol_data.lammps multi_ethanol_parm.lammps multi_ethanol.pdb \
  -t ethanol.prmtop -c 10 --charges 0 --verbose
```
If you include `--keep-temp`, the converter will also write `pairs.txt`, `bonds.txt`, `angles.txt`, and `dihedrals.txt`, each grouped by topology (handy for debugging mixed systems).
`example_lammps_input.lmp` includes `parm.lammps` and reads `data.lammps`; adjust those paths and the file names in example_lammps_input.lmp if you rename outputs.

### Python API with LAMMPS execution

```python
from amber_to_lammps import amber2lammps, validate_files
import subprocess

data_file = 'data.lammps'
param_file = 'parm.lammps'
pdb_file = 'epon.pdb'
topologies = ['epon.prmtop']
counts = [1]
charges = [0]

# Optional validation
validate_files(topologies, counts, pdb_file)

# Convert
amber2lammps(
    data_file=data_file,
    param_file=param_file,
    topologies=topologies,
    molecule_counts=counts,
    pdb_file=pdb_file,
    charges_target=charges,
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
    topologies=['epon.prmtop'],
    molecule_counts=[1],
    pdb_file='epon.pdb',
    charges_target=[0],  # Target neutral charge
    buffer=5.0,
    verbose=True,
    keep_temp=True  # Keep temporary files for inspection
)



# Example 2: Batch processing multiple molecules 
from amber_to_lammps import amber2lammps, validate_files

molecules = ['ethanol', 'benzene', 'aspirin']

for mol in molecules:
    validate_files([f'{mol}.prmtop'], [1], f'{mol}.pdb')
    amber2lammps(
        data_file=f'{mol}.data',
        param_file=f'{mol}.parm',
        topologies=[f'{mol}.prmtop'],
        molecule_counts=[1],
        pdb_file=f'{mol}.pdb',
        charges_target=[0],
        verbose=True,
        keep_temp=False
    )
```

## Multiple Copies and Mixed Molecular System Workflow

For systems containing multiple molecule types (e.g., drug + solvent mixtures), use the following workflow.

### Mixed Molecular System (Multiple Topologies)

#### 1. Prepare Individual Molecules
```bash
# Generate MOL2 files with proper bond orders (Open Babel + Antechamber)
obabel aspirin.pdb -opdb -O aspirin_obabel.pdb -h --gen3d
obabel aspirin_obabel.pdb -omol2 -O aspirin_obabel.mol2 -h
antechamber -fi mol2 -fo mol2 -i aspirin_obabel.mol2 -o aspirin_final.mol2 -c bcc

# Repeat for benzene and ethanol
obabel benzene.pdb -opdb -O benzene_obabel.pdb -h --gen3d
obabel benzene_obabel.pdb -omol2 -O benzene_obabel.mol2 -h
antechamber -fi mol2 -fo mol2 -i benzene_obabel.mol2 -o benzene_final.mol2 -c bcc

obabel ethanol.pdb -opdb -O ethanol_obabel.pdb -h --gen3d
obabel ethanol_obabel.pdb -omol2 -O ethanol_obabel.mol2 -h  
antechamber -fi mol2 -fo mol2 -i ethanol_obabel.mol2 -o ethanol_final.mol2 -c bcc
```

#### 2. Generate AMBER Topologies
```bash
# Create tleap input files and run for each molecule
tleap -f tleap_aspirin.in
tleap -f tleap_benzene.in
tleap -f tleap_ethanol.in
```

#### 3. Build Mixed System with PackMol

```bash
# Create PackMol input file
cat > packmol_mix.in << 'EOF'
tolerance 2.0
filetype pdb
output mixed_system.pdb

# Aspirin (1 molecule)
structure aspirin.pdb
  number 1
  inside box 0.0 0.0 0.0 10.0 10.0 10.0
end structure

# Benzene (2 molecules)
structure benzene.pdb
  number 2
  inside box 0.0 0.0 0.0 10.0 10.0 10.0
end structure

# Ethanol (5 molecules)
structure ethanol.pdb
  number 5
  inside box 0.0 0.0 0.0 10.0 10.0 10.0
end structure
EOF

# Run PackMol
packmol < packmol_mix.in

# The resulting mixed_system.pdb is a single combined PDB containing all molecules
# in the same order and counts you will pass to -t/-c. Keep residue/atom ordering intact.
```

#### 4. Convert to LAMMPS

```bash
# Convert the mixed molecular system to LAMMPS format
python amber_to_lammps.py mixed_data.lammps mixed_parm.lammps mixed_system.pdb \
  -t aspirin.prmtop benzene.prmtop ethanol.prmtop -c 1 2 5 --charges 0 0 0 --verbose
```

#### Python API

```python
from amber_to_lammps import amber2lammps, validate_files
import subprocess

# Validate all topology files
validate_files(['aspirin.prmtop', 'benzene.prmtop', 'ethanol.prmtop'], [1, 2, 5], 'mixed_system.pdb')

# Convert mixed system with multiple molecule types
amber2lammps(
    data_file='mixed_data.lammps',
    param_file='mixed_parm.lammps', 
    topologies=['aspirin.prmtop', 'benzene.prmtop', 'ethanol.prmtop'],
    molecule_counts=[1, 2, 5],  # 1 aspirin, 2 benzene, 5 ethanol
    pdb_file='mixed_system.pdb',
    charges_target=[0, 0, 0],  # All neutral molecules
    verbose=True,
    keep_temp=False
)

# Run LAMMPS with mixed system
subprocess.run("lmp < test_mixed_system.in", shell=True)
```

### Multiple Copies (Single Topology)

#### 1. Prepare Individual Molecules
```bash
# Generate MOL2 files with proper bond orders (Open Babel + Antechamber)
obabel ethanol.pdb -opdb -O ethanol_obabel.pdb -h --gen3d
obabel ethanol_obabel.pdb -omol2 -O ethanol_obabel.mol2 -h
antechamber -fi mol2 -fo mol2 -i ethanol_obabel.mol2 -o ethanol_final.mol2 -c bcc
```

#### 2. Generate AMBER Topologies
```bash
# Create tleap input files and run for the molecule
tleap -f tleap_ethanol.in
```

#### 3. Build Multi-Copy System with PackMol

```bash
# Create PackMol input file for multiple copies (packmol_multi_ethanol.in)
cat > packmol_multi_ethanol.in << 'EOF'
tolerance 2.0
filetype pdb
output multi_ethanol.pdb

# Ethanol (10 molecules)
structure ethanol.pdb
  number 10
  inside box 0.0 0.0 0.0 15.0 15.0 15.0
end structure
EOF

# Run PackMol
packmol < packmol_multi_ethanol.in

# The resulting multi_ethanol.pdb is a single combined PDB containing all molecules
# in the same order and counts you will pass to -t/-c. Keep residue/atom ordering intact.
```

#### 4. Convert to LAMMPS

```bash
# Convert the multi-copy molecular system to LAMMPS format
python amber_to_lammps.py multi_ethanol_data.lammps multi_ethanol_parm.lammps multi_ethanol.pdb \
  -t ethanol.prmtop -c 10 --charges 0 --verbose
```

#### Python API 

```python
from amber_to_lammps import amber2lammps, validate_files

# Validate single topology with multiple copies
validate_files(['ethanol.prmtop'], [10], 'multi_ethanol.pdb')

# Convert multiple copies of the same molecule
amber2lammps(
    data_file='multi_ethanol_data.lammps',
    param_file='multi_ethanol_parm.lammps',
    topologies=['ethanol.prmtop'],
    molecule_counts=[10],  # 10 ethanol molecules
    pdb_file='multi_ethanol.pdb',
    charges_target=[0],
    verbose=True,
    keep_temp=False
)
```

```python
# Use appropriate LAMMPS input for mixed systems
subprocess.run("lmp < test_mixed_system.in", shell=True)
```

## Validation with InterMol

Conversion results have been cross-checked with **InterMol**.

Install InterMol: https://github.com/shirtsgroup/InterMol

**Usage (InterMol):**
```bash
python convert.py --amb_in epon.prmtop epon.crd --lammps
```
InterMol expects an AMBER coordinate file (`.crd`); `amber_to_lammps.py` uses the PDB for coordinates.

**Generated Files:**
- `epon_converted.input` (LAMMPS input file)
- `epon_converted.lmp` (LAMMPS data file)


```bash
lmp < epon_converted.input
```

**Energy Comparison Results**

*Output from InterMol and AMBER2LAMMPS (see tutorial Step 3 of CLI with LAMMPS execution)*

| Energy Component | InterMol | AMBER2LAMMPS | Difference |
|------------------|----------|---------------|------------|
| E_bond           | 2.2879665 | 2.2879665     | 0.0000 |
| E_angle          | 7.0137437 | 7.0137055     | 0.0000382 |
| E_dihed          | 8.0302669 | 8.0302812     | -0.0000143 |
| E_impro          | 0.0000    | 0.0000        | 0.0000 |
| E_pair           | -11.148889 | -11.320542    | 0.171653 |
| E_vdwl           | 8.5497495 | 8.3173103     | 0.2324392 |
| E_coul           | 97.869973 | 107.25122     | -9.381247 |
| E_long           | -117.56861 | -126.88907    | 9.32046 |
| E_tail           | -0.0042213943 | -0.23666167  | 0.2324403 |
| **PotEng**       | **6.1830883** | **6.0114113** | **0.171677** |


## Contributing

When making changes:
1. Create a branch: `git checkout -b feature-name`
2. Implement the change
3. Test with varied inputs
4. Commit and push

## Acknowledgements

Thanks to Dr. Axel Kohlmeyer and Dr. Germain Clavier for tutorial feedback, and to Dr. Andrew Jewitt
(author of moltemplate) for discussions via the mailing list and email.
