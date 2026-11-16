# Pure GLYCAM06j-1 glucose validation

This fixture validates a neutral, standalone alpha-D-glucopyranose molecule
parameterized with GLYCAM06j-1. It reuses the topology and coordinates generated
by the documented `tleap` recipe in `validation/amber_family_mix/source/` and
places them in a 30 Angstrom orthorhombic periodic cell.

The validation regenerates the AMBER2LAMMPS files, runs AmberTools `sander` and
LAMMPS on identical coordinates and cell, and compares bond, angle, dihedral,
van der Waals, electrostatic and total potential energies plus every Cartesian
force component. Generated files are kept in the ignored `work/` directory.

## Input and validation files

The GLYCAM inputs already existed in the mixed-system fixture, so this folder
references them instead of committing a second copy:

| Purpose | File |
|---|---|
| GLYCAM06j-1 topology | `../amber_family_mix/source/glucose.prmtop` |
| Original coordinates | `../amber_family_mix/source/glucose.pdb` |
| Reproducible `tleap` recipe | `../amber_family_mix/source/tleap_glucose.in` |
| LAMMPS run-0 input | `glucose.input` |
| AMBER force input | `sander_force.in` |
| Validation driver | `validate_sander.py` |
| Tolerances and observed results | `expected_results.json` |

The validation driver creates `work/glucose.pdb` with the 30 Angstrom periodic
cell, along with the AMBER reference topology/restart and generated LAMMPS
data/parameter files.

## Run

Requirements: AmberTools 25, ParmEd, NumPy, SciPy, and LAMMPS with the
MOLECULE, EXTRA-MOLECULE, and KSPACE packages.

```bash
python validate_sander.py \
  --sander "$CONDA_PREFIX/bin/sander" \
  --lammps /path/to/lmp
```

## Converter CLI example

After running the validation command once to prepare `work/glucose.pdb`, run
the standalone GLYCAM conversion directly from this directory:

```bash
amber-to-lammps \
  work/glucose-cli.data \
  work/glucose-cli.parm \
  work/glucose.pdb \
  --topologies ../amber_family_mix/source/glucose.prmtop \
  --counts 1 \
  --charges 0 \
  --provenance work/glucose-cli.provenance.json
```

The source-checkout equivalent is `python ../../amber_to_lammps.py` followed
by the same arguments.

## Python API example

Run this from `validation/glycam_glucose/` after preparing the periodic PDB:

```python
from amber_to_lammps import amber2lammps

result = amber2lammps(
    data_file="work/glucose-api.data",
    param_file="work/glucose-api.parm",
    topologies=["../amber_family_mix/source/glucose.prmtop"],
    molecule_counts=[1],
    pdb_file="work/glucose.pdb",
    charges_target=[0.0],
    provenance_file="work/glucose-api.provenance.json",
)

print(result["summary"])
```

Both examples emit `pair_modify mix arithmetic`, GLYCAM's unscaled 1-4
weights (`special_bonds ... 1`), and the bonded and nonbonded coefficients
extracted from the committed AMBER topology.

The committed tolerances and observed reference values are recorded in
`expected_results.json`. This fixture validates the uniform GLYCAM SCEE/SCNB
profile of 1.0/1.0 independently of the heterogeneous GAFF2 + GLYCAM fixture.
The committed run has a total-energy difference of 0.00018 kcal/mol, maximum
force-component difference of 0.000530 kcal/mol/Angstrom, and force RMSE of
0.000118 kcal/mol/Angstrom.
