# GMSO AMBER 1-4 scaling audit

This folder characterizes whether GMSO preserves topology-native AMBER SCEE
and SCNB values when importing a parameterized `parmed.Structure`. It is a
version-specific interoperability audit, not a claim about every GMSO release.

The audit uses three existing repository fixtures:

- standalone GLYCAM06j-1 glucose (`SCEE=1.0`, `SCNB=1.0`);
- ff14SB peptide plus GAFF2 ethanol (`SCEE=1.2`, `SCNB=2.0`); and
- a combined GAFF2 plus GLYCAM reference containing both active profiles.

## Exact commands used

The committed `results.json` was produced on macOS arm64 with the official
GMSO 0.18.0 Conda package. A temporary prefix kept the test isolated from the
developer's existing GMSO checkout:

```bash
cd /path/to/AMBER2LAMMPS
conda create -y -p /private/tmp/gmso-0.18-conda \
  -c conda-forge gmso=0.18.0 parmed
conda install -y -p /private/tmp/gmso-0.18-conda \
  -c conda-forge mbuild 'setuptools<81'
/private/tmp/gmso-0.18-conda/bin/python \
  validation/gmso_scaling/audit_gmso_scaling.py \
  --output /private/tmp/gmso-0.18-results.json
```

The extra packages are needed because importing `gmso.external` in the 0.18.0
Conda environment otherwise fails first because `mbuild` is absent and then
because mBuild 0.10.9 imports `pkg_resources`. The audit contains local
compatibility aliases for that mBuild release's use of removed Python/NumPy
aliases; it does not alter GMSO or mBuild.

To rerun in any working environment containing GMSO, ParmEd, and NumPy:

```bash
cd AMBER2LAMMPS
python validation/gmso_scaling/audit_gmso_scaling.py
```

To save a new report:

```bash
python validation/gmso_scaling/audit_gmso_scaling.py \
  --output gmso-scaling-results.json
```

Always report the GMSO version printed in the JSON. Do not compare results
across versions without rerunning the audit.

## Observed behavior

With GMSO 0.18.0, all three ParmEd inputs were imported with the same default
global arrays: `[0.0, 0.0, 0.5]` for Lennard-Jones and electrostatics. The
topologies had no molecule-specific scaling factors. Consequently:

- GLYCAM's 1.0/1.0 weights were not retained;
- the ff14SB/GAFF2 Coulomb weight of `1/1.2` was not retained; and
- the two profiles in the mixed reference were collapsed to one global value.

The same result was observed in a preliminary run with GMSO 0.14.0, but the
committed JSON records only the requested latest-release audit. This test does
not establish that every GMSO version behaves this way. It also
does not imply that GMSO cannot store scaling values: it tests the specific
ParmEd import route used here. Separately, a LAMMPS data file cannot contain
the run-level `special_bonds` command, so a complete LAMMPS workflow must set
those controls in an input script or an equivalent explicit representation.
