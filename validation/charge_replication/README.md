# Charge-replication benchmark

This benchmark demonstrates a numerical consequence of preserving a small
residual charge when a molecular template is replicated. It uses the committed
53-atom GAFF2 epon topology and the committed InterMol LAMMPS data file from
the energy cross-check. It does not treat charge preservation as a conversion
error: preserving the source topology is the correct choice for a direct
engine-to-engine comparison. Instead, it documents why AMBER2LAMMPS also
offers an explicit per-template charge target.

Run the benchmark from the repository root:

```bash
python validation/charge_replication/benchmark.py \
  --output validation/charge_replication/results.json
```

The script independently measures four values:

1. the full-precision net charge in `epon.prmtop`;
2. the charge written in InterMol's committed LAMMPS `Atoms` section;
3. AMBER2LAMMPS output when its target is the raw topology charge; and
4. AMBER2LAMMPS output when its target is zero.

It then multiplies each measured per-template value by 1, 10, 100 and 1,000
copies. The result is simple but operationally important: a residual
\(\delta q\) becomes \(N\delta q\) after replicating the template \(N\) times.

| Copies | Raw topology (e) | InterMol output (e) | AMBER2LAMMPS raw target (e) | AMBER2LAMMPS neutral target (e) |
|---:|---:|---:|---:|---:|
| 1 | -0.002996 | -0.002996 | -0.002996 | 0.000000 |
| 10 | -0.029960 | -0.029960 | -0.029960 | 0.000000 |
| 100 | -0.299600 | -0.299600 | -0.299600 | 0.000000 |
| 1,000 | -2.995996 | -2.996000 | -2.995996 | -0.000001 |

AMBER2LAMMPS removes the residual by adding approximately
`+0.0000565282 e` to every atom. This is a user-selected transformation, not
an assumption: validation comparisons request the raw topology charge so that
both converters use identical charges. Large requested corrections generate a
warning, and every correction is recorded in provenance. The remaining
`-0.000001 e` at 1,000 copies comes from the decimal precision used to write
individual charges to the LAMMPS data file, not from the original residual.

The machine-readable measurements are committed in [`results.json`](results.json).
