# Charge-replication benchmark

This benchmark demonstrates the arithmetic consequence of preserving a small
residual charge when a molecular template is copied. It uses the committed
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

It then multiplies each measured per-template value by 1, 10, 100 and 1,000.
Only the one-template systems are converted; the larger totals are arithmetic
projections, not measurements from converted replicated boxes. The result is
simple but operationally important: a residual \(\delta q\) becomes
\(N\delta q\) for \(N\) unchanged copies.

| Copies | Raw topology (e) | InterMol output (e) | AMBER2LAMMPS raw target (e) | AMBER2LAMMPS neutral target (e) |
|---:|---:|---:|---:|---:|
| 1 | -0.002996 | -0.002996 | -0.002996 | 0.000000 |
| 10 | -0.029960 | -0.029960 | -0.029960 | 0.000000 |
| 100 | -0.299600 | -0.299600 | -0.299600 | 0.000000 |
| 1,000 | -2.995996 | -2.996000 | -2.995996 | -0.000001 |

When zero is supplied as the required per-template target, AMBER2LAMMPS removes
the residual by adding approximately `+0.0000565282 e` to every atom. This is a
user-selected transformation, not an assumption: validation comparisons supply
the raw topology charge so that both converters use identical charges. The
converter warns when the requested correction exceeds `0.01 e` per molecule.
The decision is recorded automatically in MCP provenance and in CLI provenance
when `--provenance` is enabled. The remaining
`-0.000001 e` at 1,000 copies comes from the decimal precision used to write
individual charges to the LAMMPS data file, not from the original residual.

InterMol's preservation of the source residual is a legitimate conversion
policy and matches what AMBER would simulate from the same topology. No
automatic InterMol normalization option was found in the version examined.

The machine-readable measurements are committed in [`results.json`](results.json).
