# Agent readiness

AMBER2LAMMPS uses **agent-ready** in a narrow, testable sense. The converter
does not contain a language model and does not autonomously select a force
field or simulation protocol. It exposes a structured, bounded interface that
an external software agent can use without loading complete molecular files or
LAMMPS logs into its context.

## Interface properties

- The MCP `convert` tool accepts typed conversion inputs and returns a compact
  summary by default. Full provenance is available only with `detail="full"`.
- Every conversion records input and output hashes, software versions,
  force-field scaling decisions, warnings, and the LAMMPS execution contract
  on disk without requiring another tool call.
- Compact warnings use stable codes and numeric fields. Output paths are
  relative to one confined workspace.
- `check_lammps_run0` verifies output hashes, invokes only the executable
  configured by the server operator, enforces a timeout, and returns compact
  energy terms or one error. The complete log remains on disk.
- Existing outputs require explicit overwrite permission. Conversion and
  validation calls are serialized to prevent cross-request output mixing.
- Unsupported scientific cases fail before output rather than inviting an
  agent to guess a replacement model.

## Evaluation required for publication

The JOSS manuscript should report a reproducible comparison between:

1. an agent workflow that reads generated data, parameter, provenance, and log
   files; and
2. the default MCP workflow that consumes only structured summaries.

For each validation fixture, record atom count, artifact bytes, estimated
context tokens using one documented tokenizer, MCP response tokens, whether
the workflow completed without opening artifacts, and validation outcome.
Include malformed-input, modified-hash, timeout, overwrite, and concurrent-call
tests. Report context reduction as an observed engineering result, not as an
improvement in scientific accuracy.

## Claim boundary

Agent readiness improves automation, traceability, context efficiency, and
operational safety. Scientific correctness continues to come from explicit
force-field handling, regression tests, and numerical comparison with AMBER
`sander`.
