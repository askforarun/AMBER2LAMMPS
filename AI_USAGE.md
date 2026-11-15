# AI usage disclosure draft

Generative AI has assisted development of AMBER2LAMMPS and must be disclosed in
any JOSS submission. This draft records the assistance known from the current
development work; the authors must add any other tools or sessions before
submission.

## Tools and scope

- **OpenAI Codex (GPT-5 family; exact hosted model revision not exposed):**
  assisted with code review, implementation and refactoring suggestions, test
  scaffolding, CI configuration, validation scripts, documentation drafting,
  release preparation and repository-readiness checks.

AI suggestions were used in converter code, regression tests, GitHub Actions,
the GAFF2/GLYCAM `sander` validation bundle, and project documentation.

## Human responsibility

The human author reviewed and edited the AI-assisted outputs, ran the automated
tests, executed the independent `sander` and LAMMPS comparisons, inspected the
reported numerical differences, and retained responsibility for scientific
design decisions, release actions, accuracy, originality and licensing.

## Before submission

The submitting authors must verify this disclosure against their complete use
of AI tools, add product/model versions where known, and explicitly approve the
final wording. AI must not be used to conduct conversational interactions with
JOSS editors or reviewers except where JOSS policy permits translation.
