# JOSS readiness evidence

This document maps repository evidence to JOSS pre-review expectations. It is
not a claim that the project is ready to submit; unresolved gates remain
explicit.

| Expectation | Evidence | Status |
|---|---|---|
| OSI-approved open source license | `LICENSE` (MIT) | Met |
| Public source, issues and contribution path | GitHub repository; `CONTRIBUTING.md` | Met |
| Installable language-standard package | `pyproject.toml`; PyPI releases | Met |
| Tests and continuous integration | Python matrix, packaging, whitespace and required `sander` validation jobs | Met |
| Scientific correctness procedures | `validation/amber_family_mix/`; `validation/glycam_glucose/`; `validation/amber_protein_water_ions/`; `docs/support-matrix.md` | Met for standalone GLYCAM06j-1 glucose, GAFF2 + GLYCAM, and ff14SB + TIP3P + NaCl fixtures |
| Agent-ready claim | `amber_to_lammps_mcp.py`; MCP tests; `docs/agent-readiness.md` | Interface and safety properties implemented; reproducible context-size benchmark still required before submission |
| Clearly bounded research scope | `docs/support-matrix.md`; `docs/scientific-design.md` | Met |
| Release history and changelog | `CHANGELOG.md`; annotated `v1.0.0` and `v1.1.0` releases | Met |
| Contribution, governance and support expectations | `CONTRIBUTING.md`; `GOVERNANCE.md` | Met |
| More than six months of public iterative development | Public Git history beginning September 2025 | Met by repository age; editors make the final assessment |
| Demonstrated research impact/use | `docs/research-use.md` | **Blocked: real use must be documented** |
| Feature-complete claimed scope | Support matrix | Validated for the specific systems listed, including a non-covalent peptide/GAFF2 ligand fixture; larger folded proteins and protein/carbohydrate mixtures remain planned and must not be claimed |
| JOSS manuscript and bibliography | `paper/paper.md`, `paper/paper.bib` | Draft present; research-impact and AI-disclosure gates remain |
| Accurate citation metadata | `CITATION.cff` | Version, date and author-supplied ORCID recorded |
| AI usage disclosure | `AI_USAGE.md` | Draft present; authors must confirm completeness before submission |
| Independent-user installation test | Research-use record | **Missing** |

## Submission gate

Do not submit until the research-use, manuscript, citation metadata and
independent-user rows are resolved. Additional protein/solvent/ion validation
is required only if the paper claims those systems as supported; alternatively,
the paper must use the narrower validated scope shown in the support matrix.
