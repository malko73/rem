# Papers and manuscript status

The `papers/` directory contains both the current preprint and earlier REM manuscript generations. The files are preserved for research-history and code-to-paper traceability, so filenames alone are **not** sufficient to determine which document is current.

## Current published preprint

**Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization**  
Version: **1.0**  
DOI: [10.5281/zenodo.21900701](https://zenodo.org/records/21900701)

Repository working snapshots:

- [`REM_paper_v0_1.md`](REM_paper_v0_1.md)
- [`REM_paper_v0_1.tex`](REM_paper_v0_1.tex)
- figures in [`figures/`](figures/)

The `v0_1` filename is a retained working filename from the drafting pipeline. **It does not indicate the publication version.** The Markdown working snapshot also still carries a pre-publication `release candidate` status label; that label is a drafting artifact, not the current publication state. The version/DOI metadata and the Zenodo record identify the released preprint as version 1.0. **Zenodo is the authoritative published artifact.**

## Current canonical specification

The canonical human-readable specification is maintained at the repository root:

- [`../REM_spec_v2_2.md`](../REM_spec_v2_2.md)

The corresponding TeX source is:

- [`REM_spec_v2_2.tex`](REM_spec_v2_2.tex)

Older specification files (`REM_spec_v1_1.tex`, `REM_spec_v2_0.tex`, root `REM_spec_v2_1.md`) are retained for historical comparison and are not canonical when they conflict with v2.2.

## Earlier REM manuscripts

The following files document the development of the research program and should be read as historical lineage rather than as the present canonical statement:

- `REM1_source.tex`
- `REM2.tex`
- `REM3.tex`
- `REM5_source.tex`
- `REM_lambda_v2_source.tex`
- `REM_lambda_v5_source.tex`

REM_lambda v5 has its own Zenodo record, DOI [10.5281/zenodo.21880505](https://zenodo.org/records/21880505), and remains the record referenced by the repository-level `CITATION.cff`. The current structural-selection preprint is a separate publication record with DOI 10.5281/zenodo.21900701.

Before citing numerical claims from an earlier manuscript, check [`../CHANGELOG.md`](../CHANGELOG.md) and [`../CURRENT_STATUS.md`](../CURRENT_STATUS.md), because some earlier objectives and numerical interpretations were later falsified, corrected, withdrawn, or re-scoped.

## Review artifacts

The [`reviews/`](reviews/) directory contains internal AI-assisted review and revision material, including files produced with OpenAI- and Gemini-based review workflows.

These files are retained for process transparency, but they are **not independent external peer review** and should not be described as journal referee reports or third-party validation.

`response_to_reviewers.md` likewise records the internal revision process; it does not imply that the preprint has completed external journal peer review.

## Which document should I read?

| Goal | Recommended document |
|---|---|
| Understand the current paper-level result | [`REM_paper_v0_1.md`](REM_paper_v0_1.md) or the Zenodo v1.0 preprint |
| Check the exact current REM definition | [`../REM_spec_v2_2.md`](../REM_spec_v2_2.md) |
| Verify the numerical evidence and caveats | [`../analysis/PHASE_D_FINAL_REPORT.md`](../analysis/PHASE_D_FINAL_REPORT.md) |
| Understand corrections and withdrawn claims | [`../CHANGELOG.md`](../CHANGELOG.md) |
| Reproduce current calculations | [`../REPRODUCING.md`](../REPRODUCING.md) |
| Trace the conceptual development of REM | Earlier REM1/2/3/5/REM_lambda manuscripts |

## Authority rule

For publication metadata and the released paper artifact, **Zenodo is authoritative**. For the current mathematical definition inside this repository, **REM Spec v2.2 is authoritative**. For the tested numerical evidence, **the Phase D validation report and its mapped committed outputs are authoritative within the repository**.
