# Numerical output archive

`analysis_output/` stores committed machine-readable outputs, tables, and selected figures generated during several stages of REM validation. It is an **evidence archive**, not a directory in which every file represents a current canonical claim.

## Current primary evidence

For the present Spec v2.2 paper claims, start with these output families:

- `d2_1_haar_singularity_audit.json` — normalized-objective singularity audit.
- `d2_2a_unnormalized_instantaneous.json` — denominator-free instantaneous functional checks.
- `d2_2b_tau_*.json`, `d2_2b_gates.json` — finite-time consistency and gate checks.
- `d1r_canonical_revalidation.json` — Hamiltonian-response revalidation.
- `d2r_state_revalidation.json` — state-response revalidation.
- `d3_timescale.json` — fixed-representative finite-time crossover analysis.
- `d4_*.json`, `d4_finite_size.json` — finite-size persistence analyses and controls.
- `e1_tauc_mechanism.json` — short-time inter-basin explanation of the crossover scale.
- `e1_sensitivity.json` — supporting sensitivity analysis.

The authoritative interpretation of those files is the consolidated [`../analysis/PHASE_D_FINAL_REPORT.md`](../analysis/PHASE_D_FINAL_REPORT.md), including its Phase E addendum.

## Earlier outputs

Files prefixed by `b*`, `c*`, `gamma_F*`, `m2_vs_var*`, `optimiser_compare*`, `early_stop*`, and `robustness*` preserve earlier experiments, proxy tests, optimizer studies, and the falsification path that preceded the current canonical formulation.

They remain useful for auditability, but older output values should **not** be promoted to current REM claims without checking:

- [`../CURRENT_STATUS.md`](../CURRENT_STATUS.md)
- [`../CHANGELOG.md`](../CHANGELOG.md)
- [`../REM_spec_v2_2.md`](../REM_spec_v2_2.md)

## Why outputs are committed

The archive allows a third party to:

1. inspect the numerical values used in reports without immediately rerunning expensive optimizations;
2. compare a fresh reproduction against the recorded result;
3. audit negative results and superseded approaches rather than seeing only successful runs.

A committed output is evidence that a recorded run produced those data; it is not, by itself, proof of generality or independent replication.

## Reproduction map

For commands and the current claim-to-script/output mapping, see [`../REPRODUCING.md`](../REPRODUCING.md). For the exact frozen protocol and commit mapping, see [Phase D §9](../analysis/PHASE_D_FINAL_REPORT.md#9-reproducibility).
