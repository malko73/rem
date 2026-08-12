# Source-code roles

The `src/` directory contains both reusable quotient-space infrastructure and earlier REM numerical scaffolds. File age does not by itself determine whether a component is used by the current validation chain.

## Quotient-space infrastructure

- `quotient_geometry.py` — quotient/tensor-factorization geometry utilities used by the validation analyses.
- `quotient_optimizer.py` — optimization utilities used by multiple analysis scripts.

Current D-series scripts under `../analysis/` import this infrastructure. The canonical equations and interpretation remain defined by [`../REM_spec_v2_2.md`](../REM_spec_v2_2.md), not by comments in an individual implementation file.

## Earlier numerical scaffolds

- `rem4_numerical.py` — earlier REM4 asymmetric-XY exact-diagonalization / structural-selection scaffold.
- `rem3.py` — earlier exploratory REM simulation.

These files are retained for historical reproducibility. In particular, the root `../reproduce.sh` executes `rem4_numerical.py`; that path does **not** reproduce the current Spec v2.2 Phase D / Phase E validation chain.

For current claim-to-script mapping, use [`../REPRODUCING.md`](../REPRODUCING.md) and [`../analysis/README.md`](../analysis/README.md).

## Authority rule

Implementation code is evidence-bearing machinery, not the specification itself. If implementation comments, earlier scripts, or historical manuscripts conflict with the current definition, use this order:

1. [`../REM_spec_v2_2.md`](../REM_spec_v2_2.md) — canonical formulation;
2. [`../analysis/PHASE_D_FINAL_REPORT.md`](../analysis/PHASE_D_FINAL_REPORT.md) — validated protocol/results;
3. current D/E analysis scripts and committed outputs;
4. earlier REM3/REM4 scaffolds — historical context.
