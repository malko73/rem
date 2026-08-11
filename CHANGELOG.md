# CHANGELOG — REM reproducibility package

## 2026-08-11 — REM_lambda v5 correction (publication package)

This entry freezes the v5 correction package. It supersedes the numerical
claims of REM_lambda v4 (Zenodo 10.5281/zenodo.21427776) on the points
listed below.

### Two independent corrections to v4

1. **Continuous-optimiser unitarity artifact.** The optimiser built
   `U = exp(i Σ θ_k G_k)` from anti-Hermitian generators `G_k`, which is a
   positive Hermitian matrix, not unitary (`||U†U − I|| ≈ 0.30`). Fixed to
   `U = exp(Σ θ_k G_k)` (`||U†U − I|| ≈ 2e-16`). All continuous-optimisation
   numbers in v4 (e.g. `Φ = 1.05`, "1.95× improvement", "180/180 trials beat
   contiguous") are artifacts and are retracted. Post-fix numbers (30 seeds):
   N=3 SGD 0.581±0.112 (60% success), Adam 0.654±0.191 (67%);
   N=4 Adam 1.752±0.031 (100%); N=5 Adam 1.480±0.028 (97%).
   Quotient dimension corrected 44 → 45 (U(1) kernel of the tensor-product
   embedding).
2. **Closed-system second moment is not an open-system stability proxy.**
   Open-system validation (Lindblad, identical fixed environment for both
   cuts) shows: C0 — under uniform pure dephasing `Γ₁^exact = Γ₂^exact =
   2.000000` exactly despite `C_H` differing by ~10×; C1 — under non-uniform
   noise the decoherence-rate ordering follows the environment and flips
   under mirroring (R_Γ: 0.586 vs 1.497), and under amplitude damping /
   mixed noise R_Γ = 1.000; C2 — norm-type Liouvillian boundary costs do not
   track Γ_F. `C_H = <H_∂F²>` is therefore re-scoped to a
   **closed-system structural surrogate**, and `λ* ≈ 0.165` is re-scoped to
   a structural crossover under that surrogate (not an open-system physical
   crossover).

### New framework introduced

- General open-system form:
  `Φ(F; λ, ρ, L, τ) = I_ρ(F) − λ C_dyn(F; ρ, L, τ)`, with `[λ] = T`
  (characteristic timescale).
- First operational realization:
  `C_Γ^(0)(F; ρ, L) = Γ_F^exact(0)`, a **signed** structural decoherence
  rate (negative values occur: 1.25% of 400 random configurations).
- REM Spec v2.0 (frozen 2026-08-11) documents the three-layer dynamical
  functional, the signed naming, and the dimensional analysis.

### Package manifest (frozen)

| Item | Location | Status |
|------|----------|--------|
| REM_lambda v5 source | `papers/REM_lambda_v5_source.tex` | committed |
| REM_lambda v5 PDF | CI artifact (`latex-build` workflow) | 6 pages |
| REM Spec v2.0 source | `papers/REM_spec_v2_0.tex` | committed |
| REM Spec v2.0 PDF | CI artifact (`latex-build` workflow) | 7 pages |
| Reproduction code | `src/`, `analysis/`, `tests/` | commit `a9fb272` |
| Test suite | `pytest -q` | 50 passed |
| LaTeX build gate | `.github/workflows/latex.yml` | green |
| ToyBox canonical sources | `~/Desktop/ToyBox/REM/REM_lambda/REM_lambda_v5.tex`, `~/Desktop/ToyBox/REM/REM_spec_v2_0.tex` | local |

### Zenodo v5 (pending)

Upload as a corrected version with a new DOI; the Description must state
explicitly that this is a correction of v4: unitarity fix and retraction of
old continuous-optimisation numbers, quotient dimension 44→45, `C_H`
re-scoped to closed-system surrogate, open-system canonical framework
introduced, `λ* ≈ 0.165` interpretation changed.
