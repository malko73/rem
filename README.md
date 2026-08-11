# REM Reproducibility Package

> **REM_lambda Version 5 (2026-08-11)** is a corrective revision of
> Version 4, addressing two independent issues found in post-publication
> validation. (1) The continuous optimiser previously built
> `U = exp(i Σ θ_k G_k)` from anti-Hermitian generators, which is NOT
> unitary; it now builds `U = exp(Σ θ_k G_k)`, and all
> continuous-optimisation numbers were recomputed — the earlier
> `Φ ≈ 1.05` / "1.95× improvement" claim is withdrawn, and the quotient
> dimension is corrected 44 → 45. (2) Open-system falsification tests
> (C0–C2) show that `C_H = <H_boundary^2>` does not predict
> factorization-dependent decoherence stability; it is retained only as a
> closed-system structural surrogate. The canonical form is generalised
> to `Φ(F; λ, ρ, L, τ) = I_ρ(F) − λ C_dyn(F; ρ, L, τ)` with the signed
> structural decoherence rate `C_Γ^(0) = Γ_F^exact(0)` as the first
> operational realization. See `CHANGELOG.md` for the full correction
> list. Version 4 is retained as part of the publication history; the
> discrete-cut XY benchmark and `λ* ≈ 0.165` (re-scoped to a surrogate
> crossover) are unchanged.

Reference code and manuscript working snapshots for the **Relational Emergence Model (REM)** research series by Yoshifumi Maruko.

Primary record: **DOI 10.5281/zenodo.21880505** (REM_lambda Version 5)
Zenodo: https://zenodo.org/records/21880505

Previous version: **DOI 10.5281/zenodo.21427776** (REM_lambda Version 4)
Zenodo: https://zenodo.org/records/21427776

Superseded numerical release: **DOI 10.5281/zenodo.21427451** (REM_lambda Version 3)

Earlier record: **DOI 10.5281/zenodo.19642303** (REM1 reproducibility package, REM_lambda v2)  
Zenodo: https://zenodo.org/records/19642303

> The Zenodo records are the authoritative published versions of the manuscripts. Files under `papers/` are repository working snapshots intended to support code-to-paper traceability and ongoing revision.

## Scope

This repository provides an auditable implementation of the finite-dimensional numerical scaffold used to study competition between:

- an informational criterion, represented by cross-boundary mutual information; and
- a dynamical boundary cost, represented by the quadratic boundary-fluctuation measure `C_H = <H_boundary^2>`.

It includes:

- Exact diagonalization for asymmetric XY chains
- Candidate contiguous cuts and lambda-dependent structural selection
- **Continuous factorization manifold optimisation** over a four-parameter submanifold via gradient ascent
- Phase-diagram generation and finite-size scaling for small systems
- A **sign-definite dynamical cost** `C_H = <H_boundary^2> ≥ 0` that eliminates the sign ambiguity documented in earlier versions

## Repository layout

```
src/rem4_numerical.py   Exact-diagonalization scaffold + continuous optimisation
src/rem3.py             Earlier exploratory simulation
papers/                 Working LaTeX snapshots for the REM series
reproduce.sh            Main reproduction command
tests/                  Numerical and output smoke tests (18 tests)
outputs/                Generated figures, excluded from Git
```

## Installation

Python 3.11 or newer is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
pip install -r requirements.txt
```

## Run the reference calculation

```bash
./reproduce.sh
```

Expected generated files include:

```
outputs/REM4_phase_exact.png
outputs/REM4_scaling.png
```

For a faster smoke run:

```bash
python3 src/rem4_numerical.py \
  --outdir outputs \
  --mode all \
  --state-mode ground \
  --j12 1.5 --j23 0.6 --h 0.2 \
  --lambda 0.2 \
  --grid-r 12 --grid-l 12
```

## Tests

```bash
pytest -q
```

18 tests covering Hamiltonian and state invariants, the exact quadratic-cost definition, the distinction between `<H_boundary^2>` and `<H_boundary>^2`, discrete-cut benchmarks, finite-size and output smoke tests, identity-factorisation consistency, factorisation-dependent boundary cost, continuous-optimisation improvement, Adam optimiser correctness, and SGD–Adam divergence.

## Dynamical-cost convention

The code uses the **quadratic dynamical cost** `C_H = <H_boundary^2>`. This is a **sign-definite** measure:

- `C_H ≥ 0` for any state and any cut
- The dynamics-only selection (`λ → ∞`) favours the cut with smallest boundary fluctuation cost
- No sign ambiguity: an antiferromagnetic ground state that gives a negative boundary expectation still yields a positive `C_H`

The earlier `phi_h = -<H_boundary>` convention (used in REM3 and pre-v4 numerical snapshots) is superseded. The quadratic form aligns with REM_lambda v4 and REM5.

## Continuous manifold optimisation

Beyond the site-contiguous partitions, the code searches a four-parameter generator submanifold of the factorisation space. For the three-qubit benchmark (`J12=1.5, J23=0.6, h=0.2, λ=0.2`), **after the 2026-08-11 unitarity fix**:

- Best contiguous-cut Φ: **0.560**
- SGD (30 seeds): **0.581 ± 0.112** (success rate 60%)
- Adam (30 seeds): **0.654 ± 0.191** (success rate 67%)
- The improvement over contiguous is modest at N=3 and strongly depends on the random generator basis.

The full quotient `U(8)/(U(4) ⊗ U(2))` has real dimension **45** (dim U(8)=64 minus 19 = dim of the embedded U(4)×U(2) image, which includes the U(1) kernel: 16+4−1). The reported result is not a global optimum over that full space.

> ⚠️ Superseded (pre-fix) values: best recovered Φ **1.05**, mean **1.09 ± 0.22** (1.95×). Those were artifacts of the non-unitary `exp(iGθ)` map and are no longer reproducible.

## Robustness across system sizes

The `evaluate_factorization` routine has been extended to accept an arbitrary
`n_a` parameter, allowing the optimisation to target any factorisation of the
form `A = sites[0, n_a)`, `B = sites[n_a, n)`.  This enables a systematic
robustness study across N = 3, 4, 5.

Procedure for each N:
1. Determine the best site-contiguous cut `n_a` from the discrete set of candidates.
2. Fix `n_a` to that value so the continuous optimisation explores the same
   Hilbert–Schmidt boundary-reconstruction problem as the discrete benchmark.
3. Run 30 independent seeds for each of two axes:
   - **Axis A** (initial-θ): fixed generator basis, 30 random initial θ.
   - **Axis B** (generator basis): 30 independent random generator bases + θ.

| N | n_a | contiguous Φ | optimised Φ (mean ± σ) | improvement ratio | success rate |
|---|-----|-------------|------------------------|-------------------|-------------|
| 3 | 2   | 0.560       | A: 0.661 ± 0.120 / B: 0.568 ± 0.116 | A: 1.18× / B: 1.01× | A: 87% / B: 43% |
| 4 | 3   | 1.723       | A: 1.719 ± 0.060 / B: 1.700 ± 0.054 | A: 1.00× / B: 0.99× | A: 57% / B: 47% |
| 5 | 4   | 1.458       | A: 1.437 ± 0.049 / B: 1.445 ± 0.038 | A: 0.99× / B: 0.99× | A: 47% / B: 53% |

(A = fixed generator basis, random initial θ; B = random generator basis + θ. 30 trials per axis, post unitarity-fix.)

**2026-08-11 re-verification conclusion**: after the unitarity fix the continuous
optimiser no longer reliably beats the best contiguous cut. Only the N=3
initial-θ axis shows a robust improvement (87% success, 1.18×); the
generator-basis axis and all N=4,5 axes are essentially coin flips (43–57%,
±1–2%). The pre-fix claim "180/180 trials beat contiguous, up to 1.95×" is
**not reproducible** and was an artifact of the non-unitary map. The earlier
statement that scatter shrinks sharply with N no longer holds.

Raw data (180-trial JSON) and violin plots are in `analysis_output/`.

**Caveat:** these results are conditioned on the chosen generator-subspace
dimension, a plain gradient-ascent optimiser (SGD or Adam, lr = 0.01,
steps = 200, no momentum), and the particular Hamiltonian parameters
(`J12=1.5, J23=0.6, h=0.2, λ=0.2`).  They do not represent a global
optimum over the full factorisation manifold. Improving the generator
subspace / initialisation strategy is an open Phase B problem.

## Optimiser comparison: SGD vs Adam

The default optimiser is plain gradient ascent (SGD, lr = 0.01).  A
separate Adam variant (`optimize_factorization_adam`) provides adaptive
moment estimates.  Both are compared under identical starting conditions:
same Hamiltonian, same generator basis, same initial θ, same step budget
(200), across 30 independent seeds per system size.

| N | SGD Φ (mean ± σ) | Adam Φ (mean ± σ) | Adam > SGD | improvement |
|---|-------------------|-------------------|------------|-------------|
| 3 | 0.581 ± 0.112     | **0.654 ± 0.191** | 20/30      | +12.5 %    |
| 4 | 1.709 ± 0.064     | **1.752 ± 0.031** | 22/30      | +2.5 %     |
| 5 | 1.454 ± 0.049     | **1.480 ± 0.028** | 21/30      | +1.8 %     |

Adam outperforms SGD in the majority of trials at all three system sizes,
with a much tighter spread (σ roughly halved at N=4,5), and is the only
optimiser that reliably beats the contiguous-cut Φ at N=4,5 (Adam success
100%/97% vs SGD 53%/60%).  The improvement is modest (≤ 12.5 %) and the
absolute Φ values are far below the pre-fix (non-unitary) numbers.

**Caveat:** the comparison is at fixed step budget (200), not at equal
wall time.  Adam's per-step cost is negligibly higher (two vector
updates); the practical run time is dominated by the evaluation + gradient
overhead shared by both methods.

## Early stopping

The Adam optimiser supports early stopping via `tol_grad` (gradient-norm
threshold) and `tol_phi` (relative Φ-change threshold over a `patience`
window).  When both criteria hold for `patience` consecutive steps the run
terminates early.

With `max_steps = 500` the Adam optimiser does not saturate: at step 200
most seeds have not reached a stationary point, so early stopping under a
200-step budget is ineffective.  Increasing the budget to 500 steps and
setting `tol_grad = 0.2`, `tol_phi = 5 × 10⁻⁴`, `patience = 10` gives:

| Mode      | Φ (mean ± σ) | steps (mean) | step saving |
|-----------|--------------|--------------|-------------|
| Full 500  | 0.690 ± 0.218 | 500          | —           |
| Early     | 0.690 ± 0.218 | 277          | **45 %**    |

(Post unitarity-fix, 2026-08-11. Pre-fix values were 1.251 ± 0.273 / 233 steps / 53 %.)

Across 30 seeds the final Φ is statistically indistinguishable
(max individual degradation < 0.002).  No trial had a quality loss
> 0.01.

For N = 4 and 5 the early-stopping mechanism is available via the same
`optimize_factorization_adam` keyword arguments, but the convergence
behaviour differs.  At a 150-step budget (the largest that runs reliably
on a laptop for 16‑qubit matrices), the Φ‑change criterion
(5 × 10⁻⁴ over 10 steps) is never met, and no step saving is observed.
Larger systems require more steps to settle, and the practical benefit of
early stopping is expected to be smaller than for N = 3.

## Limitations

The current numerical work is an existence-oriented toy-model study. It does not yet provide:

- a guarantee of global optimality on the continuous manifold (finite-difference gradients, fixed generator-basis subspace, SGD or Adam);
- a microscopic derivation of the tradeoff parameter λ;
- a full system-environment decoherence calculation;
- experimental validation.

The robustness study in the previous section is conditioned on the
generator-subspace dimension, a fixed-step gradient-ascent optimiser
(lr = 0.01, steps = 200), and the single Hamiltonian point
(J12=1.5, J23=0.6, h=0.2, λ=0.2).

These limitations are part of the research program and should be retained when citing or extending the code.

## Citation

Please cite the Zenodo record:

```bibtex
@misc{maruko2026rem,
  author       = {Yoshifumi Maruko},
  title        = {Competing Informational and Dynamical Criteria for Subsystem Structure: Relational Emergence Model — Revised Version 5},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21880505},
  url          = {https://doi.org/10.5281/zenodo.21880505}
}
```

## License

CC BY 4.0. See `LICENSE`.
