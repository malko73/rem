# REM Reproducibility Package

> **REM_lambda Version 4 (2026-07-18)** corrects the numerical implementation
> distributed with Version 3. The dynamical cost now uses
> `C_H = <H_boundary^2>` and continuous optimisation reconstructs the
> factorisation-dependent boundary interaction. Version 4 supersedes the
> numerical values reported in Version 3; the theoretical REM functional is unchanged.

Reference code and manuscript working snapshots for the **Relational Emergence Model (REM)** research series by Yoshifumi Maruko.

Primary record: **DOI 10.5281/zenodo.21427776** (REM_lambda Version 4)
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

Beyond the site-contiguous partitions, the code searches a four-parameter generator submanifold of the factorisation space. For the three-qubit benchmark (`J12=1.5, J23=0.6, h=0.2, λ=0.2`):

- Best contiguous-cut Φ: **0.560**
- Best recovered factorisation Φ: **1.05**
- Mean Φ over 30 random seeds: **1.09 ± 0.22** (1.95× improvement)
- No seed failed to improve over the contiguous best.

The full quotient `U(8)/(U(4) ⊗ U(2))` has real dimension 44. The reported result is not a global optimum over that full space.

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
| 3 | 2   | 0.560       | 1.09 ± 0.22            | 1.95×             | 100% (60/60) |
| 4 | 3   | 1.723       | 1.88 ± 0.03            | 1.09×             | 100% (60/60) |
| 5 | 4   | 1.458       | 1.62 ± 0.03            | 1.11×             | 100% (60/60) |

Across all 180 trials the continuous optimisation always beat the best
contiguous-cut Φ.  The improvement is large when the discrete best is low
(N=3: Φ~0.56 → ~1.09) and settles to a modest but consistent ~1.1× for
N=4,5.  Scatter shrinks sharply with N: the standard deviation drops from
22 % (N=3) to under 2 % (N=4,5), indicating that the random subspace
dimension we use (n_params = 4 for N=3, 16 for N=4, 20 for N=5) captures
a narrower fraction of the full manifold as the system grows.

Raw data (180-trial JSON) and violin plots are in `analysis_output/`.

**Caveat:** these results are conditioned on the chosen generator-subspace
dimension, a plain gradient-ascent optimiser (SGD, lr = 0.01, steps = 200,
no momentum), and the particular Hamiltonian parameters
(`J12=1.5, J23=0.6, h=0.2, λ=0.2`).  They do not represent a global
optimum over the full factorisation manifold.

## Optimiser comparison: SGD vs Adam

The default optimiser is plain gradient ascent (SGD, lr = 0.01).  A
separate Adam variant (`optimize_factorization_adam`) provides adaptive
moment estimates.  Both are compared under identical starting conditions:
same Hamiltonian, same generator basis, same initial θ, same step budget
(200), across 30 independent seeds per system size.

| N | SGD Φ (mean ± σ) | Adam Φ (mean ± σ) | Adam > SGD | improvement |
|---|-------------------|-------------------|------------|-------------|
| 3 | 1.083 ± 0.235     | **1.232 ± 0.269** | 30/30      | +13.8 %    |
| 4 | 1.883 ± 0.034     | **2.069 ± 0.036** | 30/30      | +9.9 %     |
| 5 | 1.624 ± 0.026     | **1.861 ± 0.056** | 30/30      | +14.6 %    |

Adam outperforms SGD in **all 90 trials** across all three system sizes.
The improvement is largest for N=5 (+14.6 %) where the manifold dimension
is highest and the adaptive step sizes provide the greatest benefit.
Adam's final gradient norm is typically 2–3× smaller than SGD's,
indicating tighter convergence at the same step count.

Both optimisers maintain a 100 % success rate — every trial beats the
best contiguous-cut Φ.

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
| Full 500  | 1.251 ± 0.273 | 500          | —           |
| Early     | 1.249 ± 0.272 | 233          | **53 %**    |

Across 30 seeds the final Φ is statistically indistinguishable
(max individual degradation < 0.002).  No trial had a quality loss
> 0.01.  The same mechanism is available for N = 4, 5 via the
`optimize_factorization_adam` keyword arguments.

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
  title        = {Variational Generation of Relational Structure in the Relational Emergence Model: Quadratic Dynamical Cost and Regime Competition},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21427776},
  url          = {https://doi.org/10.5281/zenodo.21427776}
}
```

## License

CC BY 4.0. See `LICENSE`.
