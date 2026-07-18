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
tests/                  Numerical and output smoke tests (16 tests)
outputs/                Generated figures, excluded from Git
```

## Installation

Python 3.11 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
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
python src/rem4_numerical.py \
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

16 tests covering Hamiltonian and state invariants, the exact quadratic-cost definition, the distinction between `<H_boundary^2>` and `<H_boundary>^2`, discrete-cut benchmarks, finite-size and output smoke tests, identity-factorisation consistency, factorisation-dependent boundary cost, and continuous-optimisation improvement.

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

The full quotient `U(8)/(U(4) ⊗ U(2))` has real dimension 45. The reported result is not a global optimum over that full space.

## Limitations

The current numerical work is an existence-oriented toy-model study. It does not yet provide:

- a guarantee of global optimality on the continuous manifold (finite-difference gradients, fixed generator basis);
- a microscopic derivation of the tradeoff parameter λ;
- a full system-environment decoherence calculation;
- experimental validation.

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
