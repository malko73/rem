# REM Reproducibility Package

Reference code and manuscript working snapshots for the **Relational Emergence Model (REM)** research series by Yoshifumi Maruko.

Primary record: **DOI 10.5281/zenodo.21427451** (REM_lambda Version 3)  
Zenodo: https://zenodo.org/records/21427451

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
- **Continuous factorization manifold optimisation** over `U(N)/(U(n_A) × U(n_B))` via gradient ascent
- Phase-diagram generation and finite-size scaling for small systems
- A **sign-definite dynamical cost** `C_H = <H^2> ≥ 0` (REM_lambda_v2 convention) that eliminates the sign ambiguity documented in earlier versions

## Repository layout

```
src/rem4_numerical.py   Exact-diagonalization scaffold + continuous optimisation
src/rem3.py             Earlier exploratory simulation
papers/                 Working LaTeX snapshots for the REM series
reproduce.sh            Main reproduction command
tests/                  Numerical and output smoke tests (9 tests)
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

9 tests covering: Hamiltonian hermiticity, ground-state normalisation, mutual-information non-negativity, dynamical-cost non-negativity (`C_H ≥ 0`), crossover-λ* benchmark, info/std divergence at λ=0, continuous optimisation Φ improvement, finite-size scan smoke, and phase-scan output existence.

## Dynamical-cost convention

The code uses the **quadratic dynamical cost** `C_H = <H_boundary^2>` (REM_lambda_v2 convention). This is a **sign-definite** measure:

- `C_H ≥ 0` for any state and any cut
- The dynamics-only selection (`λ → ∞`) favours the cut with smallest boundary fluctuation cost
- No sign ambiguity: an antiferromagnetic ground state that gives a negative boundary expectation still yields a positive `C_H`

The earlier `phi_h = -<H_boundary>` convention (used in REM3, REM4 pre-v2) is superseded. The quadratic form aligns with REM_lambda_v2 and REM5.

## Continuous manifold optimisation

Beyond the site-contiguous partitions, the code can search over the full factorisation manifold `U(N)/(U(n_A) × U(n_B))` via gradient ascent on anti-Hermitian generators. For the 3-qubit benchmark (`J12=1.5, J23=0.6, h=0.2, λ=0.2`):

- Best contiguous-cut Φ: **0.68**
- Continuous manifold Φ: **1.32**

This indicates that the true optimal subsystem decomposition can lie outside any site-contiguous partition.

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
  title        = {Relational Emergence Model: A Generative Extension of Relational Quantum Mechanics},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19642303},
  url          = {https://doi.org/10.5281/zenodo.19642303}
}
```

## License

CC BY 4.0. See `LICENSE`.
