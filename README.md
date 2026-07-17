# REM Reproducibility Package

Reference code and manuscript sources for the **Relational Emergence Model (REM)** research series by Yoshifumi Maruko.

Primary record: **DOI 10.5281/zenodo.19642303**  
Zenodo: https://zenodo.org/records/19642303

## Scope

This repository provides an auditable implementation of the finite-dimensional numerical scaffold used to study competition between:

- an informational criterion, represented by cross-boundary mutual information; and
- a dynamical boundary criterion, represented by a boundary-Hamiltonian proxy.

It includes exact diagonalization for asymmetric XY chains, candidate contiguous cuts, lambda-dependent structural selection, phase-diagram generation, and finite-size checks for small systems.

## Repository layout

```text
src/rem4_numerical.py   Exact-diagonalization scaffold
src/rem3.py             Earlier exploratory simulation
papers/                 LaTeX sources for the REM series
reproduce.sh            Main reproduction command
tests/                  Numerical and output smoke tests
outputs/                 Generated figures, excluded from Git
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

```text
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

## Reproducibility status and known discrepancy

The code runs deterministically for the ground-state workflow and generates the phase and finite-size figures. However, the current implementation does **not** numerically reproduce every benchmark value printed in the REM4 manuscript table for `J12=1.5`, `J23=0.6`, and `h=0.2`.

In particular, the present code uses:

```python
phi_h = -<psi|H_boundary|psi>
```

For an antiferromagnetic XY ground state, the boundary expectation can be negative, making `phi_h` positive and potentially favoring the stronger bond. This conflicts with the manuscript's stated interpretation that the dynamics-only proxy should favor weaker cross-boundary coupling.

This discrepancy is intentionally documented rather than hidden. It should be resolved by fixing and justifying one convention across the manuscript and code, for example:

- `-abs(<H_boundary>)`;
- a non-negative quadratic cost such as `<H_boundary^2>`; or
- a decoherence-rate proxy derived from an explicit open-system model.

Until that convention is fixed, the repository supports reproduction of the **computational workflow and existence-search scaffold**, not a claim of exact reproduction of all published numerical values.

## Scientific claim boundary

The current numerical work is an existence-oriented toy-model study. It does not yet provide:

- optimization over the full continuous factorization manifold;
- a microscopic derivation of the tradeoff parameter lambda;
- a full system-environment decoherence calculation; or
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
