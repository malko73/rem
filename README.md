# Relational Emergence Model (REM)

REM asks a specific structural question in quantum theory: **if the tensor-product structure (TPS) used to define subsystems is treated as a variable, which factorization should be preferred?** The current framework ranks candidate factorizations by balancing informational articulation against a dynamical functional. It is a structural-selection framework; it does **not** modify quantum dynamics.

> **Current canonical specification:** [REM Spec v2.2](REM_spec_v2_2.md)  
> **Current research status:** [CURRENT_STATUS.md](CURRENT_STATUS.md)  
> **Primary numerical record:** [Phase D validation report](analysis/PHASE_D_FINAL_REPORT.md)  
> **Current preprint:** *Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization*, DOI [10.5281/zenodo.21900701](https://zenodo.org/records/21900701)

## How to read this repository

For a first visit, use this order:

1. **This README** — five-minute orientation.
2. **[CURRENT_STATUS.md](CURRENT_STATUS.md)** — what is canonical, supported, falsified, re-scoped, and still open.
3. **[Current preprint](papers/REM_paper_v0_1.md)** — the paper-level argument and reported results. The historical filename is retained for traceability; the manuscript front matter is version 1.0.
4. **[REPRODUCING.md](REPRODUCING.md)** — how to inspect committed evidence or rerun the current validation chain.
5. **[REM Spec v2.2](REM_spec_v2_2.md)** and **[Phase D report](analysis/PHASE_D_FINAL_REPORT.md)** — exact definitions, protocols, gates, outputs, and limitations.

If an older manuscript, script, or roadmap conflicts with the documents above, the **Spec v2.2 + current-status + Phase D/Phase E record takes precedence** for the present formulation.

## Core formulation

The current canonical instantaneous functional is

$$
J_{\mathrm{dyn}}^{(0)}(F;\rho,\mathcal L)
= -\operatorname{Re}\langle Q_F\rho, Q_F\mathcal L(\rho)\rangle_{\mathrm{HS}},
$$

with the structural objective

$$
\Phi(F;\lambda)=I_\rho(F)-\lambda J_{\mathrm{dyn}}^{(0)}(F).
$$

The finite-time extension $J_{\mathrm{dyn}}^{(\tau)}$ is defined in [Spec v2.2 §§8–10](REM_spec_v2_2.md).

## What the current numerical record supports

Under the frozen asymmetric-XY / pure-dephasing protocol documented in the Phase D report:

1. **The formerly normalized instantaneous rate is singular as an unrestricted global TPS objective.** Optimization can drive the coherence denominator toward zero and produce a near-product collapse rather than a physical optimum.
2. **The denominator-free $J_{\mathrm{dyn}}^{(0)}$ removes that singularity in the tested unrestricted optimizations.** The reported product-collapse and seed-instability failure modes disappear in the D2.2 tests.
3. **The selected TPS responds to the tested Hamiltonian and state choices.** This is supported by the D1-R / D2-R cross-evaluation analyses, with documented near-degeneracies and varying separation strength.
4. **Competing TPS basins exhibit a finite-time objective crossing.** The measured crossover values lie at approximately $\tau_c=0.018$–$0.022$ for the reported $N=3$–$5$ cases, including the documented full-quotient $N=5$ closure run. The Phase E addendum further shows that the leading inter-basin expansion $\tau_c^{(1)}=-\Delta\Phi_0/\Delta\Phi_1$ explains the measured scale within about 5–11%; the quadratic fit is a consistency check, not the derivation.

A negative result is also part of the record: the observed crossover scale is **not** simply the global Liouvillian relaxation timescale. See [Phase D §§6–10](analysis/PHASE_D_FINAL_REPORT.md).

## What is not established

The repository does **not** currently establish:

- a universal result for arbitrary Hamiltonians, states, environments, or TPS dimensions;
- a thermodynamic-limit result or scaling law beyond the tested $N=3$–$5$ range;
- an experimentally confirmed laboratory realization;
- independent third-party replication;
- a claim that the finite-dimensional crossover is a thermodynamic phase transition.

The intended wording is **finite-time structural crossover** or **first-order-like basin transition in the tested realization**, not phase transition.

## Repository map

| Path | Role | Status |
|---|---|---|
| `CURRENT_STATUS.md` | Current claims, boundaries, and evidence map | **Start here** |
| `REM_spec_v2_2.md` | Canonical definitions and validation gates | **Canonical** |
| `analysis/PHASE_D_FINAL_REPORT.md` | Consolidated D0–D4 evidence plus Phase E addendum | **Primary evidence** |
| `analysis/` | Validation and exploratory scripts | See [analysis index](analysis/README.md) |
| `analysis_output/` | Committed machine-readable numerical outputs | Evidence archive |
| `papers/` | Current preprint, specifications, historical manuscripts, figures, review artifacts | See [papers index](papers/README.md) |
| `src/` | Earlier numerical scaffold and quotient-geometry utilities | Mixed current/legacy support code |
| `tests/` | Regression and numerical tests | Automated checks |
| `REM_next_steps.md` | Long-form research roadmap | Planning / historical context |
| `CHANGELOG.md` | Corrections, withdrawals, and re-scoping history | Traceability |

### Important legacy note

`reproduce.sh` and `src/rem4_numerical.py` reproduce the **earlier REM4 asymmetric-XY benchmark**, not the current Spec v2.2 Phase D validation chain. They are retained for historical reproducibility. For the current formulation, use [REPRODUCING.md](REPRODUCING.md).

The CI smoke test likewise exercises the earlier REM4 scaffold in addition to the test suite; a green CI result should not be interpreted as a full rerun of all Phase D/D4/Phase E numerical experiments.

## Active independent checks

Current follow-up work is intentionally narrow:

- [Issue #5 — D3 independent finite-time crossover reproduction](https://github.com/malko73/rem/issues/5)
- [Issue #4 — D4 finite-size persistence stress test](https://github.com/malko73/rem/issues/4)

Older roadmap Issues #1 and #2 are closed as superseded.

## Reproduce and inspect

Python 3.11+ is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pytest -q
```

The test suite is the fast integrity check. Full or partial numerical reproduction is described separately in [REPRODUCING.md](REPRODUCING.md), with the exact script/output map tied to the Phase D report.

## Publications and records

| Record | DOI | Role |
|---|---|---|
| **REM Structural Selection — paper v1.0** | [10.5281/zenodo.21900701](https://zenodo.org/records/21900701) | Current preprint: *Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization* |
| REM_lambda v5 | [10.5281/zenodo.21880505](https://zenodo.org/records/21880505) | Earlier reference-code/manuscript record; repository-level `CITATION.cff` currently points here |
| REM_lambda v4 | [10.5281/zenodo.21427776](https://zenodo.org/records/21427776) | Previous numerical record; see [CHANGELOG](CHANGELOG.md) before citing its numerical claims |
| REM1 reproducibility package | [10.5281/zenodo.19642303](https://zenodo.org/records/19642303) | Earlier research record |

Zenodo records are the authoritative published artifacts. Repository manuscript files are working snapshots used for code-to-paper traceability and ongoing research organization.

## Review-material note

Files under `papers/reviews/` are **internal AI-assisted review and revision artifacts**. They are preserved for process transparency but must not be interpreted as independent external peer review. See [papers/README.md](papers/README.md).

## Citation

For the current structural-selection preprint, cite DOI **10.5281/zenodo.21900701**. The repository-level [CITATION.cff](CITATION.cff) presently identifies the REM_lambda v5 reproducibility record (DOI **10.5281/zenodo.21880505**); the two records therefore serve different citation targets.
