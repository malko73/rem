# Relational Emergence Model (REM)

> **Current canonical REM: Spec v2.2 (2026-08-12).**  
> REM is a framework for selecting a tensor-product structure (TPS) by balancing mutual information against a dynamical functional. It does **not** modify quantum dynamics.  
> **The canonical objective is denominator-free and well-posed on the unrestricted TPS domain.** The formerly normalized rate is retained only as a fixed-TPS diagnostic after a documented falsification.  
> **What is established numerically:** a finite-time crossover between competing TPS basins is observed for the frozen protocol and persists in the tested (N=3\ldots5) range. It is not a thermodynamic-limit claim or an experimental confirmation.

**Start here:** [current status](CURRENT_STATUS.md) · [canonical specification](REM_spec_v2_2.md) · [Phase D validation report](analysis/PHASE_D_FINAL_REPORT.md) · [correction history](CHANGELOG.md)

## Five-minute map

| Question | Short answer | Evidence |
|---|---|---|
| What is selected? | A bipartite tensor-product structure (F), not a new physical dynamics. | [Spec §1–2](REM_spec_v2_2.md) |
| What is the current objective? | (Phi(F;\lambda)=I_\rho(F)-\lambda J_{\rm dyn}^{(0)}(F)), with (J_{\rm dyn}^{(0)}=-\operatorname{Re}\langle Q_F\rho,Q_F\mathcal L(\rho)\rangle_{\rm HS}). | [Spec §4–5](REM_spec_v2_2.md) |
| What was ruled out? | The normalized (Gamma_F^{\rm exact}(0)) diverges when unrestricted optimization drives initial coherence toward zero; it is not a global TPS objective. | [Spec §7](REM_spec_v2_2.md), [D2.1 report](analysis/PHASE_D_FINAL_REPORT.md#3-the-falsification-normalized-γ_f-is-singular-d21) |
| What does finite time add? | The operational extension (J_{\rm dyn}^{(\tau)}) has the correct (	au\to0) limit and can switch the preferred TPS basin at finite (	au). | [Spec §8–10](REM_spec_v2_2.md) |
| What remains open? | Independent replication, broader models/environments, larger systems, and experimental discriminability. | [current status](CURRENT_STATUS.md) |

## Current claims and boundaries

The repository supports the following bounded claims under its frozen numerical protocol:

1. The normalized instantaneous rate is singular as a global variational objective over unrestricted TPSs.
2. The unnormalized (J_{\rm dyn}^{(0)}) avoids that denominator singularity in the tested optimizations.
3. The optimized TPS responds to the tested Hamiltonian and state choices.
4. Competing TPS basins show a finite-time objective crossing, with (	au_c\approx0.018\)–(0.022) for the reported (N=3\)–(5) studies.

These results are numerical evidence, not a proof for arbitrary Hamiltonians, environments, sizes, or experiments. In particular, the (N=3\)–(5) results do not establish a scaling law or thermodynamic limit.

## Reproduce and inspect

Python 3.11+ is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pytest -q
```

The repository contains the frozen analysis scripts in `analysis/`, their machine-readable outputs in `analysis_output/`, and the consolidated results in [the Phase D report](analysis/PHASE_D_FINAL_REPORT.md). See [CURRENT_STATUS.md](CURRENT_STATUS.md) for the evidence-to-claim map before interpreting any individual result.

## Records and documentation

| Resource | Role |
|---|---|
| [CURRENT_STATUS.md](CURRENT_STATUS.md) | What is canonical, established, falsified, and still open |
| [REM Spec v2.2](REM_spec_v2_2.md) | Canonical definitions, domain, and validation gates |
| [Phase D validation report](analysis/PHASE_D_FINAL_REPORT.md) | Protocol, numerical evidence, caveats, and reproducibility map |
| [CHANGELOG.md](CHANGELOG.md) | Corrections and retractions, including the v4 unitarity artifact |
| [REM next steps](REM_next_steps.md) | Detailed research roadmap and archived planning context |
| [CITATION.cff](CITATION.cff) | Citation metadata |

## Publications and records

| Record | DOI | Content |
|---|---|---|
| **REM Structural Selection — paper v1.0** | [10.5281/zenodo.21900701](https://zenodo.org/records/21900701) | Preprint: *Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization* |
| REM_lambda v5 | [10.5281/zenodo.21880505](https://zenodo.org/records/21880505) | Reference code and manuscript record |
| REM_lambda v4 | [10.5281/zenodo.21427776](https://zenodo.org/records/21427776) | Previous numerical record; see [corrections](CHANGELOG.md) |
| REM1 reproducibility package | [10.5281/zenodo.19642303](https://zenodo.org/records/19642303) | Earlier record |

Zenodo records are the published versions. Repository files are working snapshots and reproducibility materials; their relationship to the published record is described in [CURRENT_STATUS.md](CURRENT_STATUS.md).
