# Current status — REM

**Updated:** 2026-08-13  
**Canonical specification:** [REM Spec v2.2](REM_spec_v2_2.md)  
**Primary numerical record:** [Phase D validation report](analysis/PHASE_D_FINAL_REPORT.md), including the Phase E addendum  
**Current preprint:** [Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization](papers/REM_paper_v0_1.md), DOI [10.5281/zenodo.21900701](https://zenodo.org/records/21900701)

This page is the authoritative repository-level map for interpreting the present REM formulation. It separates canonical definitions from numerical evidence, historical results, corrections, and open replication work.

## 1. Status labels used in this repository

- **Canonical** — the definition currently used for REM.
- **Supported in the tested protocol** — numerical evidence exists under the explicitly frozen finite-system protocol; this is not a universal proof.
- **Falsified / prohibited** — a proposed formulation failed its intended role and must not be used for that role.
- **Re-scoped** — a quantity remains useful, but only under a narrower interpretation than previously claimed.
- **Historical** — retained for traceability but not authoritative for the current formulation.
- **Open** — requires independent replication, extension, or experiment.

## 2. Current canonical formulation

REM treats the tensor-product structure $F$ as a variable selected under informational and dynamical constraints. It is a structural-selection framework; it does not propose a modified quantum dynamics.

The current instantaneous dynamical functional is

$$
J_{\mathrm{dyn}}^{(0)}(F;\rho,\mathcal L)
= -\operatorname{Re}\langle Q_F\rho,Q_F\mathcal L(\rho)\rangle_{\mathrm{HS}},
$$

and the variational objective is

$$
\Phi(F;\lambda)=I_\rho(F)-\lambda J_{\mathrm{dyn}}^{(0)}(F).
$$

The full definitions, sign convention, factorization domain, dimensional convention, and finite-time extension $J_{\mathrm{dyn}}^{(\tau)}$ are frozen in [Spec v2.2](REM_spec_v2_2.md).

## 3. Current evidence map

| Result | Status | Primary evidence |
|---|---|---|
| Normalized $\Gamma_F^{\mathrm{exact}}(0)$ becomes singular under unrestricted global TPS optimization | **Falsified as a global objective** | [Phase D §3](analysis/PHASE_D_FINAL_REPORT.md#3-the-falsification-normalized-γ_f-is-singular-d21), [Spec §7](REM_spec_v2_2.md) |
| Denominator-free $J_{\mathrm{dyn}}^{(0)}$ avoids the normalized-rate singularity, product collapse, and seed instability in the tested runs | **Supported in tested protocol** | [Phase D §4](analysis/PHASE_D_FINAL_REPORT.md#4-j_dyn0-well-posedness-d22-ab) |
| Selected TPS responds to tested Hamiltonian families and state classes | **Supported in tested protocol, with documented near-degeneracies** | [Phase D §5](analysis/PHASE_D_FINAL_REPORT.md#5-structural-response-d1-r-d2-r) |
| Finite-time objective crossing between competing TPS basins exists in the reported cases | **Supported in tested protocol** | [Phase D §6](analysis/PHASE_D_FINAL_REPORT.md#6-finite-time-structural-crossover-d3-commit-8c1a236) |
| Crossover persists across the tested $N=3$–$5$ finite systems, including the reduced-protocol full-quotient $N=5$ closure run | **Finite-size persistence only; no scaling-law claim** | [Phase D §7](analysis/PHASE_D_FINAL_REPORT.md#7-finite-size-persistence-d4-commit-78eb500) |
| Leading inter-basin expansion $\tau_c^{(1)}=-\Delta\Phi_0/\Delta\Phi_1$ explains the measured crossover scale within about 5–11% | **Supported mechanism result** | [Phase D §10](analysis/PHASE_D_FINAL_REPORT.md#10-phase-e-addendum--τ_c-mechanism-2026-08-12-commit-08beebb) |
| Quadratic approximation reduces the crossover residual to about 0.0–0.5% | **Consistency check, not the derivation** | [Phase D §10](analysis/PHASE_D_FINAL_REPORT.md#10-phase-e-addendum--τ_c-mechanism-2026-08-12-commit-08beebb) |
| $\tau_c$ is set by the global Liouvillian relaxation timescale | **Rejected for the tested case** | [Phase D §§6,10](analysis/PHASE_D_FINAL_REPORT.md) |

The committed scripts and raw JSON outputs are mapped in [Phase D §9](analysis/PHASE_D_FINAL_REPORT.md#9-reproducibility) and summarized for new readers in [analysis/README.md](analysis/README.md).

## 4. What was falsified, corrected, or re-scoped

| Item | Current interpretation | Source |
|---|---|---|
| Normalized $\Gamma_F^{\mathrm{exact}}(0)$ as a global TPS objective | **Falsified / prohibited** on unrestricted TPS optimization; retained only as a fixed-$F$ local diagnostic | [Spec §7](REM_spec_v2_2.md) |
| $C_H=\langle H_{\partial F}^2\rangle$ as an open-system stability predictor | **Re-scoped** to a closed-system structural surrogate | [CHANGELOG](CHANGELOG.md) |
| v4 continuous-optimization improvement values | **Retracted** after correction of the non-unitary parameterization | [CHANGELOG](CHANGELOG.md) |
| Historical $\lambda^*\approx0.165$ result | Structural crossover under the former closed-system surrogate; **not** an open-system physical crossover | [CHANGELOG](CHANGELOG.md) |

Falsification and correction are part of the research record rather than hidden history. Older specifications, releases, and scripts remain available for traceability but are not canonical when they conflict with Spec v2.2.

## 5. Limits of the current claims

The repository does **not** yet establish:

- a result for arbitrary Hamiltonians, states, environments, or tensor-product dimensions;
- a thermodynamic-limit or universal scaling law; the reported size range is $N=3$–$5$;
- experimental discriminability or a confirmed laboratory realization;
- independent third-party replication;
- that the finite-dimensional sharp crossover is a thermodynamic phase transition;
- that every trade-off scale in REM has a unique microscopic derivation from an environment.

The preferred language for the present finite-dimensional result is **finite-time structural crossover** or **first-order-like basin transition in the tested realization**.

## 6. Independent checks currently open

The active follow-up issues are deliberately narrower than the earlier roadmaps:

- [#5 — independently reproduce the finite-time TPS crossover](https://github.com/malko73/rem/issues/5)
- [#4 — stress-test finite-size persistence of the TPS crossover](https://github.com/malko73/rem/issues/4)

Issues #1 and #2 are closed as superseded historical roadmaps.

## 7. How to assess or reproduce a claim

1. Identify the exact definition in [Spec v2.2](REM_spec_v2_2.md).
2. Locate the corresponding protocol, gate, and limitation in [Phase D](analysis/PHASE_D_FINAL_REPORT.md).
3. Use [REPRODUCING.md](REPRODUCING.md) to distinguish a fast integrity check from a full numerical rerun.
4. Inspect committed raw outputs in `analysis_output/` before rerunning expensive calculations.
5. For replication or challenge work, record environment, commit, seeds, quotient coverage, and any deviation from the frozen protocol.

## 8. Document roles

| Document | Use it for |
|---|---|
| [README](README.md) | Five-minute orientation |
| [This status page](CURRENT_STATUS.md) | Current claims, boundaries, and evidence map |
| [Spec v2.2](REM_spec_v2_2.md) | Canonical theory and definitions |
| [Phase D report](analysis/PHASE_D_FINAL_REPORT.md) | Numerical evidence and exact protocol |
| [REPRODUCING.md](REPRODUCING.md) | Reproduction paths and command map |
| [analysis/README.md](analysis/README.md) | Current vs historical analysis scripts |
| [papers/README.md](papers/README.md) | Current preprint vs historical manuscript map |
| [CHANGELOG](CHANGELOG.md) | Correction, withdrawal, and re-scoping history |
| [REM next steps](REM_next_steps.md) | Long-form planning and historical roadmap material |

## 9. Publication and review interpretation

The Zenodo record DOI [10.5281/zenodo.21900701](https://zenodo.org/records/21900701) is the published v1.0 preprint artifact for the current structural-selection paper. Repository manuscript files are working snapshots for traceability.

Files under `papers/reviews/` are internal AI-assisted review and revision artifacts. They are not independent external peer review and should not be represented as such.
