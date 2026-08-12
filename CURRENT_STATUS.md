# Current status — REM

**Updated:** 2026-08-13  
**Canonical specification:** [REM Spec v2.2](REM_spec_v2_2.md)  
**Primary numerical record in this repository:** [Phase D validation report](analysis/PHASE_D_FINAL_REPORT.md)

This page is the entry point for interpreting the repository. It distinguishes the current canonical formulation, results supported by the checked numerical record, results that were falsified or re-scoped, and work that remains open.

## 1. Current canonical formulation

REM treats the tensor-product structure (F) as a variable selected under informational and dynamical constraints. It is a structural-selection framework; it does not propose a modified quantum dynamics.

The current canonical variational functional is

[
J_{\rm dyn}^{(0)}(F;\rho,\mathcal L)
= -\operatorname{Re}\langle Q_F\rho,Q_F\mathcal L(\rho)\rangle_{\rm HS},
\qquad
\Phi(F;\lambda)=I_\rho(F)-\lambda J_{\rm dyn}^{(0)}(F).
]

The full definitions, sign convention, factorization domain, and dimensional convention are frozen in [Spec v2.2](REM_spec_v2_2.md).

## 2. What the numerical record supports

Under the frozen asymmetric-XY / pure-dephasing protocol reported in [Phase D](analysis/PHASE_D_FINAL_REPORT.md):

- The normalized rate (Gamma_F^{\rm exact}(0)) is singular as an unrestricted global TPS objective. The observed failure mode is a collapse toward a product-state TPS, not a physical optimum.
- The denominator-free (J_{\rm dyn}^{(0)}) is well-posed in the reported D2.2-A/B tests: the reported singularity, product collapse, and seed instability are absent.
- The selected TPS responds to the tested Hamiltonian families and state classes; the report makes the strength and near-degeneracies explicit.
- The finite-time extension (J_{\rm dyn}^{(\tau)}) agrees with (J_{\rm dyn}^{(0)}) as (	au\to0) and exhibits a finite-(	au) crossover between competing TPS basins for the reported cases.
- That crossover is reported across (N=3)–(5), including a reduced-protocol full-quotient closure at (N=5).

The underlying scripts, raw JSON outputs, and commit-level reproduction map are listed in [Phase D §9](analysis/PHASE_D_FINAL_REPORT.md#9-reproducibility).

## 3. What was falsified, corrected, or re-scoped

| Item | Current interpretation | Source |
|---|---|---|
| Normalized (Gamma_F^{\rm exact}(0)) as global objective | **Falsified / prohibited** on unrestricted TPS optimization; retained only as a fixed-(F) local diagnostic. | [Spec §7](REM_spec_v2_2.md) |
| (C_H=\langle H_{\partial F}^2\rangle) as open-system stability predictor | **Re-scoped** to a closed-system structural surrogate. | [CHANGELOG](CHANGELOG.md) |
| v4 continuous-optimization improvement values | **Retracted** after the non-unitary parameterization was corrected. | [CHANGELOG](CHANGELOG.md) |
| (lambda^*\approx0.165) | A structural crossover under the former closed-system surrogate; not an open-system physical crossover. | [CHANGELOG](CHANGELOG.md) |

Falsification and correction are part of the research record, not hidden history. Older specifications and releases remain available for traceability but are not canonical.

## 4. Limits of the current claims

The repository does **not** yet establish:

- a result for arbitrary Hamiltonians, states, environments, or tensor-product dimensions;
- a thermodynamic-limit or universal scaling law (the reported size range is (N=3)–(5));
- a direct physical derivation of every trade-off scale from an environment;
- experimental discriminability or a confirmed laboratory realization;
- independent third-party replication.

The phase report also records a negative result: the observed crossover scale is not simply the global Liouvillian relaxation timescale.

## 5. How to assess or reproduce a claim

1. Read the corresponding section of [Spec v2.2](REM_spec_v2_2.md) to identify the exact definition.
2. Use [Phase D](analysis/PHASE_D_FINAL_REPORT.md) to locate the frozen protocol, gate, script, output, and limitation.
3. Run the test suite with `pytest -q`.
4. For a replication or challenge, use the current D3/D4 tracking issues rather than the archived pre-v2.2 roadmaps.

## 6. Document roles

| Document | Use it for |
|---|---|
| [README](README.md) | Five-minute orientation |
| [This status page](CURRENT_STATUS.md) | Current claims, boundaries, and evidence map |
| [Spec v2.2](REM_spec_v2_2.md) | Canonical theory and definitions |
| [Phase D report](analysis/PHASE_D_FINAL_REPORT.md) | Numerical evidence and exact protocol |
| [CHANGELOG](CHANGELOG.md) | Correction and retraction history |
| [REM next steps](REM_next_steps.md) | Detailed historical and future planning material |
