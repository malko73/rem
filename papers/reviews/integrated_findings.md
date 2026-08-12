# Independent Review — Integrated Findings

**Paper**: REM v0.1 (tag `paper-v0.1-review-candidate`)
**Reviews**: A = Gemini 3.6 flash (math/quantum-info lens, OpenAI-compat API),
B = OpenAI gpt-5.4-mini (numerical/referee lens, API)
**Date**: 2026-08-12
**Recommendations**: A = **Minor revision**; B = **Major revision** (1 Critical)

---

## Consolidated findings (deduplicated, severity merged)

| ID | Severity | Source | Finding | Fix |
|---|---|---|---|---|
| F1 | **Critical** | B 2.6 | "derived, not fitted" lacks: ΔΦ₂ values, uncertainty propagation, F_A/F_B sensitivity, τ-grid resolution check | numerical supplement + wording |
| F2 | Major | B 1.1 / A m1 | J₁ notation ambiguous (coefficient vs derivative); rotated-frame notation inconsistent across §2.2/§2.3/§8 | explicit J₁ definition; unify ρ_U notation |
| F3 | Major | B 1.2 | "smoothly" from Cauchy–Schwarz overclaims (Q_FL(ρ) boundedness unstated); "well-posed" needs compactness/continuity argument | add finite-dim boundedness note + compact-quotient existence argument; soften "well-posed" to "non-singular in the tested protocol" with formal note |
| F4 | Major | B 1.3 | "denominator can vanish on any open neighborhood" is a topological overreach | replace with the transitive-orbit argument (U(d) acts transitively on pure states) |
| F5 | Major | B 2.2 | "well-posed" is a mathematical claim, not established by numerics | soften; add existence-on-compact-quotient note; keep numerical evidence as is |
| F6 | Major | B 2.4 | crossover evidence under-specified: 19-point grid not shown; no seed spread/error bars for ΔΦ; A/B classification not independently verifiable; basin stability thin | add full grid + seed-spread error bars (data exists); clarify classification; add basin-stability data |
| F7 | Major | B 2.5 / A m3 | τ_c(5) protocol-sensitive (0.0173 vs 0.0217); 24% spread; "approximately flat" descriptive only | acknowledge sensitivity explicitly; N=5 2-seed note in table caption; frame strictly as finite-N observation |
| F8 | Major | B 3.1 | "finite-size persistence" reads as universality hint | tighten wording (strictly finite-N observation) |
| F9 | Major | B 3.3 | basin identity thin; "0.87 from both reps" undermines basin interpretation | add within-basin d_F stability data; acknowledge the limitation |
| F10 | Major | B 4.1 | replacement functional "numerically convenient" — physical justification ad hoc | add physical argument: J_dyn^(0) = −½ dC²/dt is the instantaneous decay rate of the residual norm (the normalization is what creates the singularity) |
| F11 | Major | B 5.1 | abstract/conclusion overstate ("demonstrate numerically", "battery of tests", "ill-posed") | soften wording |
| F12 | Major | B 5.2 | reproducibility: stopping criteria, best-of-seed vs mean, CIs, exact grid, F_A/F_B procedure, d_F definition, full-quotient parameterization | add a Reproducibility subsection (most details exist in the Validation Report) |
| F13 | Minor | A m1 | §2.2 notation should use rotated frame ρ_U | merge with F2 |
| F14 | Minor | A m2 | ground/mixed near-isometric note | add one sentence in §5 |
| F15 | Minor | A m3 | N=5 2 seeds in table caption | merge with F7 |
| F16 | Minor | B 2.1 | Claim 1 evidence narrow (one state, one size) | acknowledge: counterexample suffices; add "for the tested Haar/N=3 setting" |
| F17 | Minor | B 2.3 | Claim 3 wording too universal | "in the tested cases" phrasing |
| F18 | Minor | B 3.2 | "same state admits different structures" only shown for Haar/thermal | tighten to shown cases |
| F19 | Minor | B 4.2 | Liouvillian comparison only vs slowest mode; subtle spectral explanations not excluded | add the explicit caveat |
| F20 | Minor | B 5.3 | differentiation from prior REM work underdeveloped | add a short "relation to prior REM results" paragraph |

## Numerical supplement (F1) plan

1. ΔΦ₂ values per system (from the quadratic fits already computed in QC-3).
2. Uncertainty propagation: fit coefficient std for ΔΦ₀, ΔΦ₁, ΔΦ₂; τ_c^(1,2) error via first-order propagation.
3. Representative sensitivity: recompute τ_c^(1) with alternative F_B (τ_b = 0.2) and alternative seeds; report the spread.
4. τ-grid resolution: measured τ_c is from linear interpolation between grid points; the interpolation error is bounded by the local grid spacing (report Δτ near τ_c per system).

## Disposition

- **Agree and fix (in-scope, no new physics)**: F1–F12, F13–F20.
- **Disagree / counterargument**: none outright; several B findings (F3/F5) are addressed by adding the compact-quotient existence argument, which strengthens rather than weakens the paper.
- **Not actionable without new experiments** (per master: no new experiments): none — all fixes are wording/notation/numerical-supplement level.
