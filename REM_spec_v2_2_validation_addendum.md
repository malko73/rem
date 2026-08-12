# REM Spec v2.2 Validation Addendum — Phase D Results

**Version**: 1.0
**Date**: 2026-08-12
**Status**: Independent companion document. The Spec v2.2 text itself is
**frozen and unchanged** — this addendum maps the Phase D numerical evidence
onto the Spec's sections. Separation of theory specification from post-hoc
numerical evidence is deliberate: later results are not folded back into the
definitional document.

**Companion**: `analysis/PHASE_D_FINAL_REPORT.md` (full D0–D4 validation
report). Spec: `REM_spec_v2_2.md` / `papers/REM_spec_v2_2.tex`.

---

## 1. Evidence map (Spec section ↔ Phase D result)

| Spec v2.2 section | Content | Phase D evidence |
|---|---|---|
| §4 | Canonical J_dyn^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩ | D2.2-A: well-posed (G1–G6); no blow-up as C_F²→0 (min C_F² = 0.096, max\|γ_ref\| = 1.66) |
| §5 | Φ = I − λJ_dyn^(0) | D1-R / D2-R: non-trivial optima with cross-eval diagonal dominance |
| §6 | [J] = T⁻¹, [λ] = T | consistent with all runs (λ=0.2) |
| §7 | Γ_F demoted: fixed-(F) local diagnostic only | D2.1: singular on unrestricted TPS (corr −0.999, Φ→44) |
| §8 | J_dyn^(τ) operational extension | D2.2-B: same-basin at τ→0, crossover at larger τ (G7) |
| §9 | lim_{τ→0} J_dyn^(τ) = J_dyn^(0) | D2.2-B G7-1: dev ≈ 0.96·τ·\|J₀\|, ≤2.5×10⁻³ at τ=10⁻³ |
| §10 | finite-τ TPS crossover admissible | D3: τ_c = 0.018 (Haar), first-order-like in realization; D4: persists N=4,5 |

## 2. Validation summary (D0–D4)

| Phase | Result | Verdict |
|---|---|---|
| D1 (Γ-based) | 5/5 non-trivial F* | preliminary (superseded functional) |
| D2 (Γ-based) | Haar anomaly Φ=49.1 | re-scoped → D2.1 |
| D2.1 | normalization singularity | **falsification** (Claim 1) |
| D2.2-A | J_dyn^(0): G1–G6 | **well-posed** (Claim 2) |
| D2.2-B | J_dyn^(τ): G7 PASS (documented finite-τ crossover) | extension consistent |
| D1-R | Hamiltonian response, mean gap 0.113 | **Claim 3** |
| D2-R | state response, mean gap 0.159 | **Claim 3** |
| D3 | τ_c = 0.018 (Haar); τ_c·Δ_L = 0.0119 (O(1) rejected) | **Claim 4** + negative result |
| D4 | τ_c flat 0.017–0.022 across N=3–5 | **finite-size persistence** |

## 3. Statements currently supported vs not yet

Supported:
- "The normalized local decay rate cannot serve as a global variational
  functional over unrestricted tensor factorizations" (D2.1).
- "J_dyn^(0) is well-posed on the unrestricted TPS domain; singularity,
  product collapse, and seed instability are absent in the tested protocol"
  (D2.2-A).
- "F* responds to Hamiltonian and state" (D1-R, D2-R cross-eval).
- "Relational structure can depend systematically on the dynamical
  observation scale, with the presence and location of crossovers being
  state-dependent" (D3; ground/mixed show no crossover at τ ≤ 1).
- "The crossover scale is not simply set by the global Liouvillian
  relaxation time" (D3-C: τ_c·Δ_L = 0.0119).
- "Competing TPS basins with a finite-τ objective crossing persist at
  N = 4,5 with τ_c ≈ 0.017–0.022" (D4).

Supported with qualification:
- "Crossover persists in the tested N=5 optimization subspace" (200/945
  dims). **Resolved 2026-08-12**: the N=5 full-quotient closure run
  (945-dim, 2 seeds × 100 steps) gives τ_c(5, full) = 0.0217, consistent
  with the flat trend — **"full TPS optimization at N=5 confirms
  persistence" is now supported**.
- "τ_c ≈ 0.018–0.022 is size-independent" — flat within N=3–5 (spread
  0.178); a finite-size trend, not a scaling law.

Not claimed:
- Scaling law / thermodynamic limit (N=3–5 only).
- Phase transition (finite dimension; "sharp structural crossover" /
  "first-order-like basin transition" only).
- τ_c ≈ τ_L (rejected: τ_c·Δ_L = 0.0119).

## 4. Pending items at the time of writing

1. ~~**N=5 full quotient closure** (945-dim)~~ **DONE 2026-08-12**:
   τ_c(5, full) = 0.0217, d_F(F_A,F_B) = 0.319 — persistence confirmed
   under full TPS optimization. (The subspace value 0.0173 was a
   restricted-domain underestimate; the full value sits in the flat band.)
2. **Thermal τ_c precision** (currently bracketed (0.3, 1.0)) — lower
   priority; can be folded into Phase E.
3. **Phase E — τ_c mechanism**: test τ_c ≈ −ΔΦ₀/ΔΦ₁ against the observed
   τ_c(N). This is the next research phase, not a validation gap.

## 5. Reproducibility

Scripts, JSON outputs, and commits are listed in the Phase D report §9.
Suite: 74 passed + 1 xfailed. The addendum and the report are the
"how it was verified" documents; Spec v2.2 remains the "what is defined"
document.
