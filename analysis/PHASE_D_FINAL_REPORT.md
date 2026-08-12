# Phase D — Validation Report

**Relational Emergence Model (REM), Spec v2.2 canonical functional**

**Date**: 2026-08-12
**Status**: Final — numerical validation phase closed (D0–D4)
**Repository**: github.com/malko73/rem, branch `open`
**Suite**: 74 passed + 1 xfailed (root `conftest.py` makes plain `pytest` work)

---

## 1. Executive summary

This report consolidates the Phase D numerical validation of the REM
Spec v2.2 canonical variational functional. The central narrative is a
falsification-driven arc, not a checklist of passes:

```
normalized functional falsified (D2.1)
        ↓
J_dyn^(0) proposed and shown well-posed (D2.2-A/B, G1–G7)
        ↓
Hamiltonian response (D1-R)
        ↓
state response (D2-R)
        ↓
finite-time structural crossover (D3)
        ↓
finite-size persistence (D4)
```

The value of this arc is that **the canonical functional was falsified
mid-course and the corrected version survived a battery of tests** — the
research process itself, not just the final pass, is the result.

**Four paper claims** (Section 8):
1. The normalized rate Γ_F^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩/|Q_Fρ|² is **singular**
   as an unrestricted global TPS variational objective.
2. The unnormalized functional J_dyn^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩ is
   **well-posed**: singularity, product collapse, and seed instability
   disappear.
3. F* = F*(ρ, L, λ) **responds** to both Hamiltonian and state (cross-
   evaluation matrix diagonal dominance).
4. Competing TPS basins with a finite-τ objective crossing exist, with
   **τ_c ≈ 0.017–0.022 persisting across N = 3–5**.

---

## 2. Frozen protocol

All D-series runs use the following frozen protocol unless stated otherwise:

- **Hamiltonian**: asymmetric XY chain, J₁₂=1.5, J₂₃=0.6 (…0.6 for N>3), h=0.2
- **Environment**: pure dephasing γ = (0.5, 1.0, 2.0), periodic extension for N>3
- **Objective**: J_dyn^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩ (canonical, denominator-free);
  J_dyn^(τ) = −(C_F²(τ)−C_F²(0))/(2τ) (finite-time extension)
- **Optimizer**: Adam, 200 steps (N=3), lr = 0.01, finite-difference gradients
- **Seeds**: SEED0 = 20260813 (D2-series) / 20260812 (D1-series); 6 seeds (N=3),
  4 seeds (N=4), 4 seeds (N=5, 200-dim subspace)
- **λ = 0.2**, cut n_A = 2 (N=3: 2|1; N=4: 2|2 and 1|3; N=5: 2|3)
- **Crossover definition** (D3/D4): ΔΦ(τ) = Φ_τ(F_B) − Φ_τ(F_A) at **fixed
  basin representatives**; τ_c solves ΔΦ(τ_c) = 0. Never defined by optimizer
  basin hops.

Key infrastructure notes:
- `tps_distance` generalized to (d_a, d_b) for N ≥ 4 (commit 78eb500).
- Per-seed basin labels via d_F are **not meaningful** on the quotient (many
  distinct local maxima sit ~0.87 from any two fixed reps); value-based
  classification (closeness to Φ_τ(F_A) vs Φ_τ(F_B)) is the correct tool.

---

## 3. The falsification: normalized Γ_F is singular (D2.1)

The normalized instantaneous rate

```
Γ_F^(0) = −Re⟨Q_Fρ, Q_FL(ρ)⟩ / |Q_Fρ|²
```

was the v2.1 canonical cost candidate. On unrestricted TPS optimization
(full quotient, Haar states) the optimizer drives |Q_Fρ|² → 0, which sends
Γ_F → ±∞ and Φ → +∞ along a non-physical direction.

**Evidence** (commit `0ad8840`):
- corr(Φ, log₁₀ C_F(0)²) = **−0.999** (near-perfect anti-correlation)
- C_F(0)² → 6×10⁻⁵ drives Γ → −222, Φ → 44
- Best solution Schmidt spectrum p = [0.9997, 0.0003] (near-product collapse)

Since every pure state can be approached by a product-state TPS, the
denominator can vanish. The normalized rate is therefore well-posed only as
a **fixed-(F) local diagnostic**; it is **prohibited as an unrestricted
global TPS variational objective** (Spec v2.2 §7).

> **Claim 1 (negative result)**: the normalized local decay rate cannot
> serve as a global variational functional over unrestricted tensor
> factorizations.

---

## 4. J_dyn^(0): well-posedness (D2.2-A/B)

The canonical functional (Spec v2.2 §4–5):

```
J_dyn^(0)(F;ρ,L) = −Re⟨Q_Fρ, Q_F L(ρ)⟩          [J] = T⁻¹
Φ(F;λ)           = I_ρ(F) − λ J_dyn^(0)(F)      [λ] = T
```

By Cauchy–Schwarz, |J_dyn^(0)| ≤ |Q_Fρ|₂ |Q_FL(ρ)|₂, so J_dyn^(0) → 0
smoothly as |Q_Fρ|² → 0: **no denominator blow-up**.

### D2.2-A (commit `24db549`) — G1–G6 PASS

Frozen D2 protocol, 4 states × 6 seeds, only the dynamical term changed
(Γ → J_dyn^(0)):

| State | best Φ (normalized Γ, D2) | best Φ (J_dyn^(0)) | std |
|---|---|---|---|
| ground | 1.7993 | **1.8995** | 0.0007 |
| haar | 49.10 ⚠️ | **1.8594** | 0.0090 |
| mixed | 1.8485 | **1.9845** | 0.0051 |
| thermal | 1.8499 | **1.9465** | 0.0058 |

Haar 30-seed: best 1.8595, std 0.0077 (vs 25.71 / 5.96 under Γ). min C_F² =
0.0956 (vs 4×10⁻⁴ under Γ — the optimizer never approaches the singular
region); Schmidt p_max 0.51–0.58 (no product collapse); perturbation
continuity max|ΔΦ| = 7.7×10⁻⁵; singularity sweep (C_F²→0 direction) all
finite.

### D2.2-B (commit `2fce54d`) — G7 PASS with documented finite-τ crossover

J_dyn^(τ) = −(C_F²(τ)−C_F²(0))/(2τ), τ sweep {10⁻³ … 10⁻¹}, frozen protocol.
- **G7-1**: exact τ→0 limit (dev ≈ 0.96·τ·|J₀|, ≤2.5×10⁻³ at τ=10⁻³) —
  the derivative identity dC_F²/dt|₀ = 2Re⟨Q_Fρ, Q_FL(ρ)⟩ is verified.
- **G7-2**: same basin for all states at τ ≤ 3×10⁻³; ground/mixed/thermal at
  all τ (gaps ≤ 0.005); **Haar basin transition at larger τ** — the seed of D3.
- **G7-3**: state ranking preserved at all τ.
- **G7-4**: no singularity / product collapse / seed instability at any τ.

> **Claim 2**: with the unnormalized functional, singularity, product
> collapse, and seed instability disappear; J_dyn^(0) is well-posed as the
> canonical variational functional, with J_dyn^(τ) as its operational
> finite-time extension / consistency diagnostic.

---

## 5. Structural response (D1-R, D2-R)

### D1-R — Hamiltonian response (commit `871a104`)

Frozen old-D1 protocol (SEED0 = 20260812), only the dynamical term changed
(Γ → J_dyn^(0)). 5 Hamiltonian families, 6 seeds each.

| Family | best Φ | std |
|---|---|---|
| asymmetric_XY | 1.8994 | 0.0007 |
| transverse_field_ising | 1.8246 | 0.0010 |
| heisenberg_xxz | 1.8695 | 0.0002 |
| xyz | 1.8740 | 0.0006 |
| random_local | 1.8481 | 0.0007 |

All gates R1–R6 PASS. **R5 (key)**: cross-evaluation matrix Φ_{H_i}(F*_j)
shows diagonal dominance — mean gap = 0.113 (asym 0.068 / ising 0.092 /
heisenberg 0.106 / xyz 0.025 / random 0.274). Classification:
**family-dependent response with varying separation strength** (xyz gap is
small; do not claim all families strongly separated). d_F matrix: physically
related families closer (asym↔random 0.412, ising↔xyz 0.559; others
0.87–0.88).

### D2-R — state response (commit `64e5b03`)

Frozen old-D2 protocol (SEED0 = 20260813), 4 states × 6 seeds.

| State | best Φ | std |
|---|---|---|
| ground | 1.8995 | 0.0007 |
| haar | 1.8594 | 0.0090 |
| mixed | 1.9845 | 0.0051 |
| thermal | 1.9465 | 0.0058 |

All gates R1–R6 PASS. **R5**: cross-eval mean gap = 0.159 (ground 0.000 /
haar 0.332 / mixed 0.058 / thermal 0.247). **Precise wording**: state-
dependent response is established overall, although ground and mixed contain
nearly degenerate optima under the ground-state objective. d_F: ground↔mixed
0.343 (close), haar/thermal strongly separated (0.87–0.89); mixed vs thermal
have close Φ but d_F 0.868 (different structures). Not the same attractor.

> **Claim 3**: F* = F*(ρ, L, λ) responds to both Hamiltonian and state —
> not a state-independent common solution. Evidence: cross-evaluation
> matrices with diagonal dominance (objective response) plus geometric
> response (d_F).

---

## 6. Finite-time structural crossover (D3, commit `8c1a236`)

**Question**: does relational structure depend systematically on the
dynamical observation scale? Central target: Haar state.

### D3-A — τ_c precision identification

ΔΦ(τ) = Φ_τ(F_B) − Φ_τ(F_A) at fixed reps (F_A = D2.2-A Haar optimum, F_B =
τ=0.1 Haar optimum), fine grid τ ∈ [0.003, 0.1] (19 points, 6 seeds):

| τ | 0.003 | 0.01 | 0.016 | **0.018** | 0.02 | 0.025 | 0.03 | 0.05 | 0.1 |
|---|---|---|---|---|---|---|---|---|---|
| ΔΦ | −0.021 | −0.011 | −0.003 | **0.000** | +0.003 | +0.009 | +0.016 | +0.039 | +0.079 |

**τ_c = 0.0180** (linear interpolation of the sign change). ΔΦ is smooth and
monotone. The optimizer follows sharply: at τ=0.018 all 6 seeds sit at the A
value, at τ=0.02 all 6 at the B value (value-based classification; std
collapses 0.005 → 0.0009).

**Method corrections established here**:
1. d_F-based per-seed basin labels are meaningless (many local maxima ~0.87
   from both reps); value-based classification is correct.
2. Optimizer basin hops are NOT the physical crossover — the objective-level
   ΔΦ crossing is (exactly the separation the master anticipated).

### D3-B — transition nature

**First-order-like basin transition (sharp structural crossover) in the
realization; continuous in the objective.** Adjacent d_F jumps O(1) at
τ=0.02 and τ=0.05 (best solution hops between distinct maxima, reaching F_B
exactly at τ=0.1, d_F=0.000); ΔΦ itself is smooth/monotone. Finite
dimension: "sharp structural crossover" / "first-order-like basin
transition", NOT "phase transition".

### D3-C — Liouvillian timescale: negative result

Δ_L = min_{μ≠0}|Re μ| = 0.658, τ_L = 1.519. **τ_c·Δ_L = 0.0119 — the O(1)
conjecture is rejected**: τ_c is far shorter than the Liouvillian gap
timescale. τ_c is set by the inter-basin competition scale (the finite-time
O(τ) correction overcoming basin B's instantaneous disadvantage), not by the
slowest relaxation mode. (Thermal's τ_c ∈ (0.3, 1.0) gives ratio ~0.2–0.66 —
state-dependent, not concentrated.)

> **Constraint**: the crossover scale is not simply set by the global
> Liouvillian relaxation time.

### D3-7 — state dependence of crossover presence

- ground: **no crossover** up to τ=1.0 (d_F drift 0.11, ΔΦ ≈ 0)
- mixed: **no crossover** up to τ=1.0 (d_F 0.12, ΔΦ ≈ 0)
- **thermal: crossover with τ_c ∈ (0.3, 1.0)** (ΔΦ −0.011 → +0.001,
  d_F(F_A,F_B) = 0.879)

**Conclusion (precise wording)**: relational structure **can** depend
systematically on the dynamical observation scale, with the presence and
location of crossovers being state-dependent.

---

## 7. Finite-size persistence (D4, commit `78eb500`)

**Question**: does the D3 crossover persist as the Hilbert-space size
increases? N=3,4,5 — establishes **finite-size persistence / finite-size
trend**, NOT a scaling law or thermodynamic limit.

### D4-A — canonical well-posedness vs N (J_dyn^(0) only)

| N | cut | quotient | state | best Φ | std | Ĩ | min C_F² | p_max |
|---|---|---|---|---|---|---|---|---|
| 3 | 2\|1 | 45 (full) | ground | 1.8995 | 0.0007 | 1.000 | 0.500 | 0.51 |
| 3 | 2\|1 | 45 (full) | haar | 1.8594 | 0.0090 | 0.982 | 0.487 | 0.58 |
| 4 | 2\|2 | 225 (full) | ground | 3.6923 | 0.0008 | 0.999 | 0.749 | 0.27 |
| 4 | 2\|2 | 225 (full) | haar | 3.7041 | 0.0008 | 0.997 | 0.748 | 0.28 |
| 5 | 2\|3 | 200/945 (subspace) | ground | 3.6713 | 0.0006 | 1.000 | 0.750 | 0.26 |
| 5 | 2\|3 | 200/945 (subspace) | haar | 3.6346 | 0.0022 | 0.998 | 0.749 | 0.28 |

Ĩ = I/(2 log₂ d_min) ≈ 0.98–1.00 at all N (near-maximal informational
articulation). Well-posedness (no denominator) holds at all sizes.

### D4-B — crossover persistence

| N | cut | d_F(F_A,F_B) | **τ_c(N)** |
|---|---|---|---|
| 3 | 2\|1 (full 45) | 0.886 | **0.0180** |
| 4 | 2\|2 (full 225) | 0.968 | **0.0206** |
| 4 | 1\|3 (full 189) | 0.994 | **0.0224** |
| 5 | 2\|3 (subspace 200/945) | 0.992 | 0.0173 |
| 5 | **2\|3 (full 945)** | **0.319** | **0.0217** |

The crossover **exists at N=4 and N=5** (ΔΦ sign flip within τ ∈ [10⁻³, 1.0]).
**N=5 full-quotient closure confirms persistence**: τ_c(5, full) = 0.0217,
consistent with the flat trend (the subspace value 0.0173 was a
restricted-domain underestimate). Note d_F(F_A,F_B) differs between the
subspace (0.992) and the full quotient (0.319): the full quotient finds
competing basins that are closer in structure, but the ΔΦ crossing and the
timescale are the same. **τ_c(N) trend: FLAT** (spread 0.178 including the
full-quotient value). The N=4 balanced (2|2) vs asymmetric (1|3) control
shows τ_c nearly identical (0.021 vs 0.022) → **bipartition shape has little
effect at fixed size**, strengthening the N=5 (asymmetric 2|3) reading.

> **Claim 4**: competing TPS basins with a finite-τ objective crossing exist
> at N=4,5 with an essentially size-independent τ_c ≈ 0.018–0.022 — **full
> TPS optimization at N=5 confirms persistence** (closure run, 2026-08-12).

### D4-C — cross-size comparison

Raw values are preserved; normalized indicators added (canonical NOT
re-renormalized, the denominator is never reintroduced): Ĩ = I/(2 log₂ d_min),
J̃ = J_dyn^(0)/γ_mean.

### D4 limitations (explicit)

- N=3–5 range only; no scaling-law claim.
- All crossovers are in a finite-dimensional system; "sharp structural
  crossover" language only.
- The N=5 full-quotient run used 2 seeds × 100 steps (documented protocol
  reduction for the 945-dim quotient; std 0.0002 — tightly converged).

---

## 8. Paper claims (candidate framing)

| # | Claim | Primary evidence |
|---|---|---|
| 1 | Normalized rate Γ_F singular on unrestricted TPS | D2.1 (corr −0.999, Φ→44, product collapse) |
| 2 | J_dyn^(0) well-posed; singularity/collapse/instability gone | D2.2-A (G1–G6), D2.2-B (G7) |
| 3 | F* = F*(ρ, L, λ) responds to Hamiltonian and state | D1-R / D2-R cross-eval matrices |
| 4 | Finite-time structural crossover, τ_c ≈ 0.018–0.022, persists N=3–5 (full-quotient confirmed) | D3 (ΔΦ crossing, τ_c=0.018), D4 (τ_c flat; N=5 full closure) |

Plus the **negative result**: τ_c·Δ_L = 0.0119 — the crossover scale is not
set by the global Liouvillian relaxation time; τ_c is governed by the
objective balance between competing structural basins (→ Phase E).

---

## 9. Reproducibility

| Phase | Scripts | Outputs | Commit |
|---|---|---|---|
| D0/D1 | `d0_protocol.py`, `d1_generality.py`, `d1_hamiltonians.py` | d0_protocol.json, d1_hamiltonian_generality.json | a091f72 |
| D2 | `d2_state_generality.py` | d2_state_generality.json | a091f72 |
| D2.1 | `d2_1_haar_singularity_audit.py` | d2_1_haar_singularity_audit.json | 0ad8840 |
| D2.2-A | `d2_2a_unnormalized_instantaneous.py` | d2_2a_unnormalized_instantaneous.json | 24db549 |
| D2.2-B | `d2_2b_finite_time.py`, `d2_2b_gates.py` | d2_2b_tau_*.json, d2_2b_gates.json | 2fce54d |
| D1-R | `d1r_canonical_revalidation.py` | d1r_canonical_revalidation.json | 871a104 |
| D2-R | `d2r_state_revalidation.py` | d2r_state_revalidation.json | 64e5b03 |
| D3 | `d3_timescale.py` | d3_timescale.json | 8c1a236 |
| D4 | `d4_nscaling.py`, `d4_gates.py` | d4_*.json, d4_finite_size.json | 78eb500 |

Suite: 74 passed + 1 xfailed. Roadmap: `REM_next_steps.md` (D0–D4 status).
Spec: `REM_spec_v2_2.md` + `papers/REM_spec_v2_2.tex` (frozen, built by
`tools/gen_spec.py`, CI latex gate green). This report is independent of the
Spec (theory specification ≠ post-hoc numerical evidence).

---

## 10. Phase E outlook — τ_c mechanism

The natural next phase: **why does τ_c ≈ 0.02?** With the finite-time
expansion J_τ(F) = J₀(F) + τJ₁(F) + O(τ²), two basins A,B give

```
ΔΦ(τ) ≈ ΔΦ₀ + τ ΔΦ₁   ⇒   τ_c ≈ −ΔΦ₀/ΔΦ₁
```

which can be tested against τ_c(3..5) = 0.0180, 0.0206, 0.0224, 0.0173. If
the formula reproduces the trend, the crossover timescale becomes a derived
quantity rather than a numerical observation. This is the priority for the
next research cycle (after paper submission / Zenodo update).
