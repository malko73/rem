# Response to Independent Reviewers — REM Paper v0.1 (review candidate)

**Reviews**: A = Gemini 3.6 flash (math/quantum-info), B = OpenAI gpt-5.4-mini
(referee). Four rounds (r1–r4) with fixes between rounds.
**Current state (r4)**: A = **Minor revision**; B = **Major revision**
(5 scope/wording Majors, no Critical).
**QC**: 14/14 PASS maintained after every round; suite 74 passed + 1 xfailed.

---

## Round trajectory

| Round | Review A | Review B | Fixes |
|---|---|---|---|
| r1 | Minor revision (1 Major, 3 Minor) | Major revision (1 Critical, 11 Major) | r2: sensitivity/uncertainty supplement (ΔΦ₂, errors, reps, grid); wording/notation; compact-quotient argument; physical reading of J_dyn^(0); reproducibility §9.4 |
| r2 | Minor revision (qualification on τ_c^(2)) | Major revision (no Critical) | r3: "derived to first order" reframing (τ_c^(1) derived, τ_c^(2) fit); conditional scope of F_B; claim scoping; monotonicity; explicit 2nd-derivative derivation; D1-R full matrix (Table S3) |
| r3 | Minor revision (1 Major: conditional scope) | Major revision (5 Majors) | r4: J₁ coefficient convention; Q_F operator-space duality; "well-posed in the tested protocol"; "no statement about N→∞" |
| r4 | **Minor revision** | **Major revision** (5 scope/wording Majors, no Critical) | — |

## Reviewer A — status and remaining items

**Recommendation (r4): Minor revision.** The mathematical audit is clean
(derivations exact; factor-of-2 and signs verified; the executive question
answered YES). Remaining items are Minor: Q_F explicit form (added r4),
grid-resolution in main text (added r3/r4), state-degeneracy reiteration
(added r3). No blockers.

## Reviewer B — remaining Major items (r4) and disposition

1. **J₁/Q_F/finite-time conventions** — addressed in r4 (J₁ = coefficient of
   τ through the 1/(2τ) prefactor, stated explicitly; Q_F duality:
   operator-space projector vs Schmidt-basis construction, coincident in
   the rotated frame). **Closed.**
2. **Claim 2 "well-posed" exceeds evidence** — the claim is now scoped to
   the tested protocol in the claim statement and the conclusion; the
   heading reads "removes the singularity"; the formal existence argument
   (compact quotient + continuity) is stated. The word "well-posed"
   remains as the project's established claim identifier (Spec v2.2
   terminology). **Partly a terminology judgment call.**
3. **Claim 4b "persistence" risks scaling/universality** — scoped
   repeatedly: "finite-size persistence / finite-size trend", "no
   statement about N→∞", "not a scaling law and not a thermodynamic
   limit". **Closed in wording; the term itself is the master's
   prescribed claim name.**
4. **Claim 4c "derived" conditional** — the paper now says "derived to
   first order … at the fixed basin pair (F_A,F_B)"; F_B's finite-τ
   optimization is stated in the claim, the Fig. 4 caption, and the
   conclusion. **Closed.**
5. **Abstract/conclusion overstate generality** — softened across r2–r4
   ("show numerically", "in the tested families and states",
   "approximately flat over the tested finite-size range", "passed the
   validation checks reported here"). **Mostly closed; residual is
   stylistic.**

## Final disposition (author)

- **Review A**: Minor revision / converged.
- **Review B**: Major revision label remains, but **no unresolved Critical or
  new technical defect**.
- **Residual concerns**: claim-scope / editorial conservatism.
- **Author disposition**: addressed where scientifically justified; further
  wording attenuation **rejected as non-substantive**.
- **Status sentence**: independent adversarial review completed; all
  identified technical defects resolved, with residual disagreements
  limited to claim-scope wording.
- `well-posed` is retained with the qualifier **"in the tested protocol"**;
  `finite-size persistence` is retained, scoped to the **tested finite
  systems N = 3–5**. These scope boundaries make the scientific claims
  sufficiently precise.

## Assessment

- Every **technical** finding from both reviewers has been addressed and
  is verifiable in the manuscript (QC 14/14).
- Review A has converged to Minor revision.
- Review B's remaining Majors are **scope/wording judgment calls** on a
  finite-N numerical study. Across four rounds it re-raises the same five
  themes even after concrete fixes, and its "Major revision" standard
  appears calibrated to theorem-grade claims (e.g., "prove global
  well-posedness", "prove the optimizer generically finds singular
  directions"). The stochasticity of LLM reviewing also means severity
  can vary run to run.
- The manuscript itself now contains the bounded claims, the error
  budget, and the reproducibility trail that a referee needs.

## Decision

**r5 cancelled by the author (2026-08-12).** The current revision is
frozen as the `paper-v1.0` candidate after the final release QC. This
document and the review logs (review_A_gemini_r1..r3, review_B_openai_r1..r3,
round-4 outputs, integrated_findings.md) are kept permanently as the
review record.
