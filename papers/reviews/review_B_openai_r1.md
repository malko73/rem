## Reviewer B report

### Overall assessment
The manuscript presents an interesting numerical story, but the evidence as written does **not** fully support the strength of the claims. The main issues are: (i) several derivations are either incomplete or internally inconsistent with the stated formulas; (ii) the numerical evidence is often suggestive but not sufficient to justify broad language such as “well-posed,” “derived, not fitted,” “finite-size persistence,” or “state-dependent crossover” without stronger error budgets and reproducibility details; and (iii) the paper repeatedly overextends finite-\(N\), finite-seed, and fixed-representative results into broader claims.

### Bottom line on the most important question
**If I ignore the authors’ interpretation and inspect only the definitions, derivations, numerical evidence, and figures, Claims 1–4 do not all follow cleanly from the presented evidence.**  
- Claim 1: **mostly supported**, but the “on the unrestricted quotient” and “\(+\infty\) along non-physical near-product directions” language is stronger than the evidence shown.
- Claim 2: **partially supported**, but “well-posed” is too strong without a formal existence/compactness/continuity argument and more complete optimization diagnostics.
- Claim 3: **supported at the level of numerical nontriviality**, but the scope is limited to the tested Hamiltonians/states and does not establish general state/Hamiltonian dependence.
- Claim 4a–4c: **the weakest part**. The existence of a crossover in the tested cases is plausible, but the “derived, not fitted” claim is not fully justified, the finite-size persistence claim is overextended, and the mechanism analysis lacks enough uncertainty propagation to make the quoted 0.0–0.5% quadratic accuracy convincing.

---

# 1. Math / definition audit

## Finding 1.1 — Inconsistent notation and derivative identities for \(J_{\mathrm{dyn}}^{(\tau)}\)  
**Severity: Major**

You define
\[
J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{C_F^2(\tau)-C_F^2(0)}{2\tau},
\qquad
C_F^2(t)=|Q_F U^\dagger e^{tL}\rho\,U|^2,
\]
and then state
\[
\frac{dC_F^2}{dt}\big|_0 = 2\mathrm{Re}\langle Q_F\rho, Q_F L(\rho)\rangle
\]
and
\[
\lim_{\tau\to 0}J_{\mathrm{dyn}}^{(\tau)} = J_{\mathrm{dyn}}^{(0)}.
\]
This is fine in principle, but later in Claim 4c you introduce
\[
C^2{}''(0)=2|Q_F Y_U|^2 + 2\mathrm{Re}\langle Q_F\rho_U, Q_F(U^\dagger L^2\rho_0 U)\rangle,
\]
and conclude
\[
J_1 = -C^2{}''(0)/4.
\]
That coefficient is not obviously consistent with the Taylor expansion of
\[
J_\tau = -\frac{C^2(\tau)-C^2(0)}{2\tau}.
\]
From the expansion
\[
C^2(\tau)=C^2(0)+\tau C'^2(0)+\frac{\tau^2}{2}C''^2(0)+\cdots,
\]
one gets
\[
J_\tau = -\frac{1}{2}C'^2(0)-\frac{\tau}{4}C''^2(0)+\cdots,
\]
so \(J_0=-\tfrac12 C'^2(0)\) and \(J_1=-\tfrac14 C''^2(0)\). That part is consistent. However, the paper does not clearly define whether \(J_1\) is the coefficient in \(J_\tau=J_0+\tau J_1+\cdots\) or the derivative \(dJ_\tau/d\tau|_0\). The notation is used both ways in prose, which is a reproducibility and auditability problem.

## Finding 1.2 — The “bounded by Cauchy–Schwarz” statement is mathematically incomplete  
**Severity: Major**

You write
\[
|J_{\mathrm{dyn}}^{(0)}| \le |Q_F\rho|\,|Q_FL(\rho)|
\]
and conclude “so \(J_{\mathrm{dyn}}^{(0)}\to 0\) smoothly as \(|Q_F\rho|^2\to 0\); there is no denominator to blow up.”

The inequality is correct, but the conclusion “smoothly” is stronger than what is shown. The bound only implies \(J_{\mathrm{dyn}}^{(0)}\to 0\) if \(Q_FL(\rho)\) remains bounded along the relevant path. That is likely true for finite-dimensional \(L\), but it is not stated as a lemma. Also, “bounded by Cauchy–Schwarz” does not by itself establish global well-posedness of the optimization problem.

## Finding 1.3 — The singularity claim for the normalized functional is plausible but not fully proven as stated  
**Severity: Major**

The normalized candidate is
\[
\Gamma_F^{(0)} = -\frac{\mathrm{Re}\langle Q_F\rho, Q_FL(\rho)\rangle}{|Q_F\rho|^2}.
\]
The paper claims the denominator can vanish on the unrestricted quotient and that the optimizer drives \(|Q_F\rho|^2\to 0\). The numerical evidence in D2.1 supports this behavior for the Haar state, but the statement “Since every pure state can be approached by a product-state TPS, the denominator can vanish on any open neighborhood of the quotient” is too sweeping. The quotient is finite-dimensional and the set of TPSs is not literally an open neighborhood in the topological sense used here. This is a conceptual overreach.

---

# 2. Numerical evidence audit

## Finding 2.1 — Claim 1 is supported numerically, but the evidence is narrow  
**Severity: Minor**

The digest gives:
- \(\mathrm{corr}(\Phi,\log_{10} C_F^2) = -0.998955\),
- \(\Phi_{\max}=25.7133\),
- \(\Gamma_{\min}=-128.5251\),
- \(C_{F,0}^2{}_{\min}=4.093466\times 10^{-4}\),
- \(p_{\max}=0.999795\).

These numbers do support the claim that the normalized objective can be gamed by collapsing toward near-product directions in the tested Haar/\(N=3\) setting. However, the paper’s wording “on the unrestricted quotient” is broader than the evidence. You only show one state class and one system size in the main falsification figure/caption. That is enough for a counterexample to the normalized functional as a global objective, but not enough to characterize the full quotient behavior in general.

## Finding 2.2 — Claim 2 is supported as a numerical comparison, but “well-posed” is not established  
**Severity: Major**

The table in Sec. 4 reports:
- ground: \(\Phi^*=1.8993\) under \(\Gamma\), \(1.8995\) under \(J^{(0)}\), seed std \(0.0007\),
- haar: \(49.10\) (singular) under \(\Gamma\), \(1.8594\) under \(J^{(0)}\), seed std \(0.0090\),
- mixed: \(1.8485\to 1.9845\),
- thermal: \(1.8499\to 1.9465\).

Also:
- minimum \(C_F^2(0)=0.0956\),
- maximum \(|\gamma_{\mathrm{ref}}|=1.66\),
- seed stability over 30 Haar runs: best \(1.8595\), std \(0.0077\),
- continuity perturbation: \(\max|\Delta\Phi|=7.7\times10^{-5}\).

This is good evidence that the replacement functional behaves better numerically. But “well-posed” is a mathematical claim, not just a numerical one. You have not shown:
1. existence of a maximizer on the quotient,
2. continuity of \(\Phi\) on the full domain under the exact parameterization used,
3. absence of flat directions or optimizer pathologies beyond the tested seeds,
4. robustness to alternative initialization/step-size choices.

So the evidence supports “more stable and non-singular in the tested cases,” not a full well-posedness theorem.

## Finding 2.3 — Claim 3 is supported, but only as a finite test of nontrivial dependence  
**Severity: Minor**

The Hamiltonian-response table gives distinct \(\Phi^*\) values:
- asymmetric_XY \(1.8994\),
- transverse_field_ising \(1.8246\),
- heisenberg_xxz \(1.8695\),
- xyz \(1.8740\),
- random_local \(1.8481\),

with cross-evaluation mean gap \(0.113\).  
The state-response table gives:
- ground \(1.8995\),
- haar \(1.8594\),
- mixed \(1.9845\),
- thermal \(1.9465\),

with mean gap \(0.159\).

This is enough to show the selected structure is not trivially constant across the tested families/states. But the claim “\(F^*=F^*(\rho,L,\lambda)\), not a state-independent common solution” should be phrased as “in the tested cases, the optimizer returns different structures/values across states and Hamiltonians.” The current wording is too universal.

## Finding 2.4 — Claim 4a crossover is numerically plausible, but the evidence is under-specified  
**Severity: Major**

For Haar, \(N=3\), the table gives:
- \(\tau=0.003\): \(\Delta\Phi=-0.021\),
- \(0.01\): \(-0.011\),
- \(0.016\): \(-0.003\),
- \(0.018\): \(0.000\),
- \(0.02\): \(+0.003\),
- \(0.025\): \(+0.009\),
- \(0.03\): \(+0.016\),
- \(0.05\): \(+0.039\),
- \(0.1\): \(+0.079\).

This does support a sign change near \(0.018\). However:
- the caption says “19 points,” but only 9 are shown in the table;
- the method says 6 seeds per point, but the table reports only one \(\Delta\Phi\) per \(\tau\), with no error bars or seed spread;
- the “all 6 seeds sit at the \(A\)-value” / “all 6 at the \(B\)-value” classification is not independently verifiable from the table alone;
- the claim that the two basins are “distinct structural basins” relies on \(d_F(F_A,F_B)=0.886\), but no basin stability analysis is shown.

So the crossover is plausible, but the evidence is not yet sufficient for the stronger basin-language.

## Finding 2.5 — Claim 4b finite-size persistence is only weakly supported and partly contradicted by the paper’s own numbers  
**Severity: Major**

The main text states:
\[
\tau_c = 0.0180, 0.0206, 0.0224, 0.0217
\]
for \(N=3,4,4,5\), and calls this “approximately flat.”

But the supplementary history gives an earlier \(N=5\) restricted-domain value:
\[
\tau_c = 0.0173
\]
with \(d_F=0.992\), while the full-quotient closure gives \(0.0217\) with \(d_F=0.319\).

This is not a contradiction if clearly labeled as historical, but it does show that the \(N=5\) result is protocol-sensitive. That weakens the “persistence” claim. Also, the four values span \(0.0180\) to \(0.0224\), a relative spread of about \(24\%\) around the mean \(\sim 0.0207\). Calling this “approximately flat” is acceptable only as a descriptive statement, not as evidence of a robust finite-size trend.

## Finding 2.6 — Claim 4c “derived, not fitted” is not fully justified by the presented evidence  
**Severity: Critical**

This is the most serious numerical issue.

You claim:
\[
\tau_c^{(1)} = -\frac{\Delta\Phi_0}{\Delta\Phi_1}
\]
is derived, not fitted, and that the quadratic correction matches within \(0.0\)–\(0.5\%\).

But the evidence shows:
- \(N=3\): \(\Delta\Phi_0=-0.0260\), \(\Delta\Phi_1=1.521\), \(\tau_c^{(1)}=0.0171\), measured \(0.0180\), error \(5.1\%\).
- \(N=4,2|2\): \(-0.0099\), \(0.526\), \(0.0188\), measured \(0.0206\), error \(8.6\%\).
- \(N=4,1|3\): \(-0.0339\), \(1.633\), \(0.0208\), measured \(0.0224\), error \(7.2\%\).
- \(N=5,2|3\): \(-0.0219\), \(1.130\), \(0.0193\), measured \(0.0217\), error \(10.9\%\).

These are not small enough to justify “derived” in the strong sense unless the derivation is explicitly presented as a first-order approximation with known truncation error. More importantly, the quadratic correction is claimed to achieve \(0.0\)–\(0.5\%\) error, but the paper does not provide:
1. the actual \(\Delta\Phi_2\) values,
2. uncertainty propagation from \(\Delta\Phi_0,\Delta\Phi_1,\Delta\Phi_2\),
3. sensitivity to the choice of fixed representatives \(F_A,F_B\),
4. a check that the measured \(\tau_c\) itself is not resolution-limited by the \(\tau\)-grid.

Without those, “derived, not fitted” is too strong. The formula is derived algebraically, yes, but the numerical agreement is still a post hoc validation, not a proof of mechanism.

---

# 3. Logic / claim-scope audit

## Finding 3.1 — “Finite-size persistence” is overextended toward universality  
**Severity: Major**

The paper correctly says “not a thermodynamic limit” in several places, but then repeatedly uses language like “persists across \(N=3\text{–}5\)” and “approximately flat over the tested finite-size range” in a way that can be read as a universality hint. With only four data points and a protocol change at \(N=5\), this should be framed strictly as a finite-\(N\) observation.

## Finding 3.2 — State dependence is handled better than in many drafts, but still needs tighter wording  
**Severity: Minor**

The paper does acknowledge:
- Haar crossover at small \(\tau\),
- ground and mixed no crossover up to \(\tau=1\),
- thermal crossover only in \((0.3,1.0)\).

That is good. But the conclusion “the same state admits different optimal relational structures depending on the timescale” is only directly shown for Haar and thermal, not for all states. Ground and mixed are negative cases. This should be stated more carefully.

## Finding 3.3 — The distinction between objective crossing and optimizer basin hopping is conceptually good, but not fully demonstrated  
**Severity: Major**

The paper says the crossover is defined by \(\Delta\Phi(\tau)=0\), “never by optimizer basin hops.” That is a good methodological choice. However, the evidence for basin identity is thin: the classification is based on fixed representatives and value-based labeling, with no full landscape analysis. The statement that “many distinct local maxima sit \(\approx 0.87\) from both representatives” actually undermines the basin interpretation unless more structure is shown.

---

# 4. Falsifiability / physical meaning

## Finding 4.1 — Replacing the normalized functional is justified as a numerical fix, but the physical interpretation is still somewhat ad hoc  
**Severity: Major**

The move from
\[
\Gamma_F^{(0)} = -\frac{\mathrm{Re}\langle Q_F\rho, Q_FL(\rho)\rangle}{|Q_F\rho|^2}
\]
to
\[
J_{\mathrm{dyn}}^{(0)} = -\mathrm{Re}\langle Q_F\rho, Q_FL(\rho)\rangle
\]
is mathematically sensible if the normalized version is singular. But the paper should more explicitly justify why the unnormalized quantity is the correct physical observable rather than simply a numerically convenient surrogate. Right now the replacement reads as “the old one failed, so we adopted the one that behaved well.”

## Finding 4.2 — The Liouvillian-gap negative result is reasonable, but the comparison is too coarse  
**Severity: Minor**

You report:
\[
\tau_c\Delta_L = 0.0119 \quad (N=3), \qquad \tau_L/\tau_c \approx 83\text{–}136 \quad (N=3\text{–}5).
\]
This does support the statement that \(\tau_c\) is much shorter than the global relaxation time. However, the comparison is only to the slowest Liouvillian mode. That is enough to reject a naive one-gap hypothesis, but not enough to exclude more subtle spectral or mode-coupling explanations. The discussion should say exactly that.

---

# 5. Journal-level completeness

## Finding 5.1 — Abstract and conclusion overstate the strength of the evidence  
**Severity: Major**

Phrases like “demonstrate numerically,” “survived a battery of tests,” “derived, not fitted,” and “the canonical functional was shown to be ill-posed” are too strong relative to the limited finite-\(N\), finite-seed evidence and the missing error bars.

## Finding 5.2 — Missing reproducibility details for the key figures/tables  
**Severity: Major**

For a numerical-methods paper, the following are missing or insufficiently specified:
- exact optimizer stopping criteria,
- whether the reported \(\Phi^*\) values are best-of-seed or mean-of-seed,
- confidence intervals or seed-to-seed distributions for \(\Delta\Phi(\tau)\),
- the exact grid used for the “19 points” in Fig. 2,
- the exact procedure for identifying \(F_A\) and \(F_B\),
- the exact definition of the quotient distance \(d_F\),
- whether the “full quotient” parameterization is exhaustive or approximate.

## Finding 5.3 — Differentiation from existing work is underdeveloped  
**Severity: Minor**

The introduction situates the work in REM, but the paper does not clearly separate:
- what is new mathematically,
- what is new numerically,
- what is a protocol change,
- what is a reinterpretation of prior REM claims.

That is important because the paper explicitly says the specification is frozen and the evidence is post hoc.

---

# Axis-by-axis summary

1. **Math / definition audit:** **Major concerns**  
   Main issue: notation/derivative clarity and overclaiming from the Cauchy–Schwarz bound.

2. **Numerical evidence audit:** **Major concerns, with one Critical point**  
   Claim 4c is not sufficiently supported as “derived, not fitted.”

3. **Logic / claim-scope audit:** **Major concerns**  
   Finite-size persistence and state dependence are overextended beyond the tested regime.

4. **Falsifiability / physical meaning:** **Minor-to-Major concerns**  
   The replacement functional is numerically motivated, but the physical justification is not fully developed.

5. **Journal-level completeness:** **Major concerns**  
   Missing error budgets, reproducibility details, and stronger differentiation from prior REM work.

---

# Explicit answer to the most important question
**No, not all of Claims 1–4 actually follow from the presented evidence alone.**

- **Claim 1:** follows **approximately**, but the scope is narrower than the wording.
- **Claim 2:** the evidence supports improved numerical behavior, but not a full “well-posed” conclusion.
- **Claim 3:** the evidence supports nontrivial dependence on state/Hamiltonian in the tested cases.
- **Claim 4a–4c:** the crossover existence is plausible, but the mechanism and “derived, not fitted” claims are not yet adequately supported.

---

## Recommendation = Major revision