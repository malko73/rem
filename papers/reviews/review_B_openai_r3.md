## Reviewer B report

### Overall assessment
The manuscript presents an interesting numerical narrative, but the evidence is not yet strong enough to support several of the stronger claims as written. The main issues are: (i) some derivations are internally inconsistent or under-specified, especially in the finite-time expansion; (ii) several claims rely on post-hoc representative choices and fixed-basin comparisons that are not fully reproducible from the paper alone; (iii) the scope of the conclusions is repeatedly broadened beyond what the finite-size data justify.

**Bottom line:** the paper contains potentially publishable numerical observations, but the current version overstates what is established.

---

# 1) Math / definition audit

### 1.1 Canonical functional \(J_{\mathrm{dyn}}^{(0)}\)
The definition
\[
J_{\mathrm{dyn}}^{(0)}(F;\rho,L) = -\operatorname{Re}\langle Q_F\rho, Q_F L(\rho)\rangle_{\mathrm{HS}}
\]
is clear and consistent with the later use in Sec. 2.2 and Sec. 4.

**Minor:** The statement “Since \(\dim\mathcal{H}<\infty\), \(|Q_F L(\rho)|_{\mathrm{HS}}\) is bounded along any path, so \(J_{\mathrm{dyn}}^{(0)}\to 0\) as \(|Q_F\rho|^2\to 0\)” is not a derivation of continuity in the quotient variable \(U\); it is a bound in operator norm. The conclusion is plausible, but the paper should explicitly state the continuity assumptions on \(Q_F(U)\) and \(L(\rho)\) as functions of \(U\).

### 1.2 Finite-time extension \(J_{\mathrm{dyn}}^{(\tau)}\)
The definition
\[
J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{C_F^2(\tau)-C_F^2(0)}{2\tau},\qquad
C_F^2(t)=|Q_F U^\dagger e^{tL}\rho\,U|^2
\]
is fine as a finite-difference average decay rate.

However, the paper’s derivative identity and coefficient conventions need tightening.

You state:
\[
\frac{d}{dt}C_F^2(t)\big|_0 = 2\operatorname{Re}\langle Q_F\rho, Q_F L(\rho)\rangle
\]
and then
\[
\lim_{\tau\to 0} J_{\mathrm{dyn}}^{(\tau)} = J_{\mathrm{dyn}}^{(0)}.
\]
This is consistent.

But in Claim 4c you then write:
\[
C^2{}''(0) = 2|Q_F Y_U|^2 + 2\operatorname{Re}\langle Q_F\rho_U, Q_F(U^\dagger L(L(\rho))U)\rangle,
\]
and therefore
\[
J_1 = -C^2{}''(0)/4.
\]

This is **not fully justified as written**. If
\[
J_\tau = -\frac{C^2(\tau)-C^2(0)}{2\tau},
\]
then Taylor expanding
\[
C^2(\tau)=C^2(0)+\tau C^{2\prime}(0)+\frac{\tau^2}{2}C^{2\prime\prime}(0)+\cdots
\]
gives
\[
J_\tau = -\frac{1}{2}C^{2\prime}(0)-\frac{\tau}{4}C^{2\prime\prime}(0)+O(\tau^2).
\]
So if \(J_\tau = J_0+\tau J_1+\cdots\), then indeed
\[
J_0=-\frac12 C^{2\prime}(0),\qquad J_1=-\frac14 C^{2\prime\prime}(0).
\]
That part is correct.

**Major:** The paper does not show the derivation of the second derivative formula from the product rule in enough detail to verify the sign and factor conventions. Since Claim 4c depends on this, the derivation should be written out explicitly.

### 1.3 \(\tau_c^{(1)}\) derivation
The first-order formula
\[
\tau_c^{(1)}=-\frac{\Delta\Phi_0}{\Delta\Phi_1}
\]
is mathematically standard and follows from \(\Delta\Phi(\tau)=\Delta\Phi_0+\tau\Delta\Phi_1+\cdots\).

**Minor:** The paper alternates between calling \(J_1\) a “coefficient” and a derivative. This is acceptable if stated clearly, but the notation should be standardized to avoid ambiguity.

### 1.4 Internal consistency of the finite-time expansion
There is a notable inconsistency between the main text and the supplementary sensitivity note:

- Main text says the curve is smooth and monotone in Fig. 2 and that the quadratic correction is a “consistency check.”
- Later, the sensitivity note says \(\Delta\Phi_2<0\) “consistent with the non-monotone \(\Delta\Phi\) seen in D4.”

These are not the same statement. A negative quadratic coefficient does **not** by itself imply non-monotonicity on the interval shown.

**Major:** The paper needs to reconcile whether \(\Delta\Phi(\tau)\) is monotone on the plotted interval or whether it bends back at larger \(\tau\). The current wording is inconsistent.

---

# 2) Numerical evidence audit

## Claim 1: normalized functional is singular
The evidence is strong for the tested case \(N=3\), Haar state, full quotient:

- \(\mathrm{corr}(\Phi,\log_{10} C_F^2)=-0.998955\)
- \(\Phi_{\max}=25.7133\)
- \(\Gamma_{\min}=-128.5251\)
- \(C_{F,0}^2{}_{\min}=4.093466\times10^{-4}\)
- \(p_{\max}=0.999795\)

These numbers do support the claim that the normalized objective can be driven toward a near-product singular direction in the tested protocol.

**Major:** The paper overgeneralizes from one tested setting. The claim says “the mechanism is general for pure states,” but the numerical evidence only demonstrates the phenomenon for Haar states at \(N=3\). The transitivity argument shows existence of product-like directions, not that the optimizer generically finds them or that \(\Phi\to+\infty\) occurs for all pure states.

## Claim 2: unnormalized functional is well-posed
The table gives:

- ground: \(\Phi^*=1.8995\), std \(0.0007\)
- haar: \(\Phi^*=1.8594\), std \(0.0090\)
- mixed: \(\Phi^*=1.9845\), std \(0.0051\)
- thermal: \(\Phi^*=1.9465\), std \(0.0058\)

and the listed minima:
- minimum \(C_F^2(0)=0.0956\)
- maximum \(|\gamma_{\mathrm{ref}}|=1.66\)
- \(p_{\max}\in[0.51,0.58]\)

This supports the narrower statement that the optimizer did not encounter the same singular collapse in the tested runs.

**Major:** “Well-posed” is stronger than “no pathology observed in the tested protocol.” The existence of a maximizer on a compact quotient is a mathematical statement, but the paper does not prove continuity of the full objective in the parameterization used numerically, nor does it show that the optimizer is not missing sharper maxima. The evidence supports “numerically stable in the tested runs,” not a full well-posedness theorem.

**Minor:** The phrase “seed instability disappear” is too absolute. The Haar std is \(0.0090\), which is small, but not zero.

## Claim 3: structure responds to Hamiltonian and state
The cross-evaluation matrices do show nontrivial dependence:

- D1-R mean gap \(0.113\)
- D2-R mean gap \(0.159\)

and the diagonal entries are generally larger than off-diagonal entries.

This supports the claim that the selected structure is not universal across all Hamiltonians or states.

**Major:** The paper’s language “\(F^*=F^*(\rho,L,\lambda)\), not a state-independent common solution” is stronger than the evidence. The data show dependence in the tested families/states, but not that no common solution exists in some broader sense. Also, the ground/mixed near-degeneracy noted in the text weakens the claim of strong state separation.

**Minor:** The cross-evaluation matrices are informative, but the paper should report the full matrices for the Hamiltonian-response case, not only the mean gaps, if the claim is to be independently checked.

## Claim 4a: finite-time structural crossover
The D3 table and Supplementary Table S2 support a sign change in \(\Delta\Phi(\tau)\):

- \(\Delta\Phi(0.016)=-0.003\)
- \(\Delta\Phi(0.018)=-0.000\)
- \(\Delta\Phi(0.020)=+0.0027\)

with reported \(\tau_c=0.018024\) and \(d_F(F_A,F_B)=0.8858\).

This is good evidence for a crossover in the fixed-representative objective difference.

**Major:** The claim that “the selected structure depends on the dynamical observation scale” is only partially supported. The evidence shows that the objective difference between two fixed representatives crosses zero. It does **not** by itself prove that the optimizer’s selected structure changes in a physically meaningful way, because the basin classification is explicitly acknowledged to be value-based and not a full landscape analysis.

**Major:** The paper says the crossover is “first-order-like” and that adjacent-\(\tau\) \(d_F\) jumps are \(O(1)\), but the evidence for this is sparse and partly qualitative. The actual \(d_F\) values at the relevant \(\tau\) points are not tabulated in the main text.

## Claim 4b: finite-size persistence
The reported \(\tau_c\) values are:

- \(N=3\): \(0.0180\)
- \(N=4\), \(2|2\): \(0.0206\)
- \(N=4\), \(1|3\): \(0.0224\)
- \(N=5\), \(2|3\): \(0.0217\)

This does support a finite-size trend over \(N=3\)–5.

**Major:** The claim “persists as the Hilbert-space size increases” is too broad. The data cover only four systems, with one \(N=5\) point and two different \(N=4\) cuts. This is finite-size evidence, not persistence under size increase in any asymptotic sense.

**Major:** The paper mixes “full quotient” and “restricted-domain history” for \(N=5\). The main claim uses the full-quotient value \(0.0217\), but the supplementary notes that the earlier 200-dim subspace gave \(0.0173\). This protocol dependence should be emphasized more strongly because it affects reproducibility and the interpretation of “persistence.”

## Claim 4c: derived, not fitted
The table gives:

- \(N=3\): \(\tau_c^{(1)}=0.0171\), measured \(0.0180\), error \(5.1\%\)
- \(N=4\), \(2|2\): \(0.0188\) vs \(0.0206\), error \(8.6\%\)
- \(N=4\), \(1|3\): \(0.0208\) vs \(0.0224\), error \(7.2\%\)
- \(N=5\), \(2|3\): \(0.0193\) vs \(0.0217\), error \(10.9\%\)

This is reasonable evidence that the first-order formula is predictive at the 5–11% level.

**Major:** “Derived, not fitted” is only partly justified. The first-order formula is indeed derived from the expansion, but the paper also states that the representative choice \(F_B\) comes from optimization at \(\tau=0.1\), and the measured \(\tau_c\) depends on interpolation on a finite grid. The result is therefore derived from a chosen basin pair, not a fully prediction-only quantity.

**Major:** The quadratic correction is described as a “consistency check,” but the paper also reports it as reducing the residual to \(0.0\)–\(0.5\%\). That is fine, but it should not be presented as independent validation because it is fit-based.

---

# 3) Logic / claim-scope audit

**Major:** The manuscript repeatedly slides from “tested finite-size trend” to “persistence” and from “objective crossing” to “structural crossover” to “selected structure depends on timescale.” These are not equivalent. The evidence supports a crossover in the objective evaluated on fixed representatives; it does not fully establish a universal structural phase-like phenomenon.

**Major:** The state dependence is handled inconsistently. The paper says:
- Haar and thermal show crossovers;
- ground and mixed do not up to \(\tau=1\).

That is fine. But then it also says the selected structure responds systematically to both Hamiltonian and state, and that the crossover is a general feature of relational structure. The latter is too broad.

**Minor:** The statement that the normalized functional is “retained only in its proper role — fixed-\(F\) local diagnostic only” is a reasonable interpretation, but it should be framed as a recommendation from this dataset, not as a proven theorem.

---

# 4) Falsifiability / physical meaning

The replacement of the normalized functional by the unnormalized one is conceptually plausible and numerically motivated.

**Major:** The argument risks looking ad hoc unless the paper more clearly distinguishes:
1. a mathematical singularity of the normalized ratio on the quotient,
2. a numerical optimization pathology,
3. a physical interpretation of the unnormalized decay rate.

Right now the paper asserts all three at once, but the evidence directly supports only (1) and (2) in the tested case.

**Major:** The Liouvillian-gap negative result is useful, but the comparison is incomplete. The paper reports \(\tau_c\Delta_L=0.0119\) and \(\tau_L/\tau_c\approx 83\text{–}136\), which indeed shows a large separation. However, this only rules out the slowest global Liouvillian mode as the controlling scale. It does not exclude mode coupling, local spectral features, or basin-specific dynamical rates. The discussion acknowledges this, which is good, but the abstract still sounds more definitive than the evidence warrants.

---

# 5) Journal-level completeness

**Minor:** The abstract is clear and informative, but it overstates the generality of the conclusions.

**Major:** The introduction and discussion do not sufficiently differentiate this work from prior REM results beyond saying the specification is frozen and the evidence is post-hoc. The reader needs a sharper statement of what is genuinely new numerically versus what is a reanalysis of existing protocol outputs.

**Major:** Reproducibility is not yet fully adequate from the paper alone. Missing or under-specified items include:
- exact optimizer stopping criteria beyond fixed step counts;
- full per-seed outputs for the key claims;
- the complete D1-R and D2-R cross-evaluation matrices;
- the exact procedure for “value-based classification” of basins;
- the full derivation of the second derivative formula used in Claim 4c;
- the precise interpolation rule for \(\tau_c\) on the grid.

**Minor:** The figures are described well in captions, but the main text still relies on supplementary tables for essential evidence. That is acceptable, but the paper should make the main claims readable from the figures/tables alone more explicitly.

---

# MOST IMPORTANT QUESTION

**If I ignore the authors' interpretation and inspect only the definitions, derivations, numerical evidence, and figures, do Claims 1–4 actually follow from the presented evidence?**

**Answer: partially, but not fully.**

- **Claim 1:** Yes, for the tested Haar \(N=3\) protocol the evidence strongly supports singular behavior of the normalized functional.  
  But the broader “general for pure states” wording is not established.

- **Claim 2:** The evidence supports that the unnormalized functional is numerically stable and avoids the observed singular collapse in the tested runs.  
  But “well-posed” as a general statement is stronger than the evidence shown.

- **Claim 3:** Yes, the cross-evaluation tables support state- and Hamiltonian-dependence in the tested families.  
  But the claim should remain explicitly finite-protocol and finite-family.

- **Claim 4a–4c:** The evidence supports a fixed-representative objective crossing and a first-order estimate of \(\tau_c\) with 5–11% error.  
  But the stronger interpretation as a structural/timescale-dependent phenomenon is only partially supported, and the “derived, not fitted” phrasing needs qualification because representative selection and interpolation are part of the procedure.

---

# Findings summary

### Critical
- None that force immediate rejection on the basis of a mathematical contradiction alone.

### Major
1. Claim 1 overgeneralizes from a single tested Haar \(N=3\) case to “general for pure states.”
2. Claim 2 uses “well-posed” too strongly; the evidence shows numerical stability, not a full general proof.
3. Claim 4a/4b overstate structural/timescale persistence beyond the finite-size, finite-protocol data.
4. Claim 4c is “derived” only after fixing representatives and interpolation; this should be stated more carefully.
5. The finite-time expansion derivation, especially the second derivative formula, needs a more explicit derivation and consistency check.
6. The manuscript’s monotonicity/non-monotonicity statements about \(\Delta\Phi(\tau)\) are internally inconsistent.

### Minor
1. Standardize coefficient vs derivative notation for \(J_1\).
2. Provide full matrices and per-seed outputs in the main supplement.
3. Tighten wording around seed stability and continuity.
4. Clarify the exact basin-classification procedure and interpolation rule.

---

## Recommendation = Major revision