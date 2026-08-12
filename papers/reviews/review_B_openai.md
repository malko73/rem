## Reviewer B Report

### Overall assessment
The paper contains a coherent numerical narrative, but several central claims are only partially supported by the presented evidence, and some derivations/definitions are internally inconsistent or under-specified. The strongest results are the *existence of a numerical instability in the normalized functional* and the *empirical finite-time crossover in the tested protocol*. However, the manuscript repeatedly overstates scope beyond what the tables/figures/JSON digest establish, especially regarding “derived, not fitted,” finite-size persistence, and the universality of the crossover mechanism.

---

# 1) Math / definition audit

### Finding 1.1 — Finite-time expansion coefficient conventions are inconsistent
**Tag: Major**

The paper alternates between treating \(J_1\) as a coefficient and as a derivative. In Sec. 8:

- it writes \(J_\tau = J_0 + \tau J_1 + \tau^2 J_2 + O(\tau^3)\),
- then says “\(J_1\) is the first-order coefficient in the expansion (equivalently \(dJ_\tau/d\tau|_0 = J_1\)).”

That is fine if the convention is fixed. But the derivation of
\[
J_1 = -C^2{}''(0)/4
\]
depends on the exact relation between \(J_\tau\) and \(C^2(\tau)\). Since
\[
J_{\mathrm{dyn}}^{(\tau)} = -\frac{C^2(\tau)-C^2(0)}{2\tau},
\]
the Taylor expansion gives
\[
J_{\mathrm{dyn}}^{(\tau)} = -\frac{C^{2\,\prime}(0)}{2} - \frac{\tau}{4} C^{2\,\prime\prime}(0) + O(\tau^2),
\]
so the coefficient of \(\tau\) is indeed \(-C^{2\,\prime\prime}(0)/4\). This is mathematically consistent, but the manuscript should state explicitly that \(J_1\) is the coefficient of \(\tau\), not the derivative of \(J_\tau\) unless the factor is understood. As written, the notation is easy to misread.

### Finding 1.2 — The “canonical” projector \(Q_F\) is not fully defined consistently
**Tag: Major**

In Sec. 2.2, \(Q_F\) is described as “the projector onto the orthogonal complement of \(A_F\) acting on operators,” and then explicitly as
\[
Q_F(O) = O - \frac{1}{d_B}\mathrm{Tr}_B(O)\otimes I_B.
\]
This formula is only valid in the rotated frame and only for the specific TPS algebra in that frame. The paper also says the TPS is parameterized by a unitary \(U\) on the quotient and that the inner product is evaluated in the rotated frame with \(\rho_U = U^\dagger \rho U\). That is fine, but the manuscript never states clearly whether \(Q_F\) is a projector on the operator Hilbert space or a state-dependent dephasing map tied to the Schmidt basis of \(\rho_U\). The later use of “fixed Schmidt basis \(Q_F\)” in Sec. 2.3 and Sec. 8 suggests the latter.

This matters because the derivation of
\[
C_F^2(t)=|Q_F U^\dagger e^{tL}\rho\,U|^2
\]
assumes \(Q_F\) is fixed at \(t=0\), while the optimization over \(F\) elsewhere treats \(Q_F\) as a function of the factorization. The paper should separate:
1. the factorization-dependent projector \(Q_F\),
2. the fixed projector used in the finite-time expansion at a chosen basin representative.

### Finding 1.3 — The claim that \(J_{\mathrm{dyn}}^{(0)}\) is “bounded by Cauchy–Schwarz” is correct but incomplete
**Tag: Minor**

The inequality
\[
|J_{\mathrm{dyn}}^{(0)}| \le |Q_F\rho|\,|Q_FL(\rho)|
\]
is correct. But the statement “Since \(\dim\mathcal{H}<\infty\), \(|Q_FL(\rho)|\) is bounded along any path, so \(J_{\mathrm{dyn}}^{(0)}\to 0\) as \(|Q_F\rho|^2\to 0\)” is not the same as proving well-posedness of the optimization. It only shows the numerator vanishes in the singular direction. The actual well-posedness claim depends on the mutual-information term and on compactness of the quotient. That is plausible, but not fully demonstrated numerically.

### Finding 1.4 — The derivation of \(C^2{}''(0)\) is plausible but not fully justified from the text
**Tag: Major**

The paper states
\[
C^2{}''(0) = 2|Q_F Y_U|^2 + 2\mathrm{Re}\langle Q_F\rho_U, Q_F(U^\dagger L(L(\rho))U)\rangle.
\]
This is consistent with differentiating
\[
C^2(t)=\langle Q_F\rho_U(t),Q_F\rho_U(t)\rangle
\]
if \(Q_F\) is fixed and \(\dot\rho_U(t)=U^\dagger L(\rho(t))U\). However, the text does not address whether \(Q_F\) commutes with the time evolution in the rotated frame or whether the basis is frozen at \(t=0\). The manuscript later says it is fixed at the \(t=0\) Schmidt basis, which resolves this, but the derivation should state that assumption up front.

---

# 2) Numerical evidence audit

### Finding 2.1 — Claim 1 is supported numerically, but the evidence is protocol-specific
**Tag: Major**

The digest gives:
- \(\mathrm{corr}(\Phi,\log_{10} C_F^2) = -0.998955\),
- \(\Phi_{\max}=25.7133\),
- \(\Gamma_{\min}=-128.5251\),
- \(C_{F,0}^2{}_{\min}=4.093466\times10^{-4}\),
- \(p_{\max}=0.999795\).

This strongly supports the statement that the normalized functional can drive the optimizer toward near-product, singular directions in the tested Haar/\(N=3\) protocol. So the *existence* of the pathology is well supported.

However, the paper’s wording “the normalized local decay rate is singular on the unrestricted quotient” is stronger than the evidence. The data show a singularity in the tested protocol and a transitivity argument suggests the quotient contains arbitrarily close product directions for pure states. But the numerical evidence alone does not establish that the optimizer always exploits this, nor that the singularity is unavoidable for all states or all implementations.

### Finding 2.2 — Claim 2 is supported as a numerical comparison, but “well-posed” is stronger than the data justify
**Tag: Major**

The table in Sec. 4 shows under \(J_{\mathrm{dyn}}^{(0)}\):
- ground: \(\Phi^*=1.8995\), std \(0.0007\),
- haar: \(\Phi^*=1.8594\), std \(0.0090\),
- mixed: \(\Phi^*=1.9845\), std \(0.0051\),
- thermal: \(\Phi^*=1.9465\), std \(0.0058\),

with minimum \(C_F^2(0)=0.0956\), maximum \(|\gamma_{\mathrm{ref}}|=1.66\), and \(p_{\max}\in[0.51,0.58]\).

This supports “the singularity/product-collapse pathology disappears in the tested protocol.” But “a maximizer exists on the compact quotient by continuity” is a mathematical statement, not a numerical one, and the paper does not prove the continuity of the full objective in the exact quotient parameterization used in code. The evidence is consistent with well-posedness, but the claim should be phrased as “no pathology was observed in the tested protocol.”

### Finding 2.3 — Claim 3 is supported, but the diagonal-dominance argument is weaker than implied
**Tag: Major**

For D1-R, the cross-evaluation matrix in Supplementary Table S3 does show diagonal dominance in the sense that each row’s diagonal entry is the row maximum. The reported mean gap is \(0.113\). For D2-R, the state cross-evaluation matrix also shows diagonal dominance with mean gap \(0.159\).

This supports the claim that the selected structure depends on both Hamiltonian and state. However:
- the paper does not define the “mean gap” precisely enough to verify the reported \(0.113\) and \(0.159\) from the tables alone;
- the statement “physically related families are closer in the quotient” is not directly supported by the matrix and distances given;
- the phrase “not a state-independent common solution” is justified, but only for the tested families/states.

### Finding 2.4 — Claim 4a is supported, but the interpolation and basin language need tighter evidence
**Tag: Major**

The D3 table gives:
\[
\Delta\Phi(\tau)= -0.021,-0.011,-0.003,0.000,0.003,0.009,0.016,0.039,0.079
\]
at \(\tau=0.003,0.01,0.016,0.018,0.02,0.025,0.03,0.05,0.1\).

This clearly supports a sign change near \(\tau\approx 0.018\). The digest gives \(\tau_c=0.018024\), consistent with linear interpolation. So the crossing itself is supported.

But the claim that “the instantaneous optimum and the finite-time optimum belong to distinct structural basins” is less directly supported. The paper uses value-based classification and fixed representatives \(F_A\), \(F_B\), but the evidence shown is objective crossing, not a full landscape analysis. The statement that the optimizer “follows sharply” from 0 to 6 seeds between \(\tau=0.018\) and \(0.02\) is suggestive, but not sufficient to establish basin identity in a rigorous sense.

### Finding 2.5 — Claim 4b is only partially supported; the “finite-size persistence” wording overreaches
**Tag: Major**

The paper reports:
- \(N=3\): \(\tau_c=0.0180\),
- \(N=4\), \(2|2\): \(\tau_c=0.0206\),
- \(N=4\), \(1|3\): \(\tau_c=0.0224\),
- \(N=5\), full: \(\tau_c=0.0217\),
and also a restricted-domain \(N=5\) subspace value \(0.0173\).

This supports a finite-\(N\) observation of similar timescales. But the manuscript repeatedly says “persists” and “approximately flat” as if this were a robust trend. With only four systems and one of them having a protocol-dependent alternative value (\(0.0173\) vs \(0.0217\)), the evidence is too thin to support anything stronger than “the tested values are of the same order and vary modestly.”

### Finding 2.6 — Claim 4c is numerically supported, but “derived, not fitted” is only partly true
**Tag: Major**

The table in Sec. 8 gives:
- \(N=3\): \(\tau_c^{(1)}=0.0171\), measured \(0.0180\), error \(5.1\%\);
- \(N=4\), \(2|2\): \(0.0188\) vs \(0.0206\), error \(8.6\%\);
- \(N=4\), \(1|3\): \(0.0208\) vs \(0.0224\), error \(7.2\%\);
- \(N=5\): \(0.0193\) vs \(0.0217\), error \(10.9\%\).

This supports the first-order formula as a reasonable predictor. The analytic and numerical slopes also agree closely:
- \(1.521/1.526\),
- \(0.526/0.528\),
- \(1.633/1.639\),
- \(1.130/1.136\).

However, “derived, not fitted” is too strong unless the representative pair \((F_A,F_B)\) is fixed independently of the crossing analysis. The paper says \(F_A\) is the optimum at \(\tau=10^{-3}\) and \(F_B\) at \(\tau=10^{-1}\), which is acceptable, but the choice of those representatives is still part of the protocol and affects the result by up to \(\lesssim 4\times10^{-3}\). So the derivation is not purely parameter-free; it is derived conditional on a representative selection procedure.

---

# 3) Logic / claim-scope audit

### Finding 3.1 — “Finite-size persistence” is overextended toward trend/scaling language
**Tag: Major**

The paper correctly says “not a scaling law and not a thermodynamic limit” in places, but elsewhere it says:
- “approximately flat over the tested finite-size range,”
- “persistence across \(N=3\text{–}5\),”
- “size effect rather than shape artifact,”
- “the crossover scale is therefore set by the objective balance ... not by the global Liouvillian relaxation time.”

The first two are acceptable as finite-\(N\) observations. The last two are stronger and should be restricted. The data do not establish universality, asymptotic scaling, or thermodynamic persistence. They only show a finite-size trend over four cases.

### Finding 3.2 — State dependence is handled better than in many such papers, but still needs tighter wording
**Tag: Minor**

The manuscript does acknowledge that:
- Haar shows a crossover near \(0.018\),
- thermal crosses between \(0.3\) and \(1.0\),
- ground and mixed show no crossover up to \(1.0\).

That is good. But the discussion sometimes generalizes from Haar to “the selected structure depends on the dynamical observation scale” without emphasizing that this is state-dependent and not universal across all states. The paper does eventually say this, but the abstract and conclusion should be more careful.

### Finding 3.3 — The “first-order-like” language is acceptable only as an analogy
**Tag: Minor**

The paper explicitly says “not a phase transition,” which is good. Still, “first-order-like basin crossover” can be misleading because the evidence is from a finite-dimensional optimization landscape, not a thermodynamic order parameter. The analogy is fine, but should be clearly labeled as metaphorical.

---

# 4) Falsifiability / physical meaning

### Finding 4.1 — The replacement of the normalized functional is justified, but the narrative is somewhat ad hoc
**Tag: Major**

The paper argues that the normalized functional \(\Gamma_F\) is singular and therefore should be replaced by the unnormalized \(J_{\mathrm{dyn}}^{(0)}\). This is physically sensible: the unnormalized quantity is the instantaneous decay rate of the residual norm, while the normalized one can blow up as the residual vanishes.

That said, the manuscript frames this as a “falsification-then-replacement arc” and “canonical functional” without fully discussing whether the normalization might still be useful as a diagnostic on a restricted domain. The paper does mention “fixed-\(F\) local diagnostic only,” which is good, but the replacement still reads somewhat ad hoc unless the intended domain of validity is sharply defined.

### Finding 4.2 — The Liouvillian-gap negative result is handled appropriately, but the comparison is limited
**Tag: Minor**

The paper is careful to say the slowest Liouvillian mode is not the timescale controlling \(\tau_c\), and it explicitly notes that more subtle spectral or mode-coupling explanations are not excluded. That is good scientific practice.

However, the quantitative comparison is thin:
\[
\tau_c \Delta_L = 0.0119 \quad (N=3), \qquad \tau_L/\tau_c \approx 83\text{–}136.
\]
This shows a mismatch, but the paper should avoid implying that the Liouvillian gap hypothesis is broadly falsified beyond the tested protocol. It is only rejected as a simple explanation here.

---

# 5) Journal-level completeness

### Finding 5.1 — Abstract and conclusion overstate certainty relative to the evidence
**Tag: Major**

The abstract says:
- “We show that the natural normalized decay-rate functional is singular on the unrestricted quotient,”
- “The selected factorization responds systematically to both the Hamiltonian and the state,”
- “The crossover scale is therefore set by the objective balance ... not by the global Liouvillian relaxation time.”

These are stronger than the evidence warrants. The paper shows this in the tested protocols, not as a general theorem. The conclusion similarly reads too definitively.

### Finding 5.2 — The paper is reasonably complete on methods, but reproducibility is not fully sufficient from the manuscript alone
**Tag: Major**

Positives:
- explicit Hamiltonian and Lindbladian,
- seed counts,
- optimizer settings,
- interpolation rule,
- representative selection procedure,
- repository and script names.

Concerns:
- the exact definition of the quotient distance \(d_F\) is not fully operationalized in the manuscript;
- the “mean gap” and “diagonal dominance” metrics are not defined precisely enough;
- the JSON digest is helpful, but the paper still relies on supplementary tables for key claims;
- the \(N=5\) protocol changed from a 200-dim random subspace to the full 945-dim quotient, which complicates reproducibility and interpretation.

### Finding 5.3 — Differentiation from existing work is only partially convincing
**Tag: Minor**

The paper claims novelty in the falsification of the normalized rate, the replacement functional, and the finite-time crossover mechanism. That is plausible, but the relation to prior REM work is described mostly in programmatic terms. A journal reader would benefit from a sharper statement of what is new relative to Spec v2.2 and earlier REM papers.

---

# MOST IMPORTANT QUESTION

**If you ignore the authors' interpretation and inspect only the definitions, derivations, numerical evidence, and figures, do Claims 1-4 actually follow from the presented evidence?**

**Answer: Partially, but not fully.**

- **Claim 1:** Yes, in the tested Haar/\(N=3\) protocol the evidence strongly supports that the normalized functional can become singular and drive the optimizer toward near-product directions. But the claim should be scoped to the tested protocol, not stated as a general unrestricted theorem.
- **Claim 2:** Mostly yes for the tested protocol: the unnormalized functional removes the observed singularity/product-collapse/seed-instability. But “well-posed” is stronger than the evidence alone establishes.
- **Claim 3:** Yes, the cross-evaluation matrices support state- and Hamiltonian-dependence of the selected structure in the tested families/states.
- **Claim 4a–c:** The existence of a finite-time crossover and the first-order predictor are supported numerically. However, “finite-size persistence” and “derived, not fitted” are overstated unless carefully qualified. The evidence supports a finite-\(N\) crossover trend, not a universal mechanism or scaling law.

---

# Summary of findings by severity

## Critical
- None. The main conclusions are not outright refuted by the presented evidence.

## Major
1. Inconsistent/under-specified conventions for \(J_1\), \(Q_F\), and the finite-time expansion.
2. Claim 2’s “well-posed” language exceeds the numerical evidence.
3. Claim 4b overstates “finite-size persistence” and risks implying scaling/universality.
4. Claim 4c’s “derived, not fitted” is only conditionally true because representative selection and interpolation still matter.
5. Abstract/conclusion overstate generality relative to the finite tested protocols.

## Minor
1. “First-order-like” should be explicitly labeled as an analogy.
2. The Liouvillian-gap negative result is fine but should be framed more narrowly.
3. Definitions of mean gap / diagonal dominance / quotient distance should be made more operational.
4. Prior-work differentiation could be sharper.

---

## Recommendation = Major revision