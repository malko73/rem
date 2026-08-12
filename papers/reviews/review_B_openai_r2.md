## Reviewer B report

### Overall assessment
The paper presents an interesting numerical story, but the evidence is not yet strong enough to support the breadth of the claims as written. The most serious issues are: (i) several derivations are internally inconsistent or at least not fully aligned with the stated definitions; (ii) the finite-time crossover evidence is based on fixed representatives and interpolation, but the paper sometimes speaks as if it establishes a more general structural phenomenon than the data warrant; (iii) the “derived, not fitted” language overstates what is actually shown, because representative choice and grid resolution still enter materially; and (iv) the negative Liouvillian-gap result is too limited to justify the stronger mechanistic exclusion implied in the abstract and conclusion.

Below I audit the paper along the requested five axes, with findings tagged **Critical / Major / Minor**.

---

## 1) Math / definition audit

### 1.1 Normalized vs unnormalized dynamical functional
The paper defines
\[
\Gamma_F^{(0)} = -\frac{\mathrm{Re}\langle Q_F\rho, Q_FL(\rho)\rangle}{|Q_F\rho|^2}
\]
in Sec. 3, and
\[
J_{\mathrm{dyn}}^{(0)}(F;\rho,L) = -\operatorname{Re}\langle Q_F\rho, Q_F L(\rho)\rangle
\]
in Sec. 2.2.

This is mathematically coherent as a replacement, but the paper’s claim that the normalized functional is “singular on the unrestricted quotient” is only partially established by the presented evidence. The argument “\(\mathcal{U}(d)\) acts transitively on pure states” does imply that for a pure state one can reach product states under some factorization, but the paper does not prove that the optimizer over the quotient can approach arbitrarily small \(|Q_F\rho|^2\) while keeping the numerator sign favorable enough to drive \(\Phi\to+\infty\). The numerical evidence in D2.1 is strong, but the general mathematical statement is stronger than the evidence.

**Finding: Major** — The singularity mechanism is plausible and numerically supported, but the paper overstates the generality of the proof. The text should distinguish “observed in the tested protocol” from “proved on the unrestricted quotient.”

### 1.2 Cauchy–Schwarz and boundedness
The paper states:
\[
|J_{\mathrm{dyn}}^{(0)}| \le |Q_F\rho|\,|Q_FL(\rho)|
\]
and then concludes that since \(\dim\mathcal H<\infty\), \( |Q_FL(\rho)|\) is bounded along any path, so \(J_{\mathrm{dyn}}^{(0)}\to 0\) as \(|Q_F\rho|^2\to 0\).

This is fine as a continuity statement, but the phrase “vanishes smoothly in the singular direction” is not a theorem as written unless one specifies the path and the continuity of \(Q_F\) in the quotient coordinates. The paper does not show a uniform bound on \(Q_FL(\rho)\) over the quotient, only finite-dimensional boundedness in a loose sense.

**Finding: Minor** — The argument is acceptable but should be stated more carefully.

### 1.3 Finite-time extension and derivative identity
The paper defines
\[
J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{C_F^2(\tau)-C_F^2(0)}{2\tau},
\qquad
C_F^2(t)=|Q_F U^\dagger e^{tL}\rho\,U|^2.
\]
It then states
\[
\frac{d}{dt}C_F^2(t)\big|_{0}=2\mathrm{Re}\langle Q_F\rho,Q_FL(\rho)\rangle
\]
and therefore
\[
\lim_{\tau\to 0}J_{\mathrm{dyn}}^{(\tau)}=J_{\mathrm{dyn}}^{(0)}.
\]

This is consistent with the definitions, assuming \(L\) is the generator acting linearly on \(\rho\). However, in Sec. 8 the paper writes
\[
Y_U = U^\dagger L\rho_0 U
\]
which is notationally inconsistent with the earlier \(L(\rho)\). Also, the second derivative formula
\[
C^2{}''(0) = 2|Q_F Y_U|^2 + 2\mathrm{Re}\langle Q_F\rho_U, Q_F(U^\dagger L^2\rho_0 U)\rangle
\]
is plausible, but the paper does not show the intermediate steps, and the notation \(L^2\rho_0\) is ambiguous for a Lindbladian superoperator unless explicitly defined as composition \(L(L(\rho_0))\).

**Finding: Major** — The first-derivative consistency is fine, but the second-derivative derivation is underexplained and notation is inconsistent enough to hinder reproducibility.

### 1.4 \(\tau_c^{(1)}\) and \(\tau_c^{(2)}\)
The paper uses
\[
\tau_c^{(1)}=-\frac{\Delta\Phi_0}{\Delta\Phi_1}
\]
and reports values:
- \(N=3\): \(0.0171\) vs measured \(0.0180\)
- \(N=4,2|2\): \(0.0188\) vs \(0.0206\)
- \(N=4,1|3\): \(0.0208\) vs \(0.0224\)
- \(N=5,2|3\): \(0.0193\) vs \(0.0217\)

These are numerically consistent with the formula. The quadratic correction values
\(0.0180, 0.0206, 0.0224, 0.0216\) also match the reported measured values closely.

However, the paper’s derivation of \(J_1=-C^2{}''(0)/4\) is not fully justified from the finite-time definition as written. Since
\[
J_{\mathrm{dyn}}^{(\tau)} = -\frac{C^2(\tau)-C^2(0)}{2\tau},
\]
a Taylor expansion gives
\[
J_{\mathrm{dyn}}^{(\tau)} = -\frac{C'^2(0)}{2} - \frac{\tau C''^2(0)}{4} + O(\tau^2),
\]
so the coefficient convention is fine, but the paper should explicitly distinguish the coefficient \(J_1\) from the derivative \(dJ/d\tau|_0\). It does say this, but only after introducing the formula, which is easy to misread.

**Finding: Minor** — The formulas are likely correct, but the presentation is too compressed for a derivation that is central to Claim 4c.

---

## 2) Numerical evidence audit

### Claim 1: normalized functional is singular
Evidence given:
- \(\mathrm{corr}(\Phi,\log_{10} C_F^2)=-0.998955\)
- \(\Phi_{\max}=25.7133\)
- \(\Gamma_{\min}=-128.5251\)
- \(C_{F,0}^2{}_{\min}=4.093466\times10^{-4}\)
- \(p_{\max}=0.999795\)

This is strong evidence that the optimizer exploits the denominator in the tested Haar \(N=3\) case. The figure caption also states the best trial reaches \(C_F^2(0)=6\times10^{-4}\) with \(\Phi=25.7\), versus \(\Phi\approx1.86\) under the unnormalized functional.

But the paper’s claim that the normalized functional is singular “on the unrestricted quotient” is broader than the evidence. The evidence is one state family, one size, one protocol. The singularity direction sweep is mentioned but not quantified in the text.

**Finding: Major** — The numerical evidence supports a protocol-specific falsification, not a universal theorem.

### Claim 2: unnormalized functional is well-posed
Evidence given:
- Table in Sec. 4: \(\Phi^*\) under \(J^{(0)}\) is \(1.8995, 1.8594, 1.9845, 1.9465\) for ground, Haar, mixed, thermal.
- Minimum \(C_F^2(0)=0.0956\)
- Maximum \(|\gamma_{\mathrm{ref}}|=1.66\)
- Schmidt \(p_{\max}\in[0.51,0.58]\)
- 30-seed Haar run: best \(1.8595\), std \(0.0077\)

This supports the claim that the optimizer is stable in the tested protocol and does not collapse to product states. However, “well-posed on the unrestricted quotient” is stronger than what is shown. The paper shows no singularity in the tested runs, not a proof of global absence of singular directions. Also, the table mixes “\(\Phi^*\) under \(\Gamma\)” and “\(\Phi^*\) under \(J^{(0)}\)” but does not report the corresponding objective distributions or convergence diagnostics.

**Finding: Major** — The evidence supports numerical stability in the tested cases, not a global well-posedness theorem.

### Claim 3: selected structure responds to Hamiltonian and state
Evidence given:
- Five Hamiltonian families with \(\Phi^*\): 1.8994, 1.8246, 1.8695, 1.8740, 1.8481
- Cross-evaluation mean gap \(0.113\)
- State table: ground 1.8995, Haar 1.8594, mixed 1.9845, thermal 1.9465
- State cross-evaluation mean gap \(0.159\)

This does show that the optimum depends on the Hamiltonian family and on the state, at least in the tested set. The diagonal dominance of the cross-evaluation matrices is a reasonable indicator. However, the paper also notes near-degeneracies:
- ground vs mixed under the ground-state objective: \(1.900\) vs \(1.899\)
- mixed vs thermal objective closeness despite different structures

These caveats are important and should be emphasized more strongly. The claim is supported, but only as a finite-sample statement.

**Finding: Minor** — Supported, but the wording should remain explicitly finite-sample and avoid implying a stronger identifiability result.

### Claim 4a: finite-time structural crossover
Evidence given:
- \(\Delta\Phi(\tau)\) table: \(-0.021,-0.011,-0.003,0.000,0.003,0.009,0.016,0.039,0.079\)
- \(\tau_c=0.0180\)
- \(d_F(F_A,F_B)=0.8858\)
- seed flip from A to B between \(\tau=0.018\) and \(0.02\)

This supports a crossover in the fixed-representative objective difference. The sign change is clear. However, the paper also says “the optimizer follows sharply” and “first-order-like basin crossover.” That is a stronger dynamical statement than the evidence strictly supports, because the basin labels are value-based and the landscape is not fully mapped. The paper itself admits that many distinct local maxima sit \(\approx 0.87\) from both representatives.

**Finding: Major** — The objective crossing is supported; the “first-order-like basin crossover” language is more interpretive than demonstrated.

### Claim 4b: finite-size persistence
Evidence given:
- \(\tau_c = 0.0180, 0.0206, 0.0224, 0.0217\) for \(N=3,4,4,5\)
- “approximately flat” over \(N=3\text{–}5\)
- full-quotient \(N=5\) closure: \(d_F=0.319\), seed std \(2\times10^{-4}\)

This is enough to say the crossover persists in the tested finite-size range. It is not enough to infer any scaling law or universality. The paper does say it is not a thermodynamic limit, which is good. But the abstract’s phrasing “approximately flat over the tested finite-size range” is fine; the stronger implication that the scale is “set” in a size-independent way is not established.

**Finding: Minor** — Supported as a finite-\(N\) observation only.

### Claim 4c: derived, not fitted
Evidence given:
- \(\tau_c^{(1)}\) values and measured values
- relative errors 5.1%, 8.6%, 7.2%, 10.9%
- quadratic errors 0.1%, 0.0%, 0.2%, 0.5%
- \(\Delta\Phi_0,\Delta\Phi_1,\Delta\Phi_2\) listed
- sensitivity notes: representative shifts \(\lesssim 4\times10^{-3}\), grid resolution \(\pm 5\times10^{-3}\)

The derivation is mathematically straightforward once \(\Delta\Phi_0\) and \(\Delta\Phi_1\) are fixed. But “derived, not fitted” is only partially true: the crossing time is not fit to the full \(\Delta\Phi(\tau)\) curve, yet the choice of representatives \(F_A,F_B\) and the interpolation grid do affect the reported measured \(\tau_c\), and the paper itself quantifies this. So the claim is acceptable only if carefully qualified.

**Finding: Major** — The mechanism is derived from the \(t=0\) coefficients, but the paper overstates the independence from numerical choices.

---

## 3) Logic / claim-scope audit

### Finite-size persistence vs universality
The paper repeatedly says “finite-size persistence,” which is appropriate. But the abstract and conclusion still flirt with broader language: “the selected factorization responds systematically,” “the crossover scale is therefore set by the objective balance,” and “the same state admits different optimal relational structures depending on the timescale.” These are fine as finite-protocol statements, but they should not be read as universal claims.

**Finding: Major** — The scope is mostly controlled, but some phrasing still invites overgeneralization beyond \(N\le 5\) and the tested families.

### State dependence
The paper handles state dependence better than many such studies: it explicitly says Haar and thermal show crossovers, ground and mixed do not up to \(\tau=1\). That is good. However, the state ranking in Sec. 4 (“mixed > thermal > ground > haar”) is presented as preserved at every \(\tau\), while the crossover discussion shows that the ordering can change with \(\tau\) for some states. The paper should be careful to distinguish ranking of \(\Phi\) at fixed \(\tau\) from crossover existence.

**Finding: Minor** — Mostly clear, but the ranking language could be misread.

### Universality / thermodynamic scaling
The paper correctly says “not a scaling law and not a thermodynamic limit.” Good. But the discussion of “approximately flat trend across \(N\)” and “size effect rather than shape artifact” is still a bit too suggestive. With only four systems and one \(N=5\) full-quotient point, no scaling inference is justified.

**Finding: Minor** — The paper is mostly careful here.

---

## 4) Falsifiability / physical meaning

### Normalized functional replaced by unnormalized functional
This is conceptually reasonable: the normalized ratio can blow up when the residual norm vanishes, while the unnormalized inner product remains finite. The replacement is not ad hoc in the narrow mathematical sense. However, the paper should be more explicit that this is a modeling choice motivated by well-posedness, not a derivation from first principles. The phrase “natural unnormalized counterpart” is acceptable, but “canonical” may be too strong unless Spec v2.2 explicitly mandates it.

**Finding: Minor** — Justification is reasonable, but the normative language is stronger than the evidence.

### Liouvillian-gap negative result
The paper reports:
\[
\tau_c\Delta_L = 0.0119 \quad (N=3), \qquad \tau_L/\tau_c \approx 83\text{–}136 \quad (N=3\text{–}5).
\]
This does show the crossover is much shorter than the slowest Liouvillian relaxation time. But it does **not** exclude:
- mode-coupling explanations,
- subleading spectral structure,
- state-dependent effective gaps,
- basin-specific dynamical timescales.

The paper itself notes these are not excluded, which is good. Therefore the negative result is appropriately framed as a constraint, not a proof of irrelevance.

**Finding: Minor** — Properly cautious in the discussion, though the abstract still sounds a bit too definitive.

---

## 5) Journal-level completeness

### Abstract
Strong on narrative, but it overcompresses the evidence and uses several loaded phrases:
- “singular on the unrestricted quotient”
- “seed instability disappear”
- “approximately flat over the tested finite-size range”
- “derived, not fitted”
- “the \(O(1)\) hypothesis is rejected”

These are all defensible only with careful qualification. The abstract should be toned down.

**Finding: Major**

### Introduction
Good motivation and clear claim structure. It does a decent job distinguishing the new paper from the frozen spec. However, it would benefit from a more explicit statement of what is new numerically versus what is a restatement of Spec v2.2.

**Finding: Minor**

### Discussion
The discussion is the strongest part conceptually. It correctly separates the negative Liouvillian-gap result from the positive objective-balance mechanism. Still, the discussion occasionally reads as if the mechanism is established more generally than the finite data allow.

**Finding: Minor**

### Figures/tables alone
The tables and captions are mostly sufficient to follow the main numerical story:
- Fig. 1 supports Claim 1,
- Table in Sec. 4 supports Claim 2,
- Table in Sec. 5 supports Claim 3,
- Fig. 2 and Table S2 support Claim 4a,
- Fig. 3 and D4 table support Claim 4b,
- Fig. 4 and Phase E table support Claim 4c.

But the reproducibility trail is incomplete in the manuscript itself: the reader is told the JSONs exist, but the exact scripts, seeds, and representative-selection procedure are not fully specified in the main text. The paper also mixes “main text” and “supplementary” evidence in a way that makes independent verification harder than it should be.

**Finding: Major** — The figures/tables are informative, but not yet sufficient for fully independent reproduction from the paper alone.

---

## MOST IMPORTANT QUESTION

**If I ignore the authors’ interpretation and inspect only the definitions, derivations, numerical evidence, and figures, do Claims 1–4 actually follow from the presented evidence?**

**Answer: Partially, but not fully as stated.**

- **Claim 1:** The evidence strongly supports a **protocol-specific falsification** of the normalized functional in the tested Haar \(N=3\) setting. It does **not** fully establish the universal statement “singular on the unrestricted quotient” as a theorem from the presented material.
- **Claim 2:** The evidence supports that the unnormalized functional is **numerically stable and non-singular in the tested runs**. It does **not** prove global well-posedness on the unrestricted quotient.
- **Claim 3:** Yes, the evidence does support that the selected structure depends on both Hamiltonian and state in the tested families/states, though only as a finite-sample result.
- **Claim 4:** The existence of a fixed-representative objective crossing and a finite-size persistence over \(N=3\text{–}5\) is supported. The “derived, not fitted” mechanism is supported in a qualified sense, but the paper overstates how independent it is from representative choice and interpolation resolution.

So the claims mostly follow **numerically**, but several are stated more strongly than the evidence warrants.

---

## Summary of findings by severity

### Critical
- None that force immediate rejection of the numerical story as a whole.

### Major
1. **Claim 1 overgeneralizes** a protocol-specific numerical falsification into a universal statement about the unrestricted quotient.
2. **Claim 2 overstates global well-posedness**; the evidence shows stability in tested runs, not a proof of absence of singular directions.
3. **Claim 4a/4c overstate mechanism certainty**: the objective crossing is real, but “first-order-like basin crossover” and “derived, not fitted” are stronger than the evidence strictly supports.
4. **Reproducibility is incomplete in the manuscript**: the main text does not fully specify enough to reproduce the key analyses without the external JSONs/scripts.

### Minor
1. Cauchy–Schwarz and boundedness arguments are fine but should be stated more carefully.
2. Second-derivative notation in Sec. 8 is ambiguous.
3. State-response and finite-size language should be kept explicitly finite-sample.
4. Abstract wording is too strong in several places.

---

## Recommendation
**Major revision**