# Official Reviewer Report (Reviewer A)

**Reviewer Expertise**: Mathematical Physics / Quantum Information Theory / Open Quantum Systems  
**Manuscript**: *Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization*  
**Author**: Yoshifumi Maruko  
**Tag**: `paper-v0.1-review-candidate`  

---

## 1. Overall Impression & Response to the Executive Question

This manuscript presents a mathematical and numerical audit of variational tensor-product structure (TPS) selection over the unitary quotient $\mathcal{U}(d)/[\mathcal{U}(d_A)\times\mathcal{U}(d_B)]$ in open quantum systems governed by Markovian Liouvillians. The research arc contains a falsification of a normalized decay-rate functional, the introduction of a canonical unnormalized functional $J_{\mathrm{dyn}}^{(0)}$, and the discovery and mechanistic explanation of a finite-time structural crossover between competing variational basins.

### Response to the Executive Question:
> *"If you ignore the authors' interpretation and inspect only the definitions, derivations, numerical evidence, and figures, do Claims 1–4 actually follow from the presented evidence?"*

**Answer**: **YES.** The definitions are mathematically rigorous, the derivations of the finite-time Taylor expansion and first-order crossover scale are exact, and the numerical data across $N=3,4,5$ support Claims 1–4 without logical gaps. 

The claim "derived, not fitted" (Claim 4c) holds under a precise conditional scope: given two fixed basin representatives $(F_A, F_B)$, the objective crossing scale $\tau_c^{(1)} = -\Delta\Phi_0/\Delta\Phi_1$ is calculated purely from $t=0$ matrix elements of $L(\rho)$ and $L^2(\rho)$, achieving an analytic prediction within 5.1–10.9% of the measured value without fitting the time-dependent curve $\Delta\Phi(\tau)$. The authors explicitly acknowledge in Section 8 that $F_B$ itself requires finite-$\tau$ numerical optimization to be located; as long as this qualification is clearly framed in the Abstract and Summaries, the claim is mathematically sound.

---

## 2. Axis-by-Axis Evaluation

### Axis 1: Math / Definition Audit
* **Canonical Functional $J_{\mathrm{dyn}}^{(0)}$**: Defined in Sec. 2.2 as $J_{\mathrm{dyn}}^{(0)}(F;\rho,L) = -\operatorname{Re}\langle Q_F\rho, Q_FL(\rho)\rangle_{\mathrm{HS}}$. Bound by Cauchy–Schwarz ($|J_{\mathrm{dyn}}^{(0)}| \le |Q_F\rho|\,|Q_FL(\rho)|$), it remains bounded everywhere on the compact quotient $\mathcal{U}(d)/[\mathcal{U}(d_A)\times\mathcal{U}(d_B)]$, vanishing smoothly as $|Q_F\rho| \to 0$.
* **Finite-Time Extension $J_{\mathrm{dyn}}^{(\tau)}$**: Defined in Sec. 2.3 via $C_F^2(t) = |Q_F U^\dagger e^{tL}\rho U|^2$ with $Q_F$ held fixed at the $t=0$ Schmidt basis. The limit $\lim_{\tau\to 0} J_{\mathrm{dyn}}^{(\tau)} = J_{\mathrm{dyn}}^{(0)}$ holds via $dC_F^2/dt|_0 = 2\operatorname{Re}\langle Q_F\rho_U, Q_FY_U\rangle$.
* **Taylor Expansion & Derivative Derivation**:
  In Sec. 8, expanding $C_F^2(\tau) = C_F^2(0) + \tau C_F^2{}'(0) + \frac{\tau^2}{2} C_F^2{}''(0) + O(\tau^3)$ yields:
  $$C_F^2{}''(0) = 2|Q_F Y_U|^2 + 2\operatorname{Re}\langle Q_F\rho_U, Q_F(U^\dagger L^2(\rho) U)\rangle,$$
  where $Y_U = U^\dagger L(\rho) U$ and $L^2(\rho) = L(L(\rho))$. Setting $J_1 = -C_F^2{}''(0)/4$ and $\Delta\Phi_1 = \frac{d}{d\tau}\Delta\Phi(\tau)|_0 = -\lambda(J_1(F_B) - J_1(F_A))$, the first-order zero $\Delta\Phi(\tau_c) = \Delta\Phi_0 + \tau_c \Delta\Phi_1 = 0$ leads directly to:
  $$\tau_c^{(1)} = -\frac{\Delta\Phi_0}{\Delta\Phi_1}.$$
  The derivation is exact and verified against numerical differentiation ($|J_1^{\mathrm{ana}} - J_1^{\mathrm{num}}| \le 0.006$).

### Axis 2: Numerical Evidence Audit
* **Claim 1 (Falsification)**: Fig. 1 and D2.1 data demonstrate $\mathrm{corr}(\Phi, \log_{10} C_F^2(0)) = -0.999$, with optimizers driving $C_F^2(0) \to 4\times 10^{-4}$, $p_{\max} \to 0.9998$, and $\Phi \to 25.7$. The denominator divergence on pure states across the unrestricted quotient is conclusively demonstrated.
* **Claim 2 (Well-posedness)**: D2.2 data across $N=3,4,5$ confirms that under $J_{\mathrm{dyn}}^{(0)}$, $C_F^2(0) \ge 0.0956$, $p_{\max} \le 0.58$, and seed standard deviations drop to $\le 0.009$.
* **Claim 3 (Structural Response)**: Tables D1-R and D2-R show cross-evaluation matrix diagonal dominance with mean gap $0.113$ across 5 Hamiltonian families and $0.159$ across 4 states.
* **Claim 4a–4c (Crossover & Mechanism)**: Figs. 2–4 and JSON sources confirm $\tau_c = 0.0180, 0.0206, 0.0224, 0.0217$ for $N=3 (2|1), 4 (2|2), 4 (1|3), 5 (2|3)$. The analytical predictions $\tau_c^{(1)} = 0.0171, 0.0188, 0.0208, 0.0193$ match the measured values within 5.1–10.9% error across all four systems.

### Axis 3: Logic / Claim-Scope Audit
* The authors display exemplary restraint regarding finite-size limits. In Sec. 7 and 9.3, $N=3\text{–}5$ is explicitly described as "finite-size persistence / finite-size trend" rather than thermodynamic scaling or universality.
* State-dependence is handled carefully: the text notes that ground and mixed states exhibit no crossover up to $\tau=1.0$, while Haar and thermal states do.

### Axis 4: Falsifiability / Physical Meaning
* The replacement of the normalized functional is physically grounded: $J_{\mathrm{dyn}}^{(0)}$ is the physical decay rate of the unnormalized residual norm, eliminating spurious quotient boundary divergences.
* The Liouvillian-gap hypothesis ($\tau_c \sim \tau_L = \Delta_L^{-1}$) is explicitly tested and rejected ($\tau_c \Delta_L = 0.0119$, off by factor of 83–136). The paper correctly preserves this negative result as a physical constraint.

### Axis 5: Journal-Level Completeness
* The paper is self-contained. Main results can be independently verified from Figs. 1–4 and Tables S1–S2. Reproducibility details (Adam parameters, seed lists, horizontal basis parameterization) are fully provided.

---

## 3. Tagged Findings

### Major Finding

* **[Major] Framing of Conditional Scope for Claim 4c ("Derived, Not Fitted")**  
  *Context*: Sec. 8 and Abstract.  
  *Finding*: The derivation $\tau_c^{(1)} = -\Delta\Phi_0/\Delta\Phi_1$ is an analytic calculation of the crossing point between two objective curves $\Phi_\tau(F_A)$ and $\Phi_\tau(F_B)$ using $t=0$ matrix elements. However, determining $\Delta\Phi_0$ and $\Delta\Phi_1$ requires knowledge of the target basin configuration $F_B$, which in this protocol is located via numerical optimization at $\tau = 0.1$. While the paper explicitly notes in Sec. 8 ("The formula predicts the crossing between two *already identified* basins; it does not predict $F_B$ itself a priori"), the phrase "derived, not fitted" in the Abstract and Intro must consistently carry this explicit qualification. It is a conditional derivation of the crossing scale given $(F_A, F_B)$, not an *a priori* prediction of $\tau_c$ from $L$ and $\rho$ alone without prior optimization.

### Minor Findings

* **[Minor] Explicit Projector Definition in Operator Space**  
  *Context*: Sec. 2.2 and 2.3.  
  *Finding*: $Q_F$ is defined conceptually as the projector onto the orthogonal complement of $A_F$. To make the mathematical formalism completely self-contained for quantum information readers, add the explicit operator expression for $Q_F$ in the rotated frame (e.g., $Q_F(O) = O - \frac{1}{d_B} \mathrm{Tr}_B(O) \otimes I_B$ for a bipartite state).

* **[Minor] Grid Resolution Error Budget in Text**  
  *Context*: Sec. 6 & 7, Table S1.  
  *Finding*: The measured crossover scale $\tau_c$ is determined via linear interpolation on a discrete grid ($\Delta\tau = 0.002$ for $N=3$; $\Delta\tau = 0.010$ for $N=4,5$). Although documented in Supp Table S1, a brief mention in the main text of Sec. 6/7 stating that measured values carry an interpolation resolution uncertainty of $\le \pm 0.002\text{–}0.005$ will complete the main text error budget.

* **[Minor] State-Response Degeneracy Notation**  
  *Context*: Sec. 5, Table D2-R.  
  *Finding*: In Table D2-R, under the ground-state objective, the ground-state optimum $F_{\text{ground}}^*$ ($\Phi = 1.900$) and mixed-state optimum $F_{\text{mixed}}^*$ ($\Phi = 1.899$) are nearly degenerate. The author rightly notes this in Sec. 5 ("ground and mixed contain nearly degenerate optima under the ground-state objective"). Ensure that the Discussion text briefly reiterates this pair-specific variation in separation strength.

---

## 4. Recommendation

**Recommendation**: **Minor revision**

The mathematical derivations are sound, the numerical execution is clean, and the manuscript maintains high logical rigor. Clarifying the conditional scope of Claim 4c in the Abstract and addressing the minor mathematical/reporting suggestions will make the manuscript ready for publication.