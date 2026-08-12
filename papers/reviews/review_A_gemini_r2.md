# Formal Academic Review (Reviewer A)

**Reviewer Identity:** Mathematical Physicist / Quantum Information Theorist  
**Manuscript Title:** Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization  
**Author:** Yoshifumi Maruko  
**Document Version:** Frozen v0.1 (`paper-v0.1-review-candidate`) / REM Spec v2.2  

---

## Executive Summary & Rating Matrix

This manuscript investigates the mathematical well-posedness and finite-time observation dynamics of variational Tensor-Product Structure (TPS) selection over the unitary quotient $\mathcal{U}(d)/[\mathcal{U}(d_A) \times \mathcal{U}(d_B)]$ in the Relational Emergence Model (REM). The research arc presents a clear falsification-and-replacement sequence: demonstrating that a normalized local decay-rate functional $\Gamma_F^{(0)}$ is singular on the unrestricted quotient, replacing it with an unnormalized canonical functional $J_{\mathrm{dyn}}^{(0)}$, verifying structural response, and discovering/explaining a finite-time structural crossover scale $\tau_c \approx 0.02$.

### Rating on the 5 Review Axes

| Axis | Assessment | Rating |
| :--- | :--- | :--- |
| **1. Math / Definition Audit** | Derivations of $J_{\mathrm{dyn}}^{(0)}$, $J_{\mathrm{dyn}}^{(\tau)}$, $C_F^2{}''(0)$, and $\tau_c^{(1)}$ are mathematically exact. Hilbert-Schmidt inner products and frame transformations are rigorous. | **Excellent** |
| **2. Numerical Evidence Audit** | Claims 1, 2, 3, 4a, and 4b follow directly from the JSON outputs and tables. Claim 4c holds for the 1st-order derivation (5–11% error), but conflates a fitted 2nd-order polynomial term with a derived result. | **Minor Defect in Claim 4c** |
| **3. Logic / Claim-Scope Audit** | Exemplary discipline regarding finite-size scope ($N=3\text{--}5$). Disclaims thermodynamic scaling, universality, and phase transitions. State dependence is properly contextualized. | **Excellent** |
| **4. Falsifiability / Physical Meaning** | Falsification of normalized functional is mathematically proven. Rejection of the Liouvillian-gap hypothesis ($\tau_L/\tau_c \approx 83\text{--}136$) is cleanly handled as a physical negative constraint. | **Excellent** |
| **5. Journal-Level Completeness** | Figures, captions, and tables are self-contained. Methods, limitations, and differentiation from prior frozen specifications are clearly articulated. | **High Quality** |

---

## Explicit Answer to the Most Important Question

> **"If you ignore the authors' interpretation and inspect only the definitions, derivations, numerical evidence, and figures, do Claims 1–4 actually follow from the presented evidence?"**

**YES, WITH ONE QUALIFICATION REGARDING CLAIM 4c.**

1. **Claim 1 (Singularity of normalized functional):** **FOLLOWS.** Mathematically proven via the transitive action of $\mathcal{U}(d)$ on pure states (which forces $|Q_F \rho|^2 \to 0$), and numerically verified by $\operatorname{corr}(\Phi, \log_{10} C_F^2(0)) = -0.998955$, $p_{\max} \to 0.9998$, and objective blow-up $\Phi \to 25.71$.
2. **Claim 2 (Well-posedness of unnormalized functional):** **FOLLOWS.** $J_{\mathrm{dyn}}^{(0)}$ is bounded by Cauchy–Schwarz ($|J_{\mathrm{dyn}}^{(0)}| \le |Q_F\rho| |Q_FL(\rho)|$). Numerically, product collapse disappears ($p_{\max} \le 0.58$), minimum residuals remain strictly positive ($\min C_F^2 \ge 0.0956$), and seed variance remains small ($\le 0.009$).
3. **Claim 3 (Hamiltonian and state structural response):** **FOLLOWS.** Cross-evaluation matrices $\Phi_{H_i}(F^*_j)$ and $\Phi_{\rho_i}(F^*_j)$ display strict diagonal dominance across tested Hamiltonians (mean gap 0.113) and states (mean gap 0.159).
4. **Claim 4a & 4b (Finite-time crossover & finite-size persistence):** **FOLLOWS.** At fixed structural representatives $(F_A, F_B)$, $\Delta\Phi(\tau)$ exhibits a smooth, monotone zero-crossing at $\tau_c = 0.0180, 0.0206, 0.0224, 0.0217$ for $N=3, 4(2|2), 4(1|3), 5(2|3)$, demonstrating finite-size persistence across $N=3\text{--}5$.
5. **Claim 4c (Derived crossover mechanism):** **QUALIFIED FOLLOWS.** 
   - The 1st-order prediction $\tau_c^{(1)} = -\Delta\Phi_0 / \Delta\Phi_1^{\mathrm{ana}}$ is derived purely from $t=0$ operator matrix elements ($L^2$ action) at fixed representatives $(F_A, F_B)$, predicting $\tau_c$ within 5.1–10.9% relative error across all four systems. This portion is genuinely derived and non-fitted.
   - *Qualification:* The reported 0.0–0.5% agreement for the quadratic prediction $\tau_c^{(2)}$ relies on a coefficient $\Delta\Phi_2$ obtained via a least-squares polynomial fit to time-series data ($\Delta\Phi(\tau)$), as disclosed in Supplementary Table S1. This 2nd-order correction is therefore **fitted, not derived**. Furthermore, calculating $\tau_c^{(1)}$ requires prior identification of $F_B$ via optimization at $\tau=0.1$.

---

## Detailed Audit Findings by Axis

### Axis 1: Mathematical & Definition Audit

The mathematical framework is sound, precise, and internally consistent.

1. **Unnormalized Functional ($J_{\mathrm{dyn}}^{(0)}$):**
   $$J_{\mathrm{dyn}}^{(0)}(F;\rho,L) = -\operatorname{Re}\langle Q_F\rho, Q_F L(\rho)\rangle_{\mathrm{HS}}$$
   The dephasing projection superoperator $Q_F$ projects onto the orthogonal complement of the local TPS algebra $A_F$. In the Hilbert-Schmidt inner product $\langle A, B\rangle_{\mathrm{HS}} = \operatorname{Tr}(A^\dagger B)$, $Q_F$ is self-adjoint ($Q_F^\dagger = Q_F = Q_F^2$). Cauchy-Schwarz gives $|J_{\mathrm{dyn}}^{(0)}| \le |Q_F\rho| |Q_FL(\rho)|$. Since $\dim \mathcal{H} < \infty$, $L$ is a bounded operator, ensuring $J_{\mathrm{dyn}}^{(0)}$ vanishes smoothly as $|Q_F\rho| \to 0$.

2. **Finite-Time Extension ($J_{\mathrm{dyn}}^{(\tau)}$):**
   $$J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{C_F^2(\tau) - C_F^2(0)}{2\tau}, \qquad C_F^2(t) = |Q_F U^\dagger e^{tL}(\rho) U|^2$$
   Taking the limit $\tau \to 0$:
   $$\lim_{\tau\to 0} J_{\mathrm{dyn}}^{(\tau)} = -\frac{1}{2} \frac{d}{dt} C_F^2(t)\Big|_{t=0} = -\frac{1}{2} \cdot 2 \operatorname{Re}\langle Q_F \rho_U, Q_F Y_U\rangle_{\mathrm{HS}} = J_{\mathrm{dyn}}^{(0)}$$
   where $\rho_U = U^\dagger \rho U$ and $Y_U = U^\dagger L(\rho) U$. The factor of 2 and sign conventions match identically.

3. **Second Derivative and Analytic Expansion Coefficient ($\Delta\Phi_1$):**
   Differentiating $C_F^2(t) = \langle Q_F \rho_U(t), Q_F \rho_U(t)\rangle_{\mathrm{HS}}$ twice at $t=0$:
   $$\dot{\rho}_U(0) = Y_U = U^\dagger L(\rho) U, \qquad \ddot{\rho}_U(0) = U^\dagger L^2(\rho) U$$
   $$\frac{d^2 C_F^2}{dt^2}\Big|_{t=0} = 2 |Q_F Y_U|^2 + 2 \operatorname{Re}\langle Q_F \rho_U, Q_F (U^\dagger L^2(\rho) U)\rangle$$
   The Taylor expansion yields:
   $$J_{\mathrm{dyn}}^{(\tau)} = J_0 + \tau J_1 + O(\tau^2), \qquad \text{where } J_0 = -\frac{1}{2} C_F^2{}'(0), \quad J_1 = -\frac{1}{4} C_F^2{}''(0)$$
   With $\Phi_\tau = I_F - \lambda J_{\mathrm{dyn}}^{(\tau)}$, we have $\Delta\Phi(\tau) = \Delta\Phi_0 + \tau \Delta\Phi_1 + O(\tau^2)$, where $\Delta\Phi_0 = \Phi_0(F_B) - \Phi_0(F_A)$ and $\Delta\Phi_1 = -\lambda [J_1(F_B) - J_1(F_A)]$. Setting $\Delta\Phi(\tau_c) = 0$ yields:
   $$\tau_c^{(1)} = -\frac{\Delta\Phi_0}{\Delta\Phi_1}$$
   This derivation is mathematically exact and non-circular.

---

### Axis 2: Numerical Evidence Audit

The numerical results reported in the paper cross-check against the extracted JSON digest without discrepancy.

* **Claim 1 (Singularity):** $\operatorname{corr}(\Phi, \log_{10} C_F^2(0)) = -0.998955$, maximum $\Phi = 25.7133$, minimum $C_F^2(0) = 4.093\times 10^{-4}$, and maximum Schmidt component $p_{\max} = 0.999795$.
* **Claim 2 (Well-posedness):** $J_{\mathrm{dyn}}^{(0)}$ achieves minimum $C_F^2(0) = 0.0956$, $p_{\max} \in [0.51, 0.58]$, state seed standard deviations $\le 0.009$.
* **Claim 3 (Response):** Hamiltonian cross-evaluation mean gap = 0.113; State cross-evaluation mean gap = 0.159. Strict diagonal dominance is observed in all cross-evaluation matrices.
* **Claim 4a & 4b (Crossover & Persistence):**
  * $N=3, 2|1$: measured $\tau_c = 0.018024$, $d_F(F_A, F_B) = 0.8858$.
  * $N=4, 2|2$: measured $\tau_c = 0.0206$, $d_F = 0.968$.
  * $N=4, 1|3$: measured $\tau_c = 0.0224$, $d_F = 0.994$.
  * $N=5, 2|3$ (full quotient): measured $\tau_c = 0.021687$, $d_F = 0.3191$.

* **Claim 4c Audit (The "Derived, Not Fitted" Mechanics):**
  Using the analytic first-order formula $\tau_c^{(1)} = -\Delta\Phi_0 / \Delta\Phi_1^{\mathrm{ana}}$:
  * $N=3, 2|1$: $\Delta\Phi_0 = -0.0260$, $\Delta\Phi_1^{\mathrm{ana}} = +1.5259 \implies \tau_c^{(1)} = 0.01704$ (measured 0.0180, error **5.1%**).
  * $N=4, 2|2$: $\Delta\Phi_0 = -0.0099$, $\Delta\Phi_1^{\mathrm{ana}} = +0.5281 \implies \tau_c^{(1)} = 0.01875$ (measured 0.0206, error **8.6%**).
  * $N=4, 1|3$: $\Delta\Phi_0 = -0.0339$, $\Delta\Phi_1^{\mathrm{ana}} = +1.6389 \implies \tau_c^{(1)} = 0.02068$ (measured 0.0224, error **7.2%**).
  * $N=5, 2|3$: $\Delta\Phi_0 = -0.0219$, $\Delta\Phi_1^{\mathrm{ana}} = +1.1364 \implies \tau_c^{(1)} = 0.01927$ (measured 0.0217, error **10.9%**).

  The 1st-order analytical prediction matches the measured values within 5.1–10.9%. This confirms that the crossover scale is dictated by static $t=0$ operator matrix elements.

---

### Axis 3: Logic / Claim-Scope Audit

* **Scope discipline:** The paper restricts its claims for $N=3\text{--}5$ to "finite-size persistence" and "flat finite-size trend." It explicitly avoids claiming thermodynamic scaling laws, universality classes, or phase transitions.
* **State dependence:** The manuscript acknowledges that the finite-time crossover is present for Haar and thermal states, but absent up to $\tau=1.0$ for ground and mixed states. This state-selective occurrence is reported without overgeneralization.

---

### Axis 4: Falsifiability & Physical Meaning

* **Falsification of $\Gamma_F^{(0)}$:** Well-justified mathematically. For pure states, the transitive action of $\mathcal{U}(d)$ guarantees that $C_F^2(0) = 0$ is accessible on the quotient, making $\Gamma_F^{(0)}$ ill-conditioned as an unrestricted variational objective.
* **Rejection of the Liouvillian Gap Scale:** The relaxation timescale $\tau_L = \Delta_L^{-1}$ yields $\tau_L / \tau_c \approx 83\text{--}136$ and $\tau_c \Delta_L = 0.012$. The manuscript treats this negative result properly, ruling out global Liouvillian relaxation as the driver for structural crossovers.

---

### Axis 5: Journal-Level Completeness

The paper is structurally complete. The figures and captions clearly convey the core physics:
* **Fig 1:** Inverse correlation between objective $\Phi$ and log-residual norm $\log_{10} C_F^2(0)$ under $\Gamma_F$.
* **Fig 2:** Zero-crossing of $\Delta\Phi(\tau)$ at $\tau_c = 0.0180$ alongside optimizer seed flipping.
* **Fig 3:** Stability of $\tau_c(N) \approx 0.021$ across $N=3,4,5$.
* **Fig 4:** Parity plot of predicted vs. measured $\tau_c$.

---

## Detailed Findings (Tagged Major / Minor)

### **[Major Finding 1] Conflation of Derived 1st-Order Prediction with Fitted 2nd-Order Polynomial Correction**
* **Context:** Abstract, Section 8, and Conclusion state that the quadratic correction reproduces $\tau_c$ within 0.0–0.5% relative error as a "derived, not fitted" result.
* **Defect:** As noted in Supplementary Table S1 ("Fitted coefficients $\Delta\Phi_0, \Delta\Phi_1, \Delta\Phi_2$"), the quadratic coefficient $\Delta\Phi_2$ was obtained by fitting a 2nd-order polynomial to the time-dependent curve $\Delta\Phi(\tau)$, rather than evaluating an analytical 3rd-derivative matrix element ($C_F^2{}'''(0)$).
* **Impact:** Finding the root of a quadratic curve fitted to data points near a zero-crossing guarantees 0.0–0.5% agreement by construction. Claiming that the 0.0–0.5% quadratic result is "derived, not fitted" is inaccurate.
* **Required Revision:** The text must explicitly separate the **genuinely derived 1st-order prediction** ($\tau_c^{(1)}$, 5.1–10.9% error, computed analytically from $L^2$ matrix elements at $t=0$) from the **2nd-order polynomial fit** ($\tau_c^{(2)}$, 0.0–0.5% error, obtained via numerical fitting of $\Delta\Phi(\tau)$).

### **[Major Finding 2] Operational Dependency on Representative $F_B$**
* **Context:** Claim 4c asserts that the crossover scale $\tau_c$ is derived directly from $t=0$ matrix elements.
* **Defect:** Evaluating $\Delta\Phi_0$ and $\Delta\Phi_1^{\mathrm{ana}}$ requires prior knowledge of the alternative structural representative $F_B$. In the reported protocol, $F_B$ is located by running variational optimization at a large observation window ($\tau = 0.1$).
* **Impact:** The formula $\tau_c^{(1)} = -\Delta\Phi_0 / \Delta\Phi_1^{\mathrm{ana}}$ predicts the crossing point between two *already known* structural basins without fitting a time-dependent curve. However, it cannot predict $\tau_c$ *a priori* purely from $\rho$ and $L$ without first searching for $F_B$ via finite-$\tau$ optimization.
* **Required Revision:** Clarify in Section 8 and the Discussion that "derived, not fitted" refers to calculating the inter-basin crossing timescale between two known variational extrema $(F_A, F_B)$ using static operator matrix elements, rather than predicting $F_B$ *a priori*.

---

### **[Minor Finding 1] Interpretation of Quotient Distance for $N=5$ Full Quotient**
* **Context:** Section 7 and Table D4-B report that $d_F(F_A, F_B) = 0.3191$ for $N=5$ ($2|3$ cut) under the full 945-dimensional quotient optimization, whereas $N=3$ and $N=4$ yield $d_F \in [0.886, 0.994]$.
* **Impact:** While $d_F = 0.3191$ confirms non-identity, $F_A$ and $F_B$ are substantially closer in quotient space than those found at smaller system sizes or in the 200-dim subspace ($d_F = 0.992$).
* **Required Revision:** Add a brief remark in Section 7 discussing whether $d_F = 0.3191$ reflects a higher density of local extrema or a shallower basin landscape on the 945-dimensional manifold.

### **[Minor Finding 2] Notation Correction for Second-Order Operator Term**
* **Context:** In Section 8, the second-derivative formula reads $C_F^2{}''(0) = 2|Q_F Y_U|^2 + 2\operatorname{Re}\langle Q_F\rho_U, Q_F(U^\dagger L^2\rho_0 U)\rangle$.
* **Impact:** The notation $\rho_0$ in $L^2\rho_0$ introduces potential confusion with state notation at $t=0$.
* **Required Revision:** Change $L^2\rho_0$ to $L^2(\rho)$ or $L(L(\rho))$ for syntactic consistency with $L(\rho)$ in Section 2.

---

## Final Recommendation

**Recommendation: Major Revision**

*Rationale:* The core mathematical physics of the manuscript is sound, the falsification of $\Gamma_F^{(0)}$ is rigorous, and the 1st-order analytical derivation of the structural crossover scale $\tau_c^{(1)}$ (5.1–10.9% prediction error) is a genuine theoretical contribution. However, claiming 0.0–0.5% agreement for a quadratic term whose coefficient $\Delta\Phi_2$ was obtained via polynomial fitting conflates a derived analytical prediction with a numerical fit. Addressing Major Findings 1 and 2 requires revising the framing of Claim 4c and clarifying the operational role of $F_B$. Once these corrections are made, the manuscript will be suitable for publication.