# Mathematical Physics / Quantum Information Theory Review (Reviewer A)

**Manuscript:** *Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization*  
**Author:** Yoshifumi Maruko  
**Version:** v0.1 draft (`paper-v0.1-review-candidate`)

---

## Executive Summary & Direct Answer to the Most Important Question

### **MOST IMPORTANT QUESTION**
> *"If you ignore the authors' interpretation and inspect only the definitions, derivations, numerical evidence, and figures, do Claims 1–4 actually follow from the presented evidence?"*

**Answer:** **YES, with one necessary technical qualification on Claim 4c.**

1. **Claim 1 (Normalized Functional Singularity):** **Follows directly.** The mathematical argument (transitivity of $\mathcal{U}(d)$ on pure state orbits allowing $|Q_F \rho|^2 \to 0$) combined with the D2.1 numerical evidence ($\mathrm{corr}(\Phi, \log_{10} C_F^2(0)) = -0.999$, objective divergence to $\Phi = 25.71$, product collapse $p_{\max} = 0.9998$, residual minimum $C_F^2(0) = 4.09 \times 10^{-4}$) rigorously demonstrates that $\Gamma_F^{(0)}$ is ill-posed as an unrestricted global variational objective.
2. **Claim 2 (Unnormalized Functional Well-Posedness):** **Follows directly.** Eliminating the denominator yields a continuous objective on the compact quotient $\mathcal{U}(d)/[\mathcal{U}(d_A)\times\mathcal{U}(d_B)]$. Numerical verification (D2.2-A/B) confirms $C_F^2(0) \ge 0.0956$, $p_{\max} \le 0.58$, and multi-seed stability (std $\le 0.009$).
3. **Claim 3 (Structural Response to $H$ and $\rho$):** **Follows directly.** The cross-evaluation matrices (D1-R and D2-R) exhibit clear diagonal dominance (mean gaps $0.113$ across 5 Hamiltonian families and $0.159$ across 4 states).
4. **Claim 4a & 4b (Finite-Time Crossover and Persistence):** **Follows directly.** Objective curves $\Delta\Phi(\tau)$ evaluated at fixed representatives exhibit unambiguous zero-crossings at $\tau_c = 0.0180, 0.0206, 0.0224, 0.0217$ across $N=3,4,5$, with optimizer seed populations hopping discretely between basins.
5. **Claim 4c (Crossover Mechanism / "Derived, Not Fitted"):** **Follows conditionally.** The Taylor expansion of $J_\tau(F)$ around $\tau=0$ and the resulting formula $\tau_c^{(1)} = -\Delta\Phi_0 / \Delta\Phi_1$ are mathematically exact and evaluated purely from $t=0$ matrix elements. The prediction reproduces measured $\tau_c$ within 5.1–10.9%. However, **"derived" applies to predicting $\tau_c$ given the basin pair $(F_A, F_B)$**. Finding $F_B$ itself requires finite-$\tau$ numerical optimization ($\tau=0.1$). As long as $F_B$ is recognized as an input extracted from optimization at finite window length, the first-order prediction of $\tau_c$ is genuinely derived without fitting time-dependent curves.

---

## Axis-by-Axis Evaluation

### Axis 1: Math / Definition Audit

1. **Unnormalized Functional $J_{\mathrm{dyn}}^{(0)}$ and Projection Operator $Q_F$:**
   The projection $Q_F(O) = O - \frac{1}{d_B}\mathrm{Tr}_B(O)\otimes I_B$ defined in the frame rotated by $U$ is the exact Hilbert–Schmidt orthogonal complement projector onto the TPS sub-algebra $A_F = \mathcal{L}(\mathcal{H}_A) \otimes I_B$. Since $Q_F = Q_F^\dagger = Q_F^2$, the Hilbert–Schmidt inner product yields:
   $$J_{\mathrm{dyn}}^{(0)}(F;\rho,L) = -\operatorname{Re}\langle Q_F\rho, Q_FL(\rho)\rangle_{\mathrm{HS}} = -\operatorname{Re}\operatorname{Tr}\big[(Q_F\rho)^\dagger Q_FL(\rho)\big].$$
   By Cauchy–Schwarz, $|J_{\mathrm{dyn}}^{(0)}| \le \|Q_F\rho\|_{\mathrm{HS}} \|Q_FL(\rho)\|_{\mathrm{HS}}$. Since $\dim \mathcal{H} < \infty$, $\|Q_FL(\rho)\|_{\mathrm{HS}}$ is bounded everywhere. Thus $J_{\mathrm{dyn}}^{(0)} \to 0$ smoothly as $\|Q_F\rho\|_{\mathrm{HS}} \to 0$. The objective $\Phi_F = I_F - \lambda J_{\mathrm{dyn}}^{(0)}$ is continuous on the compact quotient $\mathcal{U}(d)/[\mathcal{U}(d_A)\times\mathcal{U}(d_B)]$, guaranteeing the existence of a global maximum by the Extreme Value Theorem.

2. **Finite-Time Derivative Identity and Consistency Limit:**
   The finite-time functional is defined as:
   $$J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{C_F^2(\tau) - C_F^2(0)}{2\tau}, \qquad C_F^2(t) = \|Q_F U^\dagger e^{tL}\rho U\|_{\mathrm{HS}}^2.$$
   Evaluating the derivative of $C_F^2(t)$ at $t=0$:
   $$\left.\frac{d}{dt} C_F^2(t)\right|_{t=0} = 2\operatorname{Re}\left\langle Q_F\rho_U, Q_F \left.\frac{d\rho_U}{dt}\right|_{t=0} \right\rangle_{\mathrm{HS}} = 2\operatorname{Re}\langle Q_F\rho_U, Q_F Y_U\rangle_{\mathrm{HS}},$$
   where $\rho_U = U^\dagger \rho U$ and $Y_U = U^\dagger L(\rho) U$. Taking $\lim_{\tau\to 0} J_{\mathrm{dyn}}^{(\tau)}(F)$:
   $$\lim_{\tau\to 0} J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{1}{2} \left.\frac{d}{dt}C_F^2(t)\right|_{t=0} = -\operatorname{Re}\langle Q_F\rho_U, Q_F Y_U\rangle_{\mathrm{HS}} = J_{\mathrm{dyn}}^{(0)}(F).$$
   The mathematical identity $\lim_{\tau\to 0} J_{\mathrm{dyn}}^{(\tau)} = J_{\mathrm{dyn}}^{(0)}$ is exact.

3. **Taylor Expansion and Derivation of $\tau_c^{(1)}$ and $\tau_c^{(2)}$:**
   Expanding $C_F^2(\tau)$ in a Taylor series around $\tau=0$:
   $$C_F^2(\tau) = C_F^2(0) + \tau C_F^2{}'(0) + \frac{\tau^2}{2} C_F^2{}''(0) + \frac{\tau^3}{6} C_F^2{}'''(0) + O(\tau^4).$$
   Substituting into $J_{\mathrm{dyn}}^{(\tau)}(F)$:
   $$J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{1}{2} C_F^2{}'(0) - \frac{\tau}{4} C_F^2{}''(0) - \frac{\tau^2}{12} C_F^2{}'''(0) + O(\tau^3) = J_0(F) + \tau J_1(F) + \tau^2 J_2(F) + O(\tau^3),$$
   identifying $J_0(F) = -\frac{1}{2} C_F^2{}'(0)$, $J_1(F) = -\frac{1}{4} C_F^2{}''(0)$, and $J_2(F) = -\frac{1}{12} C_F^2{}'''(0)$.
   Differentiating $C_F^2(t) = \langle Q_F \rho_U(t), Q_F \rho_U(t)\rangle$ twice at $t=0$:
   $$C_F^2{}''(0) = 2 \|Q_F Y_U\|_{\mathrm{HS}}^2 + 2\operatorname{Re}\langle Q_F\rho_U, Q_F(U^\dagger L^2(\rho) U)\rangle_{\mathrm{HS}} \implies J_1(F) = -\frac{1}{2} \|Q_F Y_U\|_{\mathrm{HS}}^2 - \frac{1}{2}\operatorname{Re}\langle Q_F\rho_U, Q_F(U^\dagger L^2(\rho) U)\rangle_{\mathrm{HS}}.$$
   For fixed representatives $(F_A, F_B)$, $\Delta\Phi(\tau) = \Phi_\tau(F_B) - \Phi_\tau(F_A) = \Delta\Phi_0 + \tau \Delta\Phi_1 + \tau^2 \Delta\Phi_2 + O(\tau^3)$, where $\Delta\Phi_0 = \Phi_0(F_B) - \Phi_0(F_A)$ and $\Delta\Phi_1 = -\lambda [J_1(F_B) - J_1(F_A)]$. Setting $\Delta\Phi(\tau_c)=0$:
   - First-order truncation: $\Delta\Phi_0 + \tau_c^{(1)} \Delta\Phi_1 = 0 \implies \tau_c^{(1)} = -\frac{\Delta\Phi_0}{\Delta\Phi_1}$.
   - Second-order (quadratic): $\Delta\Phi_0 + \tau_c^{(2)} \Delta\Phi_1 + (\tau_c^{(2)})^2 \Delta\Phi_2 = 0 \implies \tau_c^{(2)} = \frac{-\Delta\Phi_1 + \sqrt{\Delta\Phi_1^2 - 4\Delta\Phi_0 \Delta\Phi_2}}{2\Delta\Phi_2}$.

   The derivation of $\tau_c^{(1)}$ from $t=0$ matrix elements is mathematically sound and rigorous.

---

### Axis 2: Numerical Evidence Audit

1. **Singularity and Product Collapse (Claim 1, D2.1):**
   The numerical digest confirms $\mathrm{corr}(\Phi, \log_{10} C_F^2(0)) = -0.998955$, $\Phi_{\max} = 25.7133$, $\Gamma_{\min} = -128.5251$, and $C_F^2(0)_{\min} = 4.093 \times 10^{-4}$. Figure 1 clearly captures this systematic divergence along product state directions ($p_{\max} = 0.999795$).

2. **Well-Posed Replacement (Claim 2, D2.2-A/B):**
   Under $J_{\mathrm{dyn}}^{(0)}$, across 30 Haar seeds at $N=3$, $\Phi^* = 1.8594 \pm 0.0077$, $C_F^2(0)_{\min} = 0.0956$, and $p_{\max} \in [0.51, 0.58]$. Singularity and product collapse are completely absent.

3. **Structural Response (Claim 3, D1-R/D2-R):**
   - **Hamiltonian response:** Cross-evaluation matrix $\Phi_{H_i}(F_j^*)$ displays clear diagonal dominance (mean gap $0.113$).
   - **State response:** Mean gap $0.159$. The paper correctly notes the specific near-degeneracy between ground and mixed states under the ground state objective ($\Phi_{\mathrm{ground}}(F^*_{\mathrm{ground}}) = 1.900$ vs $\Phi_{\mathrm{ground}}(F^*_{\mathrm{mixed}}) = 1.899$).

4. **Crossover Scale and Trend (Claims 4a & 4b, D3/D4):**
   - $N=3$ ($2|1$): Measured $\tau_c = 0.0180$, $d_F(F_A, F_B) = 0.8858$.
   - $N=4$ ($2|2$): Measured $\tau_c = 0.0206$, $d_F(F_A, F_B) = 0.968$.
   - $N=4$ ($1|3$): Measured $\tau_c = 0.0224$, $d_F(F_A, F_B) = 0.994$.
   - $N=5$ ($2|3$, full 945-dim quotient): Measured $\tau_c = 0.0217$, $d_F(F_A, F_B) = 0.3191$.

5. **Analytical Prediction vs. Measurement (Claim 4c, Phase E):**
   Comparing first-order derived $\tau_c^{(1)} = -\Delta\Phi_0 / \Delta\Phi_1$ against measured $\tau_c$:
   - $N=3, 2|1$: $\tau_c^{(1)} = 0.0171$ vs $0.0180$ (Relative Error = $5.1\%$)
   - $N=4, 2|2$: $\tau_c^{(1)} = 0.0188$ vs $0.0206$ (Relative Error = $8.6\%$)
   - $N=4, 1|3$: $\tau_c^{(1)} = 0.0208$ vs $0.0224$ (Relative Error = $7.2\%$)
   - $N=5, 2|3$: $\tau_c^{(1)} = 0.0193$ vs $0.0217$ (Relative Error = $10.9\%$)

   The 5.1–10.9% agreement validates the first-order expansion mechanism.

---

### Axis 3: Logic / Claim-Scope Audit

The manuscript maintains commendable logical discipline regarding scope:
- **Finite-Size Scope:** The authors repeatedly and explicitly state that $N=3\text{–}5$ establishes a "finite-size trend / finite-size persistence," explicitly rejecting claims of thermodynamic scaling laws or universality.
- **Phase Transitions vs. Optimization Crossovers:** The transition is explicitly defined as a "sharp, first-order-like basin crossover in a finite-dimensional optimization problem," avoiding misuse of thermodynamic phase transition terminology.
- **State Dependence:** The manuscript accurately notes that structural crossovers are state-dependent (present in Haar and thermal states, absent in ground and mixed states up to $\tau=1.0$).

---

### Axis 4: Falsifiability / Physical Meaning

1. **Falsification of $\Gamma_F^{(0)}$:** The rejection of $\Gamma_F^{(0)}$ is backed by direct numerical proof of singularity. Replacing it with $J_{\mathrm{dyn}}^{(0)}$ is physically natural: $J_{\mathrm{dyn}}^{(0)}$ is the unnormalized rate of decay of the residual norm. The singularity of $\Gamma_F^{(0)}$ was an artifact of dividing by the vanishing norm $\|Q_F \rho\|^2$.
2. **Liouvillian Gap Rejection:** The paper highlights the failure of the global Liouvillian gap $\Delta_L$ to set the crossover scale ($\tau_c \Delta_L \approx 0.012$, off by a factor of 83–136). Retaining this negative result prevents naive physical misinterpretations.

---

### Axis 5: Journal-Level Completeness

The draft is self-contained. The figures, tables, abstract, and discussion cleanly convey the complete research arc. The distinction between theory specification (`REM Spec v2.2`) and empirical validation results (`Phase D/E Addendum`) is maintained throughout.

---

## Detailed Tagged Findings

### **Major Findings**

- **[Major Finding 1] Scope and Qualification of "Derived, Not Fitted" in Claim 4c (Axis 1 & Axis 2)**  
  *Context:* Abstract, Section 1, and Section 8.  
  *Detail:* The formula $\tau_c^{(1)} = -\Delta\Phi_0 / \Delta\Phi_1$ is called "derived to first order, not fitted." While calculating $\tau_c^{(1)}$ from $t=0$ matrix elements given $(F_A, F_B)$ is purely analytical, obtaining $F_B$ requires running non-convex optimization at a finite window length ($\tau = 0.1$). Although Section 8 transparently notes this distinction ("The formula predicts the crossing between two *already identified* basins; it does not predict $F_B$ itself a priori"), the Abstract and Introduction present the phrase "derived to first order" without immediate qualification.  
  *Required Action:* Clarify in the Abstract and Section 1 that $\tau_c^{(1)}$ is an analytical prediction of the crossover scale between two variational basins *given* the basin pair $(F_A, F_B)$, where $F_B$ is located via finite-$\tau$ numerical optimization.

- **[Major Finding 2] Optimization Landscape Degeneracy at $N=5$ (Axis 2 & Axis 3)**  
  *Context:* Section 7 (Table for D4-B) and Supplementary S1.  
  *Detail:* For $N=5$ ($2|3$), full 945-dimensional quotient optimization yields a basin distance $d_F(F_A, F_B) = 0.3191$, whereas the 200-dimensional subspace search gave $d_F = 0.9920$ (and $N=3,4$ gave $d_F \in [0.886, 0.994]$). Both yield similar crossover times ($\tau_c = 0.0217$ vs $0.0173$), but the drastically reduced $d_F$ on the full quotient indicates that in 945 dimensions, 2 seeds × 100 steps either located a nearby local maximum or that the objective landscape at $\tau=0.1$ contains a dense manifold of value-degenerate local maxima.  
  *Required Action:* Expand Section 7 / Section 9.3 to explicitly address why $d_F$ drops to $0.3191$ on the 945-dimensional quotient, and state whether this reflects local basin density or optimization budget limitations at $N=5$.

---

### **Minor Findings**

- **[Minor Finding 1] Notation for Superoperator Iteration in Section 8 (Axis 1)**  
  *Context:* Section 8, equation for $C^2{}''(0)$.  
  *Detail:* $C^2{}''(0)$ is written using $U^\dagger L(L(\rho)) U$. Later text references $L^2(\rho) := L(L(\rho))$.  
  *Required Action:* Explicitly state upon first use that $L^2(\rho) \equiv L(L(\rho))$ denotes the double application of the Liouvillian superoperator $L$.

- **[Minor Finding 2] Tabular Representation of Ground vs. Mixed State Degeneracy in Section 5 (Axis 2)**  
  *Context:* Section 5, D2-R cross-evaluation matrix table.  
  *Detail:* The text notes that evaluating the ground-state objective on $F^*_{\mathrm{mixed}}$ yields $\Phi = 1.899$ vs $\Phi = 1.900$ for $F^*_{\mathrm{ground}}$. The reported gap in the table is $0.000$ due to rounding.  
  *Required Action:* Add a footnote or table annotation in Section 5 explicitly highlighting that under the ground-state objective, $F^*_{\mathrm{ground}}$ and $F^*_{\mathrm{mixed}}$ are objective-degenerate to within $10^{-3}$, despite being separated by $d_F = 0.343$.

- **[Minor Finding 3] Visual Grid Markers in Figure 2 (Axis 5)**  
  *Context:* Figure 2 caption and main plot.  
  *Detail:* The text references a 19-point grid for $\Delta\Phi(\tau)$ (detailed in Supplementary Table S2), but Figure 2 renders $\Delta\Phi(\tau)$ as a solid curve without explicit data markers at the evaluated grid points.  
  *Required Action:* Include discrete point markers on the curve in Figure 2 to visually reflect the sampled grid points from Table S2.

---

## Formal Recommendation

**Recommendation:** **Minor Revision**

The paper is mathematically sound, physically well-motivated, and remarkably transparent regarding its falsification history, negative results, and finite-size limitations. Addressing the minor clarifications regarding the scope of "derived" in Claim 4c (Major Finding 1) and the $N=5$ optimization landscape (Major Finding 2) will make the manuscript publication-ready.