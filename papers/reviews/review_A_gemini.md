# Reviewer A: Mathematical Physics / Quantum Information Theory Audit

**Manuscript Title:** *Well-Posed Structural Selection and Finite-Time Crossovers in Variational Tensor-Product Factorization*  
**Author:** Yoshifumi Maruko  
**Version / Tag:** paper-v0.1-review-candidate (REM Spec v2.2)

---

## Overall Assessment & Summary

This paper presents a formal and numerical study of variational tensor-product structure (TPS) selection in open quantum systems within the Relational Emergence Model (REM). The research follows a clear two-stage arc: 
1. **Falsification & Well-Posedness:** Demonstrating that the normalized decay-rate functional $\Gamma_F^{(0)}$ possesses a singular denominator on the unrestricted unitary quotient $\mathcal{U}(d)/[\mathcal{U}(d_A)\times\mathcal{U}(d_B)]$, causing product collapse and objective divergence ($\Phi \to +\infty$). Replacing it with the unnormalized Hilbert–Schmidt functional $J_{\mathrm{dyn}}^{(0)} = -\operatorname{Re}\langle Q_F\rho, Q_F L(\rho)\rangle_{\mathrm{HS}}$, which is shown to be well-posed and bounded.
2. **Finite-Time Extension & Crossover:** Extending the functional to an observation window $\tau$, revealing structural crossovers between distinct basins ($d_F \approx 0.89$), verifying finite-size persistence ($N=3,4,5$), and deriving the crossover timescale analytically via small-$\tau$ expansion ($\tau_c^{(1)} = -\Delta\Phi_0/\Delta\Phi_1$).

The paper is mathematically sound, logically tight, and exceptionally honest about its numerical constraints and limitations. Below is an adversarial audit across the five required axes.

---

## Audit Across the 5 Review Axes

### 1. Math / Definition Audit
* **Canonical functional $J_{\mathrm{dyn}}^{(0)}$ and Cauchy–Schwarz bound:**  
  The Hilbert–Schmidt inner product $\langle A, B\rangle_{\mathrm{HS}} = \operatorname{Tr}(A^\dagger B)$ on operators yields $|J_{\mathrm{dyn}}^{(0)}| \le \|Q_F \rho\|_{\mathrm{HS}} \|Q_F L(\rho)\|_{\mathrm{HS}}$. As $\|Q_F \rho\| \to 0$, $J_{\mathrm{dyn}}^{(0)} \to 0$ smoothly, eliminating the singularity of $\Gamma_F^{(0)} \propto 1/\|Q_F\rho\|^2$.
* **Finite-time extension $J_{\mathrm{dyn}}^{(\tau)}$ consistency:**  
  $J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{C_F^2(\tau) - C_F^2(0)}{2\tau}$ with $C_F^2(t) = \|Q_F U^\dagger e^{tL}\rho U\|_{\mathrm{HS}}^2$.  
  Taking the $\tau \to 0$ limit:
  $$\lim_{\tau \to 0} J_{\mathrm{dyn}}^{(\tau)}(F) = -\frac{1}{2} \left.\frac{d C_F^2(t)}{dt}\right|_{t=0} = -\operatorname{Re}\langle Q_F U^\dagger \rho U, Q_F U^\dagger L(\rho) U\rangle_{\mathrm{HS}} = J_{\mathrm{dyn}}^{(0)}(F).$$
* **Taylor Expansion & Derivative Derivation:**  
  Let $A(t) = Q_F U^\dagger e^{tL}\rho U$. Expanding around $t=0$:
  $$A(t) = Q_F \rho_U + t Q_F Y_U + \frac{t^2}{2} Q_F Z_U + O(t^3),$$
  where $\rho_U = U^\dagger \rho U$, $Y_U = U^\dagger L(\rho) U$, and $Z_U = U^\dagger L^2(\rho) U$.  
  Evaluating $C_F^2(t) = \|A(t)\|_{\mathrm{HS}}^2$:
  $$\left.\frac{d C_F^2}{dt}\right|_0 = 2 \operatorname{Re}\langle Q_F \rho_U, Q_F Y_U \rangle = -2 J_0(F),$$
  $$\left.\frac{d^2 C_F^2}{dt^2}\right|_0 = 2 \|Q_F Y_U\|_{\mathrm{HS}}^2 + 2 \operatorname{Re}\langle Q_F \rho_U, Q_F Z_U \rangle = -4 J_1(F).$$
  Hence, $J_1(F) = -\frac{1}{4} C_F^2{}''(0)$, which matches the paper's equation in Sec. 8 exactly.
* **Crossover Formula:**  
  For $\Delta\Phi(\tau) = \Phi_\tau(F_B) - \Phi_\tau(F_A) = \Delta\Phi_0 + \tau \Delta\Phi_1 + \tau^2 \Delta\Phi_2 + O(\tau^3)$, setting $\Delta\Phi(\tau_c) = 0$ yields:
  $$\tau_c^{(1)} = -\frac{\Delta\Phi_0}{\Delta\Phi_1}, \qquad \tau_c^{(2)} = \frac{-\Delta\Phi_1 + \sqrt{\Delta\Phi_1^2 - 4\Delta\Phi_0\Delta\Phi_2}}{2\Delta\Phi_2}.$$
  This derivation is exact and free of algebraic errors.

### 2. Numerical Evidence Audit
* **D2.1 (Claim 1):** $\operatorname{corr}(\Phi, \log_{10} C_F^2(0)) = -0.998955$, $\Phi_{\max} = 25.71$, $C_{F0}^2 \to 4.09 \times 10^{-4}$, $p_{\max} = 0.999795$. Singular divergence to near-product states under $\Gamma$ is verified.
* **D2.2 (Claim 2):** Under $J_{\mathrm{dyn}}^{(0)}$, $\Phi^*$ values stabilize ($\approx 1.86 - 1.98$ for $N=3$), minimum $C_F^2 \ge 0.0956$, $p_{\max} \le 0.58$, seed std $\le 0.009$.
* **D1-R / D2-R (Claim 3):** Cross-evaluation matrix shows clear diagonal dominance (mean gap $0.113$ for Hamiltonians, $0.159$ for states).
* **D3 / D4 (Claims 4a, 4b):** Monotone objective cross $\Delta\Phi(\tau_c) = 0$ at $\tau_c = 0.0180$ ($N=3$), $0.0206$ ($N=4, 2|2$), $0.0224$ ($N=4, 1|3$), $0.0217$ ($N=5, 2|3$).
* **Phase E (Claim 4c):** First-order prediction $\tau_c^{(1)}$ reproduces measured values within 5.1–10.9%; second-order correction $\tau_c^{(2)}$ reduces error to 0.0–0.5%.

### 3. Logic / Claim-Scope Audit
* **Scope Control:** The authors carefully refrain from asserting thermodynamic scaling or universality. The findings are properly designated as "finite-size persistence / finite-size trend" ($N=3\text{–}5$).
* **State Dependence:** The paper correctly notes that the crossover is state-dependent (present for Haar and thermal; absent for ground and mixed up to $\tau=1.0$).
* **"Derived, not fitted":** The term is mathematically justified regarding the calculation of $\tau_c$: once the two structural basins $F_A$ and $F_B$ are identified, $\tau_c$ is predicted using only initial-time ($t=0$) operator matrix elements ($I_F$, $Q_F\rho$, $Q_F L\rho$, $Q_F L^2\rho$), without fitting curves to finite-$\tau$ data. (See Major Finding M1 for required clarification).

### 4. Falsifiability / Physical Meaning
* **Functional Falsification:** Adopting the unnormalized functional is not ad hoc; it directly repairs the ill-posed domain topology by removing an unconstrained denominator that approaches zero on open neighborhoods of pure product states.
* **Liouvillian Gap Rejection:** Refuting the hypothesis $\tau_c \sim \Delta_L^{-1}$ ($\tau_c \Delta_L \approx 0.012$, off by 83–136$\times$) is handled with exemplary physical clarity. It proves that the crossover is an objective-competition phenomenon rather than a global relaxation timescale.

### 5. Journal-Level Completeness
* Abstract, Intro, Formalism, Discussion, and Limitations are complete, balanced, and readable.
* Figures 1–4 cleanly communicate the main findings, and the separation between the frozen specification (`REM Spec v2.2`) and post-hoc numerical data (`Phase D/E Results`) is maintained.

---

## Tagged Findings (Critical / Major / Minor)

### Major Findings

* **[Major] M1. Precise scope clarification of "derived, not fitted" (Claim 4c)**  
  *Context:* Sec. 1, Sec. 8, and Abstract state that $\tau_c$ is "derived, not fitted".  
  *Audit:* The calculation of $\tau_c^{(1)} = -\Delta\Phi_0 / \Delta\Phi_1$ is indeed an analytical prediction derived from $t=0$ local operator derivatives ($I_F, Q_F\rho_U, Q_FY_U, Q_FZ_U$). However, finding the representative configuration $F_B$ requires variational optimization sampled at a finite time ($\tau = 0.1$).  
  *Requirement:* To prevent readers from misinterpreting "derived" as implying $F_B$ was predicted *a priori* without numerical optimization, the text in Sec. 1 and Sec. 8 must explicitly state: *"The crossover timescale $\tau_c$ is derived analytically from initial-time ($t=0$) objective derivatives at fixed basin configurations $(F_A, F_B)$, rather than fitted to the time-dependent curve $\Delta\Phi(\tau)$."*

### Minor Findings

* **[Minor] m1. Unified notation for rotated frame in Sec. 2.2 vs Sec. 2.3 & 8**  
  In Sec. 2.2, $J_{\mathrm{dyn}}^{(0)}$ is written as $-\operatorname{Re}\langle Q_F\rho, Q_FL(\rho)\rangle$, omitting the unitary rotation frame $U$. In Sec. 2.3 and Sec. 8, it is written precisely using $\rho_U = U^\dagger \rho U$ and $Y_U = U^\dagger L(\rho) U$. Sec. 2.2 should adopt $Q_F \rho_U$ or $Q_F(U^\dagger \rho U)$ for absolute formal consistency.

* **[Minor] m2. Ground vs. Mixed state cross-evaluation gap in Sec. 5**  
  In the D2-R cross-evaluation matrix (Sec. 5), evaluating $F^*_{\mathrm{mixed}}$ under the ground-state objective gives $\Phi = 1.899$, compared to $1.900$ for $F^*_{\mathrm{ground}}$ (a tiny gap of $0.001$). The text already acknowledges this ("contain nearly degenerate optima under the ground-state objective"), but adding a brief sentence noting that ground and mixed ground-state factorizations are near-isometric on the quotient would strengthen the geometric discussion.

* **[Minor] m3. Sample size confirmation for $N=5$ Full Closure**  
  The $N=5$ full closure used 2 seeds $\times$ 100 steps (documented in Sec. 9.3 and Supplementary S3). While the seed standard deviation is extremely low ($2 \times 10^{-4}$), adding a brief note in the main text table caption (Table in Sec. 7) explicitly marking $N=5$ as $N_{\mathrm{seeds}}=2$ will ensure full protocol transparency.

---

## Answering the Key Question

> **"If you ignore the authors' interpretation and inspect only the definitions, derivations, numerical evidence, and figures, do Claims 1-4 actually follow from the presented evidence?"**

**YES.** Stripping away all narrative interpretation:
1. **Claim 1:** Follows directly from Fig. 1 and D2.1 numerical data ($\mathrm{corr} = -0.999$, $C_F^2 \to 10^{-4}$, $\Phi \to 25.7$, $p_{\max} = 0.9998$).
2. **Claim 2:** Follows directly from Table 2 and D2.2 data ($C_F^2 \ge 0.0956$, $p_{\max} \le 0.58$, seed std $\le 0.009$).
3. **Claim 3:** Follows directly from the cross-evaluation matrices in Sec. 5 (diagonal dominance with mean gaps $0.113$ and $0.159$).
4. **Claim 4a–4c:** Follows directly from Figs. 2, 3, 4 and Phase E data ($\Delta\Phi(\tau_c)=0$ crossings at $\tau_c \approx 0.018\text{–}0.022$, persistence across $N=3\text{–}5$, and analytic Taylor expansion agreement within $5.1\text{–}10.9\%$ at 1st order and $0.0\text{–}0.5\%$ at 2nd order).

---

## Final Recommendation

**Recommendation:** **Minor revision**

The manuscript is mathematically sound, numerically verified, and rigorously argued. Addressing Major Finding M1 (clarification of the scope of "derived" for $\tau_c$) and Minor Findings m1–m3 will render the paper fully ready for publication.