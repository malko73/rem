# REM Specification v2.1

**Version**: 2.1  
**Date**: 2026-08-12  
**Status**: Frozen  
**Supersedes**: v2.0 (2026-08-11)

## Summary of Changes from v2.0

Phase C3.3b revealed a critical distinction: the instantaneous decoherence rate C_dyn^(0) and the finite-time decoherence functional C_dyn^(τ) are **not approximations of each other**, but **distinct structure-selection functionals operating at different timescales**.

Key findings:
- C_dyn^(0) = Γ_F^exact(0) captures only instantaneous decoherence
- C_dyn^(τ) = -(1/τ)log(C_F(τ)/C_F(0)) captures finite-time coherence preservation
- Structures selected at τ=0 and τ>0 can be fundamentally different (d_F ≈ 0.8)
- τ is an **independent physical parameter**, not a technical detail

This elevates the central equation to explicitly include τ as a variable.

---

## 1. Central Equation

The REM objective functional is:

$$
\boxed{
\Phi(F; \lambda, \rho, \mathcal{L}, \tau) = I_\rho(F) - \lambda C_{\text{dyn}}^{(\tau)}(F; \rho, \mathcal{L})
}
$$

where:
- $F$: factorization structure (unitary rotation U defining subsystem decomposition)
- $\lambda$: trade-off parameter between information and dynamics
- $\rho$: quantum state (ground state, thermal state, or general mixed state)
- $\mathcal{L}$: Lindbladian generator (open-system dynamics)
- $\tau$: observation/coarse-graining timescale (independent parameter)
- $I_\rho(F)$: mutual information of factorization F with respect to state ρ
- $C_{\text{dyn}}^{(\tau)}(F; \rho, \mathcal{L})$: dynamical stability functional at timescale τ

The optimal factorization is:

$$
\boxed{
F^* = F^*(\rho, \mathcal{L}, \lambda, \tau)
}
$$

**Key insight**: The factorization structure depends not only on the system and environment, but also on the **timescale at which structure is evaluated**.

---

## 2. Dynamical Stability Functionals

### 2.1 Finite-Time Functional (Primary)

$$
C_{\text{dyn}}^{(\tau)}(F; \rho, \mathcal{L}) = -\frac{1}{\tau} \log \frac{C_F(\tau)}{C_F(0)}
$$

where:
- $C_F(t) = \|\rho_F(t) - \mathcal{D}_F[\rho_F(t)]\|_F^2$
- $\rho_F(t) = U^\dagger \rho(t) U$ (state in factorization frame)
- $\rho(t) = e^{t\mathcal{L}}[\rho]$ (time evolution under Lindbladian)
- $\mathcal{D}_F$: dephasing projection onto Schmidt basis of F

This functional measures the average rate of coherence loss over the time interval [0, τ].

### 2.2 Instantaneous Limit (Secondary)

$$
C_{\text{dyn}}^{(0)}(F; \rho, \mathcal{L}) = \lim_{\tau \to 0} C_{\text{dyn}}^{(\tau)}(F; \rho, \mathcal{L}) = \Gamma_F^{\text{exact}}(0)
$$

where:

$$
\Gamma_F^{\text{exact}}(0) = -\frac{2 \text{Re} \langle Q_F \rho, Q_F \mathcal{L}(\rho) \rangle_{\text{HS}}}{\|Q_F \rho\|_F^2}
$$

with $Q_F = I - \mathcal{D}_F$.

**Important**: C_dyn^(0) captures only instantaneous decoherence. It does **not** predict finite-time structure selection (verified in Phase C3.3b, d_F ≈ 0.8).

### 2.3 Relationship

C_dyn^(0) and C_dyn^(τ) are related by:

$$
C_{\text{dyn}}^{(0)} = \lim_{\tau \to 0} C_{\text{dyn}}^{(\tau)}
$$

However, this is a **limit relationship**, not an approximation. For finite τ, the two functionals can select fundamentally different structures.

---

## 3. Factorization Structure

A factorization F is defined by a unitary U ∈ U(2^n) that rotates the computational basis into a tensor product structure:

$$
\mathcal{H} = \mathcal{H}_A \otimes \mathcal{H}_B
$$

The space of factorizations is the quotient:

$$
\mathcal{F} = U(2^n) / \text{Stab}(A|B)
$$

where Stab(A|B) is the stabilizer subgroup preserving the bipartition.

For n=3 with cut A|B = {0,1}|{2}, dim(ℱ) = 45.

---

## 4. Optimization

The optimal factorization is found by:

$$
F^* = \arg\max_{F \in \mathcal{F}} \Phi(F; \lambda, \rho, \mathcal{L}, \tau)
$$

Optimization is performed on the 45-dimensional quotient manifold using:
- Adam optimizer with finite-difference gradients
- Multiple random initializations (matched-seed protocol)
- Gauge-invariant distance metric d_F for comparing factorizations

---

## 5. Physical Interpretation

### 5.1 Timescale Dependence

The factorization F* depends on the observation timescale τ:

- **τ → 0**: Selects structures that minimize instantaneous decoherence rate
- **τ > 0**: Selects structures that optimize finite-time coherence preservation
- **τ → ∞**: Selects structures that maximize asymptotic coherence (if it exists)

Phase C3.3b demonstrated that **structures selected at different τ can be fundamentally different** (d_F ≈ 0.8). This means:

> "Relation" in REM includes not only the environment (ρ, L) but also the **timescale τ** at which the relation is evaluated.

### 5.2 Environmental Role

The environment (ρ, L) plays a dual role:
1. **Defines the dynamics**: L determines how coherence decays
2. **Sets natural timescales**: Liouvillian gap, decoherence rates, etc.

However, the choice of τ is **not uniquely determined by the environment**. It represents the observer's coarse-graining scale or the timescale of interest.

### 5.3 Comparison with v2.0

In v2.0, the framework was:

$$
\text{environment} \rightarrow \Gamma_F^{(0)} \rightarrow F^*
$$

In v2.1, the framework is:

$$
(\rho, \mathcal{L}, \tau) \rightarrow C_{\text{dyn}}^{(\tau)} \rightarrow F^*
$$

This is a significant conceptual shift: **REM now recognizes that "relational structure" is timescale-dependent**.

---

## 6. Gates (Validation Criteria)

### Gate 1: Non-trivial Structure Existence ✅ CONDITIONAL PASS
- Quotient-correct 45D optimization finds F* with Φ* significantly above contiguous baselines
- Verified by matched-seed control and Hessian analysis
- **Conditional**: Local maximum not rigorously proven (finite-difference Hessian unreliable)

### Gate 2: Open-System Dynamical Functional ✅ PASS
- C_Γ^(0) = Γ_F^exact(0) is a well-defined, gauge-invariant operational cost
- Validated through C0-C2.5 (well-posedness, gauge invariance, frame consistency)

### Gate 3a: λ Determination (τ fixed) 🔲 OPEN
- Given (ρ, L, τ), can λ be determined independently of factorization selection?
- Requires identifying a physical principle or measurement protocol for λ

### Gate 3b: τ Determination 🔲 OPEN (Stronger)
- Can τ be determined from (ρ, L) alone?
- Candidate: Liouvillian gap τ_L = |Re μ₁|⁻¹
- Phase C3.1-C3.2 showed τ_L does not move factorization structure (d_F < 0.1)
- But τ_L does not predict finite-time structure selection (C3.3b FAIL)

**Note**: Gate 3 is split into 3a and 3b because τ may represent an observer timescale rather than an environment property. Gate 3a (λ determination with τ given) is weaker and more achievable. Gate 3b (τ determination from environment) is stronger and may not be necessary.

### Gate 4: Predictive Power 🔲 OPEN
- Can F*(ρ, L, λ, τ) predict structure in unseen environments?
- Requires validation on hold-out test cases

### Gate 5: Generality 🔲 OPEN
- Does the framework generalize beyond 3-qubit XY chain?
- Requires testing on multiple Hamiltonian families (Phase D)

---

## 7. Benchmark Results (Phase A-C)

### 7.1 System
- 3-qubit XY chain with asymmetric couplings
- J₁₂ = 1.5, J₂₃ = 0.6, h = 0.2
- Pure dephasing environment: γ = (0.5, 1.0, 2.0)
- Temperature: T = 0 (ground state)

### 7.2 Key Findings

**Phase A-B (Closed-System Surrogate)**:
- C_H^closed = ⟨H_∂F²⟩ identifies non-trivial factorizations
- Common optimal basin across environments (B2-B)
- Hessian analysis reveals saddle-point landscape (B3.5)

**Phase C (Open-System Canonical)**:
- C_Γ^(0) is well-defined and gauge-invariant
- **C_Γ^(0) ≠ C_dyn^(τ)**: Different timescale functionals
- Finite-time prediction fails (d_F ≈ 0.8)
- τ is an independent physical parameter

### 7.3 Numerical Values

For benchmark system with λ = 0.2, τ = 0.1:
- Φ* ≈ 1.56 (A2 environment)
- d_F(F*_predicted, F*_observed) ≈ 0.80
- Γ_F^exact(0) ≈ 2.62 (B2 best structure)
- C_dyn^(0.1) ≈ 2.38 (B2 best structure)

---

## 8. Implementation Notes

### 8.1 Coherence Norm Definition

The coherence norm uses the **squared Frobenius norm**:

$$
C_F(t) = \|\rho_F(t) - \mathcal{D}_F[\rho_F(t)]\|_F^2
$$

This introduces a **factor of 2** in the instantaneous rate:

$$
\Gamma_F^{\text{exact}}(0) = -\frac{2 \text{Re} \langle Q_F \rho, Q_F \mathcal{L}(\rho) \rangle_{\text{HS}}}{\|Q_F \rho\|_F^2}
$$

(Corrected in Phase C3.3c)

### 8.2 Finite-Difference Hessian

The finite-difference Hessian on the 45D quotient manifold is **unreliable**:
- Spurious positive eigenvalues (saddle-point artifacts)
- Direct directional curvature tests show negative curvature
- **Recommendation**: Do not use Hessian eigenvalues for local-maximum verification

### 8.3 Gauge-Invariant Distance

Factorizations are compared using the projector-based distance:

$$
d_F(F_1, F_2) = \|P_{\mathcal{A}_{F_1}} - P_{\mathcal{A}_{F_2}}\|_F
$$

where $P_{\mathcal{A}_F}$ is the projector onto the local operator algebra of factorization F.

This distance is invariant under local gauge transformations V_A ⊗ V_B.

---

## 9. Open Questions

1. **Physical meaning of τ**: Is τ an observer timescale, or does the environment select a natural τ?
2. **λ determination**: Can λ be derived from microscopic principles?
3. **Timescale crossover**: Does F*(τ) show continuous deformation or sharp transitions?
4. **Generality**: Do the Phase A-C findings hold for other Hamiltonians and system sizes?

---

## 10. Version History

- **v1.0** (2026-04-25): Initial formulation with C_H^closed
- **v1.1** (2026-07-18): Added monotonic tradeoff theorem
- **v2.0** (2026-08-11): Open-system canonical form with C_Γ^(0)
- **v2.1** (2026-08-12): Timescale-dependent formulation with C_dyn^(τ)

---

## 11. References

- Phase A-C analysis scripts: `analysis/b*.py`, `analysis/c*.py`
- Numerical results: `analysis_output/b*.json`, `analysis_output/c*.json`
- Final reports: `analysis/PHASE_B_FINAL_REPORT.md`, `analysis/PHASE_C_FINAL_REPORT.md`
- Zenodo record: DOI 10.5281/zenodo.21880505 (REM_lambda v5)
