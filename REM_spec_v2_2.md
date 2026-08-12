# REM Specification v2.2

**Version**: 2.2  
**Date**: 2026-08-12  
**Status**: Frozen (canonical candidate; D1/D2 re-verification, D3, D4 open)  
**Supersedes**: v2.1 (2026-08-12), v2.0 (2026-08-11)

## Summary of Changes from v2.1

v2.1 made the normalized instantaneous rate $\Gamma_F^{\mathrm{exact}}(0)$ an
operational realization of the dynamical cost and elevated a finite-time
log-ratio functional to primary status. **Phase D2.1 falsified the normalized
rate as a global variational objective**: on unrestricted TPS optimization
the denominator $|Q_F\rho|^2$ can be driven to zero, sending
$\Gamma_F \to \pm\infty$ and $\Phi \to +\infty$ (a pure state can always be
approached by a product-state TPS). The v2.2 canonical functional is
therefore the **denominator-free (unnormalized) instantaneous functional**

$$
\boxed{\,J_{\rm dyn}^{(0)}(F;\rho,\mathcal L) =
-\operatorname{Re}\langle Q_F\rho,\; Q_F\mathcal L(\rho)\rangle_{\mathrm{HS}}\,}
$$

with

$$
\boxed{\,\Phi(F;\lambda) = I_\rho(F) - \lambda\, J_{\rm dyn}^{(0)}(F)\,}
$$

as the canonical variational functional. The normalized rate is demoted to
a fixed-$(F)$ local diagnostic (Section 7). The finite-time extension
$J_{\rm dyn}^{(\tau)}$ is repositioned as an **operational finite-time
extension / consistency diagnostic**, with the exact $\tau\to 0$ limit
$\lim_{\tau\to0}J_{\rm dyn}^{(\tau)}=J_{\rm dyn}^{(0)}$ as the main
consistency condition (Section 9). Phase D2.2-A/B validated this structure:
G1–G6 PASS (D2.2-A) and G7 **PASS with documented finite-$\tau$ crossover**
(D2.2-B): for the Haar state, $\tau \ge 10^{-2}$ induces a genuine basin
transition of $F^*$ — a new result that defines the D3 phase.

Additional corrections relative to v2.1:

- v2.1's $\Gamma_F^{\mathrm{exact}}$ carried a factor of 2 (from the squared
  norm convention). The C0 pre-registration defines the rate as
  $-\mathrm{d}\log C_F/\mathrm{d}t$ with $C_F$ the Frobenius *norm*, giving
  **no factor of 2**; the numerical-derivative test and the master's D2.2
  formula agree. Reverted in commit `fb08a3f`.
- v2.1's finite-time functional was the log-ratio
  $-(1/\tau)\log(C_F(\tau)/C_F(0))$. v2.2 freezes the difference ratio
  $J_{\rm dyn}^{(\tau)} = -(C_F^2(\tau)-C_F^2(0))/(2\tau)$, which has the
  exact derivative-based $\tau\to0$ limit above.

---

## 1. Scope / Status

REM treats the tensor-product structure of Hilbert space as a variable
structure selected by competing informational and dynamical constraints.
REM is a structural selection framework, not a modification of quantum
dynamics.

This specification freezes the **canonical variational functional** of
open-system REM after the D2.1 falsification:

- **Canonical**: $J_{\rm dyn}^{(0)}(F;\rho,\mathcal L)$ (Section 4) and
  $\Phi(F;\lambda) = I_\rho - \lambda J_{\rm dyn}^{(0)}$ (Section 5).
- **Prohibited as global objective**: the normalized rate
  $\Gamma_F^{\mathrm{exact}}(0)$ (Section 7).
- **Operational extension**: $J_{\rm dyn}^{(\tau)}$ (Section 8) with the
  $\tau\to0$ consistency relation (Section 9) and the admissible finite-$\tau$
  TPS crossover (Section 10).

Status: **frozen candidate**. Remaining validation gates: D1/D2 re-verification
under the v2.2 canonical functional, D3 (timescale / $\tau$-crossover), D4
(N-scaling) — Section 12.

## 2. Factorization domain (F)

A factorization $F$ is defined by a unitary $U \in U(2^n)$ that rotates the
computational basis into a bipartite tensor-product structure:

$$
\mathcal H \simeq \mathcal H_A \otimes \mathcal H_B,\qquad
n_A + n_B = n .
$$

The space of factorizations is the quotient

$$
\mathcal F = U(2^n)\big/ \mathrm{Stab}(A|B),
$$

where $\mathrm{Stab}(A|B)$ is the stabilizer of the bipartition. The quotient
has real dimension $D^2 - d_A^2 - d_B^2 + 1$ (the $U(1)$ kernel is included);
for $n=3$, $n_A=2$: $\dim \mathcal F = 45$.

The domain is **unrestricted**: $F$ ranges over all TPS rotations, not only
contiguous site cuts. This is the setting in which the D2.1 singularity
arises (Section 7). Factorizations are compared with the gauge-invariant
distance $d_F(F_1,F_2) = \|P_{\mathcal A_{F_1}} - P_{\mathcal A_{F_2}}\|_F$
(projector onto the local operator algebra), invariant under
$V_A \otimes V_B$ gauge transformations.

## 3. Informational functional ($I_\rho(F)$)

$$
I_\rho(F) = I(A:B)_\rho = S(\rho_A) + S(\rho_B) - S(\rho_{AB}),
\qquad
S(\rho) = -\mathrm{Tr}(\rho\log\rho).
$$

$I_\rho$ is the mutual information of the (possibly mixed) state $\rho$
across the factorization $F$. Units are nats or bits depending on the
logarithm base (the D-series numerics use $\log_2$; $I \le \min(\log_2 d_A,
\log_2 d_B)\cdot 2$ for pure states). The frozen D-series protocol evaluates
$I$ on the dominant eigenvector of the rotated state.

## 4. Canonical instantaneous dynamical functional ($J_{\rm dyn}^{(0)}$)

$$
\boxed{\,
J_{\rm dyn}^{(0)}(F;\rho,\mathcal L)
= -\operatorname{Re}\langle Q_F\rho,\; Q_F\mathcal L(\rho)\rangle_{\mathrm{HS}}
\,}
$$

with

- $Q_F = I - D_F$, where $D_F$ is the dephasing projection onto the initial
  Schmidt basis of $F$: $D_F[\sigma] = \sum_k |b_k\rangle\langle b_k|
  \sigma |b_k\rangle\langle b_k|$,
- $\mathcal L$ the Lindblad Liouvillian of the fixed environment,
- $\langle A,B\rangle_{\mathrm{HS}} = \mathrm{Tr}(A^\dagger B)$.

**Name**: the **instantaneous dynamical rate functional (signed)**. The
functional is **signed** and must never be called a "cost" unconditionally:
$J_{\rm dyn}^{(0)} > 0$ when the initial coherence $C_F^2(t)=|Q_F\rho(t)|^2$
decays (penalty in $\Phi$), $J_{\rm dyn}^{(0)} < 0$ when it grows (reward in
$\Phi$). By Cauchy–Schwarz,

$$
|J_{\rm dyn}^{(0)}| \le |Q_F\rho|_2\, |Q_F\mathcal L(\rho)|_2 \le
|Q_F\rho|_2\, \|\mathcal L(\rho)\|_2,
$$

so $J_{\rm dyn}^{(0)} \to 0$ smoothly as $|Q_F\rho|^2 \to 0$: **no
denominator blow-up on the unrestricted TPS domain**. This is the defining
property that makes $J_{\rm dyn}^{(0)}$ well-posed as a global variational
functional (D2.1 vs D2.2-A).

## 5. Canonical objective ($\Phi = I - \lambda J_{\rm dyn}^{(0)}$)

$$
\boxed{\,
\Phi(F;\lambda) = I_\rho(F) - \lambda\, J_{\rm dyn}^{(0)}(F;\rho,\mathcal L)
\,}
$$

Optimal factorization:

$$
F^*(\lambda) = \arg\max_{F\in\mathcal F} \Phi(F;\lambda).
$$

**Monotonic Tradeoff Theorem** (unchanged from v2.0, Section 8a): for
$\lambda_2 > \lambda_1 \ge 0$,

$$
J_{\rm dyn}^{(0)}(F^*(\lambda_2)) \le J_{\rm dyn}^{(0)}(F^*(\lambda_1)),
\qquad
I_\rho(F^*(\lambda_2)) \le I_\rho(F^*(\lambda_1)),
$$

for a **fixed** environment $(\rho,\mathcal L)$. The theorem holds for the
signed functional with $\lambda \ge 0$ (the proof requires only
non-negativity of $\lambda$, not of the dynamical term).

## 6. Dimensional convention ($[\lambda] = T$)

$\rho$ is dimensionless; $\mathcal L$ has units $T^{-1}$; hence

$$
[J_{\rm dyn}^{(0)}] = T^{-1},
\qquad
[\lambda] = T,
\qquad
[\lambda J_{\rm dyn}^{(0)}] = 1 .
$$

$\lambda$ is the characteristic timescale that converts the instantaneous
rate into a dimensionless contribution to $\Phi$. The optional rescaled
representation $\tilde J = \tau_0 J$, $\tilde\lambda = \lambda/\tau_0$
(with $\Phi = I - \tilde\lambda \tilde J$) is mathematically equivalent and
introduces no new physics; the canonical form keeps $[\lambda]=T$ with no
reference time $\tau_0$.

## 7. Normalized $\Gamma_F$: local-diagnostic-only status

The normalized instantaneous rate

$$
\Gamma_F^{\mathrm{exact}}(0)
= -\frac{\operatorname{Re}\langle Q_F\rho,\; Q_F\mathcal L(\rho)\rangle_{\mathrm{HS}}}
{|Q_F\rho|_2^2}
$$

is **retained but demoted**:

> **fixed-$(F)$ local diagnostic only; prohibited as an unrestricted global
> TPS variational objective.**

Rationale (Phase D2.1, commit `0ad8840`): on unrestricted TPS optimization
the optimizer drives $|Q_F\rho|^2 \to 0$, which sends $\Gamma_F \to
\pm\infty$ and $\Phi \to +\infty$ along a non-physical direction. Observed
evidence:

- $\operatorname{corr}(\Phi, \log_{10} C_F(0)^2) = -0.999$ (near-perfect
  anti-correlation),
- $C_F(0)^2 \to 6\times10^{-5}$ drives $\Gamma \to -222$, $\Phi \to 44$,
- best solution Schmidt spectrum $p = [0.9997, 0.0003]$ (near-product
  collapse).

Because every pure state can be approached by a product-state TPS, the
denominator can vanish; the normalized rate is therefore well-posed only for
a **fixed** $F$ (where $|Q_F\rho|^2$ is bounded away from zero) as a local
diagnostic of structural decoherence. The v2.1 elevation of $\Gamma_F$ to an
operational canonical cost is **abandoned** on this falsification.

## 8. Finite-time extension ($J_{\rm dyn}^{(\tau)}$)

$$
\boxed{\,
J_{\rm dyn}^{(\tau)}(F;\rho,\mathcal L)
= -\frac{C_F^2(\tau) - C_F^2(0)}{2\tau}
\,}
$$

with $C_F^2(t) = |Q_F \rho(t)|_2^2$, $\rho(t) = e^{t\mathcal L}[\rho]$, and
$Q_F$ built from the **fixed** $t=0$ Schmidt basis of $F$.

**Status**: **operational finite-time extension / consistency diagnostic** —
separate from the canonical functional, not an alternative to it. It
measures the average rate of squared-coherence change over $[0,\tau]$ and
is the natural finite-window probe of structural stability.

## 9. $\tau\to0$ consistency relation

For fixed $F$ (hence fixed $Q_F$),

$$
\frac{\mathrm{d}}{\mathrm{d}t} C_F^2(t)\bigg|_{t=0}
= 2\operatorname{Re}\langle Q_F\rho,\; Q_F\mathcal L(\rho)\rangle_{\mathrm{HS}},
$$

so

$$
\boxed{\,
\lim_{\tau\to0} J_{\rm dyn}^{(\tau)}(F;\rho,\mathcal L)
= J_{\rm dyn}^{(0)}(F;\rho,\mathcal L)
\,}
$$

is the **main consistency condition of v2.2**: the finite-time extension
reproduces the canonical instantaneous functional in the $\tau\to0$ limit.
Verified numerically (D2.2-B, G7-1): deviation scales linearly in $\tau$
($|\Delta J| \approx |J_0|\,\tau$), with $|\Delta J| \le 2.5\times10^{-3}$
at $\tau = 10^{-3}$ for all four state classes.

## 10. Finite-$\tau$ TPS crossover (admissible)

The consistency relation does **not** imply that $F^*$ is $\tau$-independent
at finite $\tau$:

> **Finite $\tau$ need not preserve the instantaneous optimum; $\tau$ may
> induce genuine transitions between competing TPS basins.**

Evidence (D2.2-B, Haar state): for $\tau \le 3\times10^{-3}$ all seeds stay
in the instantaneous basin ($d_F < 0.06$); for $\tau \ge 10^{-2}$ all seeds
move to a different basin ($I$: 1.97 → 2.00, $J_{\rm dyn}^{(\tau)}$:
0.60 → 0.44, $\Phi$ up to +4.3% at $\tau = 0.1$). The shift is robust
(6/6 seeds), monotone in $\tau$, and free of pathology (no singularity, no
product collapse, no seed instability). It is read as a **genuine structural
crossover induced by finite-time coarse-graining**, not as noise and not as
a failure of the canonical functional. This is the defining phenomenon of
the D3 phase (Section 12).

## 11. Known numerical evidence (D-series, 2026-08-12)

All runs use the frozen protocol: asymmetric-XY chain ($J_{12}=1.5$,
$J_{23}=0.6$, $h=0.2$), pure dephasing $\gamma=(0.5,1.0,2.0)$, $\lambda=0.2$,
Adam 200 steps, $\mathrm{lr}=0.01$, seeds $20260813\ldots$, quotient
$\dim\mathcal F=45$, $n_A=2$. Suite: 74 passed + 1 xfailed.

| Phase | Result | Verdict |
|---|---|---|
| D1 (Hamiltonian generality, Γ-based) | 5/5 families show non-trivial $F^*$ | preliminary (to be re-run under $J_{\rm dyn}^{(0)}$) |
| D2 (state generality, Γ-based) | ground/mixed/thermal PASS-candidate; Haar anomaly ($\Phi=49.1$) | re-scoped to D2.1 |
| **D2.1** (Haar singularity audit, `0ad8840`) | normalization singularity confirmed | **falsification** → Section 7 |
| **D2.2-A** (unnormalized instantaneous, `24db549`) | G1–G6 **PASS** | $J_{\rm dyn}^{(0)}$ well-posed |
| **D2.2-B** (finite-time, `2fce54d`) | G7 **PASS with documented finite-$\tau$ crossover** | extension + Section 10 |

D2.2-A numbers (4 states × 6 seeds): best $\Phi$ = ground 1.8995, Haar
1.8594, mixed 1.9845, thermal 1.9465; Haar 30-seed best 1.8595 / std 0.0077
(vs 25.71 / 5.96 under the normalized objective); min $C_F^2 = 0.096$ (never
approaches zero); Schmidt $p_{\max} = 0.51$–$0.58$ (no collapse); small-
perturbation continuity max $|\Delta\Phi| = 7.7\times10^{-5}$.

D2.2-B numbers (τ sweep $\{10^{-3}, 3\times10^{-3}, 10^{-2}, 3\times10^{-2},
10^{-1}\}$): G7-1 exact $\tau\to0$ limit (above); G7-2 same-basin for all
states at $\tau\le3\times10^{-3}$ and for ground/mixed/thermal at all $\tau$
(cross-eval gaps $\le 0.005$), Haar crossover at $\tau\ge10^{-2}$ (Section
10); G7-3 state ranking (mixed > thermal > ground > Haar) preserved at all
$\tau$; G7-4 min $C_F^2(0) \ge 0.09$, $p_{\max} \le 0.58$, std $\le 0.009$
at all $(\mathrm{state},\tau)$; Haar 30-seed at $\tau=0.1$: best 1.9135 /
std 0.0011.

## 12. Open validation gates

1. **D1 re-verification**: Hamiltonian generality (5 families) under the
   v2.2 canonical functional $J_{\rm dyn}^{(0)}$.
2. **D2 re-verification**: state generality (ground / Haar / mixed /
   thermal) under $J_{\rm dyn}^{(0)}$ (D2.2-A is the Haar-focused audit;
   the full 4-state protocol is the re-run).
3. **D3 Timescale**: promoted from "how to determine $\tau$" to
   **$\tau$-crossover analysis**: locate the basin-transition scale
   $\tau_c$ (for Haar the switch lies between $3\times10^{-3}$ and
   $10^{-2}$), characterize its sharpness, and connect it to
   environment timescales (e.g. Liouvillian gap).
4. **D4 N-scaling**: $N=4,5,\ldots$ behavior of $F^*$, $J_{\rm dyn}^{(0)}$,
   and the crossover.
5. Existing roadmap gates 3a/3b (λ/τ determination), Gate 4 (predictive
   hold-out), Gate 5 (generality) remain OPEN.

---

## Version History

- **v1.0** (2026-04-25): initial formulation with closed-system surrogate
  $C_H^{\mathrm{closed}} = \langle H_{\partial F}^2\rangle$.
- **v1.1** (2026-07-18): monotonic tradeoff theorem; sign-convention fix.
- **v2.0** (2026-08-11): open-system canonical form
  $\Phi = I - \lambda C_{\mathrm{dyn}}$ with $[\lambda]=T$; first
  realization $C_\Gamma^{(0)} = \Gamma_F^{\mathrm{exact}}(0)$; $C_H^{\mathrm{closed}}$
  downgraded to surrogate; signed dynamical functional named.
- **v2.1** (2026-08-12): timescale-dependent formulation; finite-time
  log-ratio functional made primary; $\Gamma_F$ retained as canonical cost.
  **Superseded by v2.2** on the D2.1 falsification (Section 7): the
  normalized rate is singular on unrestricted TPS and is demoted to a
  local diagnostic; the unnormalized instantaneous functional is frozen as
  canonical (Sections 4–5); the finite-time difference ratio is frozen as
  the operational extension (Section 8).
- **v2.2** (2026-08-12): this document. Canonical
  $J_{\rm dyn}^{(0)} = -\operatorname{Re}\langle Q_F\rho, Q_F\mathcal L(\rho)\rangle$;
  $\Phi = I - \lambda J_{\rm dyn}^{(0)}$; $[J] = T^{-1}$, $[\lambda]=T$;
  G7 recorded as PASS with documented finite-$\tau$ crossover.

## References

- D-series scripts: `analysis/d0_protocol.py`, `analysis/d1_generality.py`,
  `analysis/d2_state_generality.py`, `analysis/d2_1_haar_singularity_audit.py`,
  `analysis/d2_2a_unnormalized_instantaneous.py`,
  `analysis/d2_2b_finite_time.py`, `analysis/d2_2b_gates.py`.
- Numerical results: `analysis_output/d2_*.json`, `analysis_output/d2_2b_*.json`.
- Roadmap: `REM_next_steps.md` (D0–D2.2-B status).
- Zenodo: DOI 10.5281/zenodo.21880505 (REM_lambda v5); v6 correction plan
  will include D2.1–D2.2 and Spec v2.2.
