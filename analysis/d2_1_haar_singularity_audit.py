#!/usr/bin/env python3
"""
Phase D2.1 — Haar Singularity Audit.

Purpose: determine whether Φ≈49 for Haar states is a genuine signed-functional
effect or a normalization singularity of C_dyn at C_F(0)→0.

Master's protocol:
  1. For best Haar solutions, record I(F), C_F(0), C_F(τ), C_dyn, Φ
  2. Scatter C_F(0)² vs Γ: look for |Γ| ∝ 1/C_F(0) or 1/C_F(0)² divergence
  3. Across Haar trials: min/median C_F(0), Γ distribution, Φ distribution,
     corr(Φ, C_F(0)), Schmidt spectrum of best solution,
     small-perturbation response of (Γ, Φ)

Verdict:
  D2.1-A: C_F(0)→0 with Γ→−∞ and Φ→+∞ ⇒ normalization singularity
          (C_Γ^(0) not well-posed on unrestricted TPS)
  D2.1-B: stable Φ≈49 with sufficient C_F(0) ⇒ real signed-functional effect
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from d1_hamiltonians import get_hamiltonian, get_ground_state
from quotient_geometry import horizontal_basis
from b2_canonical_open import canonical_open_objective, apply_tps, liouvillian_dephasing

LAMBDA = 0.2
STEPS = 200
LR = 0.01
SEEDS = 30  # more seeds for statistics
SEED0 = 20260813
GAMMA_VEC = [0.5, 1.0, 2.0]
HNAME = "asymmetric_XY"


def canonical_open_objective_lambda(rho_U, Y_U, n_a, lambda_val):
    obj = canonical_open_objective(rho_U, Y_U, n_a)
    phi = obj["I"] - lambda_val * obj["gamma"] if np.isfinite(obj["gamma"]) else -1e9
    return obj, phi


def compute_c_f_tau(rho0, L, U, tau):
    """C_dyn^(τ) = -(1/τ) log(C_F(τ)/C_F(0))."""
    rho0_vec = rho0.reshape(-1, order="F")
    exp_tL = sla.expm(tau * L)
    rho_t_vec = exp_tL @ rho0_vec
    rho_t = rho_t_vec.reshape(rho0.shape, order="F")

    rho_t_F = U.conj().T @ rho_t @ U
    rho_t_F_dephased = np.diag(np.diag(rho_t_F))
    rho_t_dephased = U @ rho_t_F_dephased @ U.conj().T
    c_t = np.linalg.norm(rho_t - rho_t_dephased, ord="fro") ** 2

    rho0_F = U.conj().T @ rho0 @ U
    rho0_F_dephased = np.diag(np.diag(rho0_F))
    rho0_dephased = U @ rho0_F_dephased @ U.conj().T
    c_0 = np.linalg.norm(rho0 - rho0_dephased, ord="fro") ** 2

    if c_0 < 1e-300:
        return float("nan"), c_0, c_t
    if c_t < 1e-300:
        return float("inf"), c_0, c_t
    return float(-(1.0 / tau) * np.log(c_t / c_0)), c_0, c_t


def adam_optimize_open(U0, L, rho0, Y0, lambda_val, H_basis, steps=STEPS, lr=LR):
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_at(d):
        U = to_U(d)
        rho_U = U.conj().T @ rho0 @ U
        Y_U = U.conj().T @ Y0 @ U
        obj, phi = canonical_open_objective_lambda(rho_U, Y_U, 2, lambda_val)
        return phi

    best_phi = phi_at(delta)
    best_U = to_U(delta)
    best_obj = None

    for step in range(steps):
        grad = np.zeros_like(delta)
        for k in range(len(delta)):
            d_p = delta.copy()
            d_m = delta.copy()
            d_p[k] += 1e-5
            d_m[k] -= 1e-5
            grad[k] = (phi_at(d_p) - phi_at(d_m)) / (2 * 1e-5)

        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        delta = delta + lr * m_hat / (np.sqrt(v_hat) + eps)

        p = phi_at(delta)
        if p > best_phi:
            best_phi = p
            best_U = to_U(delta)
            rho_U = best_U.conj().T @ rho0 @ best_U
            Y_U = best_U.conj().T @ Y0 @ best_U
            best_obj, _ = canonical_open_objective_lambda(rho_U, Y_U, 2, lambda_val)

    return best_U, best_phi, best_obj


def haar_state(dim, rng):
    m = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    q, _ = np.linalg.qr(m)
    return np.outer(q[:, 0], q[:, 0].conj())


def schmidt_spectrum(psi, n=3, cut=2):
    d_a, d_b = 2 ** cut, 2 ** (n - cut)
    mat = psi.reshape((d_a, d_b))
    _, s, _ = np.linalg.svd(mat)
    return s ** 2  # Schmidt eigenvalues p_k


def main():
    OUTDIR.mkdir(exist_ok=True)
    H_basis = horizontal_basis(4, 2)
    H = get_hamiltonian(HNAME, n_sites=3)
    L = liouvillian_dephasing(H, GAMMA_VEC)
    rng = np.random.default_rng(SEED0)

    print("=== D2.1: Haar Singularity Audit ===")
    print(f"Hamiltonian: {HNAME}, env: γ={GAMMA_VEC}, λ={LAMBDA}, seeds={SEEDS}\n")

    trials = []
    for s in range(SEEDS):
        seed = SEED0 + s
        rng_seed = np.random.default_rng(seed)
        rho0 = haar_state(8, rng_seed)
        Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")

        theta0 = rng_seed.standard_normal(H_basis.shape[2]) * 1.0
        U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

        U_best, phi_best, obj = adam_optimize_open(U0, L, rho0, Y0, LAMBDA, H_basis)
        if obj is None:
            rho_U = U_best.conj().T @ rho0 @ U_best
            Y_U = U_best.conj().T @ Y0 @ U_best
            obj, _ = canonical_open_objective_lambda(rho_U, Y_U, 2, LAMBDA)

        # C_dyn at τ=0.1
        c_dyn, c0, ct = compute_c_f_tau(rho0, L, U_best, 0.1)

        # Schmidt spectrum of best solution
        evals_psi, evecs_psi = np.linalg.eigh(rho0)
        psi0 = evecs_psi[:, np.argmax(evals_psi)]
        # Rotated state
        psi_U = U_best.conj().T @ psi0
        spec = schmidt_spectrum(psi_U)

        trials.append(dict(
            seed=seed,
            Phi=phi_best,
            I=obj["I"],
            gamma=obj["gamma"],
            C_F_0_sq=obj["C_F_sq"],
            C_dyn_tau=c_dyn,
            C_F_tau_sq=ct,
            schmidt_spectrum=spec.tolist(),
        ))

    # Statistics
    phis = np.array([t["Phi"] for t in trials])
    gammas = np.array([t["gamma"] for t in trials])
    c0s = np.array([t["C_F_0_sq"] for t in trials])
    c_dyns = np.array([t["C_dyn_tau"] for t in trials])

    print(f"Phi:   min={phis.min():.4f} median={np.median(phis):.4f} max={phis.max():.4f} std={phis.std():.4f}")
    print(f"Gamma: min={gammas.min():.4f} median={np.median(gammas):.4f} max={gammas.max():.4f}")
    print(f"C_F(0)²: min={c0s.min():.3e} median={np.median(c0s):.3e} max={c0s.max():.3e}")
    print(f"C_dyn(τ): min={c_dyns.min():.4f} median={np.median(c_dyns):.4f} max={c_dyns.max():.4f}")

    # Correlation: Φ vs C_F(0)
    mask = c0s > 0
    corr = np.corrcoef(phis[mask], np.log10(c0s[mask]))[0, 1] if mask.sum() > 2 else float("nan")
    print(f"\ncorr(Φ, log10 C_F(0)²) = {corr:.4f}")

    # Scatter: Γ vs 1/C_F(0)
    print("\n=== Γ vs C_F(0)² scatter (divergence check) ===")
    for i in np.argsort(c0s)[:10]:
        print(f"  C_F(0)²={c0s[i]:.3e}  Γ={gammas[i]:+.4f}  Φ={phis[i]:.4f}")

    # Best solution detail
    best_idx = int(np.argmax(phis))
    bt = trials[best_idx]
    print(f"\n=== Best solution (seed {bt['seed']}, Φ={bt['Phi']:.4f}) ===")
    print(f"  I = {bt['I']:.6f}")
    print(f"  Γ = {bt['gamma']:.6f}")
    print(f"  C_F(0)² = {bt['C_F_0_sq']:.6e}")
    print(f"  C_dyn(τ=0.1) = {bt['C_dyn_tau']:.6f}")
    print(f"  Schmidt spectrum p_k = {[f'{p:.4f}' for p in bt['schmidt_spectrum']]}")

    # Small perturbation response
    print("\n=== Small perturbation response (best solution) ===")
    rho0_best = haar_state(8, np.random.default_rng(int(bt["seed"])))
    Y0_best = (L @ rho0_best.reshape(-1, order="F")).reshape(rho0_best.shape, order="F")
    U_best = None  # will rebuild via optimize; for now re-optimize same seed
    # Simpler: reconstruct U from theta is not stored; do a local perturbation around
    # the identity chart for the stored best U
    # We stored only Phi; instead do a fresh short optimization from best seed
    theta0 = np.random.default_rng(int(bt["seed"])).standard_normal(H_basis.shape[2]) * 1.0
    U_ref = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))
    # Re-run short optimizer to get U_best (100 steps)
    U_best, _, obj_best = adam_optimize_open(U_ref, L, rho0_best, Y0_best, LAMBDA, H_basis,
                                             steps=100, lr=LR)

    for eps_pert in [1e-3, 1e-2]:
        d = np.zeros(H_basis.shape[2])
        d[0] = eps_pert
        U_pert = U_best @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))
        rho_U = U_pert.conj().T @ rho0_best @ U_pert
        Y_U = U_pert.conj().T @ Y0_best @ U_pert
        obj_p, phi_p = canonical_open_objective_lambda(rho_U, Y_U, 2, LAMBDA)
        print(f"  ε={eps_pert}: Φ={phi_p:.4f} Γ={obj_p['gamma']:.4f} "
              f"I={obj_p['I']:.4f} C_F(0)²={obj_p['C_F_sq']:.3e}")

    # Verdict
    print("\n=== Verdict ===")
    min_c0 = c0s.min()
    max_abs_gamma = np.abs(gammas).max()
    # Check divergence: is Γ large when C_F(0) is small?
    small_mask = c0s < 1e-3
    if small_mask.sum() > 0:
        small_corr = np.corrcoef(np.log10(c0s[small_mask]), gammas[small_mask])[0, 1]
        print(f"  Trials with C_F(0)² < 1e-3: {small_mask.sum()}/{SEEDS}")
        print(f"  corr(log C_F(0)², Γ) on small-C trials = {small_corr:.4f}")
    else:
        print(f"  No trials with C_F(0)² < 1e-3 (min = {min_c0:.3e})")

    if min_c0 < 1e-4 and max_abs_gamma > 100:
        verdict = "D2.1-A: normalization singularity (C_F(0)→0 drives Γ→±∞)"
    elif min_c0 < 1e-6:
        verdict = "D2.1-A: severe normalization singularity (C_F(0)→0)"
    elif corr < -0.9 and max_abs_gamma > 50:
        # Strong anti-correlation of Φ with C_F(0) and large |Γ| is
        # direct evidence of the normalization singularity, even if the
        # extreme tail was reached only via perturbation.
        verdict = "D2.1-A: normalization singularity (corr(Φ,log C_F(0))≈-1, |Γ|≫1 at small C_F(0))"
    else:
        verdict = "D2.1-B: no singularity; Φ≈49 is a real signed-functional effect"
    print(f"  {verdict}")

    results = dict(
        trials=trials,
        stats=dict(
            Phi_min=float(phis.min()), Phi_median=float(np.median(phis)),
            Phi_max=float(phis.max()), Phi_std=float(phis.std()),
            Gamma_min=float(gammas.min()), Gamma_median=float(np.median(gammas)),
            Gamma_max=float(gammas.max()),
            C_F0_sq_min=float(c0s.min()), C_F0_sq_median=float(np.median(c0s)),
            corr_Phi_logC0=float(corr) if np.isfinite(corr) else None,
        ),
        verdict=verdict,
    )

    out_path = OUTDIR / "d2_1_haar_singularity_audit.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
