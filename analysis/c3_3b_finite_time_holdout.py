#!/usr/bin/env python3
"""
Phase C3.3b — Finite-time hold-out prediction.

Purpose: test whether the factorization predicted from the initial rate
Γ_F^(0) matches the factorization selected by actual finite-time open
dynamics.

Protocol:
  For each test environment (A2, amp, mixed):
    1. Determine τ_env from the environment only (C3.1 results)
    2. Set λ_pred = τ_env
    3. Prediction side:
         F_pred = argmax_F [I(F) - λ_pred·Γ_F^(0)]
    4. Independent evaluation side:
         F_finite = argmax_F [I(F) - λ_pred·C_dyn^(τ_env)(F)]
       where C_dyn^(τ)(F) = -(1/τ) log(C_F(τ)/C_F(0))
             C_F(t) = ||ρ(t) - D_F[ρ(t)]||_F^2
             ρ(t) = exp(tL)[ρ_0]
             D_F = dephasing projection onto Schmidt basis of F
    5. Compute d_F(F_pred, F_finite)

If d_F ≈ 0, the initial-rate prediction matches finite-time dynamics.
If d_F ≫ 0, the prediction fails.

This is a non-trivial test: it checks whether the structure predicted from
the environment-determined timescale and initial rate is actually selected
by the real open-system evolution.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from rem4_numerical import xy_chain_hamiltonian, ground_state
from quotient_geometry import horizontal_basis
from b2_canonical_open import (
    canonical_open_objective,
    apply_tps,
    liouvillian_dephasing,
    liouvillian_amp_damping,
    liouvillian_mixed,
)
from b2_1_cross_eval import tps_distance
from gamma_F import schmidt_basis, dephasing_projection

TEST_ENVS = ["A2_deph_2105", "amp_damping", "mixed"]
N_SEEDS = 6
SEED_OFFSET = 20260813


def load_c31_results() -> dict:
    """Load τ_env from C3.1."""
    c31_path = OUTDIR / "c3_1_lambda_identifiability.json"
    if not c31_path.exists():
        raise FileNotFoundError(f"Missing {c31_path}")
    return json.loads(c31_path.read_text(encoding="utf-8"))


def compute_c_finite_time(U: np.ndarray, L: np.ndarray, psi0: np.ndarray,
                          tau: float, n_a: int) -> float:
    """
    Compute C_dyn^(τ)(F) = -(1/τ) log(C_F(τ)/C_F(0)).

    C_F(t) = ||ρ(t) - D_F[ρ(t)]||_F^2
    ρ(t) = exp(tL)[ρ_0]
    D_F = dephasing projection in the F frame (computational basis).
    """
    # Initial state
    rho0 = np.outer(psi0, psi0.conj())

    # Time evolution: ρ(t) = exp(tL)[ρ_0]
    rho0_vec = rho0.reshape(-1, order="F")
    exp_tL = sla.expm(tau * L)
    rho_t_vec = exp_tL @ rho0_vec
    rho_t = rho_t_vec.reshape(rho0.shape, order="F")

    # Rotate to F frame and dephase
    rho_t_F = U.conj().T @ rho_t @ U
    rho_t_F_dephased = np.diag(np.diag(rho_t_F))
    rho_t_dephased = U @ rho_t_F_dephased @ U.conj().T

    # Coherence norm: C_F(t) = ||ρ(t) - D_F[ρ(t)]||_F^2
    c_t = np.linalg.norm(rho_t - rho_t_dephased, ord="fro") ** 2

    # Initial coherence C_F(0)
    rho0_F = U.conj().T @ rho0 @ U
    rho0_F_dephased = np.diag(np.diag(rho0_F))
    rho0_dephased = U @ rho0_F_dephased @ U.conj().T
    c_0 = np.linalg.norm(rho0 - rho0_dephased, ord="fro") ** 2

    if c_0 < 1e-300:
        return float("nan")

    # C_dyn^(τ)(F) = -(1/τ) log(C_F(τ)/C_F(0))
    if c_t < 1e-300:
        return float("inf")

    c_dyn = -(1.0 / tau) * np.log(c_t / c_0)
    return float(c_dyn)


def adam_optimize_finite_time(U0: np.ndarray, L: np.ndarray, psi0: np.ndarray,
                               lambda_val: float, tau: float,
                               steps: int = 200, lr: float = 0.01) -> tuple:
    """
    Adam optimization for the finite-time objective:
      Φ_finite(F) = I(F) - λ·C_dyn^(τ)(F)

    Returns (best_U, best_Phi_finite).
    """
    H_basis = horizontal_basis(4, 2)
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_finite_at(d):
        U = to_U(d)
        rho_U, Y_U = apply_tps(psi0, L, U)
        obj = canonical_open_objective(rho_U, Y_U, 2)
        I_F = obj["I"]
        c_dyn = compute_c_finite_time(U, L, psi0, tau, 2)
        if not np.isfinite(c_dyn):
            return -1e9
        return I_F - lambda_val * c_dyn

    best_phi = phi_finite_at(delta)
    best_U = to_U(delta)

    for step in range(steps):
        # finite-difference gradient in delta-chart
        grad = np.zeros_like(delta)
        for k in range(len(delta)):
            d_p = delta.copy()
            d_m = delta.copy()
            d_p[k] += 1e-5
            d_m[k] -= 1e-5
            grad[k] = (phi_finite_at(d_p) - phi_finite_at(d_m)) / (2 * 1e-5)

        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        delta = delta + lr * m_hat / (np.sqrt(v_hat) + eps)

        p = phi_finite_at(delta)
        if p > best_phi:
            best_phi = p
            best_U = to_U(delta)

    return best_U, best_phi


def adam_optimize_initial_rate(U0: np.ndarray, L: np.ndarray, psi0: np.ndarray,
                                lambda_val: float, steps: int = 200, lr: float = 0.01) -> tuple:
    """
    Adam optimization for the initial-rate objective:
      Φ_initial(F) = I(F) - λ·Γ_F^(0)

    Returns (best_U, best_Phi_initial).
    """
    H_basis = horizontal_basis(4, 2)
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_initial_at(d):
        U = to_U(d)
        rho_U, Y_U = apply_tps(psi0, L, U)
        obj = canonical_open_objective(rho_U, Y_U, 2)
        I_F = obj["I"]
        gamma = obj["gamma"]
        if not np.isfinite(gamma):
            return -1e9
        return I_F - lambda_val * gamma

    best_phi = phi_initial_at(delta)
    best_U = to_U(delta)

    for step in range(steps):
        # finite-difference gradient in delta-chart
        grad = np.zeros_like(delta)
        for k in range(len(delta)):
            d_p = delta.copy()
            d_m = delta.copy()
            d_p[k] += 1e-5
            d_m[k] -= 1e-5
            grad[k] = (phi_initial_at(d_p) - phi_initial_at(d_m)) / (2 * 1e-5)

        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        delta = delta + lr * m_hat / (np.sqrt(v_hat) + eps)

        p = phi_initial_at(delta)
        if p > best_phi:
            best_phi = p
            best_U = to_U(delta)

    return best_U, best_phi


def main():
    c31 = load_c31_results()

    print("=== C3.3b: Finite-time hold-out prediction ===\n")
    print(f"Test environments: {TEST_ENVS}")

    # Setup
    h_total, _ = xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi0 = ground_state(h_total)

    # Environment setup
    env_configs = {
        "A2_deph_2105": liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
        "amp_damping": liouvillian_amp_damping(h_total, [1.0, 1.0, 1.0]),
        "mixed": liouvillian_mixed(h_total, [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),
    }

    results = {}

    for env in TEST_ENVS:
        print(f"\n--- Environment: {env} ---")
        L = env_configs[env]

        # Determine τ_env from environment only (use τ_L as the most stable candidate)
        tau_env = c31[env]["tau_L"]
        lambda_pred = tau_env
        print(f"  τ_env (τ_L) = {tau_env:.4f}")
        print(f"  λ_pred = {lambda_pred:.4f}\n")

        env_results = {
            "tau_env": tau_env,
            "lambda_pred": lambda_pred,
            "seeds": [],
        }

        for seed_idx in range(N_SEEDS):
            seed = SEED_OFFSET + seed_idx
            rng = np.random.default_rng(seed)

            # Generate initial U_0 (same for both predictions)
            theta0 = rng.standard_normal(45) * 1.0
            H_basis = horizontal_basis(4, 2)
            U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

            # Prediction side: F_pred from initial rate Γ_F^(0)
            U_pred, phi_pred = adam_optimize_initial_rate(
                U0, L, psi0, lambda_val=lambda_pred, steps=200, lr=0.01
            )

            # Independent evaluation side: F_finite from finite-time dynamics
            U_finite, phi_finite = adam_optimize_finite_time(
                U0, L, psi0, lambda_val=lambda_pred, tau=tau_env, steps=200, lr=0.01
            )

            # Compute gauge-invariant distance
            d_F = tps_distance(U_pred, U_finite)

            seed_result = {
                "seed": seed,
                "tau_env": tau_env,
                "lambda_pred": lambda_pred,
                "F_pred": {
                    "Phi": phi_pred,
                    "U_real": U_pred.real.tolist(),
                    "U_imag": U_pred.imag.tolist(),
                },
                "F_finite": {
                    "Phi": phi_finite,
                    "U_real": U_finite.real.tolist(),
                    "U_imag": U_finite.imag.tolist(),
                },
                "d_F": d_F,
            }
            env_results["seeds"].append(seed_result)

            print(f"  Seed {seed_idx}: d_F = {d_F:.4f}, "
                  f"Φ_pred={phi_pred:.4f}, Φ_finite={phi_finite:.4f}")

        results[env] = env_results

    # Summary
    print("\n\n=== Summary ===")
    for env in TEST_ENVS:
        if env not in results:
            continue
        d_Fs = [s["d_F"] for s in results[env]["seeds"]]
        avg_d_F = np.mean(d_Fs)
        std_d_F = np.std(d_Fs)
        print(f"\n{env}:")
        print(f"  τ_env = {results[env]['tau_env']:.4f}")
        print(f"  avg d_F(F_pred, F_finite) = {avg_d_F:.4f} ± {std_d_F:.4f}")
        if avg_d_F < 0.1:
            print(f"  → PASS: initial-rate prediction matches finite-time dynamics")
        else:
            print(f"  → FAIL: prediction does not match dynamics")

    # Save results
    out_path = OUTDIR / "c3_3b_finite_time_holdout.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
