#!/usr/bin/env python3
"""
Phase C3.3a — Matched-seed λ-sensitivity control.

Purpose: isolate the effect of λ on factorization selection by fixing
multimodality (same initial U_0, same seed, same optimizer) and comparing
F*(λ=0.2) vs F*(λ=λ_predicted) for test environments.

Protocol:
  - Test environments: A2_deph_2105, amp_damping, mixed
  - λ candidates: τ_Γ, τ_L, τ_γ (from C3.1, averaged over training envs)
  - For each environment and each λ_candidate:
    - Run 6 matched seeds with same initial U_0
    - Compare F*(λ=0.2) vs F*(λ=λ_candidate)
    - Record: Φ*, I(F*), Γ(F*), d_F, final basin, cross-evaluation

Data:
  - c3_1_lambda_identifiability.json: τ_env per environment
  - c3_3_holdout_prediction.json: λ_predicted per candidate
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
from b2_canonical_open import canonical_open_objective, apply_tps, \
    liouvillian_dephasing, liouvillian_amp_damping, liouvillian_mixed
from b2_1_cross_eval import tps_distance

TEST_ENVS = ["A2_deph_2105", "amp_damping", "mixed"]
CANDIDATES = ["τ_Γ", "τ_L", "τ_γ"]
N_SEEDS = 6
SEED_OFFSET = 20260812  # Different from B2/B3 seeds


def load_c31_results() -> dict:
    """Load τ_env from C3.1."""
    c31_path = OUTDIR / "c3_1_lambda_identifiability.json"
    if not c31_path.exists():
        raise FileNotFoundError(f"Missing {c31_path}")
    return json.loads(c31_path.read_text(encoding="utf-8"))


def load_c33_results() -> dict:
    """Load λ_predicted from C3.3."""
    c33_path = OUTDIR / "c3_3_holdout_prediction.json"
    if not c33_path.exists():
        raise FileNotFoundError(f"Missing {c33_path}")
    return json.loads(c33_path.read_text(encoding="utf-8"))


def adam_optimize(U0: np.ndarray, L: np.ndarray, psi0: np.ndarray,
                  lambda_val: float, steps: int = 200, lr: float = 0.01) -> tuple:
    """Adam optimization in the right-multiplication chart with explicit λ.

    Returns (best_U, best_Phi_at_lambda).
    """
    H_basis = horizontal_basis(4, 2)
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_at(d):
        U = to_U(d)
        rho_U, Y_U = apply_tps(psi0, L, U)
        obj = canonical_open_objective(rho_U, Y_U, 2)
        return obj["I"] - lambda_val * obj["gamma"] if np.isfinite(obj["gamma"]) else -1e9

    best_phi = phi_at(delta)
    best_U = to_U(delta)

    for step in range(steps):
        # finite-difference gradient in delta-chart
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

    return best_U, best_phi


def main():
    c31 = load_c31_results()
    c33 = load_c33_results()

    # Compute λ_predicted for each candidate (average over training envs)
    train_envs = ["uniform_dephasing", "A1_deph_0512"]
    lambda_predicted = {}
    for cand in CANDIDATES:
        tau_vals = []
        for env in train_envs:
            if env in c31:
                tau_key = {"τ_Γ": "tau_Gamma", "τ_L": "tau_L", "τ_γ": "tau_gamma"}[cand]
                tau_val = c31[env].get(tau_key)
                if tau_val is not None and np.isfinite(tau_val):
                    tau_vals.append(tau_val)
        if tau_vals:
            lambda_predicted[cand] = np.mean(tau_vals)

    print("=== C3.3a: Matched-seed λ-sensitivity control ===\n")
    print(f"Test environments: {TEST_ENVS}")
    print(f"λ candidates: {CANDIDATES}")
    print(f"λ_predicted: {lambda_predicted}\n")

    # Setup
    h_total, _ = xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi0 = ground_state(h_total)
    H_basis = horizontal_basis(4, 2)

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
        env_results = {}

        for cand in CANDIDATES:
            if cand not in lambda_predicted:
                print(f"  SKIP {cand}: no λ_predicted")
                continue

            lambda_pred = lambda_predicted[cand]
            print(f"\n  Candidate: {cand} (λ_predicted = {lambda_pred:.4f})")

            cand_results = {
                "lambda_predicted": lambda_pred,
                "seeds": [],
            }

            for seed_idx in range(N_SEEDS):
                seed = SEED_OFFSET + seed_idx
                rng = np.random.default_rng(seed)

                # Generate initial U_0 (same for both λ=0.2 and λ=λ_pred)
                theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
                U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

                # Optimize with λ=0.2
                U_opt_02, phi_02 = adam_optimize(U0, L, psi0, lambda_val=0.2, steps=200, lr=0.01)
                rho_02, Y_02 = apply_tps(psi0, L, U_opt_02)
                obj_02 = canonical_open_objective(rho_02, Y_02, 2)

                # Optimize with λ=λ_pred
                U_opt_pred, phi_pred = adam_optimize(U0, L, psi0, lambda_val=lambda_pred, steps=200, lr=0.01)
                rho_pred, Y_pred = apply_tps(psi0, L, U_opt_pred)
                obj_pred = canonical_open_objective(rho_pred, Y_pred, 2)

                # Compute gauge-invariant distance
                d_F = tps_distance(U_opt_02, U_opt_pred)

                # Phi values at their own lambda (canonical_open_objective
                # returns Phi at the module LAMBDA=0.2, so recompute here)
                phi_02_report = obj_02["I"] - 0.2 * obj_02["gamma"]
                phi_pred_report = obj_pred["I"] - lambda_pred * obj_pred["gamma"]

                seed_result = {
                    "seed": seed,
                    "lambda_02": {
                        "Phi": phi_02_report,
                        "I": obj_02["I"],
                        "gamma": obj_02["gamma"],
                        "U_real": U_opt_02.real.tolist(),
                        "U_imag": U_opt_02.imag.tolist(),
                    },
                    "lambda_pred": {
                        "Phi": phi_pred_report,
                        "I": obj_pred["I"],
                        "gamma": obj_pred["gamma"],
                        "U_real": U_opt_pred.real.tolist(),
                        "U_imag": U_opt_pred.imag.tolist(),
                    },
                    "d_F": d_F,
                }
                cand_results["seeds"].append(seed_result)

                print(f"    Seed {seed_idx}: d_F = {d_F:.4f}, "
                      f"Φ(0.2)={phi_02_report:.4f}, Φ(pred)={phi_pred_report:.4f}")

            env_results[cand] = cand_results

        results[env] = env_results

    # Summary
    print("\n\n=== Summary ===")
    for env in TEST_ENVS:
        if env not in results:
            continue
        print(f"\n{env}:")
        for cand in CANDIDATES:
            if cand not in results[env]:
                continue
            d_Fs = [s["d_F"] for s in results[env][cand]["seeds"]]
            avg_d_F = np.mean(d_Fs)
            std_d_F = np.std(d_Fs)
            print(f"  {cand}: avg d_F = {avg_d_F:.4f} ± {std_d_F:.4f}")
            if avg_d_F < 0.1:
                print(f"    → λ has NO factorization-selection power (same basin)")
            else:
                print(f"    → λ has factorization-selection power (different basins)")

    # Save results
    out_path = OUTDIR / "c3_3a_matched_seed_control.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
