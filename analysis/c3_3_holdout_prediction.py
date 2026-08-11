#!/usr/bin/env python3
"""
Phase C3.3 — Hold-out prediction test.

Purpose: verify that λ can be determined from training environments and
applied to test environments WITHOUT re-fitting, predicting F* for unseen
environments.

Protocol:
  1. Training environments: uniform_dephasing, A1_deph_0512
  2. Test environments: A2_deph_2105, amp_damping, mixed
  3. For each candidate rule (λ = c·τ_env with c=1):
     - Compute τ_env for training environments
     - Average to get λ_predicted
     - For test environments:
       a) Predict F*_predicted using λ_predicted (from B3.3 data)
       b) Compare with F*_observed (from B2 direct optimization)
       c) Compute gauge-invariant d_F
  4. If d_F < threshold, the rule generalizes.

Data:
  - c3_1_lambda_identifiability.json: τ_env per environment
  - b2_canonical_open.json: F*_observed per environment (best U*)
  - b3_3_lambda_sweep.json: F*(λ) per environment
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

from b2_1_cross_eval import tps_distance

TRAIN_ENVS = ["uniform_dephasing", "A1_deph_0512"]
TEST_ENVS = ["A2_deph_2105", "amp_damping", "mixed"]
CANDIDATES = [
    ("τ_Γ", "tau_Gamma"),
    ("τ_L", "tau_L"),
    ("τ_γ", "tau_gamma"),
]


def load_b2_results() -> dict:
    """Load F*_observed from B2."""
    b2_path = OUTDIR / "b2_canonical_open.json"
    if not b2_path.exists():
        raise FileNotFoundError(f"Missing {b2_path}")
    return json.loads(b2_path.read_text(encoding="utf-8"))


def load_c31_results() -> dict:
    """Load τ_env from C3.1."""
    c31_path = OUTDIR / "c3_1_lambda_identifiability.json"
    if not c31_path.exists():
        raise FileNotFoundError(f"Missing {c31_path}")
    return json.loads(c31_path.read_text(encoding="utf-8"))


def load_b33_results() -> dict:
    """Load F*(λ) from B3.3."""
    b33_path = OUTDIR / "b3_3_lambda_sweep.json"
    if not b33_path.exists():
        raise FileNotFoundError(f"Missing {b33_path}")
    return json.loads(b33_path.read_text(encoding="utf-8"))


def main():
    b2 = load_b2_results()
    c31 = load_c31_results()
    b33 = load_b33_results()

    print("=== C3.3: Hold-out prediction test ===\n")
    print(f"Training environments: {TRAIN_ENVS}")
    print(f"Test environments: {TEST_ENVS}\n")

    # Step 1: For each candidate rule, compute λ_predicted from training envs
    results = {}
    for display_name, key_name in CANDIDATES:
        print(f"--- Candidate rule: λ = {display_name} (c=1) ---")

        # Compute τ for training environments
        tau_train = []
        for env in TRAIN_ENVS:
            if env not in c31:
                continue
            tau_val = c31[env].get(key_name)
            if tau_val is not None and np.isfinite(tau_val):
                tau_train.append(tau_val)
                print(f"  {env}: {key_name} = {tau_val:.4f}")

        if not tau_train:
            print(f"  SKIP: no valid {key_name} in training envs\n")
            continue

        # Average to get λ_predicted
        lambda_predicted = np.mean(tau_train)
        print(f"  → λ_predicted = mean({tau_train}) = {lambda_predicted:.4f}\n")

        # Step 2: For test environments, compare F*_predicted vs F*_observed
        # Find nearest λ grid point from B3.3
        lambda_grid = [0.05, 0.1, 0.2, 0.3, 0.5]
        lambda_nearest = min(lambda_grid, key=lambda x: abs(x - lambda_predicted))
        lambda_key = str(lambda_nearest)
        print(f"  Nearest B3.3 λ grid point: {lambda_nearest} (key: '{lambda_key}')")

        # For each test environment, get F*_observed (from B2) and F*_predicted
        # (from B3.3 at lambda_nearest)
        test_results = []
        for env in TEST_ENVS:
            if env not in b2 or env not in b33:
                continue
            if lambda_key not in b33[env]:
                print(f"  SKIP {env}: λ={lambda_key} not in B3.3 data")
                continue

            # F*_observed: best U* from B2
            best_b2 = max(b2[env], key=lambda t: t["best_phi"])
            U_observed = np.array(best_b2["best_u_real"]) + 1j * np.array(best_b2["best_u_imag"])
            phi_observed = best_b2["best_phi"]

            # F*_predicted: best U* from B3.3 at lambda_nearest
            b33_data = b33[env][lambda_key]
            U_predicted = np.array(b33_data["best_U_real"]) + 1j * np.array(b33_data["best_U_imag"])
            phi_predicted = b33_data["best_Phi"]

            # Compute gauge-invariant d_F
            d_F = tps_distance(U_observed, U_predicted)

            test_results.append({
                "environment": env,
                "phi_observed": phi_observed,
                "phi_predicted": phi_predicted,
                "delta_phi": phi_observed - phi_predicted,
                "d_F": d_F,
            })

            print(f"  {env}:")
            print(f"    F*_observed: Φ = {phi_observed:.6f}")
            print(f"    F*_predicted (λ={lambda_nearest}): Φ = {phi_predicted:.6f}")
            print(f"    ΔΦ = {phi_observed - phi_predicted:+.6f}")
            print(f"    d_F (gauge-invariant TPS distance) = {d_F:.4f}")

        results[display_name] = {
            "lambda_predicted": lambda_predicted,
            "lambda_nearest": lambda_nearest,
            "test_results": test_results,
        }
        print()

    # Step 3: Verdict
    print("=== Verdict ===")
    for tau_name, res in results.items():
        if not res["test_results"]:
            continue
        avg_delta_phi = np.mean([t["delta_phi"] for t in res["test_results"]])
        avg_d_F = np.mean([t["d_F"] for t in res["test_results"]])
        print(f"  {tau_name}:")
        print(f"    λ_predicted = {res['lambda_predicted']:.4f}")
        print(f"    avg ΔΦ = {avg_delta_phi:+.6f}")
        print(f"    avg d_F = {avg_d_F:.4f}")
        # Since B3 showed F* is insensitive to λ in [0.1, 0.5], we expect
        # small ΔΦ and small d_F for all candidates
        # d_F < 0.1 means "same structural region" (from B3.3 findings)
        if abs(avg_delta_phi) < 0.05 and avg_d_F < 0.1:
            print(f"    → PASS: rule generalizes to test environments")
        else:
            print(f"    → FAIL: rule does not generalize")
        print()

    # Save results
    out_path = OUTDIR / "c3_3_holdout_prediction.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
