#!/usr/bin/env python3
"""
Compare Adam with vs without early stopping under identical conditions.

Same seed, same Hamiltonian, same generator basis, same initial θ.
Early-stopping thresholds: tol_grad=0.1, tol_phi=1e-4, patience=5.

Measures: final Φ, converged-at step, step reduction.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).parents[1] / "src" / "rem4_numerical.py"
spec = importlib.util.spec_from_file_location("rem4_numerical", MODULE_PATH)
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)

import scipy.linalg as _sla
rem4.sla = _sla


def run_trial(
    n: int,
    couplings: list[float],
    h: float,
    lambda_value: float,
    *,
    n_a: int,
    n_params: int,
    max_steps: int,
    lr: float,
    tol_grad: float,
    tol_phi: float,
    patience: int,
    seed: int,
) -> dict:
    """Run Adam with and without early stopping from identical starting point."""
    h_total, bonds = rem4.xy_chain_hamiltonian(n, couplings, h)
    _, psi = rem4.ground_state(h_total)

    # Fix generator basis and initial theta — replicate exactly
    rem4.N_RNG = np.random.default_rng(seed=seed)
    gen_shared = rem4._random_anti_hermitian(2**n, n_params)
    theta0_shared = np.random.default_rng(seed=seed).standard_normal(n_params) * 0.1

    def _adam_run(early_stop: bool) -> dict:
        theta = theta0_shared.copy()
        generators = gen_shared.copy()
        m = np.zeros_like(theta)
        v = np.zeros_like(theta)
        beta1, beta2, eps = 0.9, 0.999, 1e-8
        best_phi = -1e9
        converged_at = max_steps
        phi_window = [] if early_stop else None

        for step in range(max_steps):
            t = step + 1
            u = rem4.sla.expm(1j * np.einsum("ijk,k->ij", generators, theta))
            phi, mi, c_h = rem4.evaluate_factorization(
                psi, h_total, bonds, u, n=n, n_a=n_a, lambda_value=lambda_value)
            grad = rem4._finite_diff_grad(
                theta, psi, h_total, bonds, generators,
                n=n, n_a=n_a, lambda_value=lambda_value)

            m = beta1 * m + (1.0 - beta1) * grad
            v = beta2 * v + (1.0 - beta2) * (grad ** 2)
            m_hat = m / (1.0 - beta1 ** t)
            v_hat = v / (1.0 - beta2 ** t)
            theta += lr * m_hat / (np.sqrt(v_hat) + eps)

            if phi > best_phi:
                best_phi = phi

            if early_stop:
                gn = float(np.linalg.norm(grad))
                phi_window.append(phi)
                if len(phi_window) > patience + 1:
                    phi_window.pop(0)
                if len(phi_window) == patience + 1:
                    phi_change = abs(phi_window[-1] - phi_window[0]) / max(1.0, abs(phi_window[0]))
                    if gn < tol_grad and phi_change < tol_phi:
                        converged_at = step + 1
                        break

        return {"best_phi": float(best_phi), "converged_at": converged_at,
                "steps_used": step + 1}

    full = _adam_run(early_stop=False)
    early = _adam_run(early_stop=True)

    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=lambda_value)
    cont_phi = max(m.phi for m in metrics)

    return {
        "seed": seed,
        "n": n,
        "n_a": n_a,
        "contiguous_phi": float(cont_phi),
        "full": full,
        "early": early,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=os.path.join(
        os.path.dirname(__file__), "..", "analysis_output"))
    parser.add_argument("--trials", type=int, default=30)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--lr", type=float, default=0.01)
    parser.add_argument("--tol-grad", type=float, default=0.0,
                        help="0 = auto: 0.05 (N=3), 0.01 (N=4), 0.005 (N=5)")
    parser.add_argument("--tol-phi", type=float, default=1e-4)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--ns", type=int, nargs="+", default=[3])
    parser.add_argument("--j12", type=float, default=1.5)
    parser.add_argument("--j23", type=float, default=0.6)
    parser.add_argument("--h", type=float, default=0.2)
    parser.add_argument("--lambda", dest="lambda_value", type=float, default=0.2)
    parser.add_argument("--seed-offset", type=int, default=2000)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for n in args.ns:
        n_params = {3: 4, 4: 16, 5: 20}.get(n, 4)
        couplings = [args.j12] + [args.j23] * (n - 2)
        tol_grad = args.tol_grad if args.tol_grad > 0 else {3: 0.05, 4: 0.01, 5: 0.005}.get(n, 0.05)

        h_total, bonds = rem4.xy_chain_hamiltonian(n, couplings, args.h)
        _, psi = rem4.ground_state(h_total)
        metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=args.lambda_value)
        best_cut = max(metrics, key=lambda m: m.phi).cut
        n_a = best_cut
        cont_phi = max(m.phi for m in metrics)

        print(f"\n{'='*65}")
        print(f"  N={n}  n_a={n_a}  n_params={n_params}  trials={args.trials}")
        print(f"  Early-stopping: tol_grad={tol_grad}, tol_phi={args.tol_phi}, patience={args.patience}")
        print(f"{'='*65}")

        all_data = []
        start_total = time.time()

        for i in range(args.trials):
            seed = args.seed_offset + i
            trial = run_trial(
                n, couplings, args.h, args.lambda_value,
                n_a=n_a, n_params=n_params,
                max_steps=args.steps, lr=args.lr,
                tol_grad=args.tol_grad, tol_phi=args.tol_phi,
                patience=args.patience, seed=seed,
            )
            all_data.append(trial)

            if (i + 1) % 10 == 0:
                print(f"  [{i+1:3d}/{args.trials}]  "
                      f"full={trial['full']['best_phi']:.4f} "
                      f"early={trial['early']['best_phi']:.4f} "
                      f"({trial['early']['converged_at']:3d} steps)")

        full_phis = np.array([t["full"]["best_phi"] for t in all_data])
        early_phis = np.array([t["early"]["best_phi"] for t in all_data])
        steps_used = np.array([t["early"]["converged_at"] for t in all_data])

        step_save = 1.0 - steps_used.mean() / args.steps
        phi_diff = early_phis - full_phis
        quality_loss = np.mean(phi_diff < -0.01) * 100

        print(f"\n  ┌──────────┬──────────────┬──────────────┬──────────┬──────────┐")
        print(f"  │ Mode     │ Φ (mean±σ)   │ steps (mean) │ save     │ quality  │")
        print(f"  ├──────────┼──────────────┼──────────────┼──────────┼──────────┤")
        print(f"  │ Full     │ {full_phis.mean():.4f}±{full_phis.std():.4f}  │ {args.steps:<13}│ —        │ —        │")
        print(f"  │ Early    │ {early_phis.mean():.4f}±{early_phis.std():.4f}  │ {steps_used.mean():<8.1f}       │ {step_save*100:.0f}%       │ {quality_loss:.0f}% loss  │")
        print(f"  └──────────┴──────────────┴──────────────┴──────────┴──────────┘")
        print(f"  Φ early > Φ full: {(early_phis > full_phis).mean()*100:.0f}%")
        print(f"  Φ early < Φ full (by >0.01): {quality_loss:.0f}%")

        elapsed = time.time() - start_total
        print(f"  Total: {elapsed:.0f}s ({elapsed/args.trials:.1f}s/trial)")

        # Save
        out = {
            "params": vars(args),
            "n": n,
            "n_a": n_a,
            "contiguous_phi": cont_phi,
            "aggregate": {
                "full_phi_mean": float(full_phis.mean()),
                "full_phi_std": float(full_phis.std()),
                "early_phi_mean": float(early_phis.mean()),
                "early_phi_std": float(early_phis.std()),
                "steps_mean": float(steps_used.mean()),
                "steps_std": float(steps_used.std()),
                "steps_save_pct": float(step_save * 100),
                "quality_loss_pct": float(quality_loss),
                "early_better_frac": float((early_phis > full_phis).mean()),
            },
            "trials": [{
                "seed": t["seed"],
                "full_phi": t["full"]["best_phi"],
                "early_phi": t["early"]["best_phi"],
                "early_steps": t["early"]["converged_at"],
            } for t in all_data],
        }
        path = os.path.join(args.outdir, f"early_stop_N{n}.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"  → {path}")


if __name__ == "__main__":
    main()
