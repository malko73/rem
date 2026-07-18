#!/usr/bin/env python3
"""
Compare SGD vs Adam optimisers under identical conditions.

For each seed, both optimisers get the same:
  - generator basis
  - initial theta
  - Hamiltonian / ground state
  - evaluation budget (steps=200)

Output: aggregated table + per-seed convergence curves + violin comparison.
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
rem4.sla = _sla  # use the module's own scipy.linalg reference


def run_both(
    n: int,
    couplings: list[float],
    h: float,
    lambda_value: float,
    *,
    n_a: int,
    n_params: int,
    steps: int,
    lr: float,
    seed: int,
) -> dict:
    """Run both SGD and Adam from identical starting conditions."""
    h_total, bonds = rem4.xy_chain_hamiltonian(n, couplings, h)
    _, psi = rem4.ground_state(h_total)

    # Fix generator basis and initial theta
    rem4.N_RNG = np.random.default_rng(seed=seed)
    generators = rem4._random_anti_hermitian(2**n, n_params)
    theta0 = np.random.default_rng(seed=seed).standard_normal(n_params) * 0.1

    beta1, beta2, eps_adam = 0.9, 0.999, 1e-8

    def _run(use_adam: bool) -> dict:
        theta = theta0.copy()
        m = np.zeros_like(theta) if use_adam else None
        v = np.zeros_like(theta) if use_adam else None
        best_phi = -1e9
        best_u = None
        history = []

        for step in range(steps):
            t = step + 1
            u = rem4.sla.expm(1j * np.einsum("ijk,k->ij", generators, theta))
            phi, mi, c_h = rem4.evaluate_factorization(
                psi, h_total, bonds, u, n=n, n_a=n_a, lambda_value=lambda_value)
            grad = rem4._finite_diff_grad(
                theta, psi, h_total, bonds, generators,
                n=n, n_a=n_a, lambda_value=lambda_value)

            if use_adam:
                m = beta1 * m + (1.0 - beta1) * grad
                v = beta2 * v + (1.0 - beta2) * (grad ** 2)
                m_hat = m / (1.0 - beta1 ** t)
                v_hat = v / (1.0 - beta2 ** t)
                theta += lr * m_hat / (np.sqrt(v_hat) + eps_adam)
            else:
                theta += lr * grad

            gn = float(np.linalg.norm(grad))
            # Lightweight history: only store every 5th step + max-phi step
            if step % 5 == 0 or phi > best_phi:
                history.append({"step": step, "phi": float(phi),
                                "mi": float(mi), "c_h": float(c_h),
                                "gn": gn})
            if phi > best_phi:
                best_phi = phi
                best_u = u.copy()

        return {"best_phi": float(best_phi), "history": history}

    sgd_result = _run(use_adam=False)
    adam_result = _run(use_adam=True)

    # contiguous best
    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=lambda_value)
    cont_phi = max(m.phi for m in metrics)

    return {
        "seed": seed,
        "n": n,
        "n_a": n_a,
        "contiguous_phi": float(cont_phi),
        "sgd": sgd_result,
        "adam": adam_result,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=os.path.join(
        os.path.dirname(__file__), "..", "analysis_output"))
    parser.add_argument("--trials", type=int, default=30)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--lr", type=float, default=0.01)
    parser.add_argument("--n-params", type=int, default=0)
    parser.add_argument("--ns", type=int, nargs="+", default=[3])
    parser.add_argument("--j12", type=float, default=1.5)
    parser.add_argument("--j23", type=float, default=0.6)
    parser.add_argument("--h", type=float, default=0.2)
    parser.add_argument("--lambda", dest="lambda_value", type=float, default=0.2)
    parser.add_argument("--seed-offset", type=int, default=1000)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for n in args.ns:
        n_params = args.n_params if args.n_params > 0 else {
            3: 4, 4: 16, 5: 20
        }.get(n, 4)
        couplings = [args.j12] + [args.j23] * (n - 2)

        # Determine n_a
        h_total, bonds = rem4.xy_chain_hamiltonian(n, couplings, args.h)
        _, psi = rem4.ground_state(h_total)
        metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=args.lambda_value)
        best_cut = max(metrics, key=lambda m: m.phi).cut
        n_a = best_cut
        cont_phi = max(m.phi for m in metrics)

        print(f"\n{'='*65}")
        print(f"  N={n}  n_a={n_a}  n_params={n_params}  trials={args.trials}")
        print(f"  Contiguous Φ = {cont_phi:.4f}")
        print(f"{'='*65}")

        all_data = []
        start_total = time.time()

        for i in range(args.trials):
            seed = args.seed_offset + i
            trial = run_both(
                n, couplings, args.h, args.lambda_value,
                n_a=n_a, n_params=n_params,
                steps=args.steps, lr=args.lr,
                seed=seed,
            )
            all_data.append(trial)

            if (i + 1) % 10 == 0:
                elapsed = time.time() - start_total
                print(f"  [{i+1:3d}/{args.trials}] "
                      f"SGD Φ={trial['sgd']['best_phi']:.4f}  "
                      f"Adam Φ={trial['adam']['best_phi']:.4f}  "
                      f"({elapsed:.0f}s)")

        # Aggregate — discard full history to save memory
        def _summary(t):
            return {
                "seed": t["seed"],
                "sgd_phi": t["sgd"]["best_phi"],
                "adam_phi": t["adam"]["best_phi"],
                "contiguous_phi": t["contiguous_phi"],
            }
        summary = [_summary(t) for t in all_data]
        del all_data  # free memory

        sgd_phis = np.array([t["sgd_phi"] for t in summary])
        adam_phis = np.array([t["adam_phi"] for t in summary])
        sgd_success = np.mean(sgd_phis >= cont_phi - 0.01)
        adam_success = np.mean(adam_phis >= cont_phi - 0.01)
        sgd_better = np.mean(sgd_phis > adam_phis)
        adam_better = np.mean(adam_phis > sgd_phis)

        print(f"\n  ┌──────────┬──────────────┬──────────────┬──────────┐")
        print(f"  │ Optimiser│ Φ (mean±σ)   │ Φ range      │ success  │")
        print(f"  ├──────────┼──────────────┼──────────────┼──────────┤")
        print(f"  │ SGD      │ {sgd_phis.mean():.4f}±{sgd_phis.std():.4f}  │ [{sgd_phis.min():.4f},{sgd_phis.max():.4f}] │ {sgd_success*100:.0f}%      │")
        print(f"  │ Adam     │ {adam_phis.mean():.4f}±{adam_phis.std():.4f}  │ [{adam_phis.min():.4f},{adam_phis.max():.4f}] │ {adam_success*100:.0f}%      │")
        print(f"  └──────────┴──────────────┴──────────────┴──────────┘")
        print(f"  SGD > Adam: {sgd_better*100:.0f}%  Adam > SGD: {adam_better*100:.0f}%  tie: {(1-sgd_better-adam_better)*100:.0f}%")
        elapsed = time.time() - start_total
        print(f"  Total: {elapsed:.0f}s ({elapsed/args.trials:.1f}s/trial)")

        # Save — append per N
        full_out = os.path.join(args.outdir, f"optimiser_compare_N{n}_full.json")
        full_data = {
            "params": vars(args),
            "trials": [{
                "seed": t["seed"],
                "sgd_phi": t["sgd_phi"],
                "adam_phi": t["adam_phi"],
                "contiguous_phi": t["contiguous_phi"],
            } for t in summary],
            "aggregate": {
                "sgd_mean": float(sgd_phis.mean()),
                "sgd_std": float(sgd_phis.std()),
                "sgd_min": float(sgd_phis.min()),
                "sgd_max": float(sgd_phis.max()),
                "adam_mean": float(adam_phis.mean()),
                "adam_std": float(adam_phis.std()),
                "adam_min": float(adam_phis.min()),
                "adam_max": float(adam_phis.max()),
                "sgd_success_rate": float(sgd_success),
                "adam_success_rate": float(adam_success),
                "sgd_better_frac": float(sgd_better),
                "adam_better_frac": float(adam_better),
            },
        }
        _save(full_out, full_data)
        # No history in this mode — keep files small


def _save(path: str, data) -> None:
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  → {path}")


if __name__ == "__main__":
    main()
