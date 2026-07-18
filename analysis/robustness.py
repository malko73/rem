#!/usr/bin/env python3
"""
Robustness analysis v2 — generalised n_a support.

For each N, the script finds the best contiguous cut, sets n_a to match
that cut, and runs continuous optimisation over that factorisation manifold.

Usage:
    python analysis/robustness.py --trials 30 --ns 3 4 5
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

import scipy.linalg as sla


# ── n_params policy ─────────────────────────────────────────────────

def _pick_n_params(n: int) -> int:
    """Heuristic: use larger subspace for bigger systems."""
    if n <= 3:
        return 4
    if n == 4:
        return 16
    if n == 5:
        return 20
    return min(32, 4 * n)


# ── per-trial runner ─────────────────────────────────────────────────

def run_trial(
    n: int,
    couplings: list[float],
    h: float,
    lambda_value: float,
    *,
    n_a: int,
    n_params: int,
    steps: int,
    lr: float,
    gen_seed: int,
    theta_seed: int,
) -> dict:
    """Single optimisation run with independent generator and theta seeds."""
    h_total, bonds = rem4.xy_chain_hamiltonian(n, couplings, h)
    _, psi = rem4.ground_state(h_total)

    # contiguous best for this n_a
    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=lambda_value)
    cut_metrics = next((m for m in metrics if m.cut == n_a), None)
    contiguous_phi = cut_metrics.phi if cut_metrics else max(m.phi for m in metrics)

    # generate generator basis with gen_seed, theta with theta_seed
    rem4.N_RNG = np.random.default_rng(seed=gen_seed)
    generators = rem4._random_anti_hermitian(2**n, n_params)
    theta = np.random.default_rng(seed=theta_seed).standard_normal(n_params) * 0.1

    best_phi = -1e9
    best_u = None
    history = []

    for step in range(steps):
        u = sla.expm(1j * np.einsum("ijk,k->ij", generators, theta))
        phi, mi, c_h = rem4.evaluate_factorization(
            psi, h_total, bonds, u, n=n, n_a=n_a, lambda_value=lambda_value)
        grad = rem4._finite_diff_grad(
            theta, psi, h_total, bonds, generators,
            n=n, n_a=n_a, lambda_value=lambda_value)
        theta += lr * grad
        history.append({
            "phi": float(phi), "mi": float(mi), "c_h": float(c_h)})
        if phi > best_phi:
            best_phi = phi
            best_u = u.copy()

    return {
        "n": n,
        "n_a": n_a,
        "gen_seed": gen_seed,
        "theta_seed": theta_seed,
        "best_phi": float(best_phi),
        "contiguous_phi": float(contiguous_phi),
        "improvement": float(best_phi - contiguous_phi),
        "improvement_ratio": float(best_phi / contiguous_phi)
        if contiguous_phi > 1e-12 else float("nan"),
        "success": bool(best_phi >= contiguous_phi - 0.01),
        "history": history,
    }


def aggregate_trials(trials: list[dict]) -> dict:
    """Summarise a list of trial results."""
    phis = np.array([t["best_phi"] for t in trials])
    conts = np.array([t["contiguous_phi"] for t in trials])
    imps = np.array([t["improvement"] for t in trials])
    ratios = np.array([t["improvement_ratio"] for t in trials])
    successes = np.array([t["success"] for t in trials])

    return {
        "n": trials[0]["n"],
        "n_a": trials[0]["n_a"],
        "n_trials": len(trials),
        "contiguous_phi": float(conts[0]),
        "phi_mean": float(np.mean(phis)),
        "phi_std": float(np.std(phis)),
        "phi_min": float(np.min(phis)),
        "phi_max": float(np.max(phis)),
        "phi_median": float(np.median(phis)),
        "improvement_mean": float(np.mean(imps)),
        "improvement_std": float(np.std(imps)),
        "improvement_ratio_mean": float(np.mean(ratios)),
        "success_rate": float(np.mean(successes)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=os.path.join(
        os.path.dirname(__file__), "..", "analysis_output"))
    parser.add_argument("--trials", type=int, default=10)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--lr", type=float, default=0.01)
    parser.add_argument("--n-params", type=int, default=0,
                        help="0 = auto per N")
    parser.add_argument("--ns", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--j12", type=float, default=1.5)
    parser.add_argument("--j23", type=float, default=0.6)
    parser.add_argument("--h", type=float, default=0.2)
    parser.add_argument("--lambda", dest="lambda_value", type=float, default=0.2)
    parser.add_argument("--seed-offset", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    all_results = {}
    all_raw: list[dict] = []

    for n in args.ns:
        n_params = args.n_params if args.n_params > 0 else _pick_n_params(n)
        couplings = [args.j12] + [args.j23] * (n - 2)

        # Determine best contiguous cut to set n_a
        h_total, bonds = rem4.xy_chain_hamiltonian(n, couplings, args.h)
        _, psi = rem4.ground_state(h_total)
        metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=args.lambda_value)
        best_cut = max(metrics, key=lambda m: m.phi).cut
        n_a = best_cut
        contiguous_phi = max(m.phi for m in metrics)

        print(f"\n{'='*60}")
        print(f"  N = {n}  n_a = {n_a}  n_params = {n_params}")
        print(f"  Couplings: J12={args.j12}, others={args.j23}")
        print(f"  Best contiguous-cut Φ = {contiguous_phi:.4f} (cut {best_cut})")
        print(f"{'='*60}")

        # ── Axis A: fixed generator, vary θ ──
        print(f"\n[A] Initial-θ dependence ({args.trials} trials, fixed generator)...")
        start = time.time()
        trials_a = []
        for i in range(args.trials):
            theta_seed = args.seed_offset + i
            tr = run_trial(
                n, couplings, args.h, args.lambda_value,
                n_a=n_a, n_params=n_params,
                steps=args.steps, lr=args.lr,
                gen_seed=0,            # fixed generator basis
                theta_seed=theta_seed,
            )
            trials_a.append(tr)
        agg_a = aggregate_trials(trials_a)
        elapsed = time.time() - start
        all_results[f"N={n}_A"] = agg_a
        all_raw.extend(trials_a)
        print(f"  Φ = {agg_a['phi_mean']:.4f} ± {agg_a['phi_std']:.4f}"
              f"  [{agg_a['phi_min']:.4f}, {agg_a['phi_max']:.4f}]"
              f"  med={agg_a['phi_median']:.4f}")
        print(f"  improvement = {agg_a['improvement_mean']:.4f} ± {agg_a['improvement_std']:.4f}"
              f"  ratio = {agg_a['improvement_ratio_mean']:.4f}x")
        print(f"  success rate = {agg_a['success_rate']*100:.0f}%"
              f"  ({elapsed:.1f}s)")

        # ── Axis B: generator basis dependence ──
        print(f"\n[B] Generator basis dependence ({args.trials} trials)...")
        start = time.time()
        trials_b = []
        for i in range(args.trials):
            gen_seed = args.seed_offset + i
            tr = run_trial(
                n, couplings, args.h, args.lambda_value,
                n_a=n_a, n_params=n_params,
                steps=args.steps, lr=args.lr,
                gen_seed=gen_seed,
                theta_seed=gen_seed,   # one θ per generator
            )
            trials_b.append(tr)
        agg_b = aggregate_trials(trials_b)
        elapsed = time.time() - start
        all_results[f"N={n}_B"] = agg_b
        all_raw.extend(trials_b)
        print(f"  Φ = {agg_b['phi_mean']:.4f} ± {agg_b['phi_std']:.4f}"
              f"  [{agg_b['phi_min']:.4f}, {agg_b['phi_max']:.4f}]"
              f"  med={agg_b['phi_median']:.4f}")
        print(f"  improvement = {agg_b['improvement_mean']:.4f} ± {agg_b['improvement_std']:.4f}"
              f"  ratio = {agg_b['improvement_ratio_mean']:.4f}x")
        print(f"  success rate = {agg_b['success_rate']*100:.0f}%"
              f"  ({elapsed:.1f}s)")

    # ── Save ──
    json_path = os.path.join(args.outdir, "robustness_results.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nAggregated → {json_path}")

    raw_path = os.path.join(args.outdir, "robustness_raw.json")
    # Strip history from raw to keep file size manageable
    raw_compact = [{k: v for k, v in t.items() if k != "history"} for t in all_raw]
    with open(raw_path, "w") as f:
        json.dump(raw_compact, f, indent=2)
    print(f"Raw data    → {raw_path}")


if __name__ == "__main__":
    main()
