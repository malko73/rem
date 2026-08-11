#!/usr/bin/env python3
"""
Phase B1 — Closed-system regression on the 45-dim quotient optimizer.

Purpose (per master 2026-08-11): NOT a proof of canonical physics, but a
regression test of what survives in the closed-system surrogate
    Phi_closed = I - lambda_closed * C_H^closed
when the FULL 45-dim quotient search is used (instead of the old random
4-generator subspace).

Outputs (analysis_output/):
    b1_closed_regression.txt / .json / .png

Statistics:
    best / median / std / success rate vs contiguous best /
    distinct minima estimate / convergence distribution
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"

spec = importlib.util.spec_from_file_location("rem4_numerical", REPO / "src" / "rem4_numerical.py")
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)

spec2 = importlib.util.spec_from_file_location("quotient_optimizer", REPO / "src" / "quotient_optimizer.py")
qopt = importlib.util.module_from_spec(spec2)
assert spec2.loader is not None
sys.modules[spec2.name] = qopt
spec2.loader.exec_module(qopt)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 3
COUPLINGS = [1.5, 0.6]
H_FIELD = 0.2
LAMBDA = 0.2
STEPS = 200
LR = 0.01
N_TRIALS = 100
SEED0 = 20260811


def contiguous_baselines(psi, h_total, bonds):
    """Evaluate contiguous cuts at identity for the chosen lambda."""
    I = np.eye(2 ** N, dtype=complex)
    rows = {}
    for n_a in (1, 2):
        phi, mi, ch = rem4.evaluate_factorization(
            psi, h_total, bonds, I, n=N, n_a=n_a, lambda_value=LAMBDA)
        rows[f"cut{n_a}"] = dict(phi=float(phi), mi=float(mi), ch=float(ch))
    return rows


def run_all(n_trials=N_TRIALS, steps=STEPS, verbose=False):
    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    _, psi = rem4.ground_state(h_total)
    baselines = contiguous_baselines(psi, h_total, bonds)
    contiguous_best = max(b["phi"] for b in baselines.values())

    results = {"sgd": [], "adam": []}
    t0 = time.time()
    for trial in range(n_trials):
        seed = SEED0 + trial
        qopt.N_RNG = np.random.default_rng(seed=seed)
        for opt_name, fn in (("sgd", qopt.quotient_optimizer_sgd),
                             ("adam", qopt.quotient_optimizer_adam)):
            r = fn(psi, h_total, bonds, n=N, n_a=2, lambda_value=LAMBDA,
                   steps=steps, lr=LR, verbose=False)
            results[opt_name].append(dict(
                seed=seed,
                best_phi=float(r["best_phi"]),
                final_phi=float(r["history"][-1][0]),
                max_unitarity=float(r["unitarity_errors"].max()),
            ))
        if verbose and (trial + 1) % 20 == 0:
            elapsed = time.time() - t0
            print(f"  trial {trial+1}/{n_trials} ({elapsed:.1f}s)")
    elapsed = time.time() - t0
    return results, baselines, contiguous_best, elapsed


def summarize(results, contiguous_best):
    summary = {}
    for opt_name, runs in results.items():
        phis = np.array([r["best_phi"] for r in runs])
        unit = np.array([r["max_unitarity"] for r in runs])
        # distinct minima estimate: cluster best_phi values in 0.05 bins
        hist, edges = np.histogram(phis, bins=np.arange(np.floor(phis.min() * 20) / 20,
                                                        np.ceil(phis.max() * 20) / 20 + 0.05,
                                                        0.05))
        distinct_minima = int(np.sum(hist > 0))
        summary[opt_name] = dict(
            best=float(phis.max()),
            median=float(np.median(phis)),
            mean=float(phis.mean()),
            std=float(phis.std()),
            min=float(phis.min()),
            success_rate=float(np.mean(phis > contiguous_best)),
            distinct_minima_estimate=distinct_minima,
            max_unitarity=float(unit.max()),
        )
    return summary


def plot(results, baselines, contiguous_best, path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, opt_name in zip(axes, ("sgd", "adam")):
        phis = np.array([r["best_phi"] for r in results[opt_name]])
        ax.hist(phis, bins=24, alpha=0.75, color="#1f77b4" if opt_name == "adam" else "#ff7f0e")
        ax.axvline(contiguous_best, color="gray", ls=":", lw=1.5)
        ax.text(contiguous_best, ax.get_ylim()[1] * 0.95, "contig",
                fontsize=8, ha="right", color="gray")
        ax.set_title(f"{opt_name.upper()}  best={phis.max():.3f}  "
                     f"med={np.median(phis):.3f}  σ={phis.std():.3f}")
        ax.set_xlabel(r"$\Phi_{\mathrm{best}}$")
        ax.set_ylabel("count")
    fig.suptitle("B1 closed-system regression — 45-dim quotient (N=3, λ=0.2)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    n_trials = N_TRIALS
    if "--smoke" in sys.argv:
        n_trials = 5
    verbose = "--verbose" in sys.argv

    results, baselines, contiguous_best, elapsed = run_all(n_trials=n_trials, verbose=verbose)
    summary = summarize(results, contiguous_best)

    out_txt = OUTDIR / "b1_closed_regression.txt"
    out_json = OUTDIR / "b1_closed_regression.json"
    out_png = OUTDIR / "b1_closed_regression.png"
    OUTDIR.mkdir(exist_ok=True)

    lines = []
    lines.append("# Phase B1 — closed-system regression (2026-08-11)")
    lines.append("# 45-dim quotient optimizer, N=3 XY, lambda=0.2, "
                 f"{n_trials} trials x {STEPS} steps")
    lines.append(f"# contiguous baselines: " + ", ".join(
        f"{k}: Φ={v['phi']:.4f}" for k, v in baselines.items()))
    lines.append(f"# contiguous best: {contiguous_best:.4f}")
    lines.append(f"# wall time: {elapsed:.1f}s")
    lines.append("")
    lines.append("| optimiser | best | median | mean | std | min | success | distinct minima | max |unit|")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for opt_name in ("sgd", "adam"):
        s = summary[opt_name]
        lines.append(
            f"| {opt_name} | {s['best']:.4f} | {s['median']:.4f} | {s['mean']:.4f} | "
            f"{s['std']:.4f} | {s['min']:.4f} | {s['success_rate']:.3f} | "
            f"{s['distinct_minima_estimate']} | {s['max_unitarity']:.1e} |")
    lines.append("")
    lines.append("convergence distribution (best_phi percentiles):")
    for opt_name in ("sgd", "adam"):
        phis = np.array([r["best_phi"] for r in results[opt_name]])
        pcts = np.percentile(phis, [0, 10, 25, 50, 75, 90, 100])
        lines.append(f"  {opt_name}: " + " ".join(f"{p:.3f}" for p in pcts))

    txt = "\n".join(lines)
    out_txt.write_text(txt, encoding="utf-8")
    out_json.write_text(json.dumps(
        dict(baselines=baselines, contiguous_best=contiguous_best,
             summary=summary, results=results, wall_time=elapsed,
             n_trials=n_trials, steps=STEPS, lambda_value=LAMBDA),
        indent=2, ensure_ascii=False), encoding="utf-8")
    plot(results, baselines, contiguous_best, out_png)

    print(txt)
    print(f"\nwrote {out_txt}")
    print(f"wrote {out_json}")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
