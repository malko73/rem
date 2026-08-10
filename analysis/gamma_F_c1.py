#!/usr/bin/env python3
"""
Phase C1 — M2 / Var vs Gamma_F across 4 pre-registered environments (2026-08-11).

Extends C0 (uniform pure dephasing, case 2: Gamma_1 ~ Gamma_2) to multiple
noise channels, with all environment parameters FIXED BEFORE looking at the
results. For every environment H, rho0, and {L_i} are identical for cut 1
and cut 2 (only the factorization F changes).

Pre-registered environments:
    A1  non-uniform dephasing   gamma = (0.5, 1, 2), kappa = 0
    A2  mirrored dephasing      gamma = (2, 1, 0.5), kappa = 0
    B   amplitude damping       gamma = 0, kappa = (1, 1, 1)
    C   mixed                   gamma = (0.5,0.5,0.5), kappa = (0.5,0.5,0.5)

Measured per cut: Gamma^exact (analytic t=0, primary), Gamma^(0) (t<=0.05
log-linear), Gamma^fit (t<=0.5). Reported ratio R_Gamma = Gamma_1 / Gamma_2
vs R_M2 = M2_1 / M2_2 ~ 10.23.

Classification (pre-registered):
    C1-A  M2 supported:   across envs Gamma_2 < Gamma_1 substantially
    C1-B  M2 unsupported: across envs Gamma_1 ~ Gamma_2 reproduced
    C1-C  env-dependent:  ordering flips between environments

Outputs (analysis_output/):
    gamma_F_c1_table.txt, gamma_F_c1_results.json, gamma_F_c1_ratios.png
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"
MODULE_PATH = REPO / "src" / "rem4_numerical.py"

spec = importlib.util.spec_from_file_location("rem4_numerical", MODULE_PATH)
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)

GPATH = REPO / "analysis" / "gamma_F.py"
gspec = importlib.util.spec_from_file_location("gamma_F", GPATH)
gf = importlib.util.module_from_spec(gspec)
assert gspec.loader is not None
sys.modules[gspec.name] = gf
gspec.loader.exec_module(gf)

# ── pre-registered environments (fixed before running) ───────────────
ENVS = {
    "A1_deph_0512": dict(gamma=(0.5, 1.0, 2.0), kappa=(0.0, 0.0, 0.0)),
    "A2_deph_2105": dict(gamma=(2.0, 1.0, 0.5), kappa=(0.0, 0.0, 0.0)),
    "B_ampdamp":    dict(gamma=(0.0, 0.0, 0.0), kappa=(1.0, 1.0, 1.0)),
    "C_mixed":      dict(gamma=(0.5, 0.5, 0.5), kappa=(0.5, 0.5, 0.5)),
}
REL_TOL = 0.10   # |R_Gamma - 1| < 10%  ->  "Gamma_1 ~ Gamma_2"
A_TOL = 1.10     # R_Gamma > 1.10      ->  "Gamma_2 < Gamma_1 substantially"


def run_env(psi, h_total, t_grid, gamma_vec, kappa_vec) -> dict:
    L = gf.liouvillian_env(h_total, list(gamma_vec), list(kappa_vec), n=3)
    d1 = gf.run_cut(psi, h_total, 1.0, t_grid, L, cut=1)
    d2 = gf.run_cut(psi, h_total, 1.0, t_grid, L, cut=2)
    r_exact = d1["gamma_exact"] / d2["gamma_exact"]
    r0 = d1["gamma0"] / d2["gamma0"]
    r_fit = d1["gamma_fit"] / d2["gamma_fit"]
    return dict(
        gamma_vec=list(gamma_vec), kappa_vec=list(kappa_vec),
        cut1=dict(gamma_exact=d1["gamma_exact"], gamma0=d1["gamma0"],
                  gamma_fit=d1["gamma_fit"], c0=d1["c0"], m2=d1["m2"], var=d1["var"]),
        cut2=dict(gamma_exact=d2["gamma_exact"], gamma0=d2["gamma0"],
                  gamma_fit=d2["gamma_fit"], c0=d2["c0"], m2=d2["m2"], var=d2["var"]),
        r_gamma_exact=float(r_exact), r_gamma0=float(r0), r_gamma_fit=float(r_fit),
    )


def classify_c1(results: dict) -> dict:
    ratios = {env: r["r_gamma_exact"] for env, r in results.items()}
    all_flat = all(abs(v - 1.0) < REL_TOL for v in ratios.values())
    all_m2 = all(v > A_TOL for v in ratios.values())
    if all_flat:
        verdict = "C1-B"
        text = ("Gamma_1 ~ Gamma_2 across all environments: M2's 10x spread is not "
                "reflected in open-system stability; lambda* ~ 0.165 loses its "
                "physical-crossover basis -> move to C_open(F; E).")
    elif all_m2:
        verdict = "C1-A"
        text = ("Gamma_2 < Gamma_1 consistently across environments: M2 is supported "
                "as an open-system stability proxy; keep canonical cost.")
    else:
        verdict = "C1-C"
        text = ("Ordering is environment-dependent: the dynamical cost must be "
                "environment-dependent, C_H = C_H(F, E) (derive from the actual "
                "open-system generator).")
    return dict(verdict=verdict, text=text,
                r_gamma_exact=ratios,
                r_m2=results[list(results)[0]]["cut1"]["m2"]
                     / results[list(results)[0]]["cut2"]["m2"])


def write_table(results: dict, verdict: dict, outpath: Path) -> None:
    lines = [
        "# Phase C1 — Gamma_F across 4 pre-registered environments (2026-08-11)",
        "# H, rho0, {L_i} identical for cut1/cut2 in every environment.",
        "# Gamma^exact is the primary indicator; R_Gamma = Gamma_1/Gamma_2.",
        "",
        "| env | L_i | Gamma1^exact | Gamma2^exact | R_Gamma(exact) | R_Gamma(0.05) | R_Gamma(0.5) |",
        "|-----|-----|-------------|-------------|----------------|---------------|--------------|",
    ]
    for env, r in results.items():
        gv = ",".join(f"{g:g}" for g in r["gamma_vec"])
        kv = ",".join(f"{k:g}" for k in r["kappa_vec"])
        desc = f"gamma=({gv})" if r["kappa_vec"][0] == 0 else (
            f"kappa=({kv})" if r["gamma_vec"][0] == 0 else f"gamma=({gv}) kappa=({kv})")
        lines.append(
            f"| {env} | {desc} | {r['cut1']['gamma_exact']:.6f} | "
            f"{r['cut2']['gamma_exact']:.6f} | {r['r_gamma_exact']:.4f} | "
            f"{r['r_gamma0']:.4f} | {r['r_gamma_fit']:.4f} |"
        )
    lines += [
        "",
        f"R_M2 = M2_1/M2_2 = {verdict['r_m2']:.4f}",
        f"verdict: {verdict['verdict']} — {verdict['text']}",
    ]
    outpath.write_text("\n".join(lines), encoding="utf-8")


def plot_ratios(results: dict, verdict: dict, outpath: Path) -> None:
    envs = list(results.keys())
    r_exact = [results[e]["r_gamma_exact"] for e in envs]
    r_fit = [results[e]["r_gamma_fit"] for e in envs]
    r_m2 = verdict["r_m2"]
    x = np.arange(len(envs))
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.axhline(1.0, color="gray", ls=":", label=r"$\Gamma_1=\Gamma_2$")
    ax.axhline(r_m2, color="tab:red", ls="--",
               label=rf"$R_{{M_2}}={r_m2:.1f}$ (M2 prediction)")
    ax.plot(x, r_exact, "o-", label=r"$R_\Gamma$ exact (t=0)", lw=1.8)
    ax.plot(x, r_fit, "s--", label=r"$R_\Gamma$ fit (t$\leq$0.5)", alpha=0.7)
    ax.set_xticks(x, envs, rotation=15, fontsize=8)
    ax.set_ylabel(r"$R_\Gamma=\Gamma_1/\Gamma_2$")
    ax.set_title(f"Phase C1: {verdict['verdict']} — "
                 f"{verdict['text'][:60]}...")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    psi, h_total = gf.build_system()
    t_grid = np.linspace(0.0, gf.T_MAX, gf.N_T)

    results = {}
    for env, params in ENVS.items():
        results[env] = run_env(psi, h_total, t_grid,
                               params["gamma"], params["kappa"])
    verdict = classify_c1(results)

    write_table(results, verdict, OUTDIR / "gamma_F_c1_table.txt")
    (OUTDIR / "gamma_F_c1_results.json").write_text(
        json.dumps(dict(environments=results, verdict=verdict,
                        rel_tol=REL_TOL, a_tol=A_TOL),
                   indent=2, ensure_ascii=False), encoding="utf-8")
    plot_ratios(results, verdict, OUTDIR / "gamma_F_c1_ratios.png")

    print(f"wrote {OUTDIR / 'gamma_F_c1_table.txt'}")
    print(f"wrote {OUTDIR / 'gamma_F_c1_results.json'}")
    print(f"wrote {OUTDIR / 'gamma_F_c1_ratios.png'}")
    print("\n" + (OUTDIR / "gamma_F_c1_table.txt").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
