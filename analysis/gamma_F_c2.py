#!/usr/bin/env python3
"""
Phase C2 — Liouvillian boundary cost vs Gamma_F (2026-08-11).

For each of the 4 pre-registered C1 environments, decompose the vectorised
Liouvillian as

    L_dF = L - P_F^local L

where P_F^local is the Hilbert-Schmidt orthogonal projection onto
    S_F = { L_A (x) I_{B^2} + I_{A^2} (x) L_B }
(built from matrix units, pseudoinverse projection — independent of the
Lindblad jump-operator representation).

Three costs per cut:
    C_L^global(F) = |L_dF|_HS^2 / |L|_HS^2          (state-independent)
    C_L^rho   (F) = |L_dF(rho0)|_2^2 / |rho0|_2^2   (state-dependent)
    C_L^align (F) = -Re<Q_F rho, Q_F L_dF(rho)> / |Q_F rho|^2   (diagnostic)

Compare R_L = C_L(F1)/C_L(F2) with R_Gamma = Gamma_1/Gamma_2 across the 4
environments (Spearman rank correlation + log-ratio error).

Verdict (pre-registered):
    C2-A: C_L reproduces the environment-dependent ordering (incl. reversal)
          -> promote C_L to canonical dynamical-cost candidate #1
    C2-B: Gamma^exact robust but C_L cannot predict
          -> keep Gamma^exact as the operational canonical cost
    C2-C: initial-rate vs finite-time-rate structure selection differs
          -> make the timescale tau an essential parameter

Outputs (analysis_output/):
    gamma_F_c2_table.txt, gamma_F_c2_results.json, gamma_F_c2_compare.png
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
from scipy.stats import spearmanr

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

ENVS = {
    "A1_deph_0512": dict(gamma=(0.5, 1.0, 2.0), kappa=(0.0, 0.0, 0.0)),
    "A2_deph_2105": dict(gamma=(2.0, 1.0, 0.5), kappa=(0.0, 0.0, 0.0)),
    "B_ampdamp":    dict(gamma=(0.0, 0.0, 0.0), kappa=(1.0, 1.0, 1.0)),
    "C_mixed":      dict(gamma=(0.5, 0.5, 0.5), kappa=(0.5, 0.5, 0.5)),
}


def run_env_c2(psi, h_total, t_grid, gamma_vec, kappa_vec) -> dict:
    L = gf.liouvillian_env(h_total, list(gamma_vec), list(kappa_vec), n=3)
    d1 = gf.run_cut(psi, h_total, 1.0, t_grid, L, cut=1)
    d2 = gf.run_cut(psi, h_total, 1.0, t_grid, L, cut=2)
    rho0 = np.outer(psi, psi.conj())

    out = {"gamma_exact": [d1["gamma_exact"], d2["gamma_exact"]],
           "gamma_fit": [d1["gamma_fit"], d2["gamma_fit"]],
           "cuts": {}}
    for cut, d in ((1, d1), (2, d2)):
        l_df = gf.boundary_residual(L, 3, cut)
        basis = gf.schmidt_basis(psi, 3, cut)
        out["cuts"][str(cut)] = dict(
            c_L_global=gf.c_L_global(L, l_df),
            c_L_rho=gf.c_L_rho(l_df, rho0),
            c_L_align=gf.c_L_align(rho0, l_df, basis),
        )
    c1, c2 = out["cuts"]["1"], out["cuts"]["2"]
    out["r_L_global"] = c1["c_L_global"] / c2["c_L_global"]
    out["r_L_rho"] = c1["c_L_rho"] / c2["c_L_rho"]
    out["r_L_align"] = c1["c_L_align"] / c2["c_L_align"]
    out["r_gamma_exact"] = d1["gamma_exact"] / d2["gamma_exact"]
    out["r_gamma_fit"] = d1["gamma_fit"] / d2["gamma_fit"]
    return out


def verdict_c2(results: dict, tol_log: float = 0.6) -> dict:
    envs = list(results.keys())
    r_g = np.array([results[e]["r_gamma_exact"] for e in envs])
    r_gf = np.array([results[e]["r_gamma_fit"] for e in envs])
    r_gl = np.array([results[e]["r_L_global"] for e in envs])
    r_rl = np.array([results[e]["r_L_rho"] for e in envs])

    log_err_global = np.abs(np.log(r_gl) - np.log(r_g))
    log_err_rho = np.abs(np.log(r_rl) - np.log(r_g))
    rho_global = float(spearmanr(r_g, r_gl).statistic)
    rho_rho = float(spearmanr(r_g, r_rl).statistic)

    # C2-C check: initial-rate vs finite-time-rate STRUCTURE SELECTION differs
    # substantially (|R_fit - 1| > 0.1 and opposite side of 1 vs exact).
    init_ord = np.sign(r_g - 1.0)
    fit_ord = np.sign(r_gf - 1.0)
    subst = np.abs(r_gf - 1.0) > 0.10
    c2c = bool(np.any((init_ord != fit_ord) & subst))

    # C2-A: sign pattern of R_Gamma reproduced by a C_L variant with small log error
    best_err = min(float(log_err_global.max()), float(log_err_rho.max()))
    sign_ok_g = bool(np.all(np.sign(r_gl - 1.0) == np.sign(r_g - 1.0)))
    sign_ok_r = bool(np.all(np.sign(r_rl - 1.0) == np.sign(r_g - 1.0)))
    sign_ok = sign_ok_g or sign_ok_r

    if c2c:
        verdict = "C2-C"
        text = ("Initial-rate and finite-time-rate structure selection differ; "
                "the timescale tau becomes an essential parameter of C_dyn.")
    elif sign_ok and best_err < tol_log:
        verdict = "C2-A"
        text = ("A C_L variant reproduces the environment-dependent ordering "
                "including the A1/A2 reversal (small log-ratio error) -> promote "
                "C_L to canonical dynamical-cost candidate #1.")
    else:
        verdict = "C2-B"
        text = ("Gamma^exact is robust but the norm-type C_L costs do not track "
                "it across environments -> keep Gamma^exact as the operational "
                "canonical cost (or investigate alignment diagnostics).")

    return dict(
        verdict=verdict, text=text,
        spearman_global=rho_global, spearman_rho=rho_rho,
        max_log_err_global=float(log_err_global.max()),
        max_log_err_rho=float(log_err_rho.max()),
        sign_ok_global=sign_ok_g, sign_ok_rho=sign_ok_r, c2c=c2c,
        r_gamma_exact=r_g.tolist(), r_L_global=r_gl.tolist(), r_L_rho=r_rl.tolist(),
    )


def write_table(results: dict, verdict: dict, outpath: Path) -> None:
    lines = [
        "# Phase C2 — Liouvillian boundary cost vs Gamma_F (2026-08-11)",
        "# L_dF = L - P_F^local L; P_F^local = HS projection onto {L_A x I + I x L_B}",
        "",
        "| env | R_Gamma(exact) | R_L(global) | R_L(rho) | R_L(align) | R_Gamma(fit) |",
        "|-----|----------------|-------------|----------|------------|--------------|",
    ]
    for env, r in results.items():
        c1, c2 = r["cuts"]["1"], r["cuts"]["2"]
        lines.append(
            f"| {env} | {r['r_gamma_exact']:.4f} | {r['r_L_global']:.4f} | "
            f"{r['r_L_rho']:.4f} | {r['r_L_align']:.4f} | {r['r_gamma_fit']:.4f} |"
        )
    lines += [
        "",
        f"Spearman(global)={verdict['spearman_global']:.3f}  "
        f"Spearman(rho)={verdict['spearman_rho']:.3f}",
        f"max |log R_L - log R_G| : global={verdict['max_log_err_global']:.3f}  "
        f"rho={verdict['max_log_err_rho']:.3f}",
        f"verdict: {verdict['verdict']} — {verdict['text']}",
    ]
    outpath.write_text("\n".join(lines), encoding="utf-8")


def plot_compare(results: dict, verdict: dict, outpath: Path) -> None:
    envs = list(results.keys())
    r_g = [results[e]["r_gamma_exact"] for e in envs]
    r_gl = [results[e]["r_L_global"] for e in envs]
    r_rl = [results[e]["r_L_rho"] for e in envs]
    x = np.arange(len(envs))
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.axhline(1.0, color="gray", ls=":", label=r"$R=1$")
    ax.plot(x, r_g, "o-", label=r"$R_\Gamma$ exact", lw=1.8)
    ax.plot(x, r_gl, "s--", label="R_L global", alpha=0.85)
    ax.plot(x, r_rl, "d-.", label="R_L rho", alpha=0.85)
    ax.set_xticks(x, envs, rotation=15, fontsize=8)
    ax.set_ylabel("ratio $C(F_1)/C(F_2)$")
    ax.set_title(f"Phase C2: {verdict['verdict']} — "
                 f"Spearman G={verdict['spearman_global']:.2f}, "
                 f"rho={verdict['spearman_rho']:.2f}")
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
        results[env] = run_env_c2(psi, h_total, t_grid,
                                  params["gamma"], params["kappa"])
    verdict = verdict_c2(results)

    write_table(results, verdict, OUTDIR / "gamma_F_c2_table.txt")
    (OUTDIR / "gamma_F_c2_results.json").write_text(
        json.dumps(dict(environments=results, verdict=verdict),
                   indent=2, ensure_ascii=False), encoding="utf-8")
    plot_compare(results, verdict, OUTDIR / "gamma_F_c2_compare.png")

    print(f"wrote {OUTDIR / 'gamma_F_c2_table.txt'}")
    print(f"wrote {OUTDIR / 'gamma_F_c2_results.json'}")
    print(f"wrote {OUTDIR / 'gamma_F_c2_compare.png'}")
    print("\n" + (OUTDIR / "gamma_F_c2_table.txt").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
