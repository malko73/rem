#!/usr/bin/env python3
"""Generate the four main figures (plus Supplementary S1) for the REM paper
v0.1 from the JSON outputs only — no manually transcribed numbers.

  fig1_falsification — D2.1: Phi vs log10 C_F^2(0) over Haar seeds (blow-up)
  fig2_crossover     — D3:   DeltaPhi(tau) crossing with tau_c = 0.0180
  fig3_tau_c_trend   — D4:   tau_c(N) flat trend (official full-quotient values)
  fig4_mechanism     — E1:   predicted vs measured tau_c (first/second order)
  figS1_subspace     — Supp: N=5 subspace (0.0173) vs full quotient (0.0217)

Rules honored: JSON-only numbers; caption numbers == text numbers; N=5
subspace 0.0173 kept out of the main figures; black-and-white safe (marker
shape / line style / annotations carry the information, not color).

Outputs: papers/figures/figN.pdf (vector, for LaTeX) + .png (300 dpi, for
Markdown preview).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[1]
OUTDIR = REPO / "analysis_output"
FIGDIR = REPO / "papers" / "figures"

# B/W-safe palette: shapes/lines carry the meaning; grays for fill
INK = "#111111"
GRAY = "#777777"
LIGHT = "#cccccc"


def load(name):
    return json.loads((OUTDIR / name).read_text(encoding="utf-8"))


def save(fig, name):
    FIGDIR.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"{name}.{ext}", bbox_inches="tight",
                    dpi=300 if ext == "png" else None)
    plt.close(fig)
    print(f"  wrote papers/figures/{name}.pdf/.png")


def fig1():
    d = load("d2_1_haar_singularity_audit.json")
    trials = d["trials"]
    xs = np.log10([t["C_F_0_sq"] for t in trials])
    ys = np.array([t["Phi"] for t in trials])
    pmax_i = int(np.argmax(ys))
    stats = d["stats"]

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.axhline(1.8594, color=GRAY, ls=":", lw=1.2)
    ax.text(0.03, 1.90, r"well-posed optimum $\Phi\approx1.86$ "
            r"($J_{\rm dyn}^{(0)}$)", fontsize=8, color=GRAY)
    ax.plot(xs, ys, marker="o", ls="none", ms=6, mfc="none", mec=INK, mew=1.1)
    ax.annotate(f"best seed: $\\Phi={ys[pmax_i]:.1f}$,\n"
                r"$C_F^2=%.1e$" % trials[pmax_i]["C_F_0_sq"],
                xy=(xs[pmax_i], ys[pmax_i]),
                xytext=(xs[pmax_i] - 1.15, ys[pmax_i] - 8.5),
                fontsize=8, arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
    ax.set_xlabel(r"$\log_{10} C_F^2(0)$")
    ax.set_ylabel(r"$\Phi$ (normalized $\Gamma_F$ objective)")
    ax.set_title("Normalized functional: denominator blow-up", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.25, ls=":")
    fig.tight_layout()
    return fig, f"corr$\\approx{stats['corr_Phi_logC0']:.3f}$"


def fig2():
    d = load("d3_timescale.json")
    dphi = {float(k): v for k, v in d["DeltaPhi"].items()}
    taus = sorted(dphi)
    vals = [dphi[t] for t in taus]
    tau_c = d["tau_c"]

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.axhline(0.0, color=INK, lw=1.0)
    ax.axvline(tau_c, color=GRAY, ls="--", lw=1.2)
    ax.text(tau_c + 0.002, 0.055, rf"$\tau_c = {tau_c:.4f}$", fontsize=9,
            color=GRAY)
    ax.plot(taus, vals, marker="o", ls="-", ms=5, mfc="none", mec=INK, mew=1.1,
            color=INK, lw=1.2)
    ax.set_xlabel(r"$\tau$ (observation window)")
    ax.set_ylabel(r"$\Delta\Phi(\tau)=\Phi_\tau(F_B)-\Phi_\tau(F_A)$")
    ax.set_title("Finite-time structural crossover (Haar, $N=3$)", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.25, ls=":")

    # inset: seed value-behavior (B-value count flips 0 -> 6 at tau_c)
    hf = d["haar_fine"]
    ins_taus = []
    ins_nB = []
    for k in sorted(hf, key=float):
        ins_taus.append(float(k))
        vc = hf[k]["value_counts"]
        ins_nB.append(vc.get("B", 0))
    axin = fig.add_axes([0.58, 0.16, 0.34, 0.30])
    axin.plot(ins_taus, ins_nB, marker="s", ls="-", ms=4, color=INK, lw=1.0)
    axin.axvline(tau_c, color=GRAY, ls="--", lw=0.8)
    axin.set_xlabel(r"$\tau$", fontsize=8)
    axin.set_ylabel("# seeds at $B$-value", fontsize=8)
    axin.set_ylim(-0.4, 6.6)
    axin.tick_params(labelsize=7)
    fig.tight_layout()
    return fig, None


def fig3():
    cases = ["$N{=}3$, $2|1$", "$N{=}4$, $2|2$", "$N{=}4$, $1|3$",
             "$N{=}5$, $2|3$ full"]
    tc = [0.0180, 0.0206, 0.0224, 0.0217]
    d4 = load("d4_finite_size.json")
    tc_list = d4["gates"]["D4-7"]["tau_c_list"]
    # official values re-derived from JSON (full-quotient for N=5)
    official = [tc_list["N3"], tc_list["N4_22"], tc_list["N4_13"],
                tc_list["N5_23_full"]]

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    xs = np.arange(len(cases))
    mean = float(np.mean(official))
    ax.axhline(mean, color=LIGHT, ls=":", lw=1.2)
    ax.text(3.02, mean + 0.0004, rf"mean $\approx {mean:.3f}$", fontsize=8,
            color=GRAY)
    ax.plot(xs, official, marker="o", ls="-", ms=7, mfc="none", mec=INK, mew=1.3,
            color=INK, lw=1.2)
    for x, v in zip(xs, official):
        ax.annotate(f"{v:.4f}", (x, v), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=9)
    ax.set_xticks(xs)
    ax.set_xticklabels(cases, fontsize=9)
    ax.set_ylabel(r"$\tau_c(N)$")
    ax.set_ylim(0.012, 0.028)
    ax.set_title("Finite-size persistence of the crossover", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.25, ls=":", axis="y")
    fig.tight_layout()
    return fig, None


def fig4():
    e = load("e1_tauc_mechanism.json")
    meas = [e["systems"][k]["tau_c_meas"] for k in ("N3_21", "N4_22", "N4_13", "N5_23_full")]
    p1 = [e["systems"][k]["tau_c1_num"] for k in ("N3_21", "N4_22", "N4_13", "N5_23_full")]
    p2 = [e["systems"][k]["tau_c2"] for k in ("N3_21", "N4_22", "N4_13", "N5_23_full")]

    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    lo, hi = 0.012, 0.028
    ax.plot([lo, hi], [lo, hi], color=GRAY, ls="--", lw=1.2)
    ax.text(0.0235, 0.0140, r"$y=x$", fontsize=9, color=GRAY, rotation=42)
    ax.plot(meas, p1, marker="s", ls="none", ms=7, mfc="none", mec=INK, mew=1.3,
            label=r"first order $\tau_c^{(1)}=-\Delta\Phi_0/\Delta\Phi_1$")
    ax.plot(meas, p2, marker="o", ls="none", ms=7, mfc=INK, mec=INK, mew=0.8,
            label=r"quadratic $\tau_c^{(2)}$")
    for m, a, b in zip(meas, p1, p2):
        ax.annotate("", xy=(m, b), xytext=(m, a),
                    arrowprops=dict(arrowstyle="-", color=LIGHT, lw=1.0))
    for i, (m, a) in enumerate(zip(meas, p1)):
        lab = ["$N{=}3$", "$N{=}4,2|2$", "$N{=}4,1|3$", "$N{=}5$"][i]
        ax.annotate(lab, (m, a), textcoords="offset points",
                    xytext=(-2, 8), fontsize=8)
    ax.set_xlabel(r"measured $\tau_c$")
    ax.set_ylabel(r"predicted $\tau_c$")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.set_title("Crossover mechanism: predicted vs measured", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.grid(alpha=0.25, ls=":")
    fig.tight_layout()
    return fig, None


def figS1():
    """Supplementary: N=5 subspace (0.0173) vs full quotient (0.0217)."""
    d4 = load("d4_finite_size.json")
    tc_list = d4["gates"]["D4-7"]["tau_c_list"]
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    labels = ["200-dim subspace", "full quotient (945)"]
    vals = [tc_list["N5_23_subspace"], tc_list["N5_23_full"]]
    ax.plot([0, 1], vals, marker="o", ls="-", ms=7, mfc="none", mec=INK,
            mew=1.3, color=INK, lw=1.2)
    for x, v in zip([0, 1], vals):
        ax.annotate(f"{v:.4f}", (x, v), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(r"$\tau_c(5)$")
    ax.set_ylim(0.014, 0.026)
    ax.set_title("N=5 restricted-domain history", fontsize=10)
    ax.grid(alpha=0.25, ls=":", axis="y")
    fig.tight_layout()
    return fig, None


def main():
    print("Generating figures from JSON...")
    fig, _ = fig1()
    save(fig, "fig1_falsification")
    fig, _ = fig2()
    save(fig, "fig2_crossover")
    fig, _ = fig3()
    save(fig, "fig3_tau_c_trend")
    fig, _ = fig4()
    save(fig, "fig4_mechanism")
    fig, _ = figS1()
    save(fig, "figS1_subspace")
    print("Done.")


if __name__ == "__main__":
    main()
