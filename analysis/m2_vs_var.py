#!/usr/bin/env python3
"""
M2 vs Var(H_boundary) — Phase A item 2 falsification test (2026-08-11).

Goal: decide whether the canonical dynamical cost C_H = <H_boundary^2> (the
second moment M2) can be maintained, or must be replaced by the true variance
Var(H) = <H^2> - <H>^2. This is a NARROWING test for the Phase A spec freeze;
the final choice is made in Phase C by comparing M2 and Var against the actual
open-system decoherence rate Gamma_F.

For each candidate TPS we record simultaneously:
    M2(F)      = <H_dF^2>
    V(F)       = <H_dF^2> - <H_dF>^2
    <H_dF>     = boundary mean
    Delta(F)   = M2(F) - V(F) = <H_dF>^2

Then, for the SAME lambda sweep, we compare
    Phi_M2  = I - lambda * M2
    Phi_Var = I - lambda * V
on: F*(lambda), lambda*, topology-switching count, and ranking.

Verdict (3-way):
    A: M2 and Var are effectively identical        -> keep M2, defer change
    B: selection results differ                     -> Phase C Gamma_F check mandatory
    C: Var is degenerate / unstable                 -> supports keeping M2

Outputs (analysis_output/):
    m2_vs_var_table.txt            numeric tables
    m2_vs_var_results.json         machine-readable results + verdict
    m2_vs_var_lambda_sweep.png     selected cut vs lambda (N=3,4,5)
    m2_vs_var_phi_curves_N3.png    Phi_M2 vs Phi_Var curves (N=3)
    m2_vs_var_delta.png            Delta(F)=<H>^2 per cut
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

BENCH = dict(j12=1.5, j23=0.6, h=0.2)   # v4 benchmark (J12=1.5, J23=0.6, h=0.2)
LAM_GRID = np.linspace(0.0, 3.0, 601)
BENCH_LAM = 0.2                          # ranking comparison point


def build_system(n: int) -> tuple:
    """XY chain, ground state. N=3 uses the v4 benchmark; N>=4 a weak-center profile."""
    if n == 3:
        couplings = [BENCH["j12"], BENCH["j23"]]
    else:
        couplings = rem4.coupling_profile(n, BENCH["j12"], BENCH["j23"], profile="weak_center")
    h_total, bond_terms = rem4.xy_chain_hamiltonian(n, couplings, BENCH["h"])
    _, psi = rem4.ground_state(h_total)
    return h_total, bond_terms, psi


def per_cut_metrics(psi: np.ndarray, n: int, bond_terms: dict) -> list[dict]:
    """Record I, <H>, M2=<H^2>, Var, Delta=<H>^2 for every contiguous cut."""
    rows = []
    for cut in range(1, n):
        bond = bond_terms[(cut - 1, cut)]
        mi = rem4.mutual_information_for_cut(psi, n, cut)
        m2 = rem4.boundary_cost_squared(psi, bond)
        mean = rem4.boundary_energy(psi, bond)
        var = m2 - mean * mean
        rows.append(dict(
            cut=cut, mi=mi, mean=mean, m2=m2, var=var, delta=mean * mean,
        ))
    return rows


def sweep_selection(rows: list[dict], lam_grid: np.ndarray) -> tuple:
    """Selected cut under Phi_M2 and Phi_Var for the same lambda sweep."""
    sel_m2, sel_var = [], []
    for lam in lam_grid:
        phi_m2 = np.array([r["mi"] - lam * r["m2"] for r in rows])
        phi_var = np.array([r["mi"] - lam * r["var"] for r in rows])
        sel_m2.append(int(np.argmax(phi_m2)) + 1)
        sel_var.append(int(np.argmax(phi_var)) + 1)
    sel_m2 = np.array(sel_m2)
    sel_var = np.array(sel_var)
    agreement = float(np.mean(sel_m2 == sel_var))
    return sel_m2, sel_var, agreement


def switching_count(sel: np.ndarray) -> int:
    return int(np.sum(sel[1:] != sel[:-1]))


def first_crossover(rows: list[dict], cost_key: str) -> float | None:
    """First positive crossover lambda from the info-selected cut, using cost_key cost."""
    i_cut = max(rows, key=lambda r: r["mi"])["cut"]
    i_row = next(r for r in rows if r["cut"] == i_cut)
    lambdas = []
    for r in rows:
        if r["cut"] == i_cut:
            continue
        dphi = i_row["mi"] - r["mi"]
        dc = i_row[cost_key] - r[cost_key]
        if abs(dc) < 1e-12:
            continue
        lam = dphi / dc
        if lam > 0:
            lambdas.append(lam)
    return float(min(lambdas)) if lambdas else None


def ranking_at(rows: list[dict], lam: float, cost_key: str) -> list[int]:
    """Cuts sorted by Phi = I - lambda*cost, descending."""
    phi = {r["cut"]: r["mi"] - lam * r[cost_key] for r in rows}
    return sorted(phi, key=lambda c: phi[c], reverse=True)


def classify(rows: list[dict], agreement: float,
             lam_star_m2: float | None, lam_star_var: float | None) -> dict:
    """Three-way verdict A / B / C with supporting flags."""
    m2s = np.array([r["m2"] for r in rows])
    deltas = np.array([r["delta"] for r in rows])
    vars_ = np.array([r["var"] for r in rows])
    max_delta_ratio = float(np.max(deltas / np.maximum(m2s, 1e-300)))
    var_degenerate = bool(np.any((vars_ < 1e-12) & (m2s > 1e-8)))
    lam_star_differs = (
        lam_star_m2 is not None and lam_star_var is not None
        and abs(lam_star_m2 - lam_star_var) > 1e-6
    )
    selection_differs = agreement < 1.0

    if var_degenerate:
        verdict = "C"
        reason = ("Var(H) is degenerate (Var~0 while M2>0): the variance proxy would "
                  "erase the dynamical term, so the second moment M2 is retained.")
    elif selection_differs or lam_star_differs:
        verdict = "B"
        reason = ("M2 and Var give different structural selections (F*, lambda*, "
                  "switching, or ranking); the Phase C Gamma_F comparison is mandatory "
                  "before freezing the canonical cost.")
    elif max_delta_ratio < 1e-3:
        verdict = "A"
        reason = ("Delta = <H>^2 is negligible relative to M2 on this benchmark: M2 and "
                  "Var are operationally identical; keep M2 and defer the decision.")
    else:
        verdict = "A"
        reason = ("Delta is non-negligible numerically but the two proxies produce the "
                  "same selections on this benchmark; keep M2, but the Phase C Gamma_F "
                  "comparison will be decisive.")

    return dict(
        verdict=verdict,
        reason=reason,
        max_delta_ratio=max_delta_ratio,
        var_degenerate=var_degenerate,
        selection_differs=selection_differs,
        lam_star_differs=lam_star_differs,
        agreement=agreement,
    )


def plot_lambda_sweep(data: dict, outpath: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.6), sharey=True)
    for ax, (key, d) in zip(axes, data.items()):
        n = d["n"]
        lam = d["lam_grid"]
        ax.step(lam, d["sel_m2"], where="post", label=r"$\Phi_{M_2}=I-\lambda M_2$", lw=1.6)
        ax.step(lam, d["sel_var"], where="post", label=r"$\Phi_{\mathrm{Var}}=I-\lambda V$",
                lw=1.6, ls="--", alpha=0.75)
        ax.set_title(f"N={n}  (agreement {d['agreement']*100:.1f}%)")
        ax.set_xlabel(r"$\lambda$")
        ax.set_yticks(list(range(1, n)))
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("selected cut $F^*(\\lambda)$")
    axes[0].legend(loc="upper right", fontsize=8)
    fig.suptitle("Topology switching: same $\\lambda$ sweep under $M_2$ vs $\\mathrm{Var}(H_\\partial)$")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def plot_phi_curves_n3(d: dict, outpath: Path) -> None:
    lam = d["lam_grid"]
    rows = d["rows"]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for r in rows:
        phi_m2 = np.array([r["mi"] - l * r["m2"] for l in lam])
        phi_var = np.array([r["mi"] - l * r["var"] for l in lam])
        ax.plot(lam, phi_m2, lw=1.8, label=f"cut {r['cut']} $M_2$")
        ax.plot(lam, phi_var, lw=1.2, ls="--", alpha=0.8,
                label=f"cut {r['cut']} $\\mathrm{{Var}}$")
    for name, val in (("$M_2$", d["lam_star_m2"]), ("$\\mathrm{Var}$", d["lam_star_var"])):
        if val is not None:
            ax.axvline(val, ls=":", alpha=0.5)
            ax.text(val, ax.get_ylim()[0], f" {name} $\\lambda^*$={val:.3f}", fontsize=8, rotation=90)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\Phi$")
    ax.set_title("N=3 benchmark: $\\Phi_{M_2}$ vs $\\Phi_{\\mathrm{Var}}$ curves")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def plot_delta(data: dict, outpath: Path) -> None:
    n_list = list(data.keys())
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    width = 0.25
    all_ticks: list[float] = []
    all_labels: list[str] = []
    for i, n in enumerate(n_list):
        rows = data[n]["rows"]
        cuts = [r["cut"] for r in rows]
        x = np.arange(len(cuts)) + i * 0.35
        ax.bar(x - width / 2, [r["m2"] for r in rows], width, label=f"N={n} $M_2$", alpha=0.85)
        ax.bar(x + width / 2, [r["var"] for r in rows], width, label=f"N={n} $\\mathrm{{Var}}$", alpha=0.85)
        ax.bar(x + 1.5 * width, [r["delta"] for r in rows], width * 0.8,
               label=f"N={n} $\\Delta=\\langle H\\rangle^2$", alpha=0.6)
        all_ticks.extend(x.tolist())
        all_labels.extend([f"N={n}\ncut {c}" for c in cuts])
    # set ticks once after the loop, otherwise later groups overwrite earlier ones
    ax.set_xticks(all_ticks, all_labels, fontsize=8)
    ax.set_ylabel("cost")
    ax.set_title(r"$\Delta(F)=M_2-V=\langle H_{\partial F}\rangle^2$ per candidate TPS")
    ax.legend(fontsize=7, ncol=3)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def run_n(n: int) -> dict:
    _, bond_terms, psi = build_system(n)
    rows = per_cut_metrics(psi, n, bond_terms)
    sel_m2, sel_var, agreement = sweep_selection(rows, LAM_GRID)
    lam_star_m2 = first_crossover(rows, "m2")
    lam_star_var = first_crossover(rows, "var")
    verdict = classify(rows, agreement, lam_star_m2, lam_star_var)
    return dict(
        n=n,
        rows=rows,
        lam_grid=LAM_GRID.tolist(),
        sel_m2=sel_m2.tolist(),
        sel_var=sel_var.tolist(),
        agreement=agreement,
        switches_m2=switching_count(sel_m2),
        switches_var=switching_count(sel_var),
        lam_star_m2=lam_star_m2,
        lam_star_var=lam_star_var,
        ranking_m2=ranking_at(rows, BENCH_LAM, "m2"),
        ranking_var=ranking_at(rows, BENCH_LAM, "var"),
        verdict=verdict,
    )


def write_table(data: dict, outpath: Path) -> None:
    lines = [
        "# M2 vs Var(H_boundary) — Phase A item 2 (2026-08-11)",
        f"# System: XY chain, J12={BENCH['j12']}, J23={BENCH['j23']}, h={BENCH['h']}, ground state.",
        f"# lambda sweep: [{LAM_GRID[0]}, {LAM_GRID[-1]}] ({len(LAM_GRID)} pts); ranking at lambda={BENCH_LAM}",
        "",
    ]
    for n, d in data.items():
        lines += [f"## N={n}", ""]
        lines += ["cut | I(A:B) [bits] | <H_dF> | M2=<H^2> | Var(H) | Delta=<H>^2 | Delta/M2",
                  "----|---------------|---------|----------|---------|-------------|---------"]
        for r in d["rows"]:
            lines.append(
                f"{r['cut']} | {r['mi']:.6f} | {r['mean']:+.6f} | {r['m2']:.6f} | "
                f"{r['var']:.6f} | {r['delta']:.6f} | {r['delta']/max(r['m2'],1e-300):.3e}"
            )
        lines += [
            "",
            f"agreement F*_M2 == F*_Var : {d['agreement']*100:.2f}%  "
            f"(switches M2={d['switches_m2']}, Var={d['switches_var']})",
            f"lambda* (M2)  = {d['lam_star_m2']}",
            f"lambda* (Var) = {d['lam_star_var']}",
            f"ranking at lambda={BENCH_LAM} (M2)  : {d['ranking_m2']}",
            f"ranking at lambda={BENCH_LAM} (Var) : {d['ranking_var']}",
            "",
        ]
    outpath.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    data = {f"N{n}": run_n(n) for n in (3, 4, 5)}

    write_table(data, OUTDIR / "m2_vs_var_table.txt")
    (OUTDIR / "m2_vs_var_results.json").write_text(
        json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")

    plot_lambda_sweep(data, OUTDIR / "m2_vs_var_lambda_sweep.png")
    plot_phi_curves_n3(data["N3"], OUTDIR / "m2_vs_var_phi_curves_N3.png")
    plot_delta(data, OUTDIR / "m2_vs_var_delta.png")

    print(f"wrote {OUTDIR / 'm2_vs_var_table.txt'}")
    print(f"wrote {OUTDIR / 'm2_vs_var_results.json'}")
    print(f"wrote {OUTDIR / 'm2_vs_var_lambda_sweep.png'}")
    print(f"wrote {OUTDIR / 'm2_vs_var_phi_curves_N3.png'}")
    print(f"wrote {OUTDIR / 'm2_vs_var_delta.png'}")
    for n, d in data.items():
        v = d["verdict"]
        print(f"\nN={n}: verdict={v['verdict']}  agreement={d['agreement']*100:.2f}%  "
              f"lambda* M2={d['lam_star_m2']} Var={d['lam_star_var']}  "
              f"max Delta/M2={v['max_delta_ratio']:.3e}")
        print(f"  {v['reason']}")


if __name__ == "__main__":
    main()
