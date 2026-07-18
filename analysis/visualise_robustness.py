#!/usr/bin/env python3
"""
Visualise robustness analysis v2 results.

Produces:
  - analysis_output/robustness_table.txt     summary table
  - analysis_output/robustness_N3.png        violin N=3
  - analysis_output/robustness_N4.png        violin N=4
  - analysis_output/robustness_N5.png        violin N=5
"""

from __future__ import annotations

import json
import os
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent.parent
RESULTS_PATH = HERE / "analysis_output" / "robustness_results.json"
RAW_PATH = HERE / "analysis_output" / "robustness_raw.json"
OUTDIR = HERE / "analysis_output"

COLOR_A = "#4C72B0"
COLOR_B = "#DD8452"
COLOR_CONT = "gray"


def main() -> None:
    with open(RESULTS_PATH) as f:
        results = json.load(f)

    with open(RAW_PATH) as f:
        raw = json.load(f)

    os.makedirs(OUTDIR, exist_ok=True)

    # Extract unique N values from results keys
    ns = sorted(set(
        int(k.split("_")[0].split("=")[1]) for k in results
    ))

    lines = [
        "=" * 75,
        "  Robustness Analysis — Continuous Factorisation Optimisation",
        "  Benchmark: XY chain, J12=1.5, J23=0.6, h=0.2, λ=0.2",
        "=" * 75,
        "",
        f"  {'N':<5} {'n_a':<5} {'axis':<5}{'trials':<7} {'contiguous Φ':<13}"
        f" {'optimised Φ':<18} {'improv. ratio':<14} {'success':<8} {'trial/s':<8}",
        "  " + "-" * 75,
    ]

    for n in ns:
        for axis in ("A", "B"):
            key = f"N={n}_{axis}"
            d = results[key]
            label = {3: "N=3", 4: "N=4", 5: "N=5"}[n]
            ax_msg = "A:init-θ" if axis == "A" else "B:generator"
            n_trials = d.get("n_trials", "?")
            stats = (
                f"{d['phi_mean']:.4f} ± {d['phi_std']:.4f}"
                f"  [{d['phi_min']:.4f}, {d['phi_max']:.4f}]"
            )
            ratio = f"{d['improvement_ratio_mean']:.3f}x"
            success = f"{d['success_rate']*100:.0f}%"
            lines.append(
                f"  {n:<5} {d.get('n_a','?'):<5} {ax_msg:<5}"
                f" {n_trials:<7} {d['contiguous_phi']:<13.4f}"
                f" {stats:<18} {ratio:<14} {success:<8}"
            )

    lines.append("")
    lines.append("  Key findings:")
    lines.append("  - All trials: success if n_a == contiguous best cut")
    lines.append("  - N=3: improvement ~1.9x (contiguous=0.56)")
    lines.append("  - N=4: improvement ~1.1x (contiguous=1.72)")
    lines.append("  - N=5: improvement ~1.1x (contiguous=1.46)")
    lines.append("  - Improvement ratio decreases with N as contiguous Φ rises")
    lines.append("")

    table_txt = "\n".join(lines)
    print(table_txt)
    with open(os.path.join(OUTDIR, "robustness_table.txt"), "w") as f:
        f.write(table_txt)

    # ── Violin plots per N ──
    for n in ns:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)

        # filter raw data for this N
        raw_n = [t for t in raw if t["n"] == n]
        raw_a = [t for t in raw_n if t["gen_seed"] == 0]
        raw_b = [t for t in raw_n if t["gen_seed"] != 0]

        if not raw_a or not raw_b:
            print(f"  WARNING: insufficient raw data for N={n}, skipping plot")
            plt.close()
            continue

        phis_a = np.array([t["best_phi"] for t in raw_a])
        phis_b = np.array([t["best_phi"] for t in raw_b])
        cont_phi = raw_a[0]["contiguous_phi"]

        # A panel
        parts_a = ax1.violinplot([phis_a], positions=[1],
                                 showmeans=True, showmedians=True, widths=0.6)
        for pc in parts_a["bodies"]:
            pc.set_facecolor(COLOR_A)
            pc.set_alpha(0.6)
        parts_a["cmeans"].set_color(COLOR_A)
        parts_a["cmedians"].set_color(COLOR_A)
        ax1.axhline(cont_phi, color=COLOR_CONT, linestyle="--",
                    linewidth=1, alpha=0.7, label=f"Best contiguous ({cont_phi:.4f})")
        ax1.set_xticks([1])
        ax1.set_xticklabels(["A: initial-θ"])
        ax1.set_ylabel(r"$\Phi$")
        ax1.set_title(f"Fixed generator basis\n{len(raw_a)} random θ seeds")
        ax1.legend(fontsize=8)

        # B panel
        parts_b = ax2.violinplot([phis_b], positions=[1],
                                 showmeans=True, showmedians=True, widths=0.6)
        for pc in parts_b["bodies"]:
            pc.set_facecolor(COLOR_B)
            pc.set_alpha(0.6)
        parts_b["cmeans"].set_color(COLOR_B)
        parts_b["cmedians"].set_color(COLOR_B)
        ax2.axhline(cont_phi, color=COLOR_CONT, linestyle="--",
                    linewidth=1, alpha=0.7, label=f"Best contiguous ({cont_phi:.4f})")
        ax2.set_xticks([1])
        ax2.set_xticklabels(["B: generator basis"])
        ax2.set_title(f"{len(raw_b)} random generator bases\n(1 θ seed each)")
        ax2.legend(fontsize=8)

        fig.suptitle(
            f"N = {n}  (n_a = {raw_a[0]['n_a']})  |  "
            f"contiguous Φ = {cont_phi:.4f}  |  "
            f"A: Φ = {np.mean(phis_a):.2f}±{np.std(phis_a):.2f}  "
            f"B: Φ = {np.mean(phis_b):.2f}±{np.std(phis_b):.2f}",
            fontsize=11, y=1.03,
        )
        plt.tight_layout()
        outpath = os.path.join(OUTDIR, f"robustness_N{n}.png")
        plt.savefig(outpath, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {outpath}")

    print("\nDone.")


if __name__ == "__main__":
    main()
