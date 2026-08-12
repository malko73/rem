#!/usr/bin/env python3
"""Phase E supplement — uncertainty and sensitivity analysis for Claim 4c.

Addresses Reviewer B's Critical finding 2.6: the "derived, not fitted"
claim needs (a) the actual DeltaPhi_2 values, (b) uncertainty propagation,
(c) sensitivity to the fixed representatives F_A/F_B, (d) a tau-grid
resolution check. Uses ONLY existing JSON data (no new optimizations).

Outputs analysis_output/e1_sensitivity.json + prints a summary.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
OUTDIR = REPO / "analysis_output"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from d4_nscaling import build_system, build_states as d4_states, CUTS, LAMBDA
from d1_hamiltonians import get_hamiltonian
from b2_canonical_open import liouvillian_dephasing
from d2_2a_unnormalized_instantaneous import build_states as d2a_states
from b2_1_cross_eval import tps_distance
from e1_tauc_mechanism import rep_data, load_U, dominant_eigenvector
import scipy.linalg as sla


SYSTEMS = [
    dict(label="N3_21", N=3, cut="21",
         FA=("d2_2a_unnormalized_instantaneous.json", ["core", "haar"]),
         FB=("d2_2b_tau_t1e-1.json", ["core", "haar"]),
         alt_tau_jsons=[("d2_2b_tau_t5e-2.json", ["core", "haar"]),
                        ("d2_2b_tau_t3e-2.json", ["core", "haar"])]),
    dict(label="N4_22", N=4, cut="22",
         FA=("d4_b_n4_cut22_haar.json", ["F_A"]),
         FB=("d4_b_n4_cut22_haar.json", ["F_B"])),
    dict(label="N4_13", N=4, cut="13",
         FA=("d4_b_n4_cut13_haar.json", ["F_A"]),
         FB=("d4_b_n4_cut13_haar.json", ["F_B"])),
    dict(label="N5_23_full", N=5, cut="23",
         FA=("d4_b_n5_cut23_haar_full.json", ["F_A"]),
         FB=("d4_b_n5_cut23_haar_full.json", ["F_B"])),
]


def build_env(N):
    if N == 3:
        H = get_hamiltonian("asymmetric_XY", n_sites=3)
        L = liouvillian_dephasing(H, [0.5, 1.0, 2.0])
        rho0 = d2a_states(H)["haar"]
        return H, [0.5, 1.0, 2.0], L, rho0
    H, gvec, L = build_system(N)
    return H, gvec, L, d4_states(H, N)["haar"]


def tau_c_from(dphi_map, taus):
    """Linear-interpolation crossing; returns (tau_c, bracket_dtau)."""
    for i in range(len(taus) - 1):
        t0, t1 = taus[i], taus[i + 1]
        d0, d1 = dphi_map[t0], dphi_map[t1]
        if d0 * d1 <= 0 and d1 != d0:
            tc = t0 + (t1 - t0) * (-d0) / (d1 - d0)
            return tc, (t1 - t0)
    return None, None


def main():
    e1 = json.loads((OUTDIR / "e1_tauc_mechanism.json").read_text(encoding="utf-8"))
    out = {}
    print("=== Phase E sensitivity / uncertainty supplement ===\n")
    for syscfg in SYSTEMS:
        label = syscfg["label"]
        N, cut = syscfg["N"], syscfg["cut"]
        n_a = CUTS[N][cut][0]
        H, gvec, L, rho0 = build_env(N)
        Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        L2rho0 = (L @ Y0.reshape(-1, order="F")).reshape(rho0.shape, order="F")

        F_A = load_U(OUTDIR / syscfg["FA"][0], syscfg["FA"][1])
        F_B = load_U(OUTDIR / syscfg["FB"][0], syscfg["FB"][1])
        ra = rep_data(F_A, rho0, Y0, L2rho0, N, n_a, LAMBDA)
        rb = rep_data(F_B, rho0, Y0, L2rho0, N, n_a, LAMBDA)
        dPhi0 = (rb["I"] - ra["I"]) - LAMBDA * (rb["J0"] - ra["J0"])

        # --- quadratic fit with covariance (DeltaPhi_2 + uncertainties) ---
        dq = {float(t): v for t, v in
              e1["systems"][label]["DeltaPhi_quad"].items()}
        taus_q = np.array(sorted(dq))
        yq = np.array([dq[t] for t in taus_q])
        Aq = np.vstack([taus_q ** 2, taus_q, np.ones_like(taus_q)]).T
        coef, _, _, _ = np.linalg.lstsq(Aq, yq, rcond=None)
        c2, c1, c0 = coef  # DeltaPhi_2, DeltaPhi_1, DeltaPhi_0
        resid = yq - Aq @ coef
        dof = max(len(yq) - 3, 1)
        s2 = float(resid @ resid) / dof
        cov = s2 * np.linalg.inv(Aq.T @ Aq)
        sd = np.sqrt(np.clip(np.diag(cov), 0, None))  # [sd2, sd1, sd0]

        # tau_c^(1) = -c0/c1 with first-order error propagation
        tc1 = -c0 / c1
        stc1 = abs(tc1) * np.sqrt((sd[2] / c0) ** 2 + (sd[1] / c1) ** 2)
        # tau_c^(2): positive physical root near tc1
        disc = c1 ** 2 - 4 * c0 * c2
        tc2 = None
        if disc >= 0:
            roots = [(-c1 + np.sqrt(disc)) / (2 * c2),
                     (-c1 - np.sqrt(disc)) / (2 * c2)]
            pos = [r for r in roots if r > 0]
            if pos:
                tc2 = min(pos, key=lambda r: abs(r - tc1))

        # --- representative sensitivity: alternative F_B from existing data ---
        sens = []
        fb_json = json.loads((OUTDIR / syscfg["FB"][0]).read_text(encoding="utf-8"))
        fb_dict = fb_json
        for k in syscfg["FB"][1]:
            fb_dict = fb_dict[k]
        seeds = fb_dict.get("seeds", [])
        if len(seeds) >= 2 and "U_real" in seeds[0]:
            for s_alt in seeds[1:3]:
                U_alt = np.array(s_alt["U_real"]) + 1j * np.array(s_alt["U_imag"])
                r_alt = rep_data(U_alt, rho0, Y0, L2rho0, N, n_a, LAMBDA)
                d0_alt = (r_alt["I"] - ra["I"]) - LAMBDA * (r_alt["J0"] - ra["J0"])
                d1_alt = -LAMBDA * (r_alt["J1_ana"] - ra["J1_ana"])
                if d1_alt != 0:
                    sens.append(-d0_alt / d1_alt)
        # (ii) alternate-tau optima already stored (N=3 only)
        for fname, keys in syscfg.get("alt_tau_jsons", []):
            p = OUTDIR / fname
            if not p.exists():
                continue
            U_alt = load_U(p, keys)
            r_alt = rep_data(U_alt, rho0, Y0, L2rho0, N, n_a, LAMBDA)
            d0_alt = (r_alt["I"] - ra["I"]) - LAMBDA * (r_alt["J0"] - ra["J0"])
            d1_alt = -LAMBDA * (r_alt["J1_ana"] - ra["J1_ana"])
            if d1_alt != 0:
                sens.append(-d0_alt / d1_alt)
        # (iii) d_F of the alternatives vs F_B (basin stability)
        d_alt = []
        for i in range(min(3, len(seeds))):
            if "U_real" not in seeds[i]:
                break
            d_alt.append(tps_distance(
                np.array(seeds[i]["U_real"]) + 1j * np.array(seeds[i]["U_imag"]),
                F_B, 2 ** n_a, 2 ** (N - n_a)))

        # --- tau-grid resolution ---
        meas = e1["systems"][label]["tau_c_meas"]
        grid = json.loads((OUTDIR / ("d3_timescale.json" if N == 3 else
                                     f"d4_b_n{N}_cut{cut}_haar"
                                     f"{'_full' if label == 'N5_23_full' else ''}.json"
                                     )).read_text(encoding="utf-8"))
        # use the stored DeltaPhi map
        dphi_map = {float(k): v for k, v in grid.get("DeltaPhi", {}).items()}
        dtau_res = None
        if dphi_map:
            taus = sorted(dphi_map)
            _, dtau = tau_c_from(dphi_map, taus)
            dtau_res = dtau

        out[label] = dict(
            DeltaPhi_0=float(c0), DeltaPhi_1=float(c1), DeltaPhi_2=float(c2),
            sd_DeltaPhi_0=float(sd[2]), sd_DeltaPhi_1=float(sd[1]),
            sd_DeltaPhi_2=float(sd[0]),
            tau_c1=float(tc1), sd_tau_c1=float(stc1),
            tau_c2=float(tc2) if tc2 else None,
            measured=float(meas),
            sens_tau_c1=sens,
            d_F_alt_FB=d_alt,
            tau_grid_resolution=dtau_res,
        )
        print(f"--- {label} ---")
        print(f"  dPhi0={c0:+.4f} +- {sd[2]:.4f}, dPhi1={c1:+.4f} +- {sd[1]:.4f}, "
              f"dPhi2={c2:+.4f} +- {sd[0]:.4f}")
        print(f"  tau_c^(1)={tc1:.4f} +- {stc1:.4f} (meas {meas:.4f}); "
              f"tau_c^(2)={tc2 if tc2 is None else round(tc2,4)}")
        print(f"  representative sensitivity: {[round(v,4) for v in sens]}")
        print(f"  d_F(alt seeds, F_B): {[round(v,3) for v in d_alt]}")
        print(f"  tau-grid resolution (interpolation bracket): {dtau_res}")

    (OUTDIR / "e1_sensitivity.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nSaved to analysis_output/e1_sensitivity.json")


if __name__ == "__main__":
    main()
