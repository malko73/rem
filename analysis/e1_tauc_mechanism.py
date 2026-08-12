#!/usr/bin/env python3
"""
Phase E — tau_c mechanism: does tau_c ~ -DeltaPhi_0 / DeltaPhi_1?

Master's question: WHY is tau_c ~ 0.02? With the finite-time expansion

    J_tau(F) = J_0(F) + tau J_1(F) + O(tau^2)

and DeltaPhi(tau) = Phi_tau(F_B) - Phi_tau(F_A) ~ DeltaPhi_0 + tau DeltaPhi_1,
the leading-order prediction is

    tau_c^(1) = -DeltaPhi_0 / DeltaPhi_1

(provided DeltaPhi_0 < 0 < DeltaPhi_1, the sign structure that produces a
crossing). If this reproduces the measured tau_c for the four fixed systems,
the crossover timescale becomes a DERIVED quantity (inter-basin objective
competition) rather than a numerical observation.

Systems (fixed reps from D3/D4; measured tau_c from the objective crossing):

  N=3, 2|1              : tau_c = 0.0180   (F_A: D2.2-A haar, F_B: tau=0.1 haar)
  N=4, 2|2 (balanced)   : tau_c = 0.0206
  N=4, 1|3 (asymmetric) : tau_c = 0.0224
  N=5, 2|3 full quotient: tau_c = 0.0217   (official value per master)

Gates:
  E1  J_tau = J_0 + tau J_1 + O(tau^2) holds numerically at fixed reps
  E2  crossover-generating sign structure (DeltaPhi_0 < 0, DeltaPhi_1 > 0)
      in all 4 cases
  E3  tau_c^(1) reproduces the measured tau_c (relative error within ~20%)
  E4  flatness of tau_c(N) explained by flatness of -DeltaPhi_0/DeltaPhi_1
  E5  adding the quadratic term improves the prediction systematically
  E6  inter-basin objective competition dominates over the Liouvillian gap

Analytic first-order coefficient (fixed Schmidt basis Q_F at t=0):

  C^2(tau) = |Q_F rho_U(tau)|^2,  rho_U(tau) = U† e^{tau L} rho_0 U
  C^2'(0)  = 2 Re<Q rho_U, Q Y_U>,            Y_U = U† L rho_0 U   -> J_0 = -C^2'(0)/2
  C^2''(0) = 2 |Q Y_U|^2 + 2 Re<Q rho_U, Q (U† L^2 rho_0 U)>
  J_1      = -C^2''(0)/4
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from gamma_F import schmidt_basis, dephasing_projection, liouvillian_env
from rem4_numerical import mutual_information_for_cut
from d1_hamiltonians import get_hamiltonian
from b2_canonical_open import liouvillian_dephasing
from d4_nscaling import build_system, build_states as d4_states, CUTS, LAMBDA

LAMBDA = 0.2
SMALL_TAU = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3]       # linear-regime grid
QUAD_TAU = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 2e-2, 3e-2]  # quadratic-fit window


def dominant_eigenvector(rho):
    evals, evecs = np.linalg.eigh(rho)
    return evecs[:, np.argmax(evals)]


def build_system_N3():
    H = get_hamiltonian("asymmetric_XY", n_sites=3)
    L = liouvillian_dephasing(H, [0.5, 1.0, 2.0])
    from d2_2a_unnormalized_instantaneous import build_states as d2a_states
    rho0 = d2a_states(H)["haar"]
    return H, [0.5, 1.0, 2.0], L, rho0


def load_U(path, keys):
    data = json.loads(path.read_text(encoding="utf-8"))
    rec = data
    for k in keys:
        rec = rec[k]
    return np.array(rec["best_U_real"]) + 1j * np.array(rec["best_U_imag"])


SYSTEMS = [
    dict(label="N3_21", N=3, cut="21",
         FA=("d2_2a_unnormalized_instantaneous.json", ["core", "haar"]),
         FB=("d2_2b_tau_t1e-1.json", ["core", "haar"]),
         tau_c_meas=0.0180),
    dict(label="N4_22", N=4, cut="22",
         FA=("d4_b_n4_cut22_haar.json", ["F_A"]),
         FB=("d4_b_n4_cut22_haar.json", ["F_B"]),
         tau_c_meas=0.0206),
    dict(label="N4_13", N=4, cut="13",
         FA=("d4_b_n4_cut13_haar.json", ["F_A"]),
         FB=("d4_b_n4_cut13_haar.json", ["F_B"]),
         tau_c_meas=0.0224),
    dict(label="N5_23_full", N=5, cut="23",
         FA=("d4_b_n5_cut23_haar_full.json", ["F_A"]),
         FB=("d4_b_n5_cut23_haar_full.json", ["F_B"]),
         tau_c_meas=0.0217),
]


def rep_data(F, rho0, Y0, L2rho0, N, n_a, lam):
    """All objective ingredients at a fixed rep F (U), rotated frame."""
    rho_U = F.conj().T @ rho0 @ F
    Y_U = F.conj().T @ Y0 @ F
    L2_U = F.conj().T @ L2rho0 @ F
    evals, evecs = np.linalg.eigh(rho_U)
    psi_U = evecs[:, np.argmax(evals)]
    I = mutual_information_for_cut(psi_U, N, n_a)
    basis = schmidt_basis(psi_U, N, n_a)
    q = lambda r: r - dephasing_projection(r, basis)  # noqa: E731
    Qrho = q(rho_U)
    QY = q(Y_U)
    C0_sq = float(np.linalg.norm(Qrho, ord="fro") ** 2)
    J0 = -float(np.real(np.trace(Qrho.conj().T @ QY)))
    # analytic J_1 = -C''(0)/4 with C''(0) = 2|QY|^2 + 2 Re<Qrho, Q L2_U>
    QL2 = q(L2_U)
    C2pp = 2.0 * float(np.linalg.norm(QY, ord="fro") ** 2) \
        + 2.0 * float(np.real(np.trace(Qrho.conj().T @ QL2)))
    J1_ana = -C2pp / 4.0
    return dict(I=I, J0=J0, J1_ana=J1_ana, C0_sq=C0_sq, basis=basis,
                Qrho=Qrho, rho_U=rho_U)


def main():
    results = {}
    gates = {}
    print("=== Phase E: tau_c mechanism ===")

    for syscfg in SYSTEMS:
        label = syscfg["label"]
        N, cut = syscfg["N"], syscfg["cut"]
        n_a = CUTS[N][cut][0]
        t0 = time.time()
        print(f"\n--- {label} (N={N}, cut={cut}, n_a={n_a}) ---")

        if N == 3:
            H, gvec, L, rho0 = build_system_N3()
        else:
            H, gvec, L = build_system(N)
            rho0 = d4_states(H, N)["haar"]
        Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        L2rho0 = (L @ Y0.reshape(-1, order="F")).reshape(rho0.shape, order="F")

        F_A = load_U(OUTDIR / syscfg["FA"][0], syscfg["FA"][1])
        F_B = load_U(OUTDIR / syscfg["FB"][0], syscfg["FB"][1])

        ra = rep_data(F_A, rho0, Y0, L2rho0, N, n_a, LAMBDA)
        rb = rep_data(F_B, rho0, Y0, L2rho0, N, n_a, LAMBDA)

        # J_tau on the small- and quadratic grids for both reps
        def jtau_curve(F, basis, c0_sq, j0):
            out = {}
            for tau in SMALL_TAU + QUAD_TAU:
                rho_tau = (sla.expm(tau * L) @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
                rho_U_tau = F.conj().T @ rho_tau @ F
                q = lambda r: r - dephasing_projection(r, basis)  # noqa: E731
                c = float(np.linalg.norm(q(rho_U_tau), ord="fro") ** 2)
                out[tau] = -(c - c0_sq) / (2.0 * tau) if tau > 0 else j0
            return out

        JA = jtau_curve(F_A, ra["basis"], ra["C0_sq"], ra["J0"])
        JB = jtau_curve(F_B, rb["basis"], rb["C0_sq"], rb["J0"])

        # numerical J1: linear fit J_tau ~ J0 + tau*J1 on SMALL_TAU
        def fit_slope(j0, jtau_map):
            xs = np.array(SMALL_TAU)
            ys = np.array([jtau_map[t] for t in xs]) - j0
            A = np.vstack([xs, np.ones_like(xs)]).T
            coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
            return float(coef[0]), float(coef[1])

        J1A_num, _ = fit_slope(ra["J0"], JA)
        J1B_num, _ = fit_slope(rb["J0"], JB)

        # DeltaPhi on both grids
        dphi_small = {t: (rb["I"] - LAMBDA * JB[t]) - (ra["I"] - LAMBDA * JA[t])
                      for t in SMALL_TAU}
        dphi_quad = {t: (rb["I"] - LAMBDA * JB[t]) - (ra["I"] - LAMBDA * JA[t])
                     for t in QUAD_TAU}

        # DeltaPhi_0 and DeltaPhi_1 (numerical slope + analytic)
        dPhi0 = (rb["I"] - ra["I"]) - LAMBDA * (rb["J0"] - ra["J0"])
        xs = np.array(SMALL_TAU)
        ys = np.array([dphi_small[t] for t in xs])
        A = np.vstack([xs, np.ones_like(xs)]).T
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        dPhi1_num = float(coef[0])
        dPhi1_ana = -LAMBDA * (rb["J1_ana"] - ra["J1_ana"])

        tau_c1_num = -dPhi0 / dPhi1_num
        tau_c1_ana = -dPhi0 / dPhi1_ana
        meas = syscfg["tau_c_meas"]

        # E5: quadratic prediction on QUAD_TAU
        xq = np.array(QUAD_TAU)
        yq = np.array([dphi_quad[t] for t in xq])
        Aq = np.vstack([xq ** 2, xq, np.ones_like(xq)]).T
        cq, bq, aq = np.linalg.lstsq(Aq, yq, rcond=None)[0]
        disc = bq ** 2 - 4 * aq * cq
        tau_c2 = None
        if disc >= 0:
            roots = [(-bq + np.sqrt(disc)) / (2 * cq), (-bq - np.sqrt(disc)) / (2 * cq)]
            pos = [r for r in roots if r > 0]
            tau_c2 = min(pos, key=lambda r: abs(r - meas)) if pos else None

        # E1 linearity check: |J_tau - J0 - tau*J1_ana| at tau=1e-3
        devA = abs(JA[1e-3] - ra["J0"] - 1e-3 * ra["J1_ana"])
        devB = abs(JB[1e-3] - rb["J0"] - 1e-3 * rb["J1_ana"])

        res = dict(
            N=N, cut=cut, tau_c_meas=meas,
            F_A=dict(I=ra["I"], J0=ra["J0"], J1_ana=ra["J1_ana"], J1_num=J1A_num,
                     C0_sq=ra["C0_sq"]),
            F_B=dict(I=rb["I"], J0=rb["J0"], J1_ana=rb["J1_ana"], J1_num=J1B_num,
                     C0_sq=rb["C0_sq"]),
            DeltaPhi_0=dPhi0,
            DeltaPhi_1_num=dPhi1_num, DeltaPhi_1_ana=dPhi1_ana,
            tau_c1_num=float(tau_c1_num), tau_c1_ana=float(tau_c1_ana),
            tau_c2=float(tau_c2) if tau_c2 is not None else None,
            rel_err1=float(abs(tau_c1_num - meas) / meas),
            rel_err2=float(abs(tau_c2 - meas) / meas) if tau_c2 is not None else None,
            lin_dev_A_1e3=float(devA), lin_dev_B_1e3=float(devB),
            J_tau_A={str(t): JA[t] for t in SMALL_TAU + QUAD_TAU},
            J_tau_B={str(t): JB[t] for t in SMALL_TAU + QUAD_TAU},
            DeltaPhi_small={str(t): dphi_small[t] for t in SMALL_TAU},
            DeltaPhi_quad={str(t): dphi_quad[t] for t in QUAD_TAU},
            elapsed=time.time() - t0,
        )
        results[label] = res
        print(f"  DeltaPhi_0 = {dPhi0:+.4f}, DeltaPhi_1 = {dPhi1_num:+.4f} "
              f"(ana {dPhi1_ana:+.4f})")
        print(f"  tau_c^(1) = {tau_c1_num:.4f} (ana {tau_c1_ana:.4f}) vs measured {meas:.4f} "
              f"(rel err {abs(tau_c1_num-meas)/meas*100:.1f}%)")
        print(f"  tau_c^(2) = {tau_c2 if tau_c2 is None else round(tau_c2, 4)} "
              f"(quadratic; rel err {(abs(tau_c2-meas)/meas*100 if tau_c2 else float('nan')):.1f}%)")
        print(f"  linearity dev at 1e-3: A {devA:.2e} B {devB:.2e}")

    # ── Gates ──────────────────────────────────────────────────────────────
    e1_ok = all(r["lin_dev_A_1e3"] < 0.01 and r["lin_dev_B_1e3"] < 0.01
                for r in results.values())
    gates["E1"] = dict(pass_=bool(e1_ok),
                       note="|J_tau - J0 - tau*J1_ana| at tau=1e-3 < 0.01 for both reps, all systems")

    signs = {k: dict(dPhi0=r["DeltaPhi_0"], dPhi1=r["DeltaPhi_1_num"])
             for k, r in results.items()}
    e2_ok = all(r["DeltaPhi_0"] < 0 < r["DeltaPhi_1_num"] for r in results.values())
    gates["E2"] = dict(pass_=bool(e2_ok), signs=signs,
                       note="crossover-generating sign structure DeltaPhi_0<0<DeltaPhi_1")

    rels = {k: r["rel_err1"] for k, r in results.items()}
    e3_ok = all(v <= 0.2 for v in rels.values())
    gates["E3"] = dict(pass_=bool(e3_ok), rel_errors_1=rels,
                       note="|tau_c^(1)-meas|/meas <= 0.2 (first-order tolerance)")

    preds = [r["tau_c1_num"] for r in results.values()]
    meas_vals = [r["tau_c_meas"] for r in results.values()]
    spread_pred = (max(preds) - min(preds)) / np.mean(preds)
    spread_meas = (max(meas_vals) - min(meas_vals)) / np.mean(meas_vals)
    e4_ok = spread_pred < 0.3 and abs(spread_pred - spread_meas) < 0.2
    gates["E4"] = dict(pass_=bool(e4_ok),
                       spread_pred=float(spread_pred), spread_meas=float(spread_meas),
                       tau_c1_all=preds, tau_c_meas_all=meas_vals,
                       note="flatness of -DeltaPhi_0/DeltaPhi_1 matches flat measured tau_c")

    err1 = np.mean([r["rel_err1"] for r in results.values()])
    err2_list = [r["rel_err2"] for r in results.values() if r["rel_err2"] is not None]
    e5_ok = bool(err2_list) and np.mean(err2_list) < err1
    gates["E5"] = dict(pass_=bool(e5_ok), mean_rel_err1=float(err1),
                       mean_rel_err2=float(np.mean(err2_list)) if err2_list else None,
                       tau_c2_all={k: r["tau_c2"] for k, r in results.items()},
                       note="quadratic prediction improves the mean relative error")

    # E6: inter-basin competition vs Liouvillian gap
    e6 = {}
    e6_ok = True
    for label, r in results.items():
        N = r["N"]
        if N == 3:
            H = get_hamiltonian("asymmetric_XY", n_sites=3)
            L = liouvillian_dephasing(H, [0.5, 1.0, 2.0])
        else:
            H, _, L = build_system(N)
        evals = np.linalg.eigvals(L)
        nz = np.abs(evals.real) > 1e-8
        Delta_L = float(np.min(np.abs(evals.real[nz])))
        tau_L = 1.0 / Delta_L
        err_comp = r["rel_err1"]
        err_gap = abs(tau_L - r["tau_c_meas"]) / r["tau_c_meas"]
        e6[label] = dict(Delta_L=Delta_L, tau_L=tau_L,
                         err_competition=err_comp, err_gap=err_gap,
                         gap_better=bool(err_gap < err_comp))
        if err_gap < err_comp:
            e6_ok = False
    gates["E6"] = dict(pass_=bool(e6_ok), by_system=e6,
                       note="inter-basin competition (tau_c^(1)) predicts tau_c far "
                            "better than the Liouvillian gap timescale")

    print("\n=== Phase E gates ===")
    for g, v in gates.items():
        print(f"  {g}: {'PASS' if v.get('pass_') else 'FAIL'}")

    out = dict(gates=gates, systems=results)
    out_path = OUTDIR / "e1_tauc_mechanism.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
