#!/usr/bin/env python3
"""
Phase D1-R — Hamiltonian Generality revalidation under Spec v2.2 canonical.

Objective (Spec v2.2 canonical ONLY):
    J_dyn^(0)(F;rho,L) = -Re<Q_F rho, Q_F L(rho)>   (denominator-free, signed)
    Phi(F;lambda)      = I_rho(F) - lambda * J_dyn^(0)(F)

Frozen protocol — identical to old D1 (d1_generality.py), only the
dynamical term changes (normalized Gamma -> unnormalized J_dyn^(0)):
    - 5 Hamiltonian families (same params)
    - lambda = 0.2, Adam 200 steps, lr = 0.01
    - 6 seeds, SEED0 = 20260812 (OLD D1 seed base, NOT 20260813)
    - theta0 ~ N(0,1)*1.0, horizontal_basis(4,2), n_a = 2
    - Environment: pure dephasing gamma = (0.5, 1.0, 2.0)
    - Stopping: fixed 200 steps (same as D1)

Old D1 is retained as:
    pre-v2.2 preliminary evidence using the superseded normalized
    dynamical functional

Gates (master's D1-R criteria):
    D1-R1  5/5 families finite optimum
    D1-R2  seed-stable solutions (std < 0.1)
    D1-R3  no product collapse (p_max < 0.99)
    D1-R4  J_dyn^(0), I, Phi non-trivial (best Phi > contiguous, I > 0.1,
           |J_dyn| > 1e-3 at optimum)
    D1-R5  F* responds appropriately to Hamiltonian change: cross-eval
           diagonal dominance Phi_{H_i}(F*_i) >= Phi_{H_i}(F*_j); mean
           gap >= 0.05 => response; < 0.05 => concentration (state/function
           dominance) — both reported
    D1-R6  no D2.1-type singularity in any family (min C_F^2 > 1e-4,
           max|gamma_ref| < 50 at optima, max|Phi| < 10)
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

from d1_hamiltonians import get_hamiltonian, get_ground_state
from quotient_geometry import horizontal_basis
from b2_canonical_open import apply_tps, liouvillian_dephasing
from b2_1_cross_eval import tps_distance
from d2_2a_unnormalized_instantaneous import unnormalized_open_objective

LAMBDA = 0.2
STEPS = 200
LR = 0.01
SEEDS = 6
SEED0 = 20260812  # OLD D1 seed base (frozen)
GAMMA_VEC = [0.5, 1.0, 2.0]
N = 3
N_A = 2

HAMILTONIANS = {
    "asymmetric_XY": {},
    "transverse_field_ising": {"J": 1.0, "h": 0.5},
    "heisenberg_xxz": {"J": 1.0, "Delta": 0.5},
    "xyz": {"J_x": 1.0, "J_y": 0.8, "J_z": 0.6},
    "random_local": {"seed": 42},
}

# Gate thresholds (pre-registered)
R2_STD_CAP = 0.1
R3_PMAX_CAP = 0.99
R4_I_MIN = 0.1
R4_J_MIN = 1e-3
R5_GAP_CAP = 0.05
R6_MIN_C0 = 1e-4
R6_MAX_ABS_GR = 50.0
R6_PHI_CAP = 10.0


def j0_objective_lambda(rho_U, Y_U, n_a, lambda_val):
    obj = unnormalized_open_objective(rho_U, Y_U, n_a)
    phi = obj["I"] - lambda_val * obj["J_dyn"] if np.isfinite(obj["J_dyn"]) else -1e9
    return obj, phi


def adam_optimize_open(U0, L, psi0, lambda_val, H_basis, steps=STEPS, lr=LR):
    """Adam optimization of Phi = I - lambda*J_dyn^(0) (frozen D1 protocol)."""
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_at(d):
        U = to_U(d)
        rho_U, Y_U = apply_tps(psi0, L, U)
        obj, phi = j0_objective_lambda(rho_U, Y_U, N_A, lambda_val)
        return phi

    best_phi = phi_at(delta)
    best_U = to_U(delta)
    best_obj = None

    for step in range(steps):
        grad = np.zeros_like(delta)
        for k in range(len(delta)):
            d_p = delta.copy()
            d_m = delta.copy()
            d_p[k] += 1e-5
            d_m[k] -= 1e-5
            grad[k] = (phi_at(d_p) - phi_at(d_m)) / (2 * 1e-5)

        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        delta = delta + lr * m_hat / (np.sqrt(v_hat) + eps)

        p = phi_at(delta)
        if p > best_phi:
            best_phi = p
            best_U = to_U(delta)
            best_obj, _ = j0_objective_lambda(*apply_tps(psi0, L, best_U), N_A, lambda_val)

    if best_obj is None:
        best_obj, _ = j0_objective_lambda(*apply_tps(psi0, L, best_U), N_A, lambda_val)

    return best_U, best_phi, best_obj


def evaluate_U(U, psi0, L, lambda_val):
    rho_U, Y_U = apply_tps(psi0, L, U)
    obj, phi = j0_objective_lambda(rho_U, Y_U, N_A, lambda_val)
    return dict(Phi=phi, I=obj["I"], J_dyn=obj["J_dyn"],
                gamma_ref=obj["gamma_ref"], C_F_sq=obj["C_F_sq"])


def schmidt_spectrum(psi, n=3, cut=2):
    d_a, d_b = 2 ** cut, 2 ** (n - cut)
    mat = psi.reshape((d_a, d_b))
    _, s, _ = np.linalg.svd(mat)
    return s ** 2


def main():
    OUTDIR.mkdir(exist_ok=True)
    H_basis = horizontal_basis(4, 2)
    results = {}

    print("=== D1-R: Hamiltonian Generality (Spec v2.2 canonical) ===")
    print("J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)> ; Phi = I - lambda*J_dyn^(0)")
    print(f"Environment: pure dephasing γ={GAMMA_VEC}, λ={LAMBDA}, "
          f"Adam {STEPS} steps, lr={LR}, seeds={SEEDS}, SEED0={SEED0}\n")

    for hname, hparams in HAMILTONIANS.items():
        t0 = time.time()
        H = get_hamiltonian(hname, n_sites=N, **hparams)
        E0, psi0 = get_ground_state(H)
        L = liouvillian_dephasing(H, GAMMA_VEC)

        cont = evaluate_U(np.eye(8, dtype=complex), psi0, L, LAMBDA)
        print(f"--- {hname} ---")
        print(f"  E0 = {E0:.4f}, contiguous Φ(cut2) = {cont['Phi']:.4f}")

        seeds_data = []
        for s in range(SEEDS):
            seed = SEED0 + s
            rng = np.random.default_rng(seed)
            theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
            U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

            U_best, phi_best, obj = adam_optimize_open(U0, L, psi0, LAMBDA, H_basis)
            ev = evaluate_U(U_best, psi0, L, LAMBDA)
            psi_U = U_best.conj().T @ psi0
            p_max = float(schmidt_spectrum(psi_U).max())

            seeds_data.append(dict(
                seed=seed, Phi=ev["Phi"], I=ev["I"], J_dyn=ev["J_dyn"],
                gamma_ref=ev["gamma_ref"], C_F_sq=ev["C_F_sq"], p_max=p_max,
                U_real=U_best.real.tolist(), U_imag=U_best.imag.tolist(),
            ))

        phis = np.array([d["Phi"] for d in seeds_data])
        best = float(phis.max())
        print(f"  best Φ = {best:.4f}, median = {np.median(phis):.4f}, "
              f"std = {phis.std():.4f}  ({time.time()-t0:.1f}s)")

        results[hname] = dict(
            E0=float(E0),
            contiguous=cont,
            seeds=seeds_data,
            best_Phi=best,
            median_Phi=float(np.median(phis)),
            std_Phi=float(phis.std()),
            success=best > cont["Phi"],
        )

    # ── D1-R5: cross-evaluation matrix and pairwise d_F ───────────────────
    names = list(HAMILTONIANS)
    best_U = {}
    for hname in names:
        b = max(results[hname]["seeds"], key=lambda d: d["Phi"])
        best_U[hname] = np.array(b["U_real"]) + 1j * np.array(b["U_imag"])

    # objective of family i evaluated at family j's best F*
    cross = {}
    for i, hi in enumerate(names):
        H_i = get_hamiltonian(hi, n_sites=N, **HAMILTONIANS[hi])
        _, psi_i = get_ground_state(H_i)
        L_i = liouvillian_dephasing(H_i, GAMMA_VEC)
        cross[hi] = {}
        for hj in names:
            U_j = best_U[hj]
            ev = evaluate_U(U_j, psi_i, L_i, LAMBDA)
            cross[hi][hj] = ev["Phi"]

    # pairwise d_F among best-of-seeds F*
    dF = {}
    for i, hi in enumerate(names):
        dF[hi] = {}
        for j, hj in enumerate(names):
            dF[hi][hj] = tps_distance(best_U[hi], best_U[hj])

    # within-family spread: mean pairwise d_F among the 6 seeds' U*
    within = {}
    for hname in names:
        us = [np.array(d["U_real"]) + 1j * np.array(d["U_imag"])
              for d in results[hname]["seeds"]]
        ds = [tps_distance(us[a], us[b])
              for a in range(SEEDS) for b in range(a + 1, SEEDS)]
        within[hname] = float(np.mean(ds))

    # diagonal gaps
    gaps = {}
    for hi in names:
        diag = cross[hi][hi]
        off = max(cross[hi][hj] for hj in names if hj != hi)
        gaps[hi] = diag - off
    mean_gap = float(np.mean(list(gaps.values())))
    response_class = "response" if mean_gap >= R5_GAP_CAP else "concentration"

    # ── Gates ─────────────────────────────────────────────────────────────
    gates = {}

    # R1: finite optima 5/5
    r1_ok = all(np.isfinite(r["best_Phi"]) for r in results.values())
    gates["D1-R1"] = dict(pass_=bool(r1_ok),
                          best_by_family={h: r["best_Phi"] for h, r in results.items()})

    # R2: seed stability
    r2_ok = all(r["std_Phi"] < R2_STD_CAP for r in results.values())
    gates["D1-R2"] = dict(pass_=bool(r2_ok),
                          std_by_family={h: r["std_Phi"] for h, r in results.items()},
                          criterion=f"std < {R2_STD_CAP}")

    # R3: no product collapse
    pmax = {h: max(d["p_max"] for d in r["seeds"]) for h, r in results.items()}
    r3_ok = all(v < R3_PMAX_CAP for v in pmax.values())
    gates["D1-R3"] = dict(pass_=bool(r3_ok), p_max_best=pmax,
                          criterion=f"p_max < {R3_PMAX_CAP}")

    # R4: non-trivial J, I, Phi
    r4_detail = {}
    r4_ok = True
    for h, r in results.items():
        b = max(r["seeds"], key=lambda d: d["Phi"])
        ok = (r["success"] and b["I"] > R4_I_MIN and abs(b["J_dyn"]) > R4_J_MIN)
        r4_detail[h] = dict(best_Phi=b["Phi"], I=b["I"], J_dyn=b["J_dyn"],
                            contiguous=r["contiguous"]["Phi"], ok=bool(ok))
        r4_ok = r4_ok and ok
    gates["D1-R4"] = dict(pass_=bool(r4_ok), by_family=r4_detail,
                          criterion="best Phi > contiguous and I > 0.1 and |J_dyn| > 1e-3")

    # R5: F* responds to Hamiltonian
    gates["D1-R5"] = dict(
        pass_=bool(mean_gap >= R5_GAP_CAP),
        classification=response_class,
        mean_gap=mean_gap,
        gap_by_family=gaps,
        cross_eval=cross,
        d_F_best=dict({hi: {hj: round(dF[hi][hj], 4) for hj in names} for hi in names}),
        within_family_d_F=within,
        note=("gap >= 0.05 => F* responds to H; gap < 0.05 => families "
              "concentrate (state/function side dominates selection)"),
    )

    # R6: no D2.1-type singularity
    all_c0 = [d["C_F_sq"] for r in results.values() for d in r["seeds"]]
    all_gr = [d["gamma_ref"] for r in results.values() for d in r["seeds"]
              if np.isfinite(d["gamma_ref"])]
    all_phi = [d["Phi"] for r in results.values() for d in r["seeds"]]
    r6 = dict(
        min_C_F_sq=float(min(all_c0)),
        max_abs_gamma_ref=float(max(abs(g) for g in all_gr)) if all_gr else float("nan"),
        max_abs_Phi=float(max(abs(p) for p in all_phi)),
    )
    r6_ok = (r6["min_C_F_sq"] > R6_MIN_C0 and r6["max_abs_gamma_ref"] < R6_MAX_ABS_GR
             and r6["max_abs_Phi"] < R6_PHI_CAP)
    gates["D1-R6"] = dict(pass_=bool(r6_ok), stats=r6,
                          criterion="min C_F^2 > 1e-4, |gamma_ref| < 50, |Phi| < 10")

    print("\n=== Gate verdicts (D1-R) ===")
    for g, v in gates.items():
        status = "PASS" if v.get("pass_") else "FAIL"
        print(f"  {g}: {status}")

    print(f"\n=== D1-R5 (response) ===")
    print(f"  mean gap = {mean_gap:.4f} -> classification: {response_class}")
    print("  cross-eval matrix (rows=objective H_i, cols=F*_j):")
    print("  " + " ".join(f"{hj[:12]:>14s}" for hj in names))
    for hi in names:
        row = " ".join(f"{cross[hi][hj]:14.4f}" for hj in names)
        print(f"  {hi[:12]:>14s} {row}")
    print("  within-family d_F:", {h: round(v, 3) for h, v in within.items()})

    results_meta = dict(
        protocol=dict(
            functional="J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)> (Spec v2.2 canonical)",
            frozen="old D1 protocol: 5 families, λ=0.2, Adam 200 steps lr=0.01, "
                   "6 seeds, SEED0=20260812, dephasing γ=(0.5,1,2), horizontal_basis(4,2), n_a=2",
            old_D1_status="pre-v2.2 preliminary evidence using the superseded "
                          "normalized dynamical functional (d1_hamiltonian_generality.json)",
        ),
        families=results,
        cross_eval=cross,
        d_F_best=dF,
        within_family_d_F=within,
        gates=gates,
    )

    out_path = OUTDIR / "d1r_canonical_revalidation.json"
    out_path.write_text(json.dumps(results_meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
