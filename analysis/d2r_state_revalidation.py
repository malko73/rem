#!/usr/bin/env python3
"""
Phase D2-R — State Generality revalidation under Spec v2.2 canonical.

Objective (Spec v2.2 canonical ONLY):
    J_dyn^(0)(F;rho,L) = -Re<Q_F rho, Q_F L(rho)>   (denominator-free, signed)
    Phi(F;lambda)      = I_rho(F) - lambda * J_dyn^(0)(F)

Frozen protocol — identical to old D2 (d2_state_generality.py), only the
dynamical term changes (normalized Gamma -> unnormalized J_dyn^(0)):
    - Hamiltonian: asymmetric_XY, n_sites = 3
    - States: ground / Haar / mixed (70/30) / thermal (beta=2),
      construction identical to D2.2-A (Haar rng 20260813)
    - lambda = 0.2, Adam 200 steps, lr = 0.01
    - 6 seeds, SEED0 = 20260813 (OLD D2 seed base)
    - theta0 ~ N(0,1)*1.0, horizontal_basis(4,2), n_a = 2
    - Environment: pure dephasing gamma = (0.5, 1.0, 2.0)
    - Stopping: fixed 200 steps

Old D2 is retained as pre-v2.2 preliminary evidence (superseded
normalized dynamical functional). D2.2-A is NOT a substitute for D2-R:
D2.2-A validated J_dyn^(0) as a well-posed singularity fix; D2-R
re-certifies state generality of the frozen Spec v2.2 canonical.

Gates (master's D2-R criteria):
    D2-R1  4/4 states finite optimum
    D2-R2  seed stability (std < 0.1)
    D2-R3  no product collapse (p_max < 0.99)
    D2-R4  I, J_dyn^(0), Phi non-trivial (best > contiguous, I > 0.1,
           |J_dyn| > 1e-3)
    D2-R5  F* responds to STATE: cross-eval matrix M_ij = Phi_{rho_i}
           (F*_{rho_j}) diagonal dominance + d_F matrix (objective AND
           geometric response); mean gap >= 0.05 => response;
           < 0.05 => same-attractor concentration (reported either way)
    D2-R6  no normalized singularity recurrence (min C_F^2 > 1e-4,
           max|gamma_ref| < 50, max|Phi| < 10)
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
from b2_canonical_open import liouvillian_dephasing
from b2_1_cross_eval import tps_distance
from d2_2a_unnormalized_instantaneous import unnormalized_open_objective, build_states

LAMBDA = 0.2
STEPS = 200
LR = 0.01
SEEDS = 6
SEED0 = 20260813  # OLD D2 seed base (frozen)
GAMMA_VEC = [0.5, 1.0, 2.0]
HNAME = "asymmetric_XY"
N = 3
N_A = 2
STATES = ["ground", "haar", "mixed", "thermal"]

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


def adam_optimize_open(U0, rho0, Y0, H_basis, lambda_val, steps=STEPS, lr=LR):
    """Adam optimization of Phi = I - lambda*J_dyn^(0) (frozen D2 protocol).

    rho0: original state; Y0 = L(rho0) precomputed (pull-back rotates both).
    """
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_at(d):
        U = to_U(d)
        rho_U = U.conj().T @ rho0 @ U
        Y_U = U.conj().T @ Y0 @ U
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
            rho_U = best_U.conj().T @ rho0 @ best_U
            Y_U = best_U.conj().T @ Y0 @ best_U
            best_obj, _ = j0_objective_lambda(rho_U, Y_U, N_A, lambda_val)

    if best_obj is None:
        rho_U = best_U.conj().T @ rho0 @ best_U
        Y_U = best_U.conj().T @ Y0 @ best_U
        best_obj, _ = j0_objective_lambda(rho_U, Y_U, N_A, lambda_val)

    return best_U, best_phi, best_obj


def evaluate_U(U, rho0, Y0, lambda_val):
    rho_U = U.conj().T @ rho0 @ U
    Y_U = U.conj().T @ Y0 @ U
    obj, phi = j0_objective_lambda(rho_U, Y_U, N_A, lambda_val)
    return dict(Phi=phi, I=obj["I"], J_dyn=obj["J_dyn"],
                gamma_ref=obj["gamma_ref"], C_F_sq=obj["C_F_sq"])


def schmidt_spectrum(psi, n=3, cut=2):
    d_a, d_b = 2 ** cut, 2 ** (n - cut)
    mat = psi.reshape((d_a, d_b))
    _, s, _ = np.linalg.svd(mat)
    return s ** 2


def dominant_eigenvector(rho):
    evals, evecs = np.linalg.eigh(rho)
    return evecs[:, np.argmax(evals)]


def main():
    OUTDIR.mkdir(exist_ok=True)
    H_basis = horizontal_basis(4, 2)
    H = get_hamiltonian(HNAME, n_sites=N)
    L = liouvillian_dephasing(H, GAMMA_VEC)
    states = build_states(H)

    print("=== D2-R: State Generality (Spec v2.2 canonical) ===")
    print("J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)> ; Phi = I - lambda*J_dyn^(0)")
    print(f"Hamiltonian: {HNAME}, env: γ={GAMMA_VEC}, λ={LAMBDA}, "
          f"Adam {STEPS} steps, lr={LR}, seeds={SEEDS}, SEED0={SEED0}\n")

    Y0s = {}
    for st in STATES:
        rho0 = states[st]
        Y0s[st] = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")

    results = {}
    for st in STATES:
        t0 = time.time()
        rho0 = states[st]
        Y0 = Y0s[st]
        cont = evaluate_U(np.eye(8, dtype=complex), rho0, Y0, LAMBDA)
        print(f"--- State: {st} ---")
        print(f"  contiguous Φ(cut2) = {cont['Phi']:.4f}")

        seeds_data = []
        for s in range(SEEDS):
            seed = SEED0 + s
            rng = np.random.default_rng(seed)
            theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
            U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

            U_best, phi_best, obj = adam_optimize_open(U0, rho0, Y0, H_basis, LAMBDA)
            ev = evaluate_U(U_best, rho0, Y0, LAMBDA)
            psi0 = dominant_eigenvector(rho0)
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

        results[st] = dict(
            contiguous=cont,
            seeds=seeds_data,
            best_Phi=best,
            median_Phi=float(np.median(phis)),
            std_Phi=float(phis.std()),
            success=best > cont["Phi"],
        )

    # ── D2-R5: cross-evaluation matrix M_ij = Phi_{rho_i}(F*_{rho_j}) ─────
    best_U = {}
    for st in STATES:
        b = max(results[st]["seeds"], key=lambda d: d["Phi"])
        best_U[st] = np.array(b["U_real"]) + 1j * np.array(b["U_imag"])

    cross = {}
    for si in STATES:
        cross[si] = {}
        for sj in STATES:
            ev = evaluate_U(best_U[sj], states[si], Y0s[si], LAMBDA)
            cross[si][sj] = ev["Phi"]

    dF = {}
    for si in STATES:
        dF[si] = {}
        for sj in STATES:
            dF[si][sj] = tps_distance(best_U[si], best_U[sj])

    # within-state seed spread (mean pairwise d_F among the 6 seeds' U*)
    within = {}
    for st in STATES:
        us = [np.array(d["U_real"]) + 1j * np.array(d["U_imag"])
              for d in results[st]["seeds"]]
        ds = [tps_distance(us[a], us[b])
              for a in range(SEEDS) for b in range(a + 1, SEEDS)]
        within[st] = float(np.mean(ds))

    gaps = {}
    for si in STATES:
        diag = cross[si][si]
        off = max(cross[si][sj] for sj in STATES if sj != si)
        gaps[si] = diag - off
    mean_gap = float(np.mean(list(gaps.values())))
    response_class = "response" if mean_gap >= R5_GAP_CAP else "same-attractor"

    # ── Gates ─────────────────────────────────────────────────────────────
    gates = {}

    r1_ok = all(np.isfinite(r["best_Phi"]) for r in results.values())
    gates["D2-R1"] = dict(pass_=bool(r1_ok),
                          best_by_state={s: r["best_Phi"] for s, r in results.items()})

    r2_ok = all(r["std_Phi"] < R2_STD_CAP for r in results.values())
    gates["D2-R2"] = dict(pass_=bool(r2_ok),
                          std_by_state={s: r["std_Phi"] for s, r in results.items()},
                          criterion=f"std < {R2_STD_CAP}")

    pmax = {s: max(d["p_max"] for d in r["seeds"]) for s, r in results.items()}
    r3_ok = all(v < R3_PMAX_CAP for v in pmax.values())
    gates["D2-R3"] = dict(pass_=bool(r3_ok), p_max_best=pmax,
                          criterion=f"p_max < {R3_PMAX_CAP}")

    r4_detail = {}
    r4_ok = True
    for s, r in results.items():
        b = max(r["seeds"], key=lambda d: d["Phi"])
        ok = (r["success"] and b["I"] > R4_I_MIN and abs(b["J_dyn"]) > R4_J_MIN)
        r4_detail[s] = dict(best_Phi=b["Phi"], I=b["I"], J_dyn=b["J_dyn"],
                            contiguous=r["contiguous"]["Phi"], ok=bool(ok))
        r4_ok = r4_ok and ok
    gates["D2-R4"] = dict(pass_=bool(r4_ok), by_state=r4_detail,
                          criterion="best Phi > contiguous and I > 0.1 and |J_dyn| > 1e-3")

    gates["D2-R5"] = dict(
        pass_=bool(mean_gap >= R5_GAP_CAP),
        classification=response_class,
        mean_gap=mean_gap,
        gap_by_state=gaps,
        cross_eval=cross,
        d_F_best=dict({si: {sj: round(dF[si][sj], 4) for sj in STATES} for si in STATES}),
        within_state_d_F=within,
        note=("gap >= 0.05 => F* responds to state; gap < 0.05 => all states "
              "fall to the same attractor; d_F gives the geometric response"),
    )

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
    gates["D2-R6"] = dict(pass_=bool(r6_ok), stats=r6,
                          criterion="min C_F^2 > 1e-4, |gamma_ref| < 50, |Phi| < 10")

    print("\n=== Gate verdicts (D2-R) ===")
    for g, v in gates.items():
        status = "PASS" if v.get("pass_") else "FAIL"
        print(f"  {g}: {status}")

    print(f"\n=== D2-R5 (response) ===")
    print(f"  mean gap = {mean_gap:.4f} -> classification: {response_class}")
    print("  cross-eval matrix (rows=objective rho_i, cols=F*_j):")
    print("  " + " ".join(f"{sj:>10s}" for sj in STATES))
    for si in STATES:
        row = " ".join(f"{cross[si][sj]:10.4f}" for sj in STATES)
        print(f"  {si:>10s} {row}")
    print("  within-state d_F:", {s: round(v, 3) for s, v in within.items()})

    results_meta = dict(
        protocol=dict(
            functional="J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)> (Spec v2.2 canonical)",
            frozen="old D2 protocol: asymmetric_XY, 4 states, λ=0.2, Adam 200 steps "
                   "lr=0.01, 6 seeds, SEED0=20260813, dephasing γ=(0.5,1,2), "
                   "horizontal_basis(4,2), n_a=2",
            old_D2_status="pre-v2.2 preliminary evidence using the superseded "
                          "normalized dynamical functional (d2_state_generality.json)",
            d22a_relation=("D2.2-A validated J_dyn^(0) as a well-posed singularity fix; "
                           "D2-R re-certifies state generality of the frozen canonical"),
        ),
        states=results,
        cross_eval=cross,
        d_F_best=dF,
        within_state_d_F=within,
        gates=gates,
    )

    out_path = OUTDIR / "d2r_state_revalidation.json"
    out_path.write_text(json.dumps(results_meta, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
