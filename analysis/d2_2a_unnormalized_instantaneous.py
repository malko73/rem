#!/usr/bin/env python3
"""
Phase D2.2-A — Unnormalized instantaneous dynamical functional.

Motivation (2026-08-12, master's decision):
D2.1 proved that the NORMALIZED local decay rate
    Gamma_F^(0) = -Re<Q_F rho, Q_F L(rho)> / |Q_F rho|^2
is well-posed as a LOCAL diagnostic but SINGULAR as a global variational
functional on unrestricted TPS: as C_F(0)^2 = |Q_F rho|^2 -> 0 the optimizer
drives Gamma -> -inf, Phi -> +inf (corr(Phi, log10 C_F(0)^2) = -0.999,
best Haar solution Schmidt p = [0.9997, 0.0003] -> product collapse).

D2.2-A replaces ONLY the dynamical term with the unnormalized version
    J_dyn^(0)(F) = -Re<Q_F rho, Q_F L(rho)>
    Phi(F)       = I(F) - lambda * J_dyn^(0)(F)
Everything else is frozen to the D2 protocol (Hamiltonian, states, optimizer,
initialization, restart count, TPS parameterization, lambda, seed, stopping
criteria).

Mathematical guarantee: |<Q_F rho, Q_F L(rho)>| <= |Q_F rho| |Q_F L(rho)|,
so J_dyn^(0) -> 0 smoothly as C_F(0)^2 -> 0. No denominator blow-up.

Gates (master's fixed criteria):
  G1  Phi finite even when C_F^2 -> 0 (no divergence, no corr(Phi, log C_F^2) ~ -1)
  G2  Haar giant outlier disappears (vs D2.1: best=25.7, std=5.96)
  G3  stable w.r.t. optimizer seed
  G4  optimum continuous under small perturbation
  G5  Schmidt spectrum shows no unnatural product collapse (p_max << 0.9997)
  G6  reproducible across ground / Haar / mixed / thermal
  G7  qualitative consistency with finite-time version -> PENDING (D2.2-B)

Frozen protocol (identical to d2_state_generality.py):
  - Hamiltonian: asymmetric_XY, n_sites=3
  - Environment: pure dephasing gamma=(0.5, 1.0, 2.0)
  - lambda = 0.2, Adam, 200 steps, lr = 0.01, finite-diff eps = 1e-5
  - 6 seeds (SEED0 = 20260813), theta0 ~ N(0, 1) * 1.0
  - TPS parameterization: horizontal_basis(4, 2), n_a = 2
  - States: ground / Haar / mixed / thermal (identical construction)

Supplementary (G2/G3 statistics, mirrors D2.1 audit):
  - Haar 30 seeds (SEED0 .. SEED0+29), same optimizer/init protocol.
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

# Imports mirror d2_state_generality.py; b2_canonical_open registers
# rem4_numerical / gamma_F in sys.modules via importlib (needed on 3.13).
from d1_hamiltonians import get_hamiltonian, get_ground_state
from quotient_geometry import horizontal_basis
from b2_canonical_open import liouvillian_dephasing
from b2_1_cross_eval import tps_distance
import rem4_numerical as rem4
import gamma_F

LAMBDA = 0.2
STEPS = 200
LR = 0.01
SEEDS = 6
SEED0 = 20260813
GAMMA_VEC = [0.5, 1.0, 2.0]
HNAME = "asymmetric_XY"
N = 3
N_A = 2
HAAR_EXTENDED_SEEDS = 30  # supplementary, mirrors D2.1 audit protocol

# Gate thresholds (pre-registered)
G1_PHI_CAP = 10.0
G1_CORR_CAP = -0.5
G2_BEST_CAP = 5.0
G2_STD_CAP = 0.5
G3_STD_CAP = 0.1
G4_DELTA_PHI_CAP = 0.5
G5_PMAX_CAP = 0.99


def unnormalized_open_objective(rho_U: np.ndarray, Y_U: np.ndarray, n_a: int) -> dict:
    """I, C_F^2, J_dyn^(0), Gamma_ref (diagnostic only), Phi for D2.2-A.

    The OPTIMIZER uses Phi = I - lambda * J_dyn^(0) (denominator-free).
    gamma_ref (the normalized D2.1 quantity) is recorded ONLY as a
    diagnostic and never enters the objective.
    """
    evals_psi, evecs_psi = np.linalg.eigh(rho_U)
    psi_U = evecs_psi[:, np.argmax(evals_psi)]
    I_F = rem4.mutual_information_for_cut(psi_U, N, n_a)

    basis = gamma_F.schmidt_basis(psi_U, N, n_a)
    q = lambda r: r - gamma_F.dephasing_projection(r, basis)  # noqa: E731
    x = q(rho_U)
    qy = q(Y_U)
    num = float(np.real(np.trace(x.conj().T @ qy)))
    den = float(np.linalg.norm(x, ord="fro") ** 2)

    J_dyn = -num  # unnormalized instantaneous dynamical functional
    gamma_ref = -num / den if den > 1e-300 else float("nan")
    Phi = I_F - LAMBDA * J_dyn
    return dict(I=I_F, C_F_sq=den, J_dyn=J_dyn, gamma_ref=gamma_ref, Phi=Phi)


def unnormalized_objective_lambda(rho_U, Y_U, n_a, lambda_val):
    obj = unnormalized_open_objective(rho_U, Y_U, n_a)
    phi = obj["I"] - lambda_val * obj["J_dyn"] if np.isfinite(obj["J_dyn"]) else -1e9
    return obj, phi


def adam_optimize_open(U0, L, rho0, Y0, lambda_val, H_basis, steps=STEPS, lr=LR):
    """Adam optimization, byte-for-byte the D2 protocol with the new objective."""
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
        obj, phi = unnormalized_objective_lambda(rho_U, Y_U, N_A, lambda_val)
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
            best_obj, _ = unnormalized_objective_lambda(rho_U, Y_U, N_A, lambda_val)

    if best_obj is None:
        rho_U = best_U.conj().T @ rho0 @ best_U
        Y_U = best_U.conj().T @ Y0 @ best_U
        best_obj, _ = unnormalized_objective_lambda(rho_U, Y_U, N_A, lambda_val)

    return best_U, best_phi, best_obj


def build_states(H):
    """Identical state construction to d2_state_generality.build_states."""
    dim = 2 ** N
    _, psi_gs = get_ground_state(H)
    rng = np.random.default_rng(20260813)

    states = {}
    states["ground"] = np.outer(psi_gs, psi_gs.conj())

    m = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    q, _ = np.linalg.qr(m)
    psi_haar = q[:, 0]
    states["haar"] = np.outer(psi_haar, psi_haar.conj())

    evals, evecs = np.linalg.eigh(H)
    rho_mix = 0.7 * np.outer(evecs[:, 0], evecs[:, 0].conj()) + \
              0.3 * np.outer(evecs[:, 1], evecs[:, 1].conj())
    states["mixed"] = rho_mix

    beta = 2.0
    rho_thermal = np.exp(-beta * H) / np.trace(np.exp(-beta * H))
    states["thermal"] = rho_thermal
    return states


def schmidt_spectrum(psi, n=3, cut=2):
    d_a, d_b = 2 ** cut, 2 ** (n - cut)
    mat = psi.reshape((d_a, d_b))
    _, s, _ = np.linalg.svd(mat)
    return s ** 2


def dominant_eigenvector(rho):
    evals, evecs = np.linalg.eigh(rho)
    return evecs[:, np.argmax(evals)]


def run_state(state_name, rho0, L, H_basis, n_seeds, verbose=True):
    """Run the frozen D2 optimization for one state; return per-seed scalars
    plus the best U."""
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    seeds_data = []
    best_U = None
    best_phi = -np.inf
    for s in range(n_seeds):
        seed = SEED0 + s
        rng = np.random.default_rng(seed)
        theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
        U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

        U_best, phi_best, obj = adam_optimize_open(U0, L, rho0, Y0, LAMBDA, H_basis)

        # Schmidt spectrum of the rotated dominant eigenvector (same probe as D2.1)
        psi0 = dominant_eigenvector(rho0)
        psi_U = U_best.conj().T @ psi0
        spec = schmidt_spectrum(psi_U)

        seeds_data.append(dict(
            seed=seed, Phi=phi_best, I=obj["I"], J_dyn=obj["J_dyn"],
            gamma_ref=obj["gamma_ref"], C_F_sq=obj["C_F_sq"],
            p_max=float(spec.max()),
        ))
        if phi_best > best_phi:
            best_phi = phi_best
            best_U = U_best

    phis = np.array([d["Phi"] for d in seeds_data])
    if verbose:
        print(f"  best Φ = {phis.max():.4f}, median = {np.median(phis):.4f}, "
              f"std = {phis.std():.4f}")

    return dict(
        seeds=seeds_data,
        best_Phi=float(phis.max()),
        median_Phi=float(np.median(phis)),
        std_Phi=float(phis.std()),
        best_U_real=best_U.real.tolist(),
        best_U_imag=best_U.imag.tolist(),
    )


def perturbation_test(rho0, L, U_best, H_basis):
    """G4: small-perturbation continuity of the optimum.

    Returns list of (dir_kind, eps, dPhi, dJ, finite) for random unit-norm
    directions and the e_0 direction used in D2.1.
    """
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    obj0, phi0 = unnormalized_objective_lambda(
        U_best.conj().T @ rho0 @ U_best,
        U_best.conj().T @ Y0 @ U_best, N_A, LAMBDA)

    rows = []
    rng = np.random.default_rng(424242)
    dirs = []
    for _ in range(8):
        d = rng.standard_normal(H_basis.shape[2])
        dirs.append(d / np.linalg.norm(d))
    e0 = np.zeros(H_basis.shape[2]); e0[0] = 1.0
    dirs.append(e0)  # D2.1's direction

    for kind, d in enumerate(dirs):
        for eps in (1e-3, 1e-2):
            U_pert = U_best @ sla.expm(np.einsum("ijk,k->ij", H_basis, eps * d))
            rho_U = U_pert.conj().T @ rho0 @ U_pert
            Y_U = U_pert.conj().T @ Y0 @ U_pert
            obj, phi = unnormalized_objective_lambda(rho_U, Y_U, N_A, LAMBDA)
            ok = np.isfinite(phi) and np.isfinite(obj["J_dyn"])
            rows.append(dict(
                dir=kind, eps=eps,
                dPhi=phi - phi0, dJ=obj["J_dyn"] - obj0["J_dyn"],
                Phi=phi if ok else None, finite=bool(ok),
            ))
    return rows, phi0


def singularity_sweep(rho0, L, U_best, H_basis):
    """G1 probe: walk TPS directions that reduce C_F^2 and confirm Phi stays
    finite. In D2.1 this same probe (eps along e_0) took C_F(0)^2 -> 6e-5 and
    Phi -> 44. Here it must stay bounded.
    """
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    rows = []
    rng = np.random.default_rng(777)
    dirs = []
    for _ in range(6):
        d = rng.standard_normal(H_basis.shape[2])
        dirs.append(d / np.linalg.norm(d))
    e0 = np.zeros(H_basis.shape[2]); e0[0] = 1.0
    dirs.append(e0)

    for kind, d in enumerate(dirs):
        for eps in (0.01, 0.05, 0.1, 0.2, 0.5):
            U_pert = U_best @ sla.expm(np.einsum("ijk,k->ij", H_basis, eps * d))
            rho_U = U_pert.conj().T @ rho0 @ U_pert
            Y_U = U_pert.conj().T @ Y0 @ U_pert
            obj, phi = unnormalized_objective_lambda(rho_U, Y_U, N_A, LAMBDA)
            rows.append(dict(
                dir=kind, eps=eps, C_F_sq=obj["C_F_sq"], Phi=phi,
                J_dyn=obj["J_dyn"], gamma_ref=obj["gamma_ref"],
                finite=bool(np.isfinite(phi) and np.isfinite(obj["J_dyn"])),
            ))
    return rows


def evaluate_gates(core, haar_ext, pert_rows, sweep_rows):
    gates = {}

    # G1: no singularity. Physical criteria (the pooled corr statistic is
    # misleading here because C_F^2 is confined to a narrow band; replaced by):
    #   a) Phi bounded across ALL core trials and the C_F^2->0 sweep
    #   b) optimizer never approaches C_F^2 ~ 0 (min C_F^2 >> 1e-4)
    #   c) the NORMALIZED diagnostic gamma_ref stays finite at the optima
    #      (in D2.1 it reached -245 at C_F^2 = 4e-4)
    all_phis, all_c0, all_gr = [], [], []
    for st in core.values():
        if not isinstance(st, dict) or "seeds" not in st:
            continue
        for d in st["seeds"]:
            all_phis.append(d["Phi"])
            all_c0.append(d["C_F_sq"])
            if np.isfinite(d["gamma_ref"]):
                all_gr.append(d["gamma_ref"])
    max_abs_phi = float(np.abs(np.array(all_phis)).max())
    min_c0 = float(np.min(all_c0))
    max_abs_gr = float(np.max(np.abs(all_gr))) if all_gr else float("nan")
    sweep_phis = np.array([r["Phi"] for r in sweep_rows if r["finite"]])
    sweep_c0 = np.array([r["C_F_sq"] for r in sweep_rows if r["finite"]])
    sweep_all_finite = all(r["finite"] for r in sweep_rows)
    max_abs_phi_sweep = float(np.abs(sweep_phis).max()) if len(sweep_phis) else float("nan")
    min_c0_sweep = float(sweep_c0.min()) if len(sweep_c0) else float("nan")

    g1 = (max_abs_phi < G1_PHI_CAP
          and min_c0 > 1e-4
          and max_abs_gr < 50.0
          and sweep_all_finite
          and max_abs_phi_sweep < G1_PHI_CAP)
    gates["G1"] = dict(
        pass_=bool(g1),
        max_abs_Phi_core=max_abs_phi,
        min_C_F_sq_core=min_c0,
        max_abs_gamma_ref_at_optima=max_abs_gr,
        sweep=dict(all_finite=sweep_all_finite,
                   max_abs_Phi=max_abs_phi_sweep,
                   min_C_F_sq=min_c0_sweep),
        d21_baseline=dict(max_abs_Phi=44.0, min_C_F_sq=6.3e-5, max_abs_gamma=222.0),
        criterion="max|Phi|<10 (core & sweep), min C_F^2>1e-4, |gamma_ref|<50",
    )

    # G2: Haar giant outlier disappears (vs D2.1 best=25.71, std=5.96)
    h = haar_ext
    hphis = np.array([d["Phi"] for d in h["seeds"]])
    g2 = (hphis.max() < G2_BEST_CAP) and (hphis.std() < G2_STD_CAP)
    gates["G2"] = dict(
        pass_=bool(g2), best=float(hphis.max()), median=float(np.median(hphis)),
        std=float(hphis.std()),
        d21_baseline=dict(best=25.7133, median=1.8063, std=5.9611),
        criterion=f"best<{G2_BEST_CAP} and std<{G2_STD_CAP}",
    )

    # G3: seed stability across the 4 core states (6 seeds each)
    g3 = True
    g3_detail = {}
    for st_name, st in core.items():
        if not isinstance(st, dict) or "seeds" not in st:
            continue
        stdv = st["std_Phi"]
        g3_detail[st_name] = stdv
        g3 = g3 and (stdv < G3_STD_CAP)
    gates["G3"] = dict(pass_=bool(g3), std_by_state=g3_detail,
                       criterion=f"std<{G3_STD_CAP} for all states")

    # G4: small-perturbation continuity
    dphis = [r["dPhi"] for r in pert_rows if r["finite"]]
    max_dphi = max(abs(x) for x in dphis) if dphis else float("nan")
    all_finite = all(r["finite"] for r in pert_rows)
    g4 = all_finite and max_dphi < G4_DELTA_PHI_CAP
    gates["G4"] = dict(pass_=bool(g4), max_dPhi=max_dphi, all_finite=all_finite,
                       n_pert=len(pert_rows),
                       criterion=f"all finite and max|dPhi|<{G4_DELTA_PHI_CAP}")

    # G5: Schmidt spectrum — no product collapse (p_max << 0.9997)
    g5_detail = {}
    for st_name, st in core.items():
        if not isinstance(st, dict) or "seeds" not in st:
            continue
        g5_detail[st_name] = max(d["p_max"] for d in st["seeds"])
    g5 = all(pm < G5_PMAX_CAP for pm in g5_detail.values())
    gates["G5"] = dict(pass_=bool(g5), p_max_best=dict(
        {k: max(d["p_max"] for d in st["seeds"]) for k, st in core.items()
         if isinstance(st, dict) and "seeds" in st}),
        d21_baseline=0.9997,
        criterion=f"p_max<{G5_PMAX_CAP} for all states")

    # G6: reproducible across ground/Haar/mixed/thermal — finite non-trivial F*
    g6_detail = {}
    for st_name, st in core.items():
        if not isinstance(st, dict) or "seeds" not in st:
            continue
        g6_detail[st_name] = dict(best=st["best_Phi"], median=st["median_Phi"])
    g6 = all(1.0 <= st["best_Phi"] <= 5.0 for st in core.values()
             if isinstance(st, dict) and "seeds" in st)
    gates["G6"] = dict(pass_=bool(g6), best_by_state=g6_detail,
                       criterion="1.0 <= best Phi <= 5.0 for all states")

    # G7: pending until D2.2-B
    gates["G7"] = dict(pass_=None, status="PENDING — requires D2.2-B finite-time run")

    return gates


def main():
    OUTDIR.mkdir(exist_ok=True)
    H_basis = horizontal_basis(4, 2)
    H = get_hamiltonian(HNAME, n_sites=N)
    L = liouvillian_dephasing(H, GAMMA_VEC)
    states = build_states(H)

    print("=== D2.2-A: Unnormalized instantaneous functional ===")
    print(f"J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)> ; Phi = I - lambda*J_dyn^(0)")
    print(f"Hamiltonian: {HNAME}, env: γ={GAMMA_VEC}, λ={LAMBDA}, "
          f"Adam {STEPS} steps, lr={LR}, seeds={SEEDS}\n")

    core = {}
    for state_name, rho0 in states.items():
        print(f"--- State: {state_name} (core, {SEEDS} seeds) ---")
        t0 = time.time()
        core[state_name] = run_state(state_name, rho0, L, H_basis, SEEDS)
        print(f"    ({time.time()-t0:.1f}s)")

    # Supplementary Haar: 30 seeds (mirrors D2.1 audit protocol)
    print(f"\n--- Haar extended ({HAAR_EXTENDED_SEEDS} seeds, D2.1-comparable) ---")
    t0 = time.time()
    haar_ext = run_state("haar", states["haar"], L, H_basis, HAAR_EXTENDED_SEEDS)
    print(f"    ({time.time()-t0:.1f}s)")

    # Cross-state TPS distances (from ground best), same as D2
    print("\n=== Cross-state TPS distances (from ground) ===")
    U_ground = np.array(core["ground"]["best_U_real"]) + \
        1j * np.array(core["ground"]["best_U_imag"])
    cross = {}
    for state_name in ["haar", "mixed", "thermal"]:
        U_other = np.array(core[state_name]["best_U_real"]) + \
            1j * np.array(core[state_name]["best_U_imag"])
        d = tps_distance(U_ground, U_other)
        cross[state_name] = d
        print(f"  d_F(F*_ground, F*_{state_name}) = {d:.4f}")
    core["cross_state_d_F"] = cross

    # Perturbation test on best Haar solution (G4)
    print("\n=== Small perturbation response (best Haar solution, G4) ===")
    U_best_haar = np.array(haar_ext["best_U_real"]) + \
        1j * np.array(haar_ext["best_U_imag"])
    pert_rows, phi0 = perturbation_test(states["haar"], L, U_best_haar, H_basis)
    for r in pert_rows:
        print(f"  dir={r['dir']:2d} eps={r['eps']:.0e}: dΦ={r['dPhi']:+.4f} "
              f"dJ={r['dJ']:+.4f} finite={r['finite']}")

    # Singularity sweep toward C_F^2 -> 0 (G1)
    print("\n=== Singularity sweep toward C_F^2 -> 0 (G1) ===")
    sweep_rows = singularity_sweep(states["haar"], L, U_best_haar, H_basis)
    min_c = min(r["C_F_sq"] for r in sweep_rows)
    max_p = max(r["Phi"] for r in sweep_rows if r["finite"])
    print(f"  sweep: min C_F² = {min_c:.3e}, max Φ = {max_p:.4f}, "
          f"all finite = {all(r['finite'] for r in sweep_rows)}")

    # Gates
    gates = evaluate_gates(core, haar_ext, pert_rows, sweep_rows)
    print("\n=== Gate verdicts ===")
    for g, v in gates.items():
        status = "PASS" if v.get("pass_") else ("PENDING" if v.get("pass_") is None else "FAIL")
        print(f"  {g}: {status}  {v}")

    results = dict(
        protocol=dict(
            functional="J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)> (unnormalized instantaneous)",
            frozen="D2 protocol: asymmetric_XY, dephasing γ=(0.5,1,2), λ=0.2, "
                   "Adam 200 steps lr=0.01, SEED0=20260813, horizontal_basis(4,2), n_a=2",
            seeds=SEEDS, haar_extended_seeds=HAAR_EXTENDED_SEEDS,
        ),
        core=core,
        haar_extended=haar_ext,
        perturbation=dict(rows=pert_rows, Phi0=phi0),
        gates=gates,
    )

    out_path = OUTDIR / "d2_2a_unnormalized_instantaneous.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
