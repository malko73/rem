#!/usr/bin/env python3
"""
Phase D2.2-B — Finite-time dynamical functional (per-tau run).

Functional (master's formula):
    J_dyn^(tau)(F) = -(C_F^2(tau) - C_F^2(0)) / (2 tau)
    Phi(F)         = I(F) - lambda * J_dyn^(tau)(F)

C_F^2(t) = ||Q_F rho_U(t)||^2 with Q_F = I - D_F, D_F the dephasing
projection onto the FIXED t=0 Schmidt basis of the rotated state's dominant
eigenvector. Q_F is time-independent, so (d/dt) C_F^2 |_0 =
2 Re<Q_F rho, Q_F L(rho)> and J_dyn^(tau) -> J_dyn^(0) as tau -> 0.

Frozen protocol: identical to D2.2-A (asymmetric_XY, dephasing
gamma=(0.5,1,2), lambda=0.2, Adam 200 steps lr=0.01, SEED0=20260813,
horizontal_basis(4,2), n_a=2, 6 seeds). Only the dynamical term changes.

Efficiency: rho_tau = e^(tau L) rho0 is U-INDEPENDENT and precomputed once
per (tau, state); every objective evaluation only does 8x8 similarity
transforms (no expm in the loop).

Usage:
    python d2_2b_finite_time.py --tau 0.1 [--haar-ext] [--seeds 6]
Writes analysis_output/d2_2b_tau_<tag>.json
"""
from __future__ import annotations

import argparse
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
import rem4_numerical as rem4
import gamma_F

LAMBDA = 0.2
STEPS = 200
LR = 0.01
SEED0 = 20260813
GAMMA_VEC = [0.5, 1.0, 2.0]
HNAME = "asymmetric_XY"
N = 3
N_A = 2
TAU_GRID = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]


def tau_tag(tau: float) -> str:
    return f"t{tau:.0e}".replace("e-0", "e-").replace("e-", "e-")


def build_states(H):
    """Identical to d2_2a / d2_state_generality."""
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


def finite_time_objective(rho0, Y0, rho_tau, U, n_a, tau):
    """Compute I, C_F^2(0), C_F^2(tau), J_dyn^(tau), J_dyn^(0) ref, Phi.

    rho0: original state; Y0 = L(rho0) (vectorised, reshaped); rho_tau =
    e^(tau L) rho0 (precomputed); U: TPS unitary.
    The Schmidt basis is fixed at t=0 (dominant eigenvector of U† rho0 U).
    """
    rho_U = U.conj().T @ rho0 @ U
    evals_psi, evecs_psi = np.linalg.eigh(rho_U)
    psi_U = evecs_psi[:, np.argmax(evals_psi)]
    I_F = rem4.mutual_information_for_cut(psi_U, N, n_a)

    basis = gamma_F.schmidt_basis(psi_U, N, n_a)
    q = lambda r: r - gamma_F.dephasing_projection(r, basis)  # noqa: E731

    x = q(rho_U)
    c0_sq = float(np.linalg.norm(x, ord="fro") ** 2)

    rho_U_tau = U.conj().T @ rho_tau @ U
    xt = q(rho_U_tau)
    ctau_sq = float(np.linalg.norm(xt, ord="fro") ** 2)

    J_tau = -(ctau_sq - c0_sq) / (2.0 * tau)

    # instantaneous reference (G7-1): J_dyn^(0) = -Re<Q_F rho, Q_F L(rho)>
    Y_U = U.conj().T @ Y0 @ U
    qy = q(Y_U)
    J_0 = -float(np.real(np.trace(x.conj().T @ qy)))

    Phi = I_F - LAMBDA * J_tau
    return dict(I=I_F, C_F_sq_0=c0_sq, C_F_sq_tau=ctau_sq,
                J_dyn_tau=J_tau, J_dyn_0=J_0, Phi=Phi)


def finite_time_objective_lambda(rho0, Y0, rho_tau, U, n_a, tau, lambda_val):
    obj = finite_time_objective(rho0, Y0, rho_tau, U, n_a, tau)
    phi = obj["I"] - lambda_val * obj["J_dyn_tau"] if np.isfinite(obj["J_dyn_tau"]) else -1e9
    return obj, phi


def adam_optimize_open(U0, rho0, Y0, rho_tau, tau, H_basis, steps=STEPS, lr=LR):
    """Adam optimization with the finite-time objective (frozen protocol)."""
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_at(d):
        U = to_U(d)
        obj, phi = finite_time_objective_lambda(rho0, Y0, rho_tau, U, N_A, tau, LAMBDA)
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
            best_obj, _ = finite_time_objective_lambda(rho0, Y0, rho_tau, best_U, N_A, tau, LAMBDA)

    if best_obj is None:
        best_obj, _ = finite_time_objective_lambda(rho0, Y0, rho_tau, best_U, N_A, tau, LAMBDA)

    return best_U, best_phi, best_obj


def schmidt_spectrum(psi, n=3, cut=2):
    d_a, d_b = 2 ** cut, 2 ** (n - cut)
    mat = psi.reshape((d_a, d_b))
    _, s, _ = np.linalg.svd(mat)
    return s ** 2


def dominant_eigenvector(rho):
    evals, evecs = np.linalg.eigh(rho)
    return evecs[:, np.argmax(evals)]


def run_state(state_name, rho0, Y0, rho_tau, tau, H_basis, n_seeds):
    seeds_data = []
    best_U = None
    best_phi = -np.inf
    for s in range(n_seeds):
        seed = SEED0 + s
        rng = np.random.default_rng(seed)
        theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
        U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

        U_best, phi_best, obj = adam_optimize_open(U0, rho0, Y0, rho_tau, tau, H_basis)

        psi0 = dominant_eigenvector(rho0)
        psi_U = U_best.conj().T @ psi0
        spec = schmidt_spectrum(psi_U)

        seeds_data.append(dict(
            seed=seed, Phi=phi_best, I=obj["I"],
            J_dyn_tau=obj["J_dyn_tau"], J_dyn_0=obj["J_dyn_0"],
            C_F_sq_0=obj["C_F_sq_0"], C_F_sq_tau=obj["C_F_sq_tau"],
            p_max=float(spec.max()),
        ))
        if phi_best > best_phi:
            best_phi = phi_best
            best_U = U_best

    phis = np.array([d["Phi"] for d in seeds_data])
    return dict(
        seeds=seeds_data,
        best_Phi=float(phis.max()),
        median_Phi=float(np.median(phis)),
        std_Phi=float(phis.std()),
        best_U_real=best_U.real.tolist(),
        best_U_imag=best_U.imag.tolist(),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tau", type=float, required=True)
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--haar-ext", action="store_true",
                    help="also run Haar with 30 seeds (D2.1/D2.2-A-comparable)")
    args = ap.parse_args()

    tau = args.tau
    if tau not in TAU_GRID:
        print(f"WARNING: tau={tau} not in canonical grid {TAU_GRID}")
    OUTDIR.mkdir(exist_ok=True)
    H_basis = horizontal_basis(4, 2)
    H = get_hamiltonian(HNAME, n_sites=N)
    L = liouvillian_dephasing(H, GAMMA_VEC)
    states = build_states(H)

    print(f"=== D2.2-B: finite-time functional, τ={tau} ===")
    print(f"J_dyn^(τ) = -(C_F²(τ)-C_F²(0))/(2τ); Φ = I - λ·J_dyn^(τ)")
    print(f"λ={LAMBDA}, Adam {STEPS} steps, lr={LR}, seeds={args.seeds}\n")

    core = {}
    for state_name, rho0 in states.items():
        t0 = time.time()
        Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        rho_tau = (sla.expm(tau * L) @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        core[state_name] = run_state(state_name, rho0, Y0, rho_tau, tau, H_basis, args.seeds)
        phis = [d["Phi"] for d in core[state_name]["seeds"]]
        print(f"  {state_name}: best Φ = {max(phis):.4f}, median = {np.median(phis):.4f}, "
              f"std = {np.std(phis):.4f}  ({time.time()-t0:.1f}s)")

    haar_ext = None
    if args.haar_ext:
        t0 = time.time()
        rho0 = states["haar"]
        Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        rho_tau = (sla.expm(tau * L) @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        haar_ext = run_state("haar", rho0, Y0, rho_tau, tau, H_basis, 30)
        phis = [d["Phi"] for d in haar_ext["seeds"]]
        print(f"  haar_ext(30): best Φ = {max(phis):.4f}, median = {np.median(phis):.4f}, "
              f"std = {np.std(phis):.4f}  ({time.time()-t0:.1f}s)")

    results = dict(
        tau=tau,
        protocol=dict(
            functional="J_dyn^(τ) = -(C_F²(τ)-C_F²(0))/(2τ)",
            frozen="D2 protocol: asymmetric_XY, dephasing γ=(0.5,1,2), λ=0.2, "
                   "Adam 200 steps lr=0.01, SEED0=20260813, horizontal_basis(4,2), n_a=2",
            seeds=args.seeds,
        ),
        core=core,
        haar_extended=haar_ext,
    )

    tag = tau_tag(tau)
    out_path = OUTDIR / f"d2_2b_tau_{tag}.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
