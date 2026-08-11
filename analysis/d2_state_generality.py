#!/usr/bin/env python3
"""
Phase D2: State Generality — optimize F* for different state types.

Protocol (d0_protocol.json):
- λ = 0.2, Adam, 200 steps, lr=0.01, 6 seeds
- Environment: pure dephasing γ=(0.5, 1.0, 2.0)
- States: ground state, Haar pure state, mixed state, thermal state

Success criterion: d_F(F*_ground, F*_thermal) < 0.1
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from d1_hamiltonians import get_hamiltonian, get_ground_state
from quotient_geometry import horizontal_basis
from b2_canonical_open import canonical_open_objective, apply_tps, liouvillian_dephasing
from b2_1_cross_eval import tps_distance

LAMBDA = 0.2
STEPS = 200
LR = 0.01
SEEDS = 6
SEED0 = 20260813
GAMMA_VEC = [0.5, 1.0, 2.0]
HNAME = "asymmetric_XY"  # reference Hamiltonian for D2


def canonical_open_objective_lambda(rho_U, Y_U, n_a, lambda_val):
    obj = canonical_open_objective(rho_U, Y_U, n_a)
    phi = obj["I"] - lambda_val * obj["gamma"] if np.isfinite(obj["gamma"]) else -1e9
    return obj, phi


def adam_optimize_open(U0, L, rho0, Y0, lambda_val, H_basis, steps=STEPS, lr=LR):
    """Adam optimization for a possibly MIXED initial state.

    Y0 = L(rho0) is precomputed; the pull-back rotates both rho0 and Y0.
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
        obj, phi = canonical_open_objective_lambda(rho_U, Y_U, 2, lambda_val)
        return phi

    best_phi = phi_at(delta)
    best_U = to_U(delta)

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

    return best_U, best_phi


def evaluate_U(U, rho0, Y0, lambda_val):
    rho_U = U.conj().T @ rho0 @ U
    Y_U = U.conj().T @ Y0 @ U
    obj, phi = canonical_open_objective_lambda(rho_U, Y_U, 2, lambda_val)
    return dict(Phi=phi, I=obj["I"], gamma=obj["gamma"], C_F_sq=obj["C_F_sq"])


def build_states(H, n_sites=3):
    """Build 4 state types: ground, Haar, mixed, thermal."""
    dim = 2 ** n_sites
    E0, psi_gs = get_ground_state(H)
    rng = np.random.default_rng(20260813)

    states = {}

    # 1. Ground state (pure)
    states["ground"] = np.outer(psi_gs, psi_gs.conj())

    # 2. Haar pure state (random unitary applied to |0...0>)
    m = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    q, _ = np.linalg.qr(m)
    psi_haar = q[:, 0]
    states["haar"] = np.outer(psi_haar, psi_haar.conj())

    # 3. Mixed state: 70/30 mixture of ground and first excited
    evals, evecs = np.linalg.eigh(H)
    rho_mix = 0.7 * np.outer(evecs[:, 0], evecs[:, 0].conj()) + \
              0.3 * np.outer(evecs[:, 1], evecs[:, 1].conj())
    states["mixed"] = rho_mix

    # 4. Thermal state at T=0.5 (β = 1/T = 2)
    beta = 2.0
    rho_thermal = np.exp(-beta * H) / np.trace(np.exp(-beta * H))
    states["thermal"] = rho_thermal

    return states


def main():
    OUTDIR.mkdir(exist_ok=True)
    H_basis = horizontal_basis(4, 2)
    H = get_hamiltonian(HNAME, n_sites=3)
    L = liouvillian_dephasing(H, GAMMA_VEC)
    states = build_states(H)

    results = {}
    print("=== D2: State Generality ===")
    print(f"Hamiltonian: {HNAME}, environment: γ={GAMMA_VEC}, λ={LAMBDA}\n")

    for state_name, rho0 in states.items():
        Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        print(f"--- State: {state_name} ---")

        seeds_data = []
        for s in range(SEEDS):
            seed = SEED0 + s
            rng = np.random.default_rng(seed)
            theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
            U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

            U_best, phi_best = adam_optimize_open(U0, L, rho0, Y0, LAMBDA, H_basis)
            ev = evaluate_U(U_best, rho0, Y0, LAMBDA)
            seeds_data.append(dict(
                seed=seed, Phi=ev["Phi"], I=ev["I"], gamma=ev["gamma"],
                U_real=U_best.real.tolist(), U_imag=U_best.imag.tolist(),
            ))

        phis = [d["Phi"] for d in seeds_data]
        best = max(phis)
        best_idx = int(np.argmax(phis))
        print(f"  best Φ = {best:.4f} (seed {best_idx}), median = {np.median(phis):.4f}, "
              f"std = {np.std(phis):.4f}")

        results[state_name] = dict(
            seeds=seeds_data,
            best_Phi=best,
            median_Phi=float(np.median(phis)),
            std_Phi=float(np.std(phis)),
        )
        print()

    # Cross-state distances: d_F(F*_ground, F*_other)
    print("=== Cross-state TPS distances (from ground) ===")
    U_ground = np.array(results["ground"]["seeds"][int(np.argmax(
        [d["Phi"] for d in results["ground"]["seeds"]]))]["U_real"]) + \
        1j * np.array(results["ground"]["seeds"][int(np.argmax(
            [d["Phi"] for d in results["ground"]["seeds"]]))]["U_imag"])

    cross = {}
    for state_name in ["haar", "mixed", "thermal"]:
        best_s = int(np.argmax([d["Phi"] for d in results[state_name]["seeds"]]))
        U_other = np.array(results[state_name]["seeds"][best_s]["U_real"]) + \
            1j * np.array(results[state_name]["seeds"][best_s]["U_imag"])
        d = tps_distance(U_ground, U_other)
        cross[state_name] = d
        print(f"  d_F(F*_ground, F*_{state_name}) = {d:.4f}")

    results["cross_state_d_F"] = cross
    results["success_d2"] = cross.get("thermal", 1.0) < 0.1

    out_path = OUTDIR / "d2_state_generality.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nD2 success (d_F(F*_ground, F*_thermal) < 0.1): {results['success_d2']}")
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()
