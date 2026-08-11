#!/usr/bin/env python3
"""
Phase D1: Hamiltonian Generality — optimize F* for each Hamiltonian family.

Protocol (from d0_protocol.json):
- λ = 0.2, τ = 0.1
- Adam, 200 steps, lr=0.01
- 6 seeds per Hamiltonian
- Environment: pure dephasing γ=(0.5, 1.0, 2.0)

Objective: Φ = I_ρ(F) − λ·Γ_F^exact(0) (initial-rate canonical)
Success: Φ* > contiguous best for each family
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
from gamma_F import schmidt_basis
from b2_canonical_open import (
    canonical_open_objective, apply_tps, liouvillian_dephasing,
)
from b2_1_cross_eval import tps_distance

LAMBDA = 0.2
STEPS = 200
LR = 0.01
SEEDS = 6
SEED0 = 20260812
GAMMA_VEC = [0.5, 1.0, 2.0]  # A1 environment

HAMILTONIANS = {
    "asymmetric_XY": {},
    "transverse_field_ising": {"J": 1.0, "h": 0.5},
    "heisenberg_xxz": {"J": 1.0, "Delta": 0.5},
    "xyz": {"J_x": 1.0, "J_y": 0.8, "J_z": 0.6},
    "random_local": {"seed": 42},
}


def canonical_open_objective_lambda(rho_U, Y_U, n_a, lambda_val):
    """Canonical objective with explicit λ (module default is 0.2)."""
    obj = canonical_open_objective(rho_U, Y_U, n_a)
    phi = obj["I"] - lambda_val * obj["gamma"] if np.isfinite(obj["gamma"]) else -1e9
    return obj, phi


def adam_optimize_open(U0, L, psi0, lambda_val, H_basis, steps=STEPS, lr=LR):
    """Adam optimization of Φ = I − λΓ^exact in the right-mult chart."""
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_at(d):
        U = to_U(d)
        rho_U, Y_U = apply_tps(psi0, L, U)
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


def evaluate_U(U, psi0, L, lambda_val):
    """Full evaluation of a U: Φ, I, Γ, plus I/Γ components."""
    rho_U, Y_U = apply_tps(psi0, L, U)
    obj, phi = canonical_open_objective_lambda(rho_U, Y_U, 2, lambda_val)
    return dict(Phi=phi, I=obj["I"], gamma=obj["gamma"], C_F_sq=obj["C_F_sq"])


def contiguous_phi(psi0, L, lambda_val, cut=2):
    """Φ at the contiguous cut (identity U)."""
    U = np.eye(8, dtype=complex)
    return evaluate_U(U, psi0, L, lambda_val)


def main():
    OUTDIR.mkdir(exist_ok=True)
    H_basis = horizontal_basis(4, 2)
    results = {}

    print("=== D1: Hamiltonian Generality ===")
    print(f"Environment: pure dephasing γ={GAMMA_VEC}, λ={LAMBDA}\n")

    for hname, hparams in HAMILTONIANS.items():
        print(f"--- {hname} ---")
        H = get_hamiltonian(hname, n_sites=3, **hparams)
        E0, psi0 = get_ground_state(H)
        L = liouvillian_dephasing(H, GAMMA_VEC)

        # Contiguous baseline (cut2 = AB|C)
        cont = contiguous_phi(psi0, L, LAMBDA)
        print(f"  E0 = {E0:.4f}, contiguous Φ(cut2) = {cont['Phi']:.4f}")

        # Optimize with 6 seeds
        seeds_data = []
        for s in range(SEEDS):
            seed = SEED0 + s
            rng = np.random.default_rng(seed)
            theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
            U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))

            U_best, phi_best = adam_optimize_open(U0, L, psi0, LAMBDA, H_basis)
            ev = evaluate_U(U_best, psi0, L, LAMBDA)

            # Distance from contiguous cut2
            d_contig = tps_distance(U_best, np.eye(8, dtype=complex))

            seeds_data.append(dict(
                seed=seed, Phi=ev["Phi"], I=ev["I"], gamma=ev["gamma"],
                d_F_contiguous=d_contig,
                U_real=U_best.real.tolist(), U_imag=U_best.imag.tolist(),
            ))

        phis = [d["Phi"] for d in seeds_data]
        best = max(phis)
        best_idx = int(np.argmax(phis))
        success = best > cont["Phi"]
        print(f"  best Φ = {best:.4f} (seed {best_idx}), median = {np.median(phis):.4f}, "
              f"std = {np.std(phis):.4f}")
        print(f"  success (Φ* > contiguous) = {success}\n")

        results[hname] = dict(
            E0=float(E0),
            contiguous=cont,
            seeds=seeds_data,
            best_Phi=best,
            median_Phi=float(np.median(phis)),
            std_Phi=float(np.std(phis)),
            success=success,
        )

    # Summary
    print("\n=== Summary ===")
    print(f"{'Hamiltonian':<25s} {'contiguous':>10s} {'best Φ':>10s} {'success':>8s}")
    print("-" * 55)
    for hname in HAMILTONIANS:
        r = results[hname]
        print(f"{hname:<25s} {r['contiguous']['Phi']:10.4f} {r['best_Phi']:10.4f} "
              f"{str(r['success']):>8s}")

    out_path = OUTDIR / "d1_hamiltonian_generality.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
