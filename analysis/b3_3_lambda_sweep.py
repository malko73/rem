#!/usr/bin/env python3
"""
B3.3: Lambda sweep to test F* = F*(λ) hypothesis.

Purpose: Determine whether the optimal factorization F* changes with λ.
If F*(λ) is constant across a range of λ, this supports the hypothesis that
the environment weights candidate stability rather than rebuilding structure.

Method:
- For uniform/A1/A2 environments, run optimization at λ = 0.05, 0.1, 0.2, 0.3, 0.5
- For each λ, use warm start (best U* from λ=0.2) + 5 random starts
- Compute gauge-invariant TPS distance d_F between F*(λ_i) and F*(λ_j)
- If d_F ≈ 0 across all λ, then F* is λ-independent

Computational cost: 3 envs × 5 λ × (1 warm + 5 random) × 50 steps ≈ 900 runs ≈ 30 minutes
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"

spec = importlib.util.spec_from_file_location("rem4_numerical", REPO / "src" / "rem4_numerical.py")
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)

spec2 = importlib.util.spec_from_file_location("quotient_geometry", REPO / "src" / "quotient_geometry.py")
qgeom = importlib.util.module_from_spec(spec2)
assert spec2.loader is not None
sys.modules[spec2.name] = qgeom
spec2.loader.exec_module(qgeom)

spec3 = importlib.util.spec_from_file_location("gamma_F", REPO / "analysis" / "gamma_F.py")
gamma_F = importlib.util.module_from_spec(spec3)
assert spec3.loader is not None
sys.modules[spec3.name] = gamma_F
spec3.loader.exec_module(gamma_F)

N, N_A = 3, 2
D_A, D_B = 4, 2
DIM = 8
LAMBDA_VALUES = [0.05, 0.1, 0.2, 0.3, 0.5]
STEPS = 50
LR = 0.01
SEED0 = 20260818


def canonical_open_objective(rho_U: np.ndarray, Y_U: np.ndarray, n_a: int, lambda_val: float) -> dict:
    """Compute I, C_F², Γ, Φ for the canonical open-system objective."""
    evals_psi, evecs_psi = np.linalg.eigh(rho_U)
    psi_U = evecs_psi[:, np.argmax(evals_psi)]
    I_F = rem4.mutual_information_for_cut(psi_U, N, n_a)

    basis = gamma_F.schmidt_basis(psi_U, N, n_a)
    q = lambda r: r - gamma_F.dephasing_projection(r, basis)
    x = q(rho_U)
    qy = q(Y_U)
    num = np.real(np.trace(x.conj().T @ qy))
    den = float(np.linalg.norm(x, ord="fro") ** 2)
    if den < 1e-300:
        gamma = float("nan")
    else:
        gamma = float(-num / den)

    Phi = I_F - lambda_val * gamma if np.isfinite(gamma) else float("nan")
    return dict(I=I_F, C_F_sq=den, gamma=gamma, Phi=Phi)


def liouvillian_dephasing(h_total: np.ndarray, gamma_vec: list[float]) -> np.ndarray:
    dim = 2 ** N
    I = np.eye(dim, dtype=complex)
    L = -1j * (np.kron(I, h_total) - np.kron(h_total.conj(), I))
    for site in range(N):
        g = gamma_vec[site]
        if g > 0:
            z = gamma_F.pauli_z(N, site)
            L += (g / 2.0) * (np.kron(z, z) - np.kron(I, I))
    return L


def apply_tps(psi0: np.ndarray, L: np.ndarray, U: np.ndarray):
    rho0 = np.outer(psi0, psi0.conj())
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    rho_U = U.conj().T @ rho0 @ U
    Y_U = U.conj().T @ Y0 @ U
    return rho_U, Y_U


def optimize_from_U(U_start: np.ndarray, L: np.ndarray, psi0: np.ndarray, lambda_val: float,
                    seed_offset: int = 0) -> tuple:
    """
    Run optimization starting from U_start.

    Parameterization: U = U_start @ expm(Σ δ_k H_k), with δ the search
    coordinates. This keeps the current iterate as the identity of the
    exponential chart, so warm start is exact (δ=0 reproduces U_start).

    Returns (best_U, best_Phi, delta_opt).
    """
    H_basis = qgeom.horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(SEED0 + seed_offset)
    delta = rng.standard_normal(H_basis.shape[2]) * 0.05  # small perturbation

    def to_U(d):
        return U_start @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    U = to_U(delta)
    rho_U, Y_U = apply_tps(psi0, L, U)
    obj = canonical_open_objective(rho_U, Y_U, N_A, lambda_val)
    best_phi = obj["Phi"]
    best_u = U
    best_delta = delta.copy()

    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    for step in range(STEPS):
        grad = np.zeros_like(delta)
        for k in range(len(delta)):
            d_p = delta.copy()
            d_m = delta.copy()
            d_p[k] += 1e-5
            d_m[k] -= 1e-5
            U_p = to_U(d_p)
            U_m = to_U(d_m)
            rho_p, Y_p = apply_tps(psi0, L, U_p)
            rho_m, Y_m = apply_tps(psi0, L, U_m)
            obj_p = canonical_open_objective(rho_p, Y_p, N_A, lambda_val)
            obj_m = canonical_open_objective(rho_m, Y_m, N_A, lambda_val)
            grad[k] = (obj_p["Phi"] - obj_m["Phi"]) / (2 * 1e-5)

        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        delta = delta + LR * m_hat / (np.sqrt(v_hat) + eps)

        U_opt = to_U(delta)
        rho_opt, Y_opt = apply_tps(psi0, L, U_opt)
        obj_opt = canonical_open_objective(rho_opt, Y_opt, N_A, lambda_val)
        if np.isfinite(obj_opt["Phi"]) and obj_opt["Phi"] > best_phi:
            best_phi = obj_opt["Phi"]
            best_u = U_opt
            best_delta = delta.copy()

    return best_u, best_phi, best_delta


def tps_distance(U1: np.ndarray, U2: np.ndarray) -> float:
    """Gauge-invariant TPS distance (same as B2.1, via A-side local algebra projectors)."""
    from b2_1_cross_eval import tps_distance as _d
    return _d(U1, U2)


def main():
    OUTDIR.mkdir(exist_ok=True)

    h_total, _ = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = rem4.ground_state(h_total)

    envs = {
        "uniform_dephasing": liouvillian_dephasing(h_total, [1.0, 1.0, 1.0]),
        "A1_deph_0512": liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
        "A2_deph_2105": liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
    }

    # Load best U* from B2 for warm start
    b2_path = OUTDIR / "b2_canonical_open.json"
    with open(b2_path) as f:
        b2_data = json.load(f)

    results = {}

    for env_name, L in envs.items():
        print(f"\n=== Lambda sweep: {env_name} ===")
        best_trial = max(b2_data[env_name], key=lambda t: t["best_phi"])
        U_warm = np.array(best_trial["best_u_real"]) + 1j * np.array(best_trial["best_u_imag"])

        env_results = {}

        for lambda_val in LAMBDA_VALUES:
            print(f"\n  λ = {lambda_val}")

            # Warm start (delta=0.05 perturbation from best U* at λ=0.2)
            U_w, Phi_w, _ = optimize_from_U(U_warm, L, psi0, lambda_val, seed_offset=0)
            print(f"    Warm start: Φ = {Phi_w:.6f}")

            # 5 random starts
            random_results = []
            for i in range(5):
                rng = np.random.default_rng(SEED0 + i)
                H_basis = qgeom.horizontal_basis(D_A, D_B)
                theta = rng.standard_normal(H_basis.shape[2]) * 1.0
                U_rand = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))
                U_r, Phi_r, _ = optimize_from_U(U_rand, L, psi0, lambda_val, seed_offset=10 + i)
                random_results.append((U_r, Phi_r))
                print(f"    Random start {i+1}: Φ = {Phi_r:.6f}")

            # Best overall
            all_results = [(U_w, Phi_w)] + random_results
            best_idx = max(range(len(all_results)), key=lambda i: all_results[i][1])
            U_best, Phi_best = all_results[best_idx]

            env_results[lambda_val] = {
                "best_U_real": U_best.real.tolist(),
                "best_U_imag": U_best.imag.tolist(),
                "best_Phi": Phi_best,
                "warm_Phi": Phi_w,
                "random_Phis": [r[1] for r in random_results],
            }

        # Compute pairwise TPS distances between best U* at different λ
        distances = {}
        for i, l1 in enumerate(LAMBDA_VALUES):
            for j, l2 in enumerate(LAMBDA_VALUES):
                if i < j:
                    U1 = np.array(env_results[l1]["best_U_real"]) + 1j * np.array(env_results[l1]["best_U_imag"])
                    U2 = np.array(env_results[l2]["best_U_real"]) + 1j * np.array(env_results[l2]["best_U_imag"])
                    d = tps_distance(U1, U2)
                    distances[f"λ={l1}_vs_λ={l2}"] = d
                    print(f"\n  d_F(F*(λ={l1}), F*(λ={l2})) = {d:.6f}")

        env_results["distances"] = distances
        results[env_name] = env_results

    # Save results
    out_path = OUTDIR / "b3_3_lambda_sweep.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
