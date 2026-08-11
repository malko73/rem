#!/usr/bin/env python3
"""
Phase B2 convergence test: 20 trials × 300 steps to assess optimizer convergence.

Purpose: determine if 50 steps is insufficient (B2-C verdict). If within-environment
distance decreases significantly with 300 steps, the issue is optimizer convergence.
If it remains ~0.87, the landscape is genuinely multi-modal.

This test runs on A1 and A2 environments only (mirror pair, most important for
F* = F*(ρ, L, λ) verification).
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
LAMBDA = 0.2
STEPS = 300
LR = 0.01
N_TRIALS = 20
SEED0 = 20260816


def canonical_open_objective(rho_U: np.ndarray, Y_U: np.ndarray, n_a: int) -> dict:
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

    Phi = I_F - LAMBDA * gamma if np.isfinite(gamma) else float("nan")
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


def run_convergence_test(env_name: str, L: np.ndarray):
    """Run 20 trials × 300 steps for one environment."""
    print(f"\n=== Convergence Test: {env_name} (20 trials × {STEPS} steps) ===")
    h_total, bonds = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = rem4.ground_state(h_total)

    H_basis = qgeom.horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(SEED0)

    results = []
    t0 = time.time()
    for trial in range(N_TRIALS):
        theta = rng.standard_normal(H_basis.shape[2]) * 1.0
        U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))
        rho_U, Y_U = apply_tps(psi0, L, U)
        obj = canonical_open_objective(rho_U, Y_U, N_A)

        theta_opt = theta.copy()
        m = np.zeros_like(theta_opt)
        v = np.zeros_like(theta_opt)
        beta1, beta2, eps = 0.9, 0.999, 1e-8
        best_phi = obj["Phi"]
        best_u = U
        best_obj = obj

        for step in range(STEPS):
            grad = np.zeros_like(theta_opt)
            for k in range(len(theta_opt)):
                theta_p = theta_opt.copy()
                theta_m = theta_opt.copy()
                theta_p[k] += 1e-5
                theta_m[k] -= 1e-5
                U_p = sla.expm(np.einsum("ijk,k->ij", H_basis, theta_p))
                U_m = sla.expm(np.einsum("ijk,k->ij", H_basis, theta_m))
                rho_p, Y_p = apply_tps(psi0, L, U_p)
                rho_m, Y_m = apply_tps(psi0, L, U_m)
                obj_p = canonical_open_objective(rho_p, Y_p, N_A)
                obj_m = canonical_open_objective(rho_m, Y_m, N_A)
                grad[k] = (obj_p["Phi"] - obj_m["Phi"]) / (2 * 1e-5)

            m = beta1 * m + (1 - beta1) * grad
            v = beta2 * v + (1 - beta2) * grad ** 2
            m_hat = m / (1 - beta1 ** (step + 1))
            v_hat = v / (1 - beta2 ** (step + 1))
            theta_opt = theta_opt + LR * m_hat / (np.sqrt(v_hat) + eps)

            U_opt = sla.expm(np.einsum("ijk,k->ij", H_basis, theta_opt))
            rho_opt, Y_opt = apply_tps(psi0, L, U_opt)
            obj_opt = canonical_open_objective(rho_opt, Y_opt, N_A)
            if np.isfinite(obj_opt["Phi"]) and obj_opt["Phi"] > best_phi:
                best_phi = obj_opt["Phi"]
                best_u = U_opt
                best_obj = obj_opt

        results.append(dict(
            trial=trial,
            seed=SEED0 + trial,
            initial_phi=obj["Phi"],
            best_phi=best_phi,
            best_u_real=best_u.real.tolist(),
            best_u_imag=best_u.imag.tolist(),
            best_I=best_obj["I"],
            best_gamma=best_obj["gamma"],
            best_C_F_sq=best_obj["C_F_sq"],
            converged=best_phi > obj["Phi"],
        ))

        if (trial + 1) % 5 == 0:
            elapsed = time.time() - t0
            print(f"  trial {trial+1}/{N_TRIALS} ({elapsed:.1f}s)")

    return results


def compute_within_distance(trials: list) -> dict:
    """Compute within-environment TPS distance statistics."""
    from b2_1_cross_eval import tps_distance
    n = len(trials)
    distances = []
    for i in range(n):
        for j in range(i + 1, n):
            u_i = np.array(trials[i]["best_u_real"]) + 1j * np.array(trials[i]["best_u_imag"])
            u_j = np.array(trials[j]["best_u_real"]) + 1j * np.array(trials[j]["best_u_imag"])
            d = tps_distance(u_i, u_j)
            distances.append(d)

    distances = np.array(distances)
    return dict(
        mean=float(np.mean(distances)),
        median=float(np.median(distances)),
        std=float(np.std(distances)),
        max=float(np.max(distances)),
        min=float(np.min(distances)),
        frac_lt_01=float(np.mean(distances < 0.1)),
    )


def main():
    OUTDIR.mkdir(exist_ok=True)

    h_total, _ = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)

    envs = {
        "A1_deph_0512": liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
        "A2_deph_2105": liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
    }

    results = {}
    for env_name, L in envs.items():
        trials = run_convergence_test(env_name, L)
        stats = compute_within_distance(trials)
        results[env_name] = dict(trials=trials, within_stats=stats)

        print(f"\n{env_name}:")
        print(f"  within d_F: mean={stats['mean']:.4f}  median={stats['median']:.4f}  "
              f"std={stats['std']:.4f}  frac<0.1={stats['frac_lt_01']:.3f}")
        print(f"  best Φ: {max(t['best_phi'] for t in trials):.4f}")

    out_path = OUTDIR / "b2_convergence_test.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
