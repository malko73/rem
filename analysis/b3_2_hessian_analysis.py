#!/usr/bin/env python3
"""
B3.2: Hessian analysis at best U* for uniform/A1/A2 environments.

Purpose: Compute 45-dimensional Hessian matrix at best solutions to verify
they are local maxima (all eigenvalues ≤ 0) or identify saddle points.

Method:
- Load best U* from B2 results for uniform/A1/A2
- Compute Hessian H_ij = ∂²Φ/∂θ_i∂θ_j via finite differences
- Diagonalize H and report eigenvalue spectrum
- Classify: strict local max (all λ < 0), flat max (some λ ≈ 0), or saddle (some λ > 0)

Computational cost: ~4050 Φ evaluations per environment × 3 environments ≈ 20 minutes
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
H_EPS = 1e-4  # Finite difference step size


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


def phi_at_U(U: np.ndarray, psi0: np.ndarray, L: np.ndarray) -> float:
    """Compute Φ at given unitary U."""
    rho_U, Y_U = apply_tps(psi0, L, U)
    obj = canonical_open_objective(rho_U, Y_U, N_A)
    return obj["Phi"]


def compute_hessian(U: np.ndarray, H_basis: np.ndarray, psi0: np.ndarray, L: np.ndarray) -> np.ndarray:
    """
    Compute Hessian matrix H_ij = ∂²Φ/∂θ_i∂θ_j via finite differences.

    Perturb U directly: U' = U @ expm(ε * H_k)
    This avoids extracting theta from U.

    Diagonal: H_ii = (f(U@expm(ε*H_i)) - 2f(U) + f(U@expm(-ε*H_i))) / ε^2
    Off-diagonal: H_ij = (f(U@expm(ε*H_i)@expm(ε*H_j)) - ...) / (4ε^2)
    """
    n = H_basis.shape[2]  # 45 dimensions
    H = np.zeros((n, n))
    f0 = phi_at_U(U, psi0, L)

    print(f"  Computing diagonal elements ({n} evaluations)...")
    # Diagonal elements
    for i in range(n):
        H_i = H_basis[:, :, i]
        U_p = U @ sla.expm(H_EPS * H_i)
        U_m = U @ sla.expm(-H_EPS * H_i)
        f_p = phi_at_U(U_p, psi0, L)
        f_m = phi_at_U(U_m, psi0, L)
        H[i, i] = (f_p - 2 * f0 + f_m) / (H_EPS ** 2)

    print(f"  Computing off-diagonal elements ({n*(n-1)//2} pairs × 4 evaluations)...")
    # Off-diagonal elements (symmetric, only compute upper triangle)
    for i in range(n):
        if i % 10 == 0:
            print(f"    Row {i}/{n}")
        H_i = H_basis[:, :, i]
        for j in range(i + 1, n):
            H_j = H_basis[:, :, j]

            U_pp = U @ sla.expm(H_EPS * H_i) @ sla.expm(H_EPS * H_j)
            U_pm = U @ sla.expm(H_EPS * H_i) @ sla.expm(-H_EPS * H_j)
            U_mp = U @ sla.expm(-H_EPS * H_i) @ sla.expm(H_EPS * H_j)
            U_mm = U @ sla.expm(-H_EPS * H_i) @ sla.expm(-H_EPS * H_j)

            f_pp = phi_at_U(U_pp, psi0, L)
            f_pm = phi_at_U(U_pm, psi0, L)
            f_mp = phi_at_U(U_mp, psi0, L)
            f_mm = phi_at_U(U_mm, psi0, L)

            H[i, j] = (f_pp - f_pm - f_mp + f_mm) / (4 * H_EPS ** 2)
            H[j, i] = H[i, j]  # Symmetric

    return H


def analyze_hessian(H: np.ndarray, env_name: str) -> dict:
    """Diagonalize Hessian and classify critical point."""
    eigenvalues = np.linalg.eigvalsh(H)
    eigenvalues_sorted = np.sort(eigenvalues)[::-1]  # Descending order

    # Classification
    eps = 1e-3
    n_positive = np.sum(eigenvalues > eps)
    n_zero = np.sum(np.abs(eigenvalues) <= eps)
    n_negative = np.sum(eigenvalues < -eps)

    if n_positive == 0 and n_zero == 0:
        classification = "strict local maximum"
    elif n_positive == 0 and n_zero > 0:
        classification = "flat maximum"
    elif n_positive > 0:
        classification = "saddle point"
    else:
        classification = "unknown"

    result = {
        "environment": env_name,
        "classification": classification,
        "eigenvalues": eigenvalues_sorted.tolist(),
        "lambda_max": float(eigenvalues_sorted[0]),
        "lambda_min": float(eigenvalues_sorted[-1]),
        "n_positive": int(n_positive),
        "n_zero": int(n_zero),
        "n_negative": int(n_negative),
        "condition_number": float(abs(eigenvalues_sorted[0] / eigenvalues_sorted[-1])) if eigenvalues_sorted[-1] != 0 else float("inf"),
    }

    print(f"\n{env_name}:")
    print(f"  Classification: {classification}")
    print(f"  λ_max = {result['lambda_max']:.6e}")
    print(f"  λ_min = {result['lambda_min']:.6e}")
    print(f"  #positive (> {eps}): {n_positive}")
    print(f"  #zero (|λ| ≤ {eps}): {n_zero}")
    print(f"  #negative (< -{eps}): {n_negative}")
    print(f"  Condition number: {result['condition_number']:.2e}")

    return result


def main():
    OUTDIR.mkdir(exist_ok=True)

    # Load B2 best U* for uniform/A1/A2
    b2_path = OUTDIR / "b2_canonical_open.json"
    if not b2_path.exists():
        print("ERROR: B2 results not found. Run b2_canonical_open.py first.")
        return

    with open(b2_path) as f:
        b2_data = json.load(f)

    h_total, _ = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = rem4.ground_state(h_total)

    envs = {
        "uniform_dephasing": liouvillian_dephasing(h_total, [1.0, 1.0, 1.0]),
        "A1_deph_0512": liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
        "A2_deph_2105": liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
    }

    H_basis = qgeom.horizontal_basis(D_A, D_B)
    results = {}

    for env_name, L in envs.items():
        print(f"\n=== Hessian analysis: {env_name} ===")

        # Load best U* from B2
        best_trial = max(b2_data[env_name], key=lambda t: t["best_phi"])
        U_best = np.array(best_trial["best_u_real"]) + 1j * np.array(best_trial["best_u_imag"])

        print(f"  Best Φ = {best_trial['best_phi']:.6f}")

        # Compute Hessian
        H = compute_hessian(U_best, H_basis, psi0, L)

        # Analyze
        result = analyze_hessian(H, env_name)
        results[env_name] = result

    # Save results
    out_path = OUTDIR / "b3_2_hessian_analysis.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
