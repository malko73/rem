#!/usr/bin/env python3
"""
Phase B2 — Canonical Open REM continuous TPS optimization.

Purpose: test whether F* = F*(ρ, L, λ) by running the canonical open-system
objective Φ(F; λ, ρ, L) = I_ρ(F) - λ Γ_F^exact(0) on the 45-dim quotient
for the 5 C1 environments.

Environments (from C1):
  1. uniform dephasing (γ=1.0)
  2. A1: dephasing γ=(0.5, 1.0, 2.0)
  3. A2: dephasing γ=(2.0, 1.0, 0.5)  [mirror of A1]
  4. amplitude damping κ=1.0
  5. mixed γ=κ=0.5

For each environment, run 100 trials (50 SGD + 50 Adam) and compute:
  - best/median/mean/std/min Φ
  - success rate (vs contiguous best)
  - distinct maxima estimate
  - convergence distribution

Key test: do A1 and A2 (mirror environments) yield different F*?
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

N = 3
N_A = 2
D_A, D_B = 2 ** N_A, 2 ** (N - N_A)
DIM = 2 ** N
COUPLINGS = [1.5, 0.6]
H_FIELD = 0.2
LAMBDA = 0.2
STEPS = 200
LR = 0.01
N_TRIALS_PER_ENV = 100
SEED0 = 20260816


def canonical_open_objective(rho_U: np.ndarray, Y_U: np.ndarray, n_a: int) -> dict:
    """Compute I, C_F², Γ, Φ for the canonical open-system objective."""
    evals_psi, evecs_psi = np.linalg.eigh(rho_U)
    psi_U = evecs_psi[:, np.argmax(evals_psi)]
    I_F = rem4.mutual_information_for_cut(psi_U, N, n_a)

    # Schmidt basis of the rotated state
    basis = gamma_F.schmidt_basis(psi_U, N, n_a)
    q = lambda r: r - gamma_F.dephasing_projection(r, basis)  # noqa: E731
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
    """Build vectorised Lindblad Liouvillian L for pure dephasing."""
    dim = 2 ** N
    I = np.eye(dim, dtype=complex)
    L = -1j * (np.kron(I, h_total) - np.kron(h_total.conj(), I))
    for site in range(N):
        g = gamma_vec[site]
        if g > 0:
            z = gamma_F.pauli_z(N, site)
            L += (g / 2.0) * (np.kron(z, z) - np.kron(I, I))
    return L


def liouvillian_amp_damping(h_total: np.ndarray, kappa_vec: list[float]) -> np.ndarray:
    """Build vectorised Lindblad Liouvillian L for amplitude damping."""
    dim = 2 ** N
    I = np.eye(dim, dtype=complex)
    L = -1j * (np.kron(I, h_total) - np.kron(h_total.conj(), I))
    for site in range(N):
        k = kappa_vec[site]
        if k > 0:
            sm = gamma_F.sigma_minus(N, site)
            spsm = sm.conj().T @ sm
            L += k * (np.kron(sm.conj(), sm)
                      - 0.5 * np.kron(I, spsm)
                      - 0.5 * np.kron(spsm.conj(), I))
    return L


def liouvillian_mixed(h_total: np.ndarray, gamma_vec: list[float], kappa_vec: list[float]) -> np.ndarray:
    """Build vectorised Lindblad Liouvillian L for mixed dephasing + damping."""
    dim = 2 ** N
    I = np.eye(dim, dtype=complex)
    L = -1j * (np.kron(I, h_total) - np.kron(h_total.conj(), I))
    for site in range(N):
        g = gamma_vec[site]
        k = kappa_vec[site]
        if g > 0:
            z = gamma_F.pauli_z(N, site)
            L += (g / 2.0) * (np.kron(z, z) - np.kron(I, I))
        if k > 0:
            sm = gamma_F.sigma_minus(N, site)
            spsm = sm.conj().T @ sm
            L += k * (np.kron(sm.conj(), sm)
                      - 0.5 * np.kron(I, spsm)
                      - 0.5 * np.kron(spsm.conj(), I))
    return L


def apply_tps(psi0: np.ndarray, L: np.ndarray, U: np.ndarray):
    """Return (rho_U, Y_U) for a TPS U on the given environment."""
    rho0 = np.outer(psi0, psi0.conj())
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    rho_U = U.conj().T @ rho0 @ U
    Y_U = U.conj().T @ Y0 @ U
    return rho_U, Y_U


def run_environment(env_name: str, L: np.ndarray, n_trials: int = N_TRIALS_PER_ENV):
    """Run optimization for one environment."""
    print(f"\n=== Environment: {env_name} ===")
    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    _, psi0 = rem4.ground_state(h_total)

    H_basis = qgeom.horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(SEED0)

    results = []
    t0 = time.time()
    for trial in range(n_trials):
        theta = rng.standard_normal(H_basis.shape[2]) * 1.0
        U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))
        rho_U, Y_U = apply_tps(psi0, L, U)

        # Evaluate objective
        obj = canonical_open_objective(rho_U, Y_U, N_A)

        # Short optimization (Adam, 50 steps)
        from rem4_numerical import _finite_diff_grad
        theta_opt = theta.copy()
        m = np.zeros_like(theta_opt)
        v = np.zeros_like(theta_opt)
        beta1, beta2, eps = 0.9, 0.999, 1e-8
        best_phi = obj["Phi"]
        for step in range(50):
            # Gradient via finite differences
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

            # Adam update
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

        results.append(dict(
            trial=trial,
            initial_phi=obj["Phi"],
            best_phi=best_phi,
            C_F_sq=obj["C_F_sq"],
            gamma=obj["gamma"],
        ))

        if (trial + 1) % 20 == 0:
            elapsed = time.time() - t0
            print(f"  trial {trial+1}/{n_trials} ({elapsed:.1f}s)")

    return results


def main():
    OUTDIR.mkdir(exist_ok=True)

    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    _, psi0 = rem4.ground_state(h_total)

    # Contiguous baselines
    I = np.eye(DIM, dtype=complex)
    rho0 = np.outer(psi0, psi0.conj())
    I_cut2 = rem4.mutual_information_for_cut(psi0, N, N_A)
    contiguous_best = I_cut2  # Approximate; actual depends on environment

    # 5 environments
    envs = {
        "uniform_dephasing": liouvillian_dephasing(h_total, [1.0, 1.0, 1.0]),
        "A1_deph_0512": liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
        "A2_deph_2105": liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
        "amp_damping": liouvillian_amp_damping(h_total, [1.0, 1.0, 1.0]),
        "mixed": liouvillian_mixed(h_total, [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),
    }

    all_results = {}
    for env_name, L in envs.items():
        results = run_environment(env_name, L, n_trials=100)  # Full statistics
        all_results[env_name] = results

        # Statistics
        phis = np.array([r["best_phi"] for r in results if np.isfinite(r["best_phi"])])
        if len(phis) > 0:
            hist, edges = np.histogram(phis, bins=np.arange(np.floor(phis.min() * 20) / 20,
                                                            np.ceil(phis.max() * 20) / 20 + 0.05,
                                                            0.05))
            distinct_maxima = int(np.sum(hist > 0))
            print(f"\n{env_name}:")
            print(f"  best={phis.max():.4f}  median={np.median(phis):.4f}  "
                  f"mean={phis.mean():.4f}  std={phis.std():.4f}")
            print(f"  distinct maxima: {distinct_maxima}")
        else:
            print(f"\n{env_name}: no finite results")

    # Save results
    out_path = OUTDIR / "b2_canonical_open.json"
    out_path.write_text(json.dumps(all_results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
