"""
Phase B0 — Quotient-aware optimizer (45-dim horizontal TPS search).

Connects the SVD-based quotient geometry (src/quotient_geometry.py) to the
continuous TPS optimizer: the search runs on the 45-dimensional horizontal
basis, so the gauge directions (vertical, dim 19) are excluded by
construction — no random 4-generator subspace, no manual kernel removal.

The old random-subspace optimizers (optimize_factorization /
optimize_factorization_adam in rem4_numerical.py) remain for backward
compatibility with the v4/v5 regression baselines; B1 will re-run the
closed-system regression with this optimizer.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla

from quotient_geometry import horizontal_basis, vertical_basis, project_horizontal

N_RNG = np.random.default_rng(seed=20260811)


def quotient_optimizer_sgd(
    psi: np.ndarray,
    h_total: np.ndarray,
    bond_terms: dict,
    *,
    n: int = 3,
    n_a: int = 2,
    lambda_value: float = 0.2,
    steps: int = 200,
    lr: float = 0.01,
    verbose: bool = True,
) -> dict:
    """
    Plain gradient-ascent (SGD) on Φ over the 45-dim horizontal quotient.

    Returns
    -------
    dict with keys: best_phi, best_u, theta, history, h_basis,
    unitarity_errors (per step).
    """
    from rem4_numerical import evaluate_factorization, _finite_diff_grad

    d_a, d_b = 2 ** n_a, 2 ** (n - n_a)
    h_basis = horizontal_basis(d_a, d_b)
    n_params = h_basis.shape[2]

    theta = N_RNG.standard_normal(n_params) * 0.1

    best_phi = -1e9
    best_u = None
    history = []
    unitarity_errors = []

    for step in range(steps):
        u = sla.expm(np.einsum("ijk,k->ij", h_basis, theta))
        phi, mi, c_h = evaluate_factorization(
            psi, h_total, bond_terms, u, n=n, n_a=n_a, lambda_value=lambda_value
        )
        grad = _finite_diff_grad(
            theta, psi, h_total, bond_terms, h_basis,
            n=n, n_a=n_a, lambda_value=lambda_value,
        )
        theta = theta + lr * grad

        history.append((phi, mi, c_h))
        unitarity_errors.append(
            float(np.linalg.norm(u.conj().T @ u - np.eye(2 ** n, dtype=complex), ord=2))
        )

        if phi > best_phi:
            best_phi = phi
            best_u = u

        if verbose and step % 50 == 0:
            print(f"  [{step:3d}] Φ={phi:.6f}  MI={mi:.6f}  C_H={c_h:.6f}")

    return dict(
        best_phi=best_phi,
        best_u=best_u,
        theta=theta,
        history=np.array(history),
        h_basis=h_basis,
        unitarity_errors=np.array(unitarity_errors),
    )


def quotient_optimizer_adam(
    psi: np.ndarray,
    h_total: np.ndarray,
    bond_terms: dict,
    *,
    n: int = 3,
    n_a: int = 2,
    lambda_value: float = 0.2,
    steps: int = 200,
    lr: float = 0.01,
    beta1: float = 0.9,
    beta2: float = 0.999,
    eps: float = 1e-8,
    tol_grad: float | None = None,
    tol_phi: float | None = None,
    patience: int = 5,
    verbose: bool = True,
) -> dict:
    """
    Adam gradient-ascent on Φ over the 45-dim horizontal quotient.

    Parameters
    ----------
    tol_grad / tol_phi / patience : early stopping (same semantics as
        rem4_numerical.optimize_factorization_adam).

    Returns
    -------
    dict with keys:
        best_phi, best_u, theta, history, h_basis,
        unitarity_errors (per step), converged_at
    """
    from rem4_numerical import evaluate_factorization, _finite_diff_grad

    d_a, d_b = 2 ** n_a, 2 ** (n - n_a)
    h_basis = horizontal_basis(d_a, d_b)          # (D, D, 45)
    n_params = h_basis.shape[2]

    theta = N_RNG.standard_normal(n_params) * 0.1
    m = np.zeros_like(theta)
    v = np.zeros_like(theta)

    best_phi = -1e9
    best_u = None
    history = []
    unitarity_errors = []
    converged_at = steps

    phi_window = [] if (tol_grad is not None and tol_phi is not None) else None

    for step in range(steps):
        u = sla.expm(np.einsum("ijk,k->ij", h_basis, theta))
        phi, mi, c_h = evaluate_factorization(
            psi, h_total, bond_terms, u, n=n, n_a=n_a, lambda_value=lambda_value
        )
        grad = _finite_diff_grad(
            theta, psi, h_total, bond_terms, h_basis,
            n=n, n_a=n_a, lambda_value=lambda_value,
        )

        # Adam update
        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        theta = theta + lr * m_hat / (np.sqrt(v_hat) + eps)

        history.append((phi, mi, c_h))
        unitarity_errors.append(
            float(np.linalg.norm(u.conj().T @ u - np.eye(2 ** n, dtype=complex), ord=2))
        )

        if phi > best_phi:
            best_phi = phi
            best_u = u

        # Early stopping
        if phi_window is not None:
            phi_window.append(phi)
            if len(phi_window) > patience:
                phi_window.pop(0)
            grad_norm = float(np.linalg.norm(grad))
            if (grad_norm < tol_grad
                    and len(phi_window) == patience
                    and max(phi_window) - min(phi_window) < tol_phi):
                converged_at = step + 1
                if verbose:
                    print(f"  early stop @ step {step+1} (|grad|={grad_norm:.2e})")
                break

        if verbose and step % 50 == 0:
            print(f"  [{step:3d}] Φ={phi:.6f}  MI={mi:.6f}  C_H={c_h:.6f}")

    return dict(
        best_phi=best_phi,
        best_u=best_u,
        theta=theta,
        history=np.array(history),
        h_basis=h_basis,
        unitarity_errors=np.array(unitarity_errors),
        converged_at=converged_at,
    )
