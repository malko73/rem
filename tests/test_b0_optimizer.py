"""
Phase B0 — Optimizer connection smoke test.

Verifies that the quotient-aware optimizer:
  * uses the 45-dimensional horizontal basis
  * preserves unitarity to machine precision at every step
  * improves Φ over the identity/contiguous baseline on the XY benchmark
"""
import numpy as np
import scipy.linalg as sla
import pytest

from quotient_optimizer import quotient_optimizer_adam
from rem4_numerical import xy_chain_hamiltonian, evaluate_factorization

N, N_A = 3, 2
DIM = 2 ** N
LAMBDA = 0.2


@pytest.fixture(scope="module")
def system():
    h_total, bond_terms = xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    psi = np.linalg.eigh(h_total)[1][:, 0]
    return psi, h_total, bond_terms


def test_optimizer_uses_45_dimensions(system):
    psi, h_total, bond_terms = system
    r = quotient_optimizer_adam(
        psi, h_total, bond_terms, n=N, n_a=N_A, lambda_value=LAMBDA,
        steps=20, lr=0.01, verbose=False,
    )
    assert r["h_basis"].shape[2] == 45


def test_optimizer_preserves_unitarity(system):
    psi, h_total, bond_terms = system
    r = quotient_optimizer_adam(
        psi, h_total, bond_terms, n=N, n_a=N_A, lambda_value=LAMBDA,
        steps=20, lr=0.01, verbose=False,
    )
    assert r["unitarity_errors"].max() < 1e-10


def test_optimizer_improves_phi(system):
    psi, h_total, bond_terms = system
    # Contiguous baseline at identity (cut AB|C)
    I = np.eye(DIM, dtype=complex)
    phi_id, _, _ = evaluate_factorization(
        psi, h_total, bond_terms, I, n=N, n_a=N_A, lambda_value=LAMBDA
    )
    # Short optimizer run (deterministic seed)
    r = quotient_optimizer_adam(
        psi, h_total, bond_terms, n=N, n_a=N_A, lambda_value=LAMBDA,
        steps=100, lr=0.01, verbose=False,
    )
    assert r["best_phi"] > phi_id + 1e-6, (
        f"optimizer did not improve: Φ_opt={r['best_phi']} ≤ Φ_id={phi_id}"
    )
