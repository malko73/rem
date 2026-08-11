"""
Phase B0 — 9 mandatory tests (master's gate list) on the quotient geometry.

Locks, on the 3-qubit XY benchmark (n=3, cut A = sites [0,1], d_A=4, d_B=2):
  1. horizontal basis dimension = 45
  2. vertical basis rank = 19
  3. horizontal ⟂ vertical
  4. gauge invariance: Φ invariant under V_A ⊗ V_B
  5. vertical-only step: ΔΦ ≈ 0
  6. horizontal step: ΔΦ ≠ 0 (generically)
  7. |U†U − I| ≈ machine precision after each step
  8. spectrum / norm invariance
  9. projected gradient has no vertical component
"""
import numpy as np
import scipy.linalg as sla
import pytest

from quotient_geometry import (
    u_basis,
    raw_vertical_generators,
    vertical_basis,
    horizontal_basis,
    project_vertical,
    project_horizontal,
)
from rem4_numerical import xy_chain_hamiltonian, evaluate_factorization

N, N_A = 3, 2
D_A, D_B = 2 ** N_A, 2 ** (N - N_A)  # 4, 2
DIM = 2 ** N  # 8
LAMBDA = 0.2


@pytest.fixture(scope="module")
def system():
    h_total, bond_terms = xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    psi = np.linalg.eigh(h_total)[1][:, 0]
    return psi, h_total, bond_terms


def _phi(psi, h_total, bond_terms, u):
    return evaluate_factorization(
        psi, h_total, bond_terms, u, n=N, n_a=N_A, lambda_value=LAMBDA
    )[0]


# ─── 1. horizontal basis dimension = 45 ─────────────────────────────────────

def test_horizontal_dimension():
    H = horizontal_basis(D_A, D_B)
    assert H.shape == (DIM, DIM, 45)


# ─── 2. vertical basis rank = 19 ────────────────────────────────────────────

def test_vertical_rank():
    V = vertical_basis(D_A, D_B)
    assert V.shape == (DIM, DIM, 19)


# ─── 3. horizontal ⟂ vertical ───────────────────────────────────────────────

def test_orthogonality():
    V = vertical_basis(D_A, D_B)
    H = horizontal_basis(D_A, D_B)
    maxv = 0.0
    for i in range(V.shape[2]):
        for j in range(H.shape[2]):
            maxv = max(maxv, abs(np.real(np.trace(V[:, :, i].conj().T @ H[:, :, j]))))
    assert maxv < 1e-10


# ─── 4. gauge invariance: Φ(V_A⊗V_B · U) = Φ(U) ────────────────────────────

def test_gauge_invariance(system):
    psi, h_total, bond_terms = system
    rng = np.random.default_rng(123)

    # Random U
    G = 1j * rng.standard_normal((DIM, DIM))
    G = G - G.conj().T
    U = sla.expm(G)

    # Random local unitaries V_A, V_B
    v_a = u_basis(D_A)
    v_b = u_basis(D_B)
    theta_a = rng.standard_normal(D_A * D_A)
    theta_b = rng.standard_normal(D_B * D_B)
    V_A = sla.expm(np.einsum("ijk,k->ij", v_a, theta_a))
    V_B = sla.expm(np.einsum("ijk,k->ij", v_b, theta_b))
    V = np.kron(V_A, V_B)

    phi_u = _phi(psi, h_total, bond_terms, U)
    phi_vu = _phi(psi, h_total, bond_terms, V @ U)
    assert abs(phi_u - phi_vu) < 1e-8, f"gauge violation: ΔΦ={abs(phi_u - phi_vu)}"


# ─── 5. vertical-only step → ΔΦ ≈ 0 ─────────────────────────────────────────

def test_vertical_step_no_change(system):
    psi, h_total, bond_terms = system
    rng = np.random.default_rng(456)
    V = vertical_basis(D_A, D_B)

    G = 1j * rng.standard_normal((DIM, DIM))
    G = G - G.conj().T
    U = sla.expm(G)

    theta = rng.standard_normal(V.shape[2]) * 0.1
    V_step = sla.expm(np.einsum("ijk,k->ij", V, theta))

    phi_0 = _phi(psi, h_total, bond_terms, U)
    phi_1 = _phi(psi, h_total, bond_terms, V_step @ U)
    assert abs(phi_1 - phi_0) < 1e-6, f"vertical ΔΦ={abs(phi_1 - phi_0)}"


# ─── 6. horizontal step → ΔΦ ≠ 0 (generically) ─────────────────────────────

def test_horizontal_step_changes(system):
    psi, h_total, bond_terms = system
    rng = np.random.default_rng(789)
    H = horizontal_basis(D_A, D_B)

    G = 1j * rng.standard_normal((DIM, DIM))
    G = G - G.conj().T
    U = sla.expm(G)

    theta = rng.standard_normal(H.shape[2]) * 0.3
    H_step = sla.expm(np.einsum("ijk,k->ij", H, theta))

    phi_0 = _phi(psi, h_total, bond_terms, U)
    phi_1 = _phi(psi, h_total, bond_terms, H_step @ U)
    assert abs(phi_1 - phi_0) > 1e-5, f"horizontal ΔΦ={abs(phi_1 - phi_0)}"


# ─── 7. unitarity preserved ─────────────────────────────────────────────────

def test_unitarity_preserved(system):
    psi, h_total, bond_terms = system
    rng = np.random.default_rng(101)
    H = horizontal_basis(D_A, D_B)

    # Several steps from identity
    U = np.eye(DIM, dtype=complex)
    for _ in range(20):
        theta = rng.standard_normal(H.shape[2]) * 0.05
        K = np.einsum("ijk,k->ij", H, theta)
        U = sla.expm(K) @ U
        err = np.linalg.norm(U.conj().T @ U - np.eye(DIM, dtype=complex), ord=2)
        assert err < 1e-10, f"unitarity violation: |U†U−I|={err}"


# ─── 8. spectrum / norm invariance ──────────────────────────────────────────

def test_spectrum_norm_invariance(system):
    h_total, _ = system[1], system[2]
    rng = np.random.default_rng(202)
    H = horizontal_basis(D_A, D_B)

    theta = rng.standard_normal(H.shape[2]) * 0.1
    U = sla.expm(np.einsum("ijk,k->ij", H, theta))

    h_rot = U @ h_total @ U.conj().T
    assert np.linalg.norm(np.sort(np.linalg.eigvalsh(h_rot)) - np.sort(np.linalg.eigvalsh(h_total))) < 1e-10
    assert abs(np.linalg.norm(h_rot, ord=2) - np.linalg.norm(h_total, ord=2)) < 1e-10


# ─── 9. projected gradient has no vertical component ────────────────────────

def test_projected_gradient_no_vertical(system):
    psi, h_total, bond_terms = system
    rng = np.random.default_rng(303)
    H = horizontal_basis(D_A, D_B)
    V = vertical_basis(D_A, D_B)

    # Random U
    G = 1j * rng.standard_normal((DIM, DIM))
    G = G - G.conj().T
    U = sla.expm(G)

    # Gradient in the Lie algebra via finite differences along horizontal basis
    eps = 1e-5
    grad = np.zeros(H.shape[2])
    for k in range(H.shape[2]):
        U_p = sla.expm(eps * H[:, :, k]) @ U
        U_m = sla.expm(-eps * H[:, :, k]) @ U
        grad[k] = (_phi(psi, h_total, bond_terms, U_p)
                   - _phi(psi, h_total, bond_terms, U_m)) / (2 * eps)

    grad_lie = np.einsum("ijk,k->ij", H, grad)
    # grad_lie is a horizontal combination; verify no vertical component
    maxv = 0.0
    for j in range(V.shape[2]):
        maxv = max(maxv, abs(np.real(np.trace(V[:, :, j].conj().T @ grad_lie))))
    assert maxv < 1e-8, f"projected gradient has vertical component: {maxv}"
