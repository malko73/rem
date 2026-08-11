"""
Phase B2a — Well-posedness regression tests for C_Γ^(0) on the quotient.

Locks the fixes found by the B2a audit:
  1. Liouvillian commutator uses the correct vec_F convention
     (kron(I,H) - kron(H.conj(), I)); the old reversed-kron version was
     wrong for general (non-ground) states.
  2. Rotated-environment dephasing uses kron(z.T, z) (z^rot = U†ZU is
     Hermitian but not symmetric; vec_F(z ρ z) = (z^T ⊗ z) vec_F(ρ)).
  3. Frame consistency: pull-back (U†ρ₀U, U†Y₀U) equals lab-frame with
     rotated H and Z operators.
"""
import numpy as np
import scipy.linalg as sla
import pytest

from quotient_geometry import horizontal_basis
from gamma_F import (
    liouvillian_env,
    pauli_z,
    schmidt_basis,
    dephasing_projection,
    gamma_exact,
    build_system,
)
from rem4_numerical import xy_chain_hamiltonian, ground_state

N, N_A = 3, 2
D_A, D_B = 4, 2
DIM = 8
GAMMA_VEC = [1.0, 1.0, 1.0]


def _vectorized_L(h_total):
    """Commutator-only Liouvillian for a Hermitian h (vec_F convention)."""
    I = np.eye(DIM, dtype=complex)
    return -1j * (np.kron(I, h_total) - np.kron(h_total.conj(), I))


def test_liouvillian_commutator_matches_direct():
    """The vectorised commutator must equal -i[H, rho] for a general state."""
    h_total, _ = xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    rng = np.random.default_rng(7)
    # general (not ground) state
    psi = rng.standard_normal(DIM) + 1j * rng.standard_normal(DIM)
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())

    L = _vectorized_L(h_total)
    y_vec = (L @ rho.reshape(-1, order="F")).reshape(rho.shape, order="F")
    y_dir = -1j * (h_total @ rho - rho @ h_total)
    assert np.linalg.norm(y_vec - y_dir) < 1e-10


def test_liouvillian_env_commutator_correct():
    """liouvillian_env must match direct Lindblad for a general state."""
    h_total, _ = xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    rng = np.random.default_rng(8)
    psi = rng.standard_normal(DIM) + 1j * rng.standard_normal(DIM)
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())

    L = liouvillian_env(h_total, GAMMA_VEC, [0.0] * N, n=N)
    y_vec = (L @ rho.reshape(-1, order="F")).reshape(rho.shape, order="F")
    # direct
    I = np.eye(DIM, dtype=complex)
    y_dir = -1j * (h_total @ rho - rho @ h_total)
    for site in range(N):
        z = pauli_z(N, site)
        y_dir += (GAMMA_VEC[site] / 2.0) * (z @ rho @ z - rho)
    assert np.linalg.norm(y_vec - y_dir) < 1e-10


def test_frame_consistency_rotated_environment():
    """Pull-back equals lab-frame when the environment is rotated too."""
    h_total, _ = xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = ground_state(h_total)
    rho0 = np.outer(psi0, psi0.conj())

    H_basis = horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(20260815)
    theta = rng.standard_normal(H_basis.shape[2]) * 0.5
    U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))

    # Pull-back: Y_U = U† L(rho0) U
    L0 = liouvillian_env(h_total, GAMMA_VEC, [0.0] * N, n=N)
    Y0 = (L0 @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    Y_pull = U.conj().T @ Y0 @ U

    # Lab-frame with rotated H and Z (note kron(z.T, z) for z^rot)
    rho_rot = U.conj().T @ rho0 @ U
    h_rot = U.conj().T @ h_total @ U
    I = np.eye(DIM, dtype=complex)
    L_rot = -1j * (np.kron(I, h_rot) - np.kron(h_rot.conj(), I))
    for site in range(N):
        z = U.conj().T @ pauli_z(N, site) @ U
        L_rot += (GAMMA_VEC[site] / 2.0) * (np.kron(z.T, z) - np.kron(I, I))
    Y_lab = (L_rot @ rho_rot.reshape(-1, order="F")).reshape(rho_rot.shape, order="F")

    assert np.linalg.norm(Y_lab - Y_pull) < 1e-10


@pytest.mark.xfail(
    strict=False,
    reason=("B2a finding 2026-08-11: Γ^exact is NOT invariant under rotation of "
            "the Schmidt basis within a NON-ZERO degenerate block (p1=p2). "
            "Null-subspace completion invariance (C1.5) does NOT extend to "
            "non-zero degeneracy; D_F depends on the SVD basis choice there. "
            "Open question for the Spec: pinching-style projection onto the "
            "full degenerate block, or accept a domain restriction."),
)
def test_gamma_exact_gauge_invariance_nonzero_schmidt():
    """Γ^exact invariance under degenerate-block Schmidt rotations.

    KNOWN FAILURE (xfail): for exactly degenerate Schmidt values p₁=p₂,
    the Schmidt basis is non-unique and D_F (hence C_F² and Γ) depends on
    the SVD's basis choice. This is a real well-posedness limitation of
    C_Γ^(0) that B2a uncovered; it must be resolved before the canonical
    open-system functional is finalized.
    """
    h_total, _ = xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    L = liouvillian_env(h_total, GAMMA_VEC, [0.0] * N, n=N)

    # Bell-like state: p1 = p2 = 1/2 on cut A|B
    psi = np.zeros(DIM, dtype=complex)
    psi[0b000] = 1.0 / np.sqrt(2.0)
    psi[0b011] = 1.0 / np.sqrt(2.0)
    rho0 = np.outer(psi, psi.conj())

    basis1 = schmidt_basis(psi, N, N_A)
    g1 = gamma_exact(rho0, L, basis1)

    rng = np.random.default_rng(20260814)
    m = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
    q, _ = np.linalg.qr(m)
    basis2 = basis1.copy()
    for l in range(D_B):
        col0 = basis1[:, l].copy()
        col1 = basis1[:, D_B + l].copy()
        basis2[:, l] = q[0, 0] * col0 + q[1, 0] * col1
        basis2[:, D_B + l] = q[0, 1] * col0 + q[1, 1] * col1
    g2 = gamma_exact(rho0, L, basis2)

    assert abs(g1 - g2) < 1e-8, f"Γ depends on degenerate Schmidt basis: {g1} vs {g2}"
