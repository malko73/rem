"""Unitarity regression tests for the continuous factorization optimiser.

2026-08-11 hotfix: generators are anti-Hermitian (G_k^dagger = -G_k), so the
unitary must be U = exp(Σ θ_k G_k). The previous code used exp(i Σ θ_k G_k),
which is a positive Hermitian matrix, NOT unitary — the optimiser was not
moving on U(8). These tests lock in:
  * U^dagger U = I
  * state norm preservation: ||U psi|| = ||psi||
  * Hamiltonian spectrum preservation under U H U^dagger (unitary similarity)
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).parents[1] / "src" / "rem4_numerical.py"
spec = importlib.util.spec_from_file_location("rem4_numerical", MODULE_PATH)
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)


def _random_u_and_theta(n: int = 3, n_params: int = 4, seed: int = 1234):
    rng = np.random.default_rng(seed)
    generators = rem4._random_anti_hermitian(2**n, n_params, rng)
    theta = rng.standard_normal(n_params) * 0.1
    u = __import__("scipy.linalg", fromlist=["expm"]).expm(
        np.einsum("ijk,k->ij", generators, theta))
    return u, generators, theta


def test_generators_are_anti_hermitian():
    rng = np.random.default_rng(7)
    for dim, n_params in ((8, 4), (16, 16)):
        g = rem4._random_anti_hermitian(dim, n_params, rng)
        for k in range(n_params):
            assert np.allclose(g[..., k], -g[..., k].conj().T), \
                f"generator {k} not anti-Hermitian"


def test_u_is_unitary():
    """U = exp(sum theta_k G_k) must satisfy U^dagger U = I."""
    u, _, _ = _random_u_and_theta()
    assert np.allclose(u.conj().T @ u, np.eye(8), atol=1e-10)


def test_u_not_hermitian_for_nonzero_theta():
    """The fixed unitary must NOT be Hermitian for a non-trivial theta
    (regression: the old exp(iGθ) construction produced positive Hermitian
    matrices — unitary and Hermitian coincide only at θ=0)."""
    u, _, _ = _random_u_and_theta()
    assert not np.allclose(u, u.conj().T, atol=1e-8), \
        "U is Hermitian — suspicious: exp(iGθ) construction may be back"


def test_state_norm_preserved():
    """Unitary rotations must preserve state norm."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)
    u, _, _ = _random_u_and_theta()
    assert np.isclose(np.linalg.norm(u @ psi), np.linalg.norm(psi), rtol=1e-10)


def test_hamiltonian_spectrum_preserved():
    """U H U^dagger is a unitary similarity: spectrum must be unchanged."""
    n = 3
    h_total, _ = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    u, _, _ = _random_u_and_theta()
    spec_orig = np.linalg.eigvalsh(h_total)
    spec_rot = np.linalg.eigvalsh(u @ h_total @ u.conj().T)
    assert np.allclose(np.sort(spec_orig), np.sort(spec_rot), rtol=1e-8, atol=1e-8)


def test_identity_at_zero_theta():
    """θ=0 must give U=I exactly."""
    n = 3
    rng = np.random.default_rng(5)
    generators = rem4._random_anti_hermitian(2**n, 4, rng)
    theta = np.zeros(4)
    u = __import__("scipy.linalg", fromlist=["expm"]).expm(
        np.einsum("ijk,k->ij", generators, theta))
    assert np.allclose(u, np.eye(8), atol=1e-12)
