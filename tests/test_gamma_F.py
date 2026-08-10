"""Phase C0: Gamma_F measurement — unit tests.

Lock in the mathematical invariants the Lindblad solver and the
factorisation-specific coherence rely on:
  * trace preservation: Tr(rho(t)) = 1
  * Hermiticity preservation: rho(t)^dagger = rho(t)
  * pure dephasing conserves <Z_i> (diagonal elements untouched)
  * D_F is a projection: D_F[D_F[rho]] = D_F[rho]
  * C_F(0) > 0 (non-trivial coherence exists for both cuts)
  * Gamma^(0) and Gamma^fit are finite and mutually consistent
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

ANALYSIS_PATH = Path(__file__).parents[1] / "analysis" / "gamma_F.py"
aspec = importlib.util.spec_from_file_location("gamma_F", ANALYSIS_PATH)
gf = importlib.util.module_from_spec(aspec)
assert aspec.loader is not None
sys.modules[aspec.name] = gf
aspec.loader.exec_module(gf)


def _setup():
    psi, h_total = gf.build_system()
    L = gf.liouvillian(h_total, gf.GAMMA, n=3)
    t_grid = np.linspace(0.0, 1.0, 21)
    return psi, h_total, L, t_grid


def test_liouvillian_trace_preserving():
    _, _, L, t_grid = _setup()
    rho0 = np.outer(_setup()[0], _setup()[0].conj())
    rho_t = gf.solve_dynamics(L, rho0, t_grid)
    for i in range(0, len(t_grid), 5):
        rho = rho_t[i].reshape((8, 8), order="F")
        assert np.isclose(np.trace(rho).real, 1.0, atol=1e-8), f"t={t_grid[i]}"


def test_liouvillian_hermiticity_preserving():
    _, _, L, t_grid = _setup()
    psi, _ = gf.build_system()
    rho0 = np.outer(psi, psi.conj())
    rho_t = gf.solve_dynamics(L, rho0, t_grid)
    for i in range(0, len(t_grid), 5):
        rho = rho_t[i].reshape((8, 8), order="F")
        assert np.allclose(rho, rho.conj().T, atol=1e-8), f"t={t_grid[i]}"


def test_pure_dephasing_conserves_z_expectation_without_hamiltonian():
    """With H=0, local Z dephasing conserves <Z_i> (Z-diagonal invariant)."""
    psi, h_total = gf.build_system()
    L_deph = gf.liouvillian(np.zeros_like(h_total), gf.GAMMA, n=3)  # H=0
    t_grid = np.linspace(0.0, 1.0, 21)
    rho0 = np.outer(psi, psi.conj())
    rho_t = gf.solve_dynamics(L_deph, rho0, t_grid)
    for site in range(3):
        z = gf.pauli_z(3, site)
        z0 = np.real(np.trace(rho0 @ z))
        for i in (0, 5, 10, 20):
            rho = rho_t[i].reshape((8, 8), order="F")
            zt = np.real(np.trace(rho @ z))
            assert np.isclose(zt, z0, atol=1e-6), f"site {site}, t={t_grid[i]}"


def test_dephasing_dissipator_annihilates_z_diagonal():
    """The pure-dephasing dissipator leaves Z-diagonal matrices invariant."""
    rng = np.random.default_rng(0)
    rho_diag = np.diag(rng.random(8) + 0.1)
    rho_diag = rho_diag / np.trace(rho_diag)
    L_deph = gf.liouvillian(np.zeros((8, 8), dtype=complex), gf.GAMMA, n=3)
    out = (L_deph @ rho_diag.reshape(-1, order="F")).reshape((8, 8), order="F")
    off = out - np.diag(np.diag(out))
    assert np.allclose(off, 0.0, atol=1e-10), "dissipator mixes off-diagonals"


def test_dephasing_projection_is_idempotent():
    psi, _ = gf.build_system()
    rho = np.outer(psi, psi.conj())
    for cut in (1, 2):
        basis = gf.schmidt_basis(psi, 3, cut)
        d1 = gf.dephasing_projection(rho, basis)
        d2 = gf.dephasing_projection(d1, basis)
        assert np.allclose(d1, d2, atol=1e-10), f"cut {cut}: D_F^2 != D_F"


def test_coherence_positive_for_both_cuts():
    """C_F(0) > 0: both cuts have non-trivial coherence in the Schmidt basis."""
    psi, _ = gf.build_system()
    rho = np.outer(psi, psi.conj())
    for cut in (1, 2):
        basis = gf.schmidt_basis(psi, 3, cut)
        c0 = gf.coherence_norm(rho, basis)
        assert c0 > 1e-6, f"cut {cut}: C_F(0) = {c0} (no coherence?)"


def test_gamma_measurements_finite_and_consistent():
    """run_cut must return finite Gamma^(0) and Gamma^fit, close to each other."""
    psi, h_total = gf.build_system()
    L = gf.liouvillian(h_total, gf.GAMMA, n=3)
    t_grid = np.linspace(0.0, 1.0, 101)
    for cut in (1, 2):
        d = gf.run_cut(psi, h_total, gf.GAMMA, t_grid, L, cut)
        assert np.isfinite(d["gamma0"]), f"cut {cut}: Gamma^(0) not finite"
        assert np.isfinite(d["gamma_fit"]), f"cut {cut}: Gamma^fit not finite"
        # for the short window the two definitions should agree within 50%
        assert abs(d["gamma0"] - d["gamma_fit"]) <= 0.5 * max(
            abs(d["gamma0"]), abs(d["gamma_fit"]), 1e-9), \
            f"cut {cut}: Gamma^(0)={d['gamma0']:.4f} vs Gamma^fit={d['gamma_fit']:.4f}"


def test_classify_returns_case():
    # g1 < g2 -> case 3 (Gamma_1 < Gamma_2, falsification)
    d1 = dict(gamma0=0.5, gamma_fit=0.5, gamma_exact=0.5)
    d2 = dict(gamma0=1.0, gamma_fit=1.0, gamma_exact=1.0)
    assert gf.classify_case(d1, d2)["verdict"] == 3
    # g1 > g2 -> case 1 (Gamma_2 < Gamma_1, supports M2)
    d1 = dict(gamma0=1.0, gamma_fit=1.0, gamma_exact=1.0)
    d2 = dict(gamma0=0.5, gamma_fit=0.5, gamma_exact=0.5)
    assert gf.classify_case(d1, d2)["verdict"] == 1
    # equal -> case 2 (Var degeneracy consistent)
    d1 = dict(gamma0=0.5, gamma_fit=0.5, gamma_exact=0.5)
    d2 = dict(gamma0=0.5, gamma_fit=0.5, gamma_exact=0.5)
    assert gf.classify_case(d1, d2)["verdict"] == 2
