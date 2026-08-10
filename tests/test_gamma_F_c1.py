"""Phase C0.1/C1 — unit tests.

Covers:
  * Gamma^exact (analytic t=0) matches a numerical derivative of log C_F(t)
  * amplitude-damping Liouvillian preserves trace
  * C1 environment classification (C1-A / C1-B / C1-C)
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

GPATH = Path(__file__).parents[1] / "analysis" / "gamma_F.py"
gspec = importlib.util.spec_from_file_location("gamma_F", GPATH)
gf = importlib.util.module_from_spec(gspec)
assert gspec.loader is not None
sys.modules[gspec.name] = gf
gspec.loader.exec_module(gf)

C1PATH = Path(__file__).parents[1] / "analysis" / "gamma_F_c1.py"
cspec = importlib.util.spec_from_file_location("gamma_F_c1", C1PATH)
c1 = importlib.util.module_from_spec(cspec)
assert cspec.loader is not None
sys.modules[cspec.name] = c1
cspec.loader.exec_module(c1)


def test_gamma_exact_matches_numerical_derivative():
    """Analytic t=0 rate equals -(d/dt log C)|_0 by central difference."""
    psi, h_total = gf.build_system()
    L = gf.liouvillian(h_total, gf.GAMMA, n=3)
    rho0 = np.outer(psi, psi.conj())
    Lrho = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    h = 1e-6
    for cut in (1, 2):
        basis = gf.schmidt_basis(psi, 3, cut)
        g_exact = gf.gamma_exact(rho0, L, basis)
        qq = lambda r: r - gf.dephasing_projection(r, basis)          # noqa: E731
        c_p = np.linalg.norm(qq(rho0 + h * Lrho), ord="fro")
        c_m = np.linalg.norm(qq(rho0 - h * Lrho), ord="fro")
        dlog = (np.log(c_p) - np.log(c_m)) / (2 * h)
        assert np.isclose(g_exact, -dlog, rtol=1e-4), \
            f"cut {cut}: exact={g_exact:.6f} vs numeric={-dlog:.6f}"


def test_amplitude_damping_trace_preserving():
    psi, h_total = gf.build_system()
    L = gf.liouvillian_env(h_total, [0.0] * 3, [1.0] * 3, n=3)
    rho0 = np.outer(psi, psi.conj())
    rho_t = gf.solve_dynamics(L, rho0, np.linspace(0.0, 1.0, 11))
    for i in range(len(rho_t)):
        rho = rho_t[i].reshape((8, 8), order="F")
        assert np.isclose(np.trace(rho).real, 1.0, atol=1e-8), f"t index {i}"


def test_mixed_env_trace_preserving():
    psi, h_total = gf.build_system()
    L = gf.liouvillian_env(h_total, [0.5] * 3, [0.5] * 3, n=3)
    rho0 = np.outer(psi, psi.conj())
    rho_t = gf.solve_dynamics(L, rho0, np.linspace(0.0, 1.0, 11))
    for i in range(len(rho_t)):
        rho = rho_t[i].reshape((8, 8), order="F")
        assert np.isclose(np.trace(rho).real, 1.0, atol=1e-8), f"t index {i}"


def _fake_results(ratios):
    """Build minimal results dicts for classification tests."""
    return {f"env{i}": dict(r_gamma_exact=r, cut1=dict(m2=8.379, var=0.62),
                            cut2=dict(m2=0.819, var=0.62))
            for i, r in enumerate(ratios)}


def test_c1_classify_flat_is_B():
    v = c1.classify_c1(_fake_results([1.02, 0.99, 1.01]))
    assert v["verdict"] == "C1-B"


def test_c1_classify_all_m2_is_A():
    v = c1.classify_c1(_fake_results([1.5, 1.3, 1.4]))
    assert v["verdict"] == "C1-A"


def test_c1_classify_env_dependent_is_C():
    v = c1.classify_c1(_fake_results([1.02, 1.5, 0.99]))
    assert v["verdict"] == "C1-C"
