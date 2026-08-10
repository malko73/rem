"""Phase C2 — Liouvillian boundary cost: unit tests.

Lock in:
  * local superoperator embedding correctness (A-local L -> residual ~ 0)
  * projection idempotence (P^2 = P)
  * the three C_L costs are finite
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


def _random_op(d: int, rng: np.random.Generator) -> np.ndarray:
    return rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))


def test_embed_a_acts_only_on_a():
    """A local superoperator must act on the A index and leave B untouched."""
    d_a, d_b, dim = 2, 4, 8
    rng = np.random.default_rng(3)
    op_a = _random_op(d_a**2, rng)
    emb = gf._super_embed_A(op_a, d_a, d_b)
    # product state rho = rho_A (x) rho_B
    rho_a = (rng.standard_normal((d_a, d_a)) + 1j * rng.standard_normal((d_a, d_a)))
    rho_b = (rng.standard_normal((d_b, d_b)) + 1j * rng.standard_normal((d_b, d_b)))
    rho = np.kron(rho_a, rho_b)
    out = (emb @ rho.reshape(-1, order="F")).reshape((dim, dim), order="F")
    # NOTE: op_a uses row-major flat index a = p*d_a + q, so reshape must be
    # C-order to match (F-order would transpose the A index mapping).
    out_a = (op_a @ rho_a.reshape(-1)).reshape((d_a, d_a))
    expected = np.kron(out_a, rho_b)
    assert np.allclose(out, expected, atol=1e-10), \
        "embedding does not act as op_a (x) I_B on product states"


def test_a_local_liouvillian_has_zero_boundary_residual():
    """If L is purely A-local (in S_F), then L_dF must vanish."""
    d_a, d_b, dim, n = 2, 4, 8, 3
    rng = np.random.default_rng(5)
    op_a = _random_op(d_a**2, rng)
    l_local = gf._super_embed_A(op_a, d_a, d_b)
    for cut in (1, 2):
        # for cut 2 (A=4, B=2) a pure A-local op with d_a=4 is different
        if cut == 2:
            d_a2, d_b2 = 4, 2
            op_a2 = _random_op(d_a2**2, rng)
            l_local = gf._super_embed_A(op_a2, d_a2, d_b2)
        resid = gf.boundary_residual(l_local, n, cut)
        assert np.linalg.norm(resid) < 1e-8, f"cut {cut}: residual {np.linalg.norm(resid)}"


def test_projection_idempotent():
    """Applying the local projection twice must equal applying it once."""
    n, cut = 3, 1
    psi, h_total = gf.build_system()
    L = gf.liouvillian(h_total, gf.GAMMA, n=3)
    bmat, pinv_gram = gf.local_projection_data(n, cut)
    l_flat = L.reshape(-1)
    proj1 = bmat @ (pinv_gram @ (bmat.conj().T @ l_flat))
    proj2 = bmat @ (pinv_gram @ (bmat.conj().T @ proj1))
    assert np.allclose(proj1, proj2, atol=1e-8), "projection not idempotent"


def test_c_L_costs_finite():
    """All three C_L costs must be finite for both cuts on the benchmark env."""
    psi, h_total = gf.build_system()
    L = gf.liouvillian(h_total, gf.GAMMA, n=3)
    rho0 = np.outer(psi, psi.conj())
    for cut in (1, 2):
        l_df = gf.boundary_residual(L, 3, cut)
        basis = gf.schmidt_basis(psi, 3, cut)
        g = gf.c_L_global(L, l_df)
        r = gf.c_L_rho(l_df, rho0)
        a = gf.c_L_align(rho0, l_df, basis)
        assert np.isfinite(g) and np.isfinite(r) and np.isfinite(a), \
            f"cut {cut}: costs not finite"
        assert g >= 0 and r >= 0
