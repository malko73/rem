"""Phase A item 2: M2 vs Var(H_boundary) — unit tests.

These lock in the mathematical identities that the analysis script
(analysis/m2_vs_var.py) relies on:
  * evaluate_cuts now populates <H_dF> and Var(H_dF) = M2 - <H>^2
  * Var >= 0 (Cauchy-Schwarz), Delta = <H>^2 = M2 - Var
  * selection-agreement and classification helpers are deterministic
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

ANALYSIS_PATH = Path(__file__).parents[1] / "analysis" / "m2_vs_var.py"
aspec = importlib.util.spec_from_file_location("m2_vs_var", ANALYSIS_PATH)
m2v = importlib.util.module_from_spec(aspec)
assert aspec.loader is not None
sys.modules[aspec.name] = m2v
aspec.loader.exec_module(m2v)


def _benchmark():
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    return bonds, psi


def test_evaluate_cuts_populates_mean_and_var():
    """boundary_energy and var must no longer be NaN / default."""
    bonds, psi = _benchmark()
    metrics = rem4.evaluate_cuts(psi, 3, bonds, lambda_value=0.0)
    for m in metrics:
        assert np.isfinite(m.boundary_energy), f"cut {m.cut}: <H> not finite"
        assert np.isfinite(m.var), f"cut {m.cut}: Var not finite"


def test_var_identity():
    """Var(H) == M2 - <H>^2 for every cut, consistent with direct computation."""
    bonds, psi = _benchmark()
    metrics = rem4.evaluate_cuts(psi, 3, bonds, lambda_value=0.0)
    for m in metrics:
        bond = bonds[(m.cut - 1, m.cut)]
        m2_direct = rem4.boundary_cost_squared(psi, bond)
        mean_direct = rem4.boundary_energy(psi, bond)
        assert np.isclose(m.c_h, m2_direct)
        assert np.isclose(m.boundary_energy, mean_direct)
        assert np.isclose(m.var, m2_direct - mean_direct * mean_direct)


def test_variance_is_nonnegative():
    """Cauchy-Schwarz: <H^2> >= <H>^2  =>  Var >= 0 for all cuts, N=3..5."""
    for n in (3, 4, 5):
        _, bonds, psi = m2v.build_system(n)
        for r in m2v.per_cut_metrics(psi, n, bonds):
            assert r["var"] >= -1e-10, f"N={n} cut {r['cut']}: Var={r['var']}"


def test_delta_equals_mean_squared():
    """Delta(F) = M2 - Var = <H>^2."""
    _, bonds, psi = m2v.build_system(3)
    for r in m2v.per_cut_metrics(psi, 3, bonds):
        assert np.isclose(r["delta"], r["mean"] * r["mean"])
        assert np.isclose(r["m2"] - r["var"], r["delta"])


def test_benchmark_crossover_m2_matches_locked_value():
    """lambda* under M2 must match the locked benchmark 0.16512732363451144."""
    bonds, psi = _benchmark()
    rows = m2v.per_cut_metrics(psi, 3, bonds)
    lam = m2v.first_crossover(rows, "m2")
    assert lam is not None
    assert np.isclose(lam, 0.16512732363451144, rtol=1e-12)


def test_classify_returns_valid_verdict():
    """Classification always returns one of A/B/C with a reason."""
    bonds, psi = _benchmark()
    rows = m2v.per_cut_metrics(psi, 3, bonds)
    _, _, agreement = m2v.sweep_selection(rows, np.linspace(0.0, 3.0, 101))
    lam_m2 = m2v.first_crossover(rows, "m2")
    lam_var = m2v.first_crossover(rows, "var")
    v = m2v.classify(rows, agreement, lam_m2, lam_var)
    assert v["verdict"] in ("A", "B", "C")
    assert isinstance(v["reason"], str) and len(v["reason"]) > 0

