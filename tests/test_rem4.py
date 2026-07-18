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


def test_hamiltonian_is_hermitian():
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    assert np.allclose(h, h.conj().T)
    assert len(bonds) == 2


def test_ground_state_is_normalized():
    h, _ = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    assert np.isclose(np.linalg.norm(psi), 1.0)


def test_mutual_information_is_nonnegative():
    h, _ = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    for cut in (1, 2):
        assert rem4.mutual_information_for_cut(psi, 3, cut) >= -1e-10


def test_phase_scan_writes_output(tmp_path):
    result = rem4.make_three_qubit_phase_scan(
        outdir=str(tmp_path),
        j23=0.6,
        h=0.2,
        lambda_max=1.0,
        grid_r=8,
        grid_l=8,
        alpha=0.85,
        beta=0.35,
        state_mode="ground",
        evolve_time=0.8,
    )
    assert Path(result["phase_png"]).exists()


# ── new tests ────────────────────────────────────────────────────────

def test_dynamical_cost_is_nonnegative():
    """C_H = <H^2> must be >= 0 for any state and any cut."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    metrics = rem4.evaluate_cuts(psi, 3, bonds, lambda_value=0.0)
    for m in metrics:
        assert m.c_h >= 0.0, f"cut {m.cut}: C_H={m.c_h} < 0"
        assert np.isfinite(m.c_h)


def test_cross_over_lambda_at_benchmark():
    """
    At J12=1.5, J23=0.6, h=0.2 the crossover should be non-negative
    and finite.  (Exact value depends on the convention; this test
    checks stability, not a precise published number.)
    """
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    lam_star = rem4.first_crossover_lambda(psi, 3, bonds)
    assert lam_star is not None, "λ* should exist for this parameter set"
    assert 0.1 < lam_star < 0.3, f"λ*={lam_star} outside expected range (0.1, 0.3)"


def test_info_and_std_diverge_at_lambda_zero():
    """At λ=0, info-selected and std-selected cuts should differ."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    metrics = rem4.evaluate_cuts(psi, 3, bonds, lambda_value=0.0)
    info_cut = rem4.info_selected_cut(metrics)
    std_cut = rem4.std_selected_cut(metrics)
    assert info_cut != std_cut, (
        f"info and std cuts should differ; both={info_cut}"
    )


def test_continuous_optimization_improves_phi():
    """Continuous manifold optimization should yield Φ >= best contiguous cut Φ."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)
    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=0.2)
    best_contiguous = max(m.phi for m in metrics)

    # Use a deterministic seed for reproducibility in tests
    rem4.N_RNG = np.random.default_rng(seed=12345)
    result = rem4.optimize_factorization(
        psi, h_total, bonds, n=n, lambda_value=0.2,
        n_params=4, steps=100, lr=0.01, verbose=False,
    )
    assert result["best_phi"] >= best_contiguous - 0.01, (
        f"cont Φ={best_contiguous:.4f} but opt Φ={result['best_phi']:.4f}"
    )
    assert result["best_u"] is not None


def test_finite_size_scan_runs():
    """finite_size_scan should run without error and return a dict."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        result = rem4.finite_size_scan(
            outdir=tmp,
            ns=[3, 4],
            h=0.2,
            j_scale_left=1.5,
            j_scale_right=0.6,
            alpha=0.85,
            beta=0.35,
            state_mode="ground",
            evolve_time=0.8,
            coupling_profile_name="weak_center",
        )
    assert isinstance(result, dict)
    for n in (3, 4):
        val = result.get(n)
        assert val is not None, f"missing N={n}"
        assert np.isfinite(val) or np.isnan(val)
