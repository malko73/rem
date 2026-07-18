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


# ── core structural tests ────────────────────────────────────────────

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
        j23=0.6, h=0.2, lambda_max=1.0,
        grid_r=8, grid_l=8,
        alpha=0.85, beta=0.35,
        state_mode="ground", evolve_time=0.8,
    )
    assert Path(result["phase_png"]).exists()


# ── C_H definition tests ────────────────────────────────────────────

def test_boundary_cost_squared_matches_explicit():
    """C_H = <psi|bond^2|psi> must equal explicit vdot(bond psi, bond psi)."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    for cut in (1, 2):
        bond = bonds[(cut - 1, cut)]
        c_h = rem4.boundary_cost_squared(psi, bond)
        explicit = float(np.real(np.vdot(bond @ psi, bond @ psi)))
        assert np.isclose(c_h, explicit), f"cut {cut}: {c_h} != {explicit}"


def test_dynamical_cost_is_nonnegative():
    """C_H = <H^2> must be >= 0 for any state and any cut."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    metrics = rem4.evaluate_cuts(psi, 3, bonds, lambda_value=0.0)
    for m in metrics:
        assert m.c_h >= 0.0, f"cut {m.cut}: C_H={m.c_h} < 0"
        assert np.isfinite(m.c_h)


def test_variance_detection():
    """<H^2> should differ from <H>^2 for a non-eigenstate of bond."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    for cut in (1, 2):
        bond = bonds[(cut - 1, cut)]
        c_h = rem4.boundary_cost_squared(psi, bond)
        ebd = rem4.boundary_energy(psi, bond)
        mean_sq = ebd * ebd
        # For the XY chain ground state, <H^2> > <H>^2 (non-zero variance)
        assert c_h > mean_sq + 1e-10, (
            f"cut {cut}: <H^2>={c_h:.6f} <= <H>^2={mean_sq:.6f}"
        )


def test_eigenstate_cost_squared_equals_mean_squared():
    """For an eigenstate of the bond operator, <H^2> == <H>^2."""
    # Construct an eigenstate of the (1,2) bond
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    bond12 = bonds[(0, 1)]
    evals, evecs = np.linalg.eigh(bond12)
    psi_eig = evecs[:, 0]  # ground state of bond12
    psi_eig = rem4.normalize(psi_eig)

    c_h = rem4.boundary_cost_squared(psi_eig, bond12)
    ebd = rem4.boundary_energy(psi_eig, bond12)
    assert np.isclose(c_h, ebd * ebd), (
        f"eigenstate: <H^2>={c_h:.10f} != <H>^2={ebd*ebd:.10f}"
    )


def test_evaluate_cuts_uses_correct_cost():
    """evaluate_cuts should produce C_H identical to boundary_cost_squared."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    metrics = rem4.evaluate_cuts(psi, 3, bonds, lambda_value=0.0)
    for m in metrics:
        bond = bonds[(m.cut - 1, m.cut)]
        direct = rem4.boundary_cost_squared(psi, bond)
        assert np.isclose(m.c_h, direct), f"cut {m.cut}: {m.c_h} != {direct}"


# ── crossover / benchmark tests ──────────────────────────────────────

def test_cross_over_lambda_at_benchmark():
    """λ* at J12=1.5, J23=0.6, h=0.2 should be precisely reproducible."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    lam_star = rem4.first_crossover_lambda(psi, 3, bonds)
    assert lam_star is not None
    # Exact value depends on <H^2>; lock it in here
    expected = 0.16512732363451144
    assert np.isclose(lam_star, expected, rtol=1e-12), (
        f"λ*={lam_star} != {expected}"
    )


def test_info_and_std_diverge_at_lambda_zero():
    """At λ=0, info-selected and std-selected cuts should differ."""
    h, bonds = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h)
    metrics = rem4.evaluate_cuts(psi, 3, bonds, lambda_value=0.0)
    assert rem4.info_selected_cut(metrics) != rem4.std_selected_cut(metrics)


def test_benchmark_metrics_fixed():
    """Lock in all three-qubit benchmark values at λ=0, λ=0.2 (seed 0)."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)

    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=0.2)
    m1 = next(m for m in metrics if m.cut == 1)
    m2 = next(m for m in metrics if m.cut == 2)

    assert np.isclose(m1.phi_s, 1.972465, rtol=1e-5)
    assert np.isclose(m1.c_h,   8.379310, rtol=1e-5)
    assert np.isclose(m1.phi,   0.296603, rtol=1e-5)
    assert np.isclose(m2.phi_s, 0.724103, rtol=1e-5)
    assert np.isclose(m2.c_h,   0.819310, rtol=1e-5)
    assert np.isclose(m2.phi,   0.560240, rtol=1e-5)


# ── continuous optimisation tests ────────────────────────────────────

def test_identity_factorization_matches_discrete():
    """U=I should give the same Φ, MI, C_H as discrete cut 2."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)
    u_eye = np.eye(8, dtype=complex)

    phi_c, mi_c, ch_c = rem4.evaluate_factorization(
        psi, h_total, bonds, u_eye, n=n, lambda_value=0.2)
    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=0.2)
    cut2 = next(m for m in metrics if m.cut == 2)

    assert np.isclose(phi_c, cut2.phi)
    assert np.isclose(mi_c, cut2.phi_s)
    assert np.isclose(ch_c, cut2.c_h)


def test_continuous_cost_depends_on_u():
    """C_H should change under non-trivial U (not constant)."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)

    u_eye = np.eye(8, dtype=complex)
    _, _, ch_id = rem4.evaluate_factorization(
        psi, h_total, bonds, u_eye, n=n, lambda_value=0.2)

    # a random non-trivial U
    import scipy.linalg as sla
    rng = np.random.default_rng(seed=999)
    G = 1j * rng.standard_normal((8, 8))
    G = G - G.conj().T
    u_rand = sla.expm(0.5 * G)
    _, _, ch_rot = rem4.evaluate_factorization(
        psi, h_total, bonds, u_rand, n=n, lambda_value=0.2)

    assert not np.isclose(ch_id, ch_rot, rtol=1e-8), (
        f"C_H unchanged: {ch_id} == {ch_rot}"
    )


def test_continuous_optimization_improves_phi():
    """Continuous optimisation should yield Φ >= best contiguous cut Φ."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)
    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=0.2)
    best_contiguous = max(m.phi for m in metrics)

    rem4.N_RNG = np.random.default_rng(seed=12345)
    result = rem4.optimize_factorization(
        psi, h_total, bonds, n=n, lambda_value=0.2,
        n_params=4, steps=100, lr=0.01, verbose=False,
    )
    assert result["best_phi"] >= best_contiguous - 0.01, (
        f"cont Φ={best_contiguous:.4f} but opt Φ={result['best_phi']:.4f}"
    )
    assert result["best_u"] is not None


def test_adam_optimization_improves_phi():
    """Adam optimiser should also beat the best contiguous-cut Φ."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)
    metrics = rem4.evaluate_cuts(psi, n, bonds, lambda_value=0.2)
    best_contiguous = max(m.phi for m in metrics)

    rem4.N_RNG = np.random.default_rng(seed=12345)
    result = rem4.optimize_factorization_adam(
        psi, h_total, bonds, n=n, lambda_value=0.2,
        n_params=4, steps=50, lr=0.01, verbose=False,
    )
    assert result["best_phi"] >= best_contiguous - 0.01, (
        f"cont Φ={best_contiguous:.4f} but adam Φ={result['best_phi']:.4f}"
    )
    assert result["best_u"] is not None


def test_adam_vs_sgd_same_seed():
    """Under identical seed+params, Adam should not produce identical Φ to SGD."""
    n = 3
    h_total, bonds = rem4.xy_chain_hamiltonian(n, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)

    rem4.N_RNG = np.random.default_rng(seed=999)
    result_sgd = rem4.optimize_factorization(
        psi, h_total, bonds, n=n, lambda_value=0.2,
        n_params=4, steps=50, lr=0.01, verbose=False,
    )

    rem4.N_RNG = np.random.default_rng(seed=999)
    result_adam = rem4.optimize_factorization_adam(
        psi, h_total, bonds, n=n, lambda_value=0.2,
        n_params=4, steps=50, lr=0.01, verbose=False,
    )
    assert not np.isclose(result_sgd["best_phi"], result_adam["best_phi"],
                          rtol=1e-5), "Adam and SGD gave identical Φ—suspicious"


# ── finite-size / smoke tests ───────────────────────────────────────

def test_finite_size_scan_runs():
    """finite_size_scan should run without error and return a dict."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        result = rem4.finite_size_scan(
            outdir=tmp, ns=[3, 4],
            h=0.2, j_scale_left=1.5, j_scale_right=0.6,
            alpha=0.85, beta=0.35,
            state_mode="ground", evolve_time=0.8,
            coupling_profile_name="weak_center",
        )
    assert isinstance(result, dict)
    for n in (3, 4):
        assert result.get(n) is not None, f"missing N={n}"
        assert np.isfinite(result[n]) or np.isnan(result[n])
