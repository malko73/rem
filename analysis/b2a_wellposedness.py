#!/usr/bin/env python3
"""
B2a — Continuous-TPS well-posedness audit of C_Γ^(0).

Purpose: test whether Γ_F^exact is well-defined as a global variational
functional on the full 45-dim quotient, or whether it diverges (especially
to -∞, which would make Φ→+∞ and break the optimization).

Tests (per master 2026-08-11):
  1. product limit: construct a continuous path to a product state and
     track I(F), C_F², Γ_F^exact, Φ(F).
  2. random TPS sweep: sample ~10k points from the 45-dim quotient and
     compute statistics (min C_F², min/max Γ_F, min/max Φ, Γ vs 1/C_F).
  3. Schmidt non-zero degeneracy: check gauge independence when p₁=p₂.
  4. local gauge invariance: I, Γ, Φ invariant under V_A⊗V_B.
  5. frame consistency: lab-frame vs pull-back computation agree.

Cost reduction (per master): for each environment, compute ρ₀ and Y=L(ρ₀)
once. For each TPS U, compute ρ_U = U†ρ₀U and Y_U = U†YU (8×8 ops only),
then evaluate Γ_U on the reference (A|B) cut.

Γ_F^exact is evaluated with the same Q_F = I - D_F construction as
gamma_F.gamma_exact (dephasing projection onto the Schmidt basis of the
rotated state), which was verified gauge-invariant in C1.5.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"

spec = importlib.util.spec_from_file_location("rem4_numerical", REPO / "src" / "rem4_numerical.py")
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)

spec2 = importlib.util.spec_from_file_location("quotient_geometry", REPO / "src" / "quotient_geometry.py")
qgeom = importlib.util.module_from_spec(spec2)
assert spec2.loader is not None
sys.modules[spec2.name] = qgeom
spec2.loader.exec_module(qgeom)

spec3 = importlib.util.spec_from_file_location("gamma_F", REPO / "analysis" / "gamma_F.py")
gamma_F = importlib.util.module_from_spec(spec3)
assert spec3.loader is not None
sys.modules[spec3.name] = gamma_F
spec3.loader.exec_module(gamma_F)

N = 3
N_A = 2
D_A, D_B = 2 ** N_A, 2 ** (N - N_A)
DIM = 2 ** N
COUPLINGS = [1.5, 0.6]
H_FIELD = 0.2
LAMBDA = 0.2


def liouvillian_dephasing(h_total: np.ndarray, gamma_vec: list[float]) -> np.ndarray:
    """Build vectorised Lindblad Liouvillian L for pure dephasing.

    L(rho) = -i[H, rho] + sum_i (gamma_i/2)(Z_i rho Z_i - rho)
    Returns a 64x64 matrix acting on vectorised rho (column-major).
    """
    dim = 2 ** N
    I = np.eye(dim, dtype=complex)
    # Same correction as gamma_F.liouvillian_env (B2a): vec_F convention
    # requires kron(I, h) - kron(h.conj(), I) for the commutator term.
    L = -1j * (np.kron(I, h_total) - np.kron(h_total.conj(), I))
    for site in range(N):
        g = gamma_vec[site]
        if g > 0:
            z = gamma_F.pauli_z(N, site)
            L += (g / 2.0) * (np.kron(z, z) - np.kron(I, I))
    return L


def gamma_exact_fast(rho_U: np.ndarray, Y_U: np.ndarray, n_a: int) -> tuple[float, float]:
    """
    Compute Γ_F^exact = -Re<Q_F ρ, Q_F Y> / |Q_F ρ|² on the reference cut (A|B).

    Uses the same Q_F = I - D_F construction as gamma_F.gamma_exact
    (dephasing projection onto the Schmidt basis of the rotated state),
    verified gauge-invariant in C1.5. Y_U = U† L(ρ₀) U is supplied directly.

    Returns (gamma, C_F²) where C_F² = |Q_F ρ|².
    """
    # Extract the pure state vector from rho_U (rho_U = |ψ><ψ|)
    evals_psi, evecs_psi = np.linalg.eigh(rho_U)
    psi_U = evecs_psi[:, np.argmax(evals_psi)]

    # Schmidt basis of the rotated state on the reference cut
    basis = gamma_F.schmidt_basis(psi_U, N, n_a)

    # Q_F = I - D_F applied to rho_U and Y_U
    q = lambda r: r - gamma_F.dephasing_projection(r, basis)  # noqa: E731
    x = q(rho_U)
    qy = q(Y_U)
    num = np.real(np.trace(x.conj().T @ qy))
    den = float(np.linalg.norm(x, ord="fro") ** 2)
    if den < 1e-300:
        return float("nan"), den
    return float(-num / den), den


def apply_tps(psi0: np.ndarray, h_total: np.ndarray, gamma_vec: list[float],
              U: np.ndarray):
    """Return (rho_U, Y_U) for a TPS U on the given environment."""
    rho0 = np.outer(psi0, psi0.conj())
    L = liouvillian_dephasing(h_total, gamma_vec)
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    rho_U = U.conj().T @ rho0 @ U
    Y_U = U.conj().T @ Y0 @ U
    return rho_U, Y_U


def evaluate_objective(rho_U: np.ndarray, Y_U: np.ndarray) -> dict:
    """Return {I, C_F_sq, gamma, Phi} for the rotated state."""
    evals_psi, evecs_psi = np.linalg.eigh(rho_U)
    psi_U = evecs_psi[:, np.argmax(evals_psi)]
    I_F = rem4.mutual_information_for_cut(psi_U, N, N_A)
    gamma, C_F_sq = gamma_exact_fast(rho_U, Y_U, N_A)
    Phi = I_F - LAMBDA * gamma if np.isfinite(gamma) else np.nan
    return dict(I=I_F, C_F_sq=C_F_sq, gamma=gamma, Phi=Phi)


def test_product_limit():
    """Test 1: construct a path toward a product state and track Γ."""
    print("\n=== Test 1: Product Limit ===")
    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    _, psi0 = rem4.ground_state(h_total)
    gamma_vec = [1.0, 1.0, 1.0]

    # Random horizontal direction (fixed seed for reproducibility)
    H_basis = qgeom.horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(20260811)
    theta_max = rng.standard_normal(H_basis.shape[2]) * 3.0

    n_points = 60
    t_values = np.linspace(0, 1, n_points)
    results = []
    for t in t_values:
        theta = t * theta_max
        U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))
        rho_U, Y_U = apply_tps(psi0, h_total, gamma_vec, U)
        results.append(dict(t=t, **evaluate_objective(rho_U, Y_U)))

    print("t       | I(F)    | C_F²     | Γ        | Φ")
    print("-" * 50)
    for r in results[::10]:
        g = r["gamma"] if np.isfinite(r["gamma"]) else float("nan")
        p = r["Phi"] if np.isfinite(r["Phi"]) else float("nan")
        print(f"{r['t']:.2f}   | {r['I']:.4f} | {r['C_F_sq']:.6f} | {g:.4f} | {p:.4f}")

    cf = np.array([r["C_F_sq"] for r in results if np.isfinite(r["gamma"])])
    ga = np.array([r["gamma"] for r in results if np.isfinite(r["gamma"])])
    ph = np.array([r["Phi"] for r in results if np.isfinite(r["Phi"])])
    print(f"\nmin C_F² = {cf.min():.6e}")
    print(f"Γ range: [{ga.min():.4f}, {ga.max():.4f}]")
    print(f"Φ range: [{ph.min():.4f}, {ph.max():.4f}]")

    if cf.min() < 1e-6:
        print("NOTE: C_F² → 0 approached")
    if ga.max() > 1e3:
        print("WARNING: Γ → +∞ (divergence detected)")
    if ga.min() < -1e3:
        print("CRITICAL: Γ → -∞ (Φ → +∞, optimization broken)")
    return results


def test_random_sweep(n_samples: int = 2000):
    """Test 2: random TPS sweep from the 45-dim quotient."""
    print(f"\n=== Test 2: Random TPS Sweep ({n_samples} samples) ===")
    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    _, psi0 = rem4.ground_state(h_total)
    gamma_vec = [1.0, 1.0, 1.0]

    H_basis = qgeom.horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(20260812)

    results = []
    for i in range(n_samples):
        theta = rng.standard_normal(H_basis.shape[2]) * 1.0
        U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))
        rho_U, Y_U = apply_tps(psi0, h_total, gamma_vec, U)
        results.append(evaluate_objective(rho_U, Y_U))

    cf = np.array([r["C_F_sq"] for r in results if np.isfinite(r["gamma"])])
    ga = np.array([r["gamma"] for r in results if np.isfinite(r["gamma"])])
    ph = np.array([r["Phi"] for r in results if np.isfinite(r["Phi"])])
    print(f"min C_F² = {cf.min():.6e}")
    print(f"Γ range: [{ga.min():.4f}, {ga.max():.4f}]")
    print(f"Φ range: [{ph.min():.4f}, {ph.max():.4f}]")

    mask = cf > 1e-6
    if np.sum(mask) > 10:
        corr = np.corrcoef(ga[mask], 1.0 / cf[mask])[0, 1]
        print(f"corr(Γ, 1/C_F²) = {corr:.4f}")

    n_neg = int(np.sum(ga < 0))
    n_pos = int(np.sum(ga > 0))
    print(f"Γ: {n_neg} negative, {n_pos} positive, "
          f"{len(ga)-n_neg-n_pos} zero/nan")
    return results


def test_schmidt_degeneracy():
    """Test 3: Schmidt non-zero degeneracy (p₁=p₂) gauge independence.

    Uses a state with exactly degenerate Schmidt values p₁=p₂ on cut A|B.
    D_F projects onto the full degenerate block, so Γ must be invariant
    under rotation of the Schmidt basis within that block.
    """
    print("\n=== Test 3: Schmidt Degeneracy ===")
    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    gamma_vec = [1.0, 1.0, 1.0]

    # Bell-like state on cut A|B (A=sites[0,1], B=site[2]):
    #   |psi> = (|000> + |011>)/sqrt(2)  -> Schmidt values 1/sqrt2, 1/sqrt2 (degenerate)
    psi = np.zeros(DIM, dtype=complex)
    psi[0b000] = 1.0 / np.sqrt(2.0)
    psi[0b011] = 1.0 / np.sqrt(2.0)
    rho0 = np.outer(psi, psi.conj())
    L = liouvillian_dephasing(h_total, gamma_vec)

    basis1 = gamma_F.schmidt_basis(psi, N, N_A)
    g1 = gamma_F.gamma_exact(rho0, L, basis1)
    c1 = float(np.linalg.norm(
        rho0 - gamma_F.dephasing_projection(rho0, basis1), ord="fro") ** 2)

    # Rotate the Schmidt basis within the degenerate A-block (indices 0,1)
    rng = np.random.default_rng(20260814)
    m = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
    q, _ = np.linalg.qr(m)  # random unitary 2x2
    basis2 = basis1.copy()
    for l in range(D_B):
        col0 = basis1[:, l].copy()
        col1 = basis1[:, D_B + l].copy()
        basis2[:, l] = q[0, 0] * col0 + q[1, 0] * col1
        basis2[:, D_B + l] = q[0, 1] * col0 + q[1, 1] * col1
    g2 = gamma_F.gamma_exact(rho0, L, basis2)
    c2 = float(np.linalg.norm(
        rho0 - gamma_F.dephasing_projection(rho0, basis2), ord="fro") ** 2)

    print(f"basis1: C_F²={c1:.6e} Γ={g1:.6f}")
    print(f"basis2 (degenerate-block rotated): C_F²={c2:.6e} Γ={g2:.6f}")
    print(f"ΔΓ = {abs(g1-g2):.2e}")
    if abs(g1 - g2) < 1e-8:
        print("PASS: degenerate Schmidt block rotation leaves Γ invariant")
    else:
        print("FAIL: Γ depends on degenerate Schmidt basis choice")
    return dict(base=dict(gamma=g1, C_F_sq=c1),
                rotated=dict(gamma=g2, C_F_sq=c2))


def test_gauge_invariance():
    """Test 4: local gauge invariance (V_A⊗V_B)."""
    print("\n=== Test 4: Local Gauge Invariance ===")
    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    _, psi0 = rem4.ground_state(h_total)
    gamma_vec = [1.0, 1.0, 1.0]

    H_basis = qgeom.horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(20260813)
    theta = rng.standard_normal(H_basis.shape[2]) * 0.5
    U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))
    rho_U, Y_U = apply_tps(psi0, h_total, gamma_vec, U)
    base = evaluate_objective(rho_U, Y_U)

    # Random local gauge V_A⊗V_B
    u_a = qgeom.u_basis(D_A)
    u_b = qgeom.u_basis(D_B)
    theta_a = rng.standard_normal(D_A * D_A) * 0.3
    theta_b = rng.standard_normal(D_B * D_B) * 0.3
    V_A = sla.expm(np.einsum("ijk,k->ij", u_a, theta_a))
    V_B = sla.expm(np.einsum("ijk,k->ij", u_b, theta_b))
    V = np.kron(V_A, V_B)

    rho_V = V.conj().T @ rho_U @ V
    Y_V = V.conj().T @ Y_U @ V
    gauged = evaluate_objective(rho_V, Y_V)

    print(f"Before gauge: I={base['I']:.6f}, Γ={base['gamma']:.6f}, "
          f"Φ={base['Phi']:.6f}")
    print(f"After gauge:  I={gauged['I']:.6f}, Γ={gauged['gamma']:.6f}, "
          f"Φ={gauged['Phi']:.6f}")
    dI = abs(base["I"] - gauged["I"])
    dG = abs(base["gamma"] - gauged["gamma"])
    dP = abs(base["Phi"] - gauged["Phi"])
    print(f"ΔI = {dI:.2e}, ΔΓ = {dG:.2e}, ΔΦ = {dP:.2e}")

    if dI < 1e-8 and dG < 1e-8 and dP < 1e-8:
        print("PASS: gauge invariance confirmed")
    else:
        print("FAIL: gauge invariance violated")
    return base, gauged


def test_frame_consistency():
    """Test 5: lab-frame vs pull-back computation agree.

    The pull-back (master's prescription) rotates rho0 and Y0 = L(rho0)
    into the reference frame. The lab-frame equivalent must rotate the
    ENVIRONMENT too (dephasing operators Z_i -> U† Z_i U), otherwise the
    two computations describe different environments and disagree.
    """
    print("\n=== Test 5: Frame Consistency ===")
    h_total, bonds = rem4.xy_chain_hamiltonian(N, COUPLINGS, H_FIELD)
    _, psi0 = rem4.ground_state(h_total)
    gamma_vec = [1.0, 1.0, 1.0]

    H_basis = qgeom.horizontal_basis(D_A, D_B)
    rng = np.random.default_rng(20260815)
    theta = rng.standard_normal(H_basis.shape[2]) * 0.5
    U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))

    # Pull-back: rho_U = U† rho0 U, Y_U = U† Y0 U (Y0 = L(rho0) in lab frame)
    rho_U, Y_U = apply_tps(psi0, h_total, gamma_vec, U)
    pulled = evaluate_objective(rho_U, Y_U)

    # Lab frame with ROTATED environment: psi_rot = U† psi0, h_rot = U† h U,
    # Z_i^rot = U† Z_i U. Then rho_rot = rho_U and L_rot(rho_rot) must equal Y_U.
    dim = 2 ** N
    I = np.eye(dim, dtype=complex)
    psi_rot = U.conj().T @ psi0
    rho_rot = np.outer(psi_rot, psi_rot.conj())
    h_rot = U.conj().T @ h_total @ U
    L_rot = -1j * (np.kron(I, h_rot) - np.kron(h_rot.conj(), I))
    for site in range(N):
        z = U.conj().T @ gamma_F.pauli_z(N, site) @ U
        # NOTE: z^rot = U† Z U is Hermitian but generally NOT symmetric;
        # vec_F(z ρ z) = (z^T ⊗ z) vec_F(ρ), so use kron(z.T, z).
        L_rot += (gamma_vec[site] / 2.0) * (np.kron(z.T, z) - np.kron(I, I))
    Y_rot = (L_rot @ rho_rot.reshape(-1, order="F")).reshape(rho_rot.shape, order="F")
    lab = evaluate_objective(rho_rot, Y_rot)

    # Also verify the environment rotation directly: Y_rot should equal Y_U
    dY = float(np.linalg.norm(Y_rot - Y_U, ord="fro"))

    print(f"pull-back: I={pulled['I']:.6f}, Γ={pulled['gamma']:.6f}")
    print(f"lab-frame: I={lab['I']:.6f}, Γ={lab['gamma']:.6f}")
    print(f"ΔI = {abs(pulled['I']-lab['I']):.2e}, "
          f"ΔΓ = {abs(pulled['gamma']-lab['gamma']):.2e}, "
          f"|Y_rot - Y_U| = {dY:.2e}")
    if abs(pulled["I"] - lab["I"]) < 1e-8 and abs(pulled["gamma"] - lab["gamma"]) < 1e-8:
        print("PASS: frame consistency confirmed (with rotated environment)")
    else:
        print("FAIL: frame consistency violated")
    return pulled, lab


def main():
    OUTDIR.mkdir(exist_ok=True)

    results1 = test_product_limit()
    results2 = test_random_sweep(n_samples=2000)
    res3 = test_schmidt_degeneracy()
    res4 = test_gauge_invariance()
    res5 = test_frame_consistency()

    output = {
        "product_limit": results1,
        "random_sweep": results2,
        "schmidt_degeneracy": {
            "base": res3["base"], "rotated": res3["rotated"]},
        "gauge_invariance": {
            "base": res4[0], "gauged": res4[1]},
        "frame_consistency": {
            "pullback": res5[0], "lab": res5[1]},
        "params": dict(N=N, N_A=N_A, lambda_value=LAMBDA,
                       couplings=COUPLINGS, h_field=H_FIELD),
    }
    out_path = OUTDIR / "b2a_wellposedness.json"
    out_path.write_text(json.dumps(output, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
