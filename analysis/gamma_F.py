#!/usr/bin/env python3
"""
Phase C0 — M2 / Var vs Gamma_F falsification test (2026-08-11).

Question: does the canonical dynamical cost C_H = M2 = <H_dF^2> (or its
variance alternative Var(H_dF)) predict which contiguous factorization is
more stable under a FIXED, IDENTICAL open-system environment?

Setup (identical for cut 1 and cut 2):
    H      = 1.5(X1X2+Y1Y2) + 0.6(X2X3+Y2Y3) + 0.2 Σ_i Z_i   (3-qubit XY)
    rho_0  = ground state
    L_i    = sqrt(gamma/2) Z_i,  i = 1,2,3   (local pure dephasing, all qubits)
    gamma  = 1.0 (model input; Gamma_F scales linearly with gamma)

Master-equation (Lindblad):
    d rho/dt = -i[H, rho] + (gamma/2) Σ_i (Z_i rho Z_i - rho)
Vectorised Liouvillian (64x64) built with NumPy/SciPy only (no QuTiP).

Factorisation-specific coherence:
    D_F = dephasing projection onto the initial Schmidt basis of cut F
    C_F(t) = || rho(t) - D_F[rho(t)] ||_2   (Frobenius)

Measured quantities (pre-registered):
    Gamma_F^(0) = - d/dt log C_F(t) |_{t=0}     (initial structural decoherence rate)
    Gamma_F^fit = slope of log C_F(t) over short-time window [0, t_fit]
    check Gamma_F^(0) ≈ Gamma_F^fit

Verdict table and 3-case classification:
    case 1: Gamma_2 < Gamma_1  -> supports M2 (cut2 lower M2, lower Gamma)
    case 2: Gamma_1 ≈ Gamma_2  -> consistent with Var degeneracy (M2's
             lambda* ≈ 0.165 may be a second-moment proxy artifact)
    case 3: Gamma_1 < Gamma_2  -> contradicts both M2 and Var; rethink C_H

Outputs (analysis_output/):
    gamma_F_table.txt, gamma_F_results.json, gamma_F_decay.png
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"
MODULE_PATH = REPO / "src" / "rem4_numerical.py"

spec = importlib.util.spec_from_file_location("rem4_numerical", MODULE_PATH)
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)

# ── C0 fixed conditions ──────────────────────────────────────────────
GAMMA = 1.0                       # dephasing rate model input
T_MAX = 5.0                       # integration horizon
T_FIT = 0.5                       # short-time fit window for Gamma^fit
T_INIT = 0.05                     # ultra-short window for Gamma^(0)
N_T = 400                         # time grid points


def build_system():
    """Ground state VECTOR and Hamiltonian of the 3-qubit XY benchmark."""
    h_total, _ = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
    _, psi = rem4.ground_state(h_total)
    return psi, h_total


def pauli_z(n: int, site: int) -> np.ndarray:
    """Z on `site` (0-indexed) in the n-qubit Hilbert space."""
    z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    ops = [np.eye(2, dtype=complex)] * n
    ops[site] = z
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def sigma_minus(n: int, site: int) -> np.ndarray:
    """sigma_- = |0><1| on `site` (0-indexed) in the n-qubit Hilbert space."""
    sm = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=complex)
    ops = [np.eye(2, dtype=complex)] * n
    ops[site] = sm
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def liouvillian_env(h_total: np.ndarray, gamma_vec, kappa_vec, n: int = 3) -> np.ndarray:
    """Vectorised Lindblad Liouvillian with per-site dephasing + amplitude damping.

    L_i^deph = sqrt(gamma_i/2) Z_i   ->  (gamma_i/2)(Z rho Z - rho)
    L_i^damp = sqrt(kappa_i) sigma_-  ->  kappa (s_- rho s_+ - 1/2 {s_+ s_-, rho})
    """
    dim = 2**n
    I = np.eye(dim, dtype=complex)
    L = -1j * (np.kron(h_total, I) - np.kron(I, h_total.conj()))
    for site in range(n):
        g = gamma_vec[site]
        k = kappa_vec[site]
        if g > 0:
            z = pauli_z(n, site)
            L += (g / 2.0) * (np.kron(z, z) - np.kron(I, I))
        if k > 0:
            sm = sigma_minus(n, site)
            spsm = sm.conj().T @ sm                 # sigma_+ sigma_-
            L += k * (np.kron(sm.conj(), sm)
                      - 0.5 * np.kron(I, spsm)
                      - 0.5 * np.kron(spsm.conj(), I))
    return L


def liouvillian(h_total: np.ndarray, gamma: float, n: int = 3) -> np.ndarray:
    """Uniform pure-dephasing Liouvillian (backward-compatible wrapper)."""
    return liouvillian_env(h_total, [gamma] * n, [0.0] * n, n=n)


def solve_dynamics(L: np.ndarray, rho0: np.ndarray, t_grid: np.ndarray) -> np.ndarray:
    """Integrate the vectorised master equation; return rho(t) per grid point."""
    rho0_vec = rho0.reshape(-1, order="F")

    def rhs(_t: float, rho_vec: np.ndarray) -> np.ndarray:
        return L @ rho_vec

    sol = solve_ivp(rhs, (0.0, t_grid[-1]), rho0_vec,
                    t_eval=t_grid, method="RK45", rtol=1e-9, atol=1e-11)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    return sol.y.T  # (n_t, 64)


def schmidt_basis(psi: np.ndarray, n: int, cut: int) -> np.ndarray:
    """Full product basis |a_k>⊗|b_l> from the Schmidt decomposition at cut.

    Uses the full SVD so that U (d_a x d_a) and Vh (d_b x d_b) are COMPLETE
    orthonormal bases of both subsystems (the Schmidt vectors plus a completion).
    """
    d_a, d_b = 2**cut, 2 ** (n - cut)
    mat = psi.reshape((d_a, d_b))
    u, _, vh = np.linalg.svd(mat, full_matrices=True)
    basis = np.zeros((d_a * d_b, d_a * d_b), dtype=complex)
    for k in range(d_a):
        for l in range(d_b):
            basis[:, k * d_b + l] = np.kron(u[:, k], vh[l, :])
    return basis


def dephasing_projection(rho: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """D_F[rho] = basis diag(basis^† rho basis) basis^† (off-diagonals zeroed)."""
    b = basis.conj().T @ rho @ basis
    return basis @ (np.diag(np.diag(b.real))) @ basis.conj().T


def coherence_norm(rho: np.ndarray, basis: np.ndarray) -> float:
    """C_F(rho) = || rho - D_F[rho] ||_2 (Frobenius)."""
    return float(np.linalg.norm(rho - dephasing_projection(rho, basis), ord="fro"))


def gamma_exact(rho0: np.ndarray, L: np.ndarray, basis: np.ndarray) -> float:
    """Exact t=0 structural decoherence rate.

    With Q_F = I - D_F and X_F = Q_F rho0, C_F = |X_F|_2, the exact
    initial logarithmic derivative is

        Gamma_F^exact(0) = - Re <X_F, Q_F L(rho0)>_HS / |X_F|_2^2
    """
    q = lambda r: r - dephasing_projection(r, basis)          # noqa: E731
    x = q(rho0)
    lrho = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    qlrho = q(lrho)
    num = np.real(np.trace(x.conj().T @ qlrho))
    den = float(np.linalg.norm(x, ord="fro") ** 2)
    if den < 1e-300:
        return float("nan")
    return float(-num / den)


def log_linear_slope(t: np.ndarray, c: np.ndarray) -> float:
    """Slope of log(c) vs t by least squares (i.e. -Gamma if c ~ c0 e^{-G t})."""
    mask = c > 1e-12
    if np.sum(mask) < 2:
        return float("nan")
    x = t[mask]
    y = np.log(c[mask])
    slope, _ = np.polyfit(x, y, 1)
    return float(slope)


def run_cut(psi: np.ndarray, h_total: np.ndarray, gamma: float,
            t_grid: np.ndarray, L: np.ndarray, cut: int) -> dict:
    """Measure M2, Var, Gamma^(0), Gamma^fit for one contiguous cut."""
    n = 3
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    i2 = np.eye(2, dtype=complex)
    if cut == 1:
        bond = 1.5 * (np.kron(np.kron(x, x), i2) + np.kron(np.kron(y, y), i2))
    else:
        bond = 0.6 * (np.kron(i2, np.kron(x, x)) + np.kron(i2, np.kron(y, y)))

    m2 = rem4.boundary_cost_squared(psi, bond)
    mean = rem4.boundary_energy(psi, bond)
    var = m2 - mean * mean

    rho0 = np.outer(psi, psi.conj())
    basis = schmidt_basis(psi, n, cut)
    c0 = coherence_norm(rho0, basis)
    g_exact = gamma_exact(rho0, L, basis)

    rho_t = solve_dynamics(L, rho0, t_grid)
    c_t = np.array([coherence_norm(rho_t[i].reshape((8, 8), order="F"), basis)
                    for i in range(len(t_grid))])

    # Gamma^(0): ultra-short window
    mask0 = t_grid <= T_INIT
    gamma0 = -log_linear_slope(t_grid[mask0], c_t[mask0])
    # Gamma^fit: short-time window
    maskf = t_grid <= T_FIT
    gamma_fit = -log_linear_slope(t_grid[maskf], c_t[maskf])

    # auxiliary: <H_dF>_t decay
    hd_t = np.array([np.real(np.trace(rho_t[i].reshape((8, 8), order="F") @ bond))
                     for i in range(len(t_grid))])

    return dict(
        cut=cut, m2=float(m2), var=float(var), mean=float(mean),
        c0=float(c0), gamma_exact=float(g_exact),
        gamma0=float(gamma0), gamma_fit=float(gamma_fit),
        c_t=c_t.tolist(), hd_t=hd_t.tolist(),
    )


def classify_case(d1: dict, d2: dict, rel_tol: float = 0.05) -> dict:
    # primary indicator: exact t=0 derivative (fall back to window estimate)
    g1e, g2e = float(d1["gamma_exact"]), float(d2["gamma_exact"])
    g1 = g1e if np.isfinite(g1e) else float(d1["gamma0"])
    g2 = g2e if np.isfinite(g2e) else float(d2["gamma0"])
    if g2 < g1 * (1 - rel_tol):
        verdict = 1
        text = ("Gamma_2 < Gamma_1: the cut with lower M2 is more stable under "
                "identical dephasing -> supports the canonical M2 proxy.")
    elif g1 < g2 * (1 - rel_tol):
        verdict = 3
        text = ("Gamma_1 < Gamma_2: contradicts both M2 and Var ordering -> "
                "rethink the definition of C_H (most important falsification).")
    else:
        verdict = 2
        text = ("Gamma_1 ~ Gamma_2: consistent with Var degeneracy; the M2 "
                "crossover lambda* ~ 0.165 may be a second-moment proxy artifact.")
    return dict(verdict=verdict, text=text,
                gamma1=float(g1), gamma2=float(g2),
                gamma1_fit=float(d1["gamma_fit"]), gamma2_fit=float(d2["gamma_fit"]))


def write_table(d1: dict, d2: dict, verdict: dict, outpath: Path) -> None:
    rows = [
        "# Phase C0 — M2 / Var vs Gamma_F (identical fixed environment, pure dephasing)",
        f"# H = XY(1.5,0.6,h=0.2), rho0 = ground state, L_i = sqrt(gamma/2) Z_i, gamma={GAMMA}",
        "# Gamma^exact: analytic t=0 derivative (primary). Gamma^(0): log-linear fit over t<= "
        f"{T_INIT}; Gamma^fit over t<= {T_FIT}",
        "",
        "| cut | M2=<H^2> | Var(H) | <H> | C_F(0) | Gamma^exact | Gamma^(0) | Gamma^fit |",
        "|-----|----------|--------|-----|--------|-------------|-----------|-----------|",
    ]
    for d in (d1, d2):
        rows.append(
            f"| {d['cut']} | {d['m2']:.6f} | {d['var']:.6f} | {d['mean']:+.6f} | "
            f"{d['c0']:.6f} | {d['gamma_exact']:.6f} | {d['gamma0']:.6f} | {d['gamma_fit']:.6f} |"
        )
    rows += [
        "",
        f"case {verdict['verdict']}: {verdict['text']}",
        f"Gamma1={verdict['gamma1']:.6f}  Gamma2={verdict['gamma2']:.6f}",
    ]
    outpath.write_text("\n".join(rows), encoding="utf-8")


def plot_decay(d1: dict, d2: dict, t_grid: np.ndarray, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for d, ls, mk in ((d1, "-", "o"), (d2, "--", "s")):
        ax.semilogy(t_grid, d["c_t"], ls, markersize=2, label=f"cut {d['cut']} "
                     rf"($\Gamma^{{0}}={d['gamma0']:.3f}$, fit={d['gamma_fit']:.3f})")
        # fit line
        mask = t_grid <= T_FIT
        t = t_grid[mask]
        c0, g = d["c0"], d["gamma_fit"]
        ax.plot(t, c0 * np.exp(-g * t), ":", color="gray", alpha=0.7)
    ax.set_xlabel("time $t$")
    ax.set_ylabel(r"$C_F(t)=\|\rho(t)-\mathcal{D}_F[\rho(t)]\|_2$")
    ax.set_title("Factorisation-specific coherence decay (identical dephasing env.)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    psi, h_total = build_system()
    L = liouvillian(h_total, GAMMA, n=3)
    t_grid = np.linspace(0.0, T_MAX, N_T)

    d1 = run_cut(psi, h_total, GAMMA, t_grid, L, cut=1)
    d2 = run_cut(psi, h_total, GAMMA, t_grid, L, cut=2)
    verdict = classify_case(d1, d2)

    write_table(d1, d2, verdict, OUTDIR / "gamma_F_table.txt")
    (OUTDIR / "gamma_F_results.json").write_text(
        json.dumps(dict(cut1=d1, cut2=d2, verdict=verdict,
                        params=dict(gamma=GAMMA, t_max=T_MAX, t_fit=T_FIT,
                                    t_init=T_INIT)),
                   indent=2, ensure_ascii=False), encoding="utf-8")
    plot_decay(d1, d2, t_grid, OUTDIR / "gamma_F_decay.png")

    print(f"wrote {OUTDIR / 'gamma_F_table.txt'}")
    print(f"wrote {OUTDIR / 'gamma_F_results.json'}")
    print(f"wrote {OUTDIR / 'gamma_F_decay.png'}")
    print("\n" + (OUTDIR / "gamma_F_table.txt").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
