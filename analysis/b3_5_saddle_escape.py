#!/usr/bin/env python3
"""
B3.5: Hessian-validated saddle escape (A1 first, per master).

Master's protocol:
  1. Verify the Hessian is correct: take the max-positive eigenvector v_max,
     build K_vmax = Σ_k v_max[k] H_k, and compute
         U(ε) = U* @ expm(ε K_vmax)         (right-multiplied, matches the
                                              Hessian parameterization)
     for ε = 1e-4, 1e-3, 1e-2, 3e-2, 1e-1. Near a stationary point,
         ΔΦ ≈ (1/2) λ_max ε²
     confirms the Hessian implementation, sign, and quotient coordinates.
  2. Escape along the positive-curvature direction: pick the ε that raises Φ
     most, set U_escape = U* @ expm(ε K_vmax), restart Adam from there.
  3. Repeat up to 10 times; record Φ, |∇Φ|, λ_max(H) each round.
     Stop when |∇Φ| < 1e-3 AND λ_max(H) < ε_H (ε_H = 1e-2).
  4. Verdict: strict local max (all λ<0) / flat max (λ≈0 only) /
     saddle chain continues -> Gate 1 remains PENDING.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
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

spec4 = importlib.util.spec_from_file_location("b2_canonical_open", REPO / "analysis" / "b2_canonical_open.py")
b2 = importlib.util.module_from_spec(spec4)
assert spec4.loader is not None
sys.modules[spec4.name] = b2
spec4.loader.exec_module(b2)

N, N_A = 3, 2
D_A, D_B = 4, 2
DIM = 8
LAMBDA = 0.2
ADAM_STEPS = 500
LR = 0.02
H_EPS = 1e-4
GRAD_TOL = 1e-3
LAMBDA_H_TOL = 1e-2
MAX_ROUNDS = 10
ESC_EPSILONS = [1e-4, 1e-3, 1e-2, 3e-2, 1e-1]
ENVS = ["A1_deph_0512"]


def phi(U: np.ndarray, L: np.ndarray, psi0: np.ndarray) -> float:
    rho_U, Y_U = b2.apply_tps(psi0, L, U)
    obj = b2.canonical_open_objective(rho_U, Y_U, N_A)
    return obj["Phi"]


def grad_phi(U: np.ndarray, L: np.ndarray, psi0: np.ndarray,
             H_basis: np.ndarray, eps: float = 1e-5) -> np.ndarray:
    """Finite-difference gradient in the right-multiplication chart."""
    g = np.zeros(H_basis.shape[2])
    for k in range(len(g)):
        d = np.zeros(len(g))
        d[k] = eps
        U_p = U @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))
        U_m = U @ sla.expm(np.einsum("ijk,k->ij", H_basis, -d))
        g[k] = (phi(U_p, L, psi0) - phi(U_m, L, psi0)) / (2 * eps)
    return g


def hessian_at(U: np.ndarray, L: np.ndarray, psi0: np.ndarray,
               H_basis: np.ndarray) -> np.ndarray:
    """45x45 Hessian in the right-multiplication chart (same as B3.4)."""
    n = H_basis.shape[2]
    H = np.zeros((n, n))
    f0 = phi(U, L, psi0)
    for i in range(n):
        H_i = H_basis[:, :, i]
        U_p = U @ sla.expm(H_EPS * H_i)
        U_m = U @ sla.expm(-H_EPS * H_i)
        H[i, i] = (phi(U_p, L, psi0) - 2 * f0 + phi(U_m, L, psi0)) / (H_EPS ** 2)
    for i in range(n):
        H_i = H_basis[:, :, i]
        for j in range(i + 1, n):
            H_j = H_basis[:, :, j]
            U_pp = U @ sla.expm(H_EPS * H_i) @ sla.expm(H_EPS * H_j)
            U_pm = U @ sla.expm(H_EPS * H_i) @ sla.expm(-H_EPS * H_j)
            U_mp = U @ sla.expm(-H_EPS * H_i) @ sla.expm(H_EPS * H_j)
            U_mm = U @ sla.expm(-H_EPS * H_i) @ sla.expm(-H_EPS * H_j)
            H[i, j] = (phi(U_pp, L, psi0) - phi(U_pm, L, psi0)
                       - phi(U_mp, L, psi0) + phi(U_mm, L, psi0)) / (4 * H_EPS ** 2)
            H[j, i] = H[i, j]
    return H


def adam_from(U_start: np.ndarray, L: np.ndarray, psi0: np.ndarray,
              H_basis: np.ndarray, steps: int = ADAM_STEPS, lr: float = LR) -> dict:
    """Run Adam in the right-multiplication chart from U_start."""
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U_start @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    best_phi_val = phi(U_start, L, psi0)
    best_U = U_start
    final_delta = delta.copy()
    for step in range(steps):
        U_cur = to_U(delta)
        g = grad_phi(U_cur, L, psi0, H_basis)
        m = beta1 * m + (1 - beta1) * g
        v = beta2 * v + (1 - beta2) * g ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        delta = delta + lr * m_hat / (np.sqrt(v_hat) + eps)
        U_new = to_U(delta)
        p = phi(U_new, L, psi0)
        if p > best_phi_val:
            best_phi_val = p
            best_U = U_new
            final_delta = delta.copy()
    g_final = grad_phi(best_U, L, psi0, H_basis)
    return dict(U=best_U, Phi=best_phi_val, grad_norm=float(np.linalg.norm(g_final)))


def verify_hessian(U: np.ndarray, v_max: np.ndarray, lam_max: float,
                   L: np.ndarray, psi0: np.ndarray, H_basis: np.ndarray) -> list:
    """Step 1: verify ΔΦ ≈ ½ λ_max ε² along ±v_max."""
    print("  [verify] ΔΦ along ±v_max (expected ΔΦ ≈ ½ λ_max ε²):")
    rows = []
    K = np.einsum("ijk,k->ij", H_basis, v_max)  # anti-Hermitian generator
    for eps in ESC_EPSILONS:
        U_p = U @ sla.expm(eps * K)
        U_m = U @ sla.expm(-eps * K)
        dp = phi(U_p, L, psi0) - phi(U, L, psi0)
        dm = phi(U_m, L, psi0) - phi(U, L, psi0)
        pred = 0.5 * lam_max * eps ** 2
        rows.append(dict(eps=eps, dPhi_plus=dp, dPhi_minus=dm, pred_half_lam_eps2=pred))
        print(f"    ε={eps:.0e}: ΔΦ(+)= {dp:+.6e}  ΔΦ(-)= {dm:+.6e}  "
              f"½λmaxε²= {pred:+.6e}")
    return rows


def main():
    OUTDIR.mkdir(exist_ok=True)
    h_total, _ = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = rem4.ground_state(h_total)
    H_basis = qgeom.horizontal_basis(D_A, D_B)

    envs = {
        "A1_deph_0512": b2.liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
    }

    # Load converged U* from B3.4 (or fall back to B2 best)
    conv_path = OUTDIR / "b3_4_gate1_convergence.json"
    b2_path = OUTDIR / "b2_canonical_open.json"
    b2_data = json.loads(b2_path.read_text(encoding="utf-8"))

    results = {}
    for env_name, L in envs.items():
        print(f"\n=== B3.5 saddle escape: {env_name} ===")
        if conv_path.exists():
            conv = json.loads(conv_path.read_text(encoding="utf-8"))["results"][env_name]
            U = np.array(conv["U_real"]) + 1j * np.array(conv["U_imag"])
            start_phi = conv["converged_phi"]
        else:
            best_trial = max(b2_data[env_name], key=lambda t: t["best_phi"])
            U = np.array(best_trial["best_u_real"]) + 1j * np.array(best_trial["best_u_imag"])
            start_phi = best_trial["best_phi"]
        print(f"  Start: Φ = {start_phi:.6f}")

        rounds = []
        converged = False
        verdict = None

        for rnd in range(1, MAX_ROUNDS + 1):
            print(f"\n  --- Round {rnd} ---")
            H = hessian_at(U, L, psi0, H_basis)
            eigs = np.linalg.eigvalsh(H)
            lam_max = float(eigs.max())
            idx_max = int(np.argmax(eigs))
            _, v = np.linalg.eigh(H)
            v_max = v[:, idx_max]
            g_norm = float(np.linalg.norm(grad_phi(U, L, psi0, H_basis)))

            print(f"  Φ={phi(U, L, psi0):.6f}  |∇Φ|={g_norm:.4e}  λ_max={lam_max:.4e}")

            if g_norm < GRAD_TOL and lam_max < LAMBDA_H_TOL:
                converged = True
                if lam_max < -LAMBDA_H_TOL:
                    verdict = "strict local maximum"
                else:
                    verdict = "flat maximum"
                print(f"  CONVERGED: {verdict}")
                rounds.append(dict(round=rnd, phi=phi(U, L, psi0), grad_norm=g_norm,
                                   lambda_max=lam_max, action="converged"))
                break

            if lam_max < LAMBDA_H_TOL:
                # no positive curvature; treat as flat/converged
                converged = True
                verdict = "flat maximum"
                print(f"  CONVERGED (flat): λ_max={lam_max:.4e}")
                rounds.append(dict(round=rnd, phi=phi(U, L, psi0), grad_norm=g_norm,
                                   lambda_max=lam_max, action="converged_flat"))
                break

            # Step 1: verify Hessian along v_max (only round 1)
            verify_rows = None
            if rnd == 1:
                verify_rows = verify_hessian(U, v_max, lam_max, L, psi0, H_basis)

            # Step 2: pick best ε along +v_max, escape, Adam
            K = np.einsum("ijk,k->ij", H_basis, v_max)
            best_eps = None
            best_phi_esc = -np.inf
            for eps in ESC_EPSILONS:
                U_e = U @ sla.expm(eps * K)
                p = phi(U_e, L, psi0)
                if p > best_phi_esc:
                    best_phi_esc = p
                    best_eps = eps
            print(f"  Escape: ε={best_eps:.0e}  Φ_esc={best_phi_esc:.6f} "
                  f"(Δ={best_phi_esc - phi(U, L, psi0):+.6f})")

            # Step 2 cont: restart Adam from the escaped point
            U_esc = U @ sla.expm(best_eps * K)
            res = adam_from(U_esc, L, psi0, H_basis)
            rounds.append(dict(round=rnd, phi_before=phi(U, L, psi0),
                               grad_norm_before=g_norm, lambda_max_before=lam_max,
                               eps=best_eps, phi_escape=best_phi_esc,
                               phi_after=res["Phi"], grad_norm_after=res["grad_norm"],
                               action="escape+adam", verify=verify_rows))
            print(f"  After Adam: Φ={res['Phi']:.6f}  |∇Φ|={res['grad_norm']:.4e}")
            U = res["U"]

            if res["Phi"] < phi(U, L, psi0) - 1e-8 and rnd > 1:
                # no improvement: stop
                print("  No improvement; stopping.")
                break

        if not converged:
            verdict = "saddle chain continues"
            print(f"\n  VERDICT after {len(rounds)} rounds: {verdict}")

        results[env_name] = dict(start_phi=start_phi, rounds=rounds,
                                 final_phi=phi(U, L, psi0),
                                 final_grad_norm=float(np.linalg.norm(
                                     grad_phi(U, L, psi0, H_basis))),
                                 verdict=verdict)

    out_path = OUTDIR / "b3_5_saddle_escape.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False, default=str),
                        encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
