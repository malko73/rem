#!/usr/bin/env python3
"""
Phase B2.1 — Gauge-invariant TPS distance + cross-environment objective matrix.

Consumes the saved U* from b2_canonical_open.py (analysis_output/b2_canonical_open.json)
and computes:

  1. Gauge-invariant TPS distance d_F between factorizations:
       A_F = U (u(d_A) (x) I_B) U†   (A-side local operator algebra)
       P_A_F = operator-space projector onto the traceless part of A_F
       d_F(F1, F2) = ||P_A_F1 - P_A_F2||_F / sqrt(2 * rank(P_A))
     This is invariant under U -> U (V_A (x) V_B).

  2. Cross-environment objective matrix:
       Phi_{A1}(U*_{A1}) vs Phi_{A1}(U*_{A2})
       Phi_{A2}(U*_{A2}) vs Phi_{A2}(U*_{A1})
     If each environment prefers its own optimum, F* = F*(rho, L, lambda).

  3. B2 verdict: B2-A (environment-dependent factorization) /
                B2-B (same factorization, env-dependent score) /
                B2-C (multiple/unstable optima).
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

N, N_A = 3, 2
D_A, D_B = 4, 2
DIM = 8
LAMBDA = 0.2


def tps_projector(U: np.ndarray, d_a: int = D_A, d_b: int = D_B) -> np.ndarray:
    """Operator-space projector onto A_F = U (u(d_A) (x) I_B) U† (traceless part)."""
    dim = d_a * d_b
    G = qgeom.u_basis(d_a)                       # (d_a, d_a, d_a²)
    vecs = []
    for k in range(d_a * d_a):
        M = U @ np.kron(G[:, :, k], np.eye(d_b, dtype=complex)) @ U.conj().T
        # remove trace to make the generator traceless (u -> su)
        M = M - (np.trace(M) / dim) * np.eye(dim, dtype=complex)
        vecs.append(M.reshape(-1))
    V = np.stack(vecs, axis=1)                   # (dim², d_a²)
    gram = V.conj().T @ V
    pinv = np.linalg.pinv(gram)
    P = V @ pinv @ V.conj().T
    return P


def tps_distance(U1: np.ndarray, U2: np.ndarray) -> float:
    """Gauge-invariant TPS distance d_F(F1, F2) in [0, 1]."""
    P1 = tps_projector(U1)
    P2 = tps_projector(U2)
    rank = np.linalg.matrix_rank(P1, tol=1e-8)
    return float(np.linalg.norm(P1 - P2, ord="fro") / np.sqrt(2.0 * max(rank, 1)))


def load_results(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for env, trials in data.items():
        cleaned = []
        for t in trials:
            u = np.array(t["best_u_real"], dtype=complex) \
                + 1j * np.array(t["best_u_imag"], dtype=complex)
            cleaned.append(dict(trial=t["trial"], best_phi=t["best_phi"],
                                U=u, I=t.get("best_I", float("nan")),
                                gamma=t.get("best_gamma", float("nan"))))
        out[env] = cleaned
    return out


def pairwise_distance_matrix(trials: list) -> np.ndarray:
    """Pairwise d_F among all U* in one environment."""
    n = len(trials)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = tps_distance(trials[i]["U"], trials[j]["U"])
            D[i, j] = D[j, i] = d
    return D


def within_env_stats(D: np.ndarray) -> dict:
    n = D.shape[0]
    off = D[np.triu_indices(n, k=1)]
    return dict(mean=float(off.mean()), median=float(np.median(off)),
                std=float(off.std()), max=float(off.max()),
                frac_lt_01=float(np.mean(off < 0.1)))


def main():
    data_path = OUTDIR / "b2_canonical_open.json"
    if not data_path.exists():
        print("Run analysis/b2_canonical_open.py first (saves U*).")
        return

    env_results = load_results(data_path)
    env_names = list(env_results.keys())

    # 1. Within-environment distance
    print("=== Within-environment TPS distance (d_F) ===")
    within = {}
    for env in env_names:
        trials = env_results[env]
        D = pairwise_distance_matrix(trials)
        stats = within_env_stats(D)
        within[env] = stats
        print(f"{env}: n={len(trials)}  mean d_F={stats['mean']:.4f}  "
              f"median={stats['median']:.4f}  max={stats['max']:.4f}  "
              f"frac<0.1={stats['frac_lt_01']:.3f}")

    # 2. Cross-environment distance (best U* per env)
    print("\n=== Cross-environment d_F (best U* per env) ===")
    best_u = {env: max(env_results[env], key=lambda t: t["best_phi"])["U"]
              for env in env_names}
    cross = {}
    for i, e1 in enumerate(env_names):
        for e2 in env_names[i + 1:]:
            d = tps_distance(best_u[e1], best_u[e2])
            cross[f"{e1} vs {e2}"] = d
            print(f"d_F(F*_{e1}, F*_{e2}) = {d:.4f}")

    # 3. Cross-environment objective matrix
    print("\n=== Cross-environment objective matrix ===")
    h_total, _ = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = rem4.ground_state(h_total)
    rho0 = np.outer(psi0, psi0.conj())

    # Rebuild Liouvillians (import b2 functions)
    import b2_canonical_open as b2
    L_map = {
        "uniform_dephasing": b2.liouvillian_dephasing(h_total, [1.0, 1.0, 1.0]),
        "A1_deph_0512": b2.liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
        "A2_deph_2105": b2.liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
        "amp_damping": b2.liouvillian_amp_damping(h_total, [1.0, 1.0, 1.0]),
        "mixed": b2.liouvillian_mixed(h_total, [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),
    }

    matrix = {}
    for env_row in env_names:
        for env_col in env_names:
            U = best_u[env_col]
            rho_U = U.conj().T @ rho0 @ U
            L = L_map[env_row]
            Y0 = (L @ rho0.reshape(-1, order="F")).reshape(DIM, DIM, order="F")
            Y_U = U.conj().T @ Y0 @ U
            obj = b2.canonical_open_objective(rho_U, Y_U, N_A)
            matrix[f"Phi_{env_row}(U*_{env_col})"] = obj
            if env_row == "A1_deph_0512" or env_row == "A2_deph_2105":
                print(f"Phi_{env_row}(U*_{env_col}) = {obj['Phi']:.6f}  "
                      f"(I={obj['I']:.4f}, G={obj['gamma']:.4f})")

    # 4. Verdict
    print("\n=== B2 Verdict ===")
    # A1/A2 cross check
    phi_a1_a1 = matrix["Phi_A1_deph_0512(U*_A1_deph_0512)"]["Phi"]
    phi_a1_a2 = matrix["Phi_A1_deph_0512(U*_A2_deph_2105)"]["Phi"]
    phi_a2_a2 = matrix["Phi_A2_deph_2105(U*_A2_deph_2105)"]["Phi"]
    phi_a2_a1 = matrix["Phi_A2_deph_2105(U*_A1_deph_0512)"]["Phi"]
    d_a1a2 = cross["A1_deph_0512 vs A2_deph_2105"]
    w_a1 = within["A1_deph_0512"]["mean"]
    w_a2 = within["A2_deph_2105"]["mean"]

    print(f"A1 prefers own optimum: {phi_a1_a1:.4f} vs {phi_a1_a2:.4f} "
          f"(delta={phi_a1_a1 - phi_a1_a2:+.4f})")
    print(f"A2 prefers own optimum: {phi_a2_a2:.4f} vs {phi_a2_a1:.4f} "
          f"(delta={phi_a2_a2 - phi_a2_a1:+.4f})")
    print(f"d_F(F*_A1, F*_A2) = {d_a1a2:.4f} (within-env mean {w_a1:.4f}/{w_a2:.4f})")

    verdict = "B2-C"
    if d_a1a2 > 0.05 and phi_a1_a1 > phi_a1_a2 and phi_a2_a2 > phi_a2_a1:
        verdict = "B2-A"
    elif d_a1a2 < 0.05:
        verdict = "B2-B"
    print(f"VERDICT: {verdict}")

    # Save
    out = dict(within=within, cross=cross, matrix=matrix, verdict=verdict)
    out_path = OUTDIR / "b2_1_cross_eval.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str),
                        encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
