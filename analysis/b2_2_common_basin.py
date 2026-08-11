#!/usr/bin/env python3
"""
Phase B2.2 — Common optimal basin test for A1/A2.

Purpose: verify that the top A1/A2 solutions share a common optimal basin,
not just "the best U happened to be close".

Method (per master 2026-08-11):
  1. Collect the top 20 solutions from A1 and A2 (candidate pool S).
  2. Cross-evaluate every candidate F_i in S under BOTH environments:
       Phi_A1(F_i), Phi_A2(F_i)
  3. Rank candidates under each environment.
  4. If the same candidate group occupies the top ranks under both
     environments -> common optimal basin evidence.

Also reports the overlap of the top-k sets under both rankings.
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
DIM = 8
LAMBDA = 0.2
TOP_K = 20
ENVS = ["A1_deph_0512", "A2_deph_2105"]


def load_U_from_json(env_results: dict, env: str) -> list:
    """Load all U* (with best_phi) from saved results."""
    trials = env_results[env]
    out = []
    for t in trials:
        u = np.array(t["best_u_real"]) + 1j * np.array(t["best_u_imag"])
        out.append(dict(env=env, trial=t["trial"], best_phi=t["best_phi"], U=u))
    return out


def cross_evaluate_pool(pool: list, L_map: dict, psi0: np.ndarray) -> dict:
    """Evaluate every candidate under every environment."""
    rho0 = np.outer(psi0, psi0.conj())
    results = {}
    for env_row in ENVS:
        L = L_map[env_row]
        Y0 = (L @ rho0.reshape(-1, order="F")).reshape(DIM, DIM, order="F")
        scores = []
        for cand in pool:
            U = cand["U"]
            rho_U = U.conj().T @ rho0 @ U
            Y_U = U.conj().T @ Y0 @ U
            obj = b2.canonical_open_objective(rho_U, Y_U, N_A)
            scores.append(dict(
                env=cand["env"], trial=cand["trial"],
                Phi=obj["Phi"], I=obj["I"], gamma=obj["gamma"],
            ))
        # Rank descending by Phi
        scores.sort(key=lambda s: s["Phi"], reverse=True)
        results[env_row] = scores
    return results


def topk_overlap(rank_a1: list, rank_a2: list, k: int = 10) -> float:
    """Fraction of the top-k under A1 that also appear in the top-k under A2."""
    top_a1 = {(r["env"], r["trial"]) for r in rank_a1[:k]}
    top_a2 = {(r["env"], r["trial"]) for r in rank_a2[:k]}
    inter = top_a1 & top_a2
    return len(inter) / k


def main():
    data_path = OUTDIR / "b2_canonical_open.json"
    if not data_path.exists():
        print("Run analysis/b2_canonical_open.py first.")
        return
    data = json.loads(data_path.read_text(encoding="utf-8"))

    # Build candidate pool: top 20 from A1 and top 20 from A2
    pool = []
    for env in ENVS:
        trials = data[env]
        sorted_trials = sorted(trials, key=lambda t: t["best_phi"], reverse=True)[:TOP_K]
        for t in sorted_trials:
            u = np.array(t["best_u_real"]) + 1j * np.array(t["best_u_imag"])
            pool.append(dict(env=env, trial=t["trial"], best_phi=t["best_phi"], U=u))

    h_total, _ = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = rem4.ground_state(h_total)
    L_map = {
        "A1_deph_0512": b2.liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
        "A2_deph_2105": b2.liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
    }

    ranked = cross_evaluate_pool(pool, L_map, psi0)

    print("=== Common optimal basin test (A1/A2) ===")
    print(f"candidate pool: top {TOP_K} from A1 + top {TOP_K} from A2 = {len(pool)}")
    print("\nTop 10 under A1 environment:")
    for i, r in enumerate(ranked["A1_deph_0512"][:10]):
        print(f"  {i+1:2d}. {r['env']:4s} trial {r['trial']:3d}  "
              f"Phi_A1={r['Phi']:.6f}  I={r['I']:.4f}  G={r['gamma']:.4f}")
    print("\nTop 10 under A2 environment:")
    for i, r in enumerate(ranked["A2_deph_2105"][:10]):
        print(f"  {i+1:2d}. {r['env']:4s} trial {r['trial']:3d}  "
              f"Phi_A2={r['Phi']:.6f}  I={r['I']:.4f}  G={r['gamma']:.4f}")

    overlap10 = topk_overlap(ranked["A1_deph_0512"], ranked["A2_deph_2105"], k=10)
    overlap5 = topk_overlap(ranked["A1_deph_0512"], ranked["A2_deph_2105"], k=5)
    print(f"\ntop-5 overlap: {overlap5:.2f}")
    print(f"top-10 overlap: {overlap10:.2f}")

    # What fraction of A1's top-10 are A1-native vs A2-native?
    a1_native = sum(1 for r in ranked["A1_deph_0512"][:10] if r["env"] == "A1_deph_0512")
    a2_native = sum(1 for r in ranked["A2_deph_2105"][:10] if r["env"] == "A2_deph_2105")
    print(f"A1 top-10 native fraction: {a1_native}/10")
    print(f"A2 top-10 native fraction: {a2_native}/10")

    # Common basin conclusion
    if overlap10 >= 0.7:
        basin = "COMMON BASIN (strong evidence)"
    elif overlap10 >= 0.4:
        basin = "COMMON BASIN (moderate evidence)"
    else:
        basin = "NO common basin (rankings diverge)"
    print(f"\nConclusion: {basin}")

    out = dict(pool_size=len(pool), ranked=ranked,
               overlap_top5=overlap5, overlap_top10=overlap10,
               a1_native_top10=a1_native, a2_native_top10=a2_native,
               conclusion=basin)
    out_path = OUTDIR / "b2_2_common_basin.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str),
                        encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
