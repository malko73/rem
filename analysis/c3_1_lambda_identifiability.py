#!/usr/bin/env python3
"""
Phase C3.1 — λ Identification Audit: environment-only timescales.

Purpose: pre-define candidate λ values from the ENVIRONMENT ONLY (no
factorization input), to test whether λ can be determined independently
of F* optimization (Gate 3).

Candidates ([λ] = T):
  1. τ_Γ  = 1/Γ_ref       : reference structural decoherence time
     Γ_ref = <Γ_F> over a uniform random TPS sample (env-only average)
  2. τ_L  = 1/|Re μ₁|     : Liouvillian gap (slowest non-zero relaxation mode)
  3. τ_γ  = 1/γ̄           : mean local dephasing rate
  4. τ_κ  = 1/κ̄           : mean local damping rate
  5. τ_env = τ_L           : (alias; we keep the gap as the principal candidate)

For each of the 5 C1 environments, compute all candidate timescales.
NO factorization is used (except the env-only Γ_ref average, which is a
TPS-independent reference — it uses random U but does not optimize F).
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
D_A, D_B = 4, 2
DIM = 8
N_SAMPLE = 200  # random TPS for the env-only Γ_ref average


def liouvillian_gap(L: np.ndarray) -> float:
    """|Re μ₁| where μ₁ is the slowest non-zero eigenvalue of L.

    Zero eigenvalue (steady state) is excluded. For a Lindblad generator
    the non-zero eigenvalues have Re μ <= 0; the gap is the smallest
    |Re μ| among them.
    """
    evals = np.linalg.eigvals(L)
    reals = evals.real
    # exclude near-zero eigenvalues (steady-state / conserved)
    nonzero = reals[np.abs(reals) > 1e-8]
    if len(nonzero) == 0:
        return float("nan")
    return float(np.max(np.abs(nonzero)))  # slowest non-zero mode


def gamma_ref_env_only(L: np.ndarray, psi0: np.ndarray,
                       H_basis: np.ndarray, n_sample: int = N_SAMPLE) -> float:
    """Average Γ_F over random TPS (env-only reference, no optimization)."""
    rng = np.random.default_rng(20260812)
    rho0 = np.outer(psi0, psi0.conj())
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    gammas = []
    for _ in range(n_sample):
        theta = rng.standard_normal(H_basis.shape[2]) * 1.0
        U = sla.expm(np.einsum("ijk,k->ij", H_basis, theta))
        rho_U = U.conj().T @ rho0 @ U
        Y_U = U.conj().T @ Y0 @ U
        obj = b2.canonical_open_objective(rho_U, Y_U, N_A)
        if np.isfinite(obj["gamma"]):
            gammas.append(obj["gamma"])
    return float(np.mean(gammas))


def main():
    OUTDIR.mkdir(exist_ok=True)
    h_total, _ = rem4.xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
    _, psi0 = rem4.ground_state(h_total)
    H_basis = None
    spec2 = importlib.util.spec_from_file_location("quotient_geometry", REPO / "src" / "quotient_geometry.py")
    qgeom = importlib.util.module_from_spec(spec2)
    assert spec2.loader is not None
    sys.modules[spec2.name] = qgeom
    spec2.loader.exec_module(qgeom)
    H_basis = qgeom.horizontal_basis(D_A, D_B)

    envs = {
        "uniform_dephasing": (b2.liouvillian_dephasing(h_total, [1.0, 1.0, 1.0]),
                              [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]),
        "A1_deph_0512": (b2.liouvillian_dephasing(h_total, [0.5, 1.0, 2.0]),
                         [0.5, 1.0, 2.0], [0.0, 0.0, 0.0]),
        "A2_deph_2105": (b2.liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
                         [2.0, 1.0, 0.5], [0.0, 0.0, 0.0]),
        "amp_damping": (b2.liouvillian_amp_damping(h_total, [1.0, 1.0, 1.0]),
                        [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        "mixed": (b2.liouvillian_mixed(h_total, [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),
                  [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),
    }

    results = {}
    print("=== C3.1: environment-only λ candidates ===")
    print(f"{'env':20s} {'τ_Γ':>10s} {'τ_L(gap)':>10s} {'τ_γ':>8s} {'τ_κ':>8s} "
          f"{'Γ_ref':>8s} {'gap':>8s}")
    for env_name, (L, gv, kv) in envs.items():
        gap = liouvillian_gap(L)
        g_ref = gamma_ref_env_only(L, psi0, H_basis)
        tau_G = 1.0 / g_ref if g_ref > 0 else float("nan")
        tau_L = 1.0 / gap if gap > 0 else float("nan")
        g_mean = float(np.mean(gv)) if np.any(np.array(gv) > 0) else float("nan")
        k_mean = float(np.mean(kv)) if np.any(np.array(kv) > 0) else float("nan")
        tau_g = 1.0 / g_mean if g_mean > 0 else float("nan")
        tau_k = 1.0 / k_mean if k_mean > 0 else float("nan")

        results[env_name] = dict(
            tau_Gamma=tau_G, tau_L=tau_L, tau_gamma=tau_g, tau_kappa=tau_k,
            Gamma_ref=g_ref, liouvillian_gap=gap,
            gamma_vec=gv, kappa_vec=kv,
        )
        print(f"{env_name:20s} {tau_G:10.4f} {tau_L:10.4f} "
              f"{tau_g:8.4f} {tau_k:8.4f} {g_ref:8.4f} {gap:8.4f}")

    out_path = OUTDIR / "c3_1_lambda_identifiability.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
