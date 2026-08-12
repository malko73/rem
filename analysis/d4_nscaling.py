#!/usr/bin/env python3
"""
Phase D4 — N-scaling / finite-size persistence of the structural crossover.

Question (master's D4): does the D3 structural crossover persist as the
Hilbert-space size increases? N=3,4,5 only — establishes "finite-size
persistence / finite-size trend", NOT a scaling law or thermodynamic limit.

D4-A  canonical well-posedness vs N: optimize with J_dyn^(0) ONLY.
D4-B  crossover persistence: Haar primary. For each N extract competing
      basins F_A^(N) (optimize at small tau) and F_B^(N) (optimize at large
      tau), then measure the objective crossing
          DeltaPhi_N(tau) = Phi_{tau,N}(F_B^(N)) - Phi_{tau,N}(F_A^(N))
      and tau_c(N) from DeltaPhi_N(tau_c(N)) = 0. As in D3, tau_c is NEVER
      defined by optimizer basin hops.
D4-C  cross-size comparison: keep raw quantities, add normalized
      indicators (I~ = I/(2 log d_min), J~ = J/gamma_mean). The canonical
      definition is NOT renormalized (no denominator).

Protocol (documented reductions for large N):
  N=3: full 45-dim quotient, 200 steps, 6 seeds (frozen D2/D3 protocol)
  N=4: full quotient (cut 2|2: 225-dim; cut 1|3: 189-dim), 100 steps, 4 seeds
  N=5: full 945-dim quotient (cut 2|3), 80 steps, 3 seeds
  lambda=0.2, lr=0.01, SEED0=20260813, dephasing gamma=(0.5,1,2) periodic
  extension, asymmetric XY chain J=(1.5,0.6,...), h=0.2.

Usage:
  python d4_nscaling.py d4a --n 4 --state haar            # D4-A one run
  python d4_nscaling.py d4b --n 4 --cut 22 --tau-b 0.1    # D4-B one N (Haar)
  python d4_nscaling.py d4b --n 4 --cut 13 --tau-b 0.1    # N=4 asymmetric bonus
Writes analysis_output/d4_{a,b}_n{N}_cut{cut}_{state|haar}.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from quotient_geometry import horizontal_basis
from rem4_numerical import xy_chain_hamiltonian, ground_state, mutual_information_for_cut
from gamma_F import schmidt_basis, dephasing_projection, liouvillian_env
from b2_1_cross_eval import tps_distance

LAMBDA = 0.2
LR = 0.01
SEED0 = 20260813
COUPLINGS = {3: [1.5, 0.6], 4: [1.5, 0.6, 0.6], 5: [1.5, 0.6, 0.6, 0.6]}
GAMMA_BASE = [0.5, 1.0, 2.0]
STEPS = {3: 200, 4: 100, 5: 100}
SEEDS = {3: 6, 4: 4, 5: 4}
# N=5: full quotient is 945-dim (too slow in-session); use a random
# 200-dim horizontal subspace (~21% of the quotient) — subspace-restricted
# lower bound, following the repo's established n_params heuristic for
# large N. The DeltaPhi crossing is objective-level and unaffected.
SUBSPACE_DIM = {5: 200}
SUBSPACE_SEED = 777
# cuts: (n_a, d_A, d_B) per N; N=4 offers balanced 2|2 and asymmetric 1|3
CUTS = {
    3: {"21": (2, 4, 2)},
    4: {"22": (2, 4, 4), "13": (1, 2, 8)},
    5: {"23": (2, 4, 8)},
}


def gamma_vec(N):
    return [GAMMA_BASE[i % len(GAMMA_BASE)] for i in range(N)]


def haar_state(dim, rng):
    m = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    q, _ = np.linalg.qr(m)
    return np.outer(q[:, 0], q[:, 0].conj())


def build_system(N):
    H, _ = xy_chain_hamiltonian(N, COUPLINGS[N], 0.2)
    gvec = gamma_vec(N)
    L = liouvillian_env(H, gvec, [0.0] * N, n=N)
    return H, gvec, L


def build_states(H, N):
    """ground + haar (same construction style as the D-series)."""
    _, psi_gs = ground_state(H)
    rng = np.random.default_rng(20260813)
    return {"ground": np.outer(psi_gs, psi_gs.conj()),
            "haar": haar_state(2 ** N, rng)}


def objective(rho_U, Y_U, N, n_a, tau, lam):
    """Return (Phi, I, J) for J_dyn^(tau) (tau=0 limit handled by caller).

    J_dyn^(0) is used when tau == 'j0' (no rho_tau needed); otherwise the
    finite-time difference ratio with the FIXED t=0 Schmidt basis.
    """
    evals, evecs = np.linalg.eigh(rho_U)
    psi_U = evecs[:, np.argmax(evals)]
    I = mutual_information_for_cut(psi_U, N, n_a)
    basis = schmidt_basis(psi_U, N, n_a)
    q = lambda r: r - dephasing_projection(r, basis)  # noqa: E731
    x = q(rho_U)
    c0_sq = float(np.linalg.norm(x, ord="fro") ** 2)
    if tau == "j0":
        qy = q(Y_U)
        J = -float(np.real(np.trace(x.conj().T @ qy)))
        return I - lam * J, I, J, c0_sq
    xt = q(Y_U)  # Y_U holds rho(tau) in the finite-time path (see caller)
    ctau_sq = float(np.linalg.norm(xt, ord="fro") ** 2)
    J = -(ctau_sq - c0_sq) / (2.0 * tau)
    return I - lam * J, I, J, c0_sq


def adam_optimize(U0, rho0, Y0, rho_tau, N, n_a, tau, H_basis, steps, lr=LR):
    """Adam over Phi = I - lambda*J^(tau). tau='j0' uses J_dyn^(0).

    Y0 must be L(rho0) for 'j0', or rho(tau) for finite tau (the caller
    prepares the right second operand).
    """
    delta = np.zeros(H_basis.shape[2])
    m = np.zeros_like(delta)
    v = np.zeros_like(delta)
    beta1, beta2, eps = 0.9, 0.999, 1e-8

    def to_U(d):
        return U0 @ sla.expm(np.einsum("ijk,k->ij", H_basis, d))

    def phi_at(d):
        U = to_U(d)
        rho_U = U.conj().T @ rho0 @ U
        Y_U = U.conj().T @ Y0 @ U
        o, _, _, _ = objective(rho_U, Y_U, N, n_a, tau, LAMBDA)
        return o

    best_phi = phi_at(delta)
    best_U = to_U(delta)
    best_obj = None
    for step in range(steps):
        grad = np.zeros_like(delta)
        for k in range(len(delta)):
            d_p = delta.copy()
            d_m = delta.copy()
            d_p[k] += 1e-5
            d_m[k] -= 1e-5
            grad[k] = (phi_at(d_p) - phi_at(d_m)) / (2 * 1e-5)
        m = beta1 * m + (1 - beta1) * grad
        v = beta2 * v + (1 - beta2) * grad ** 2
        m_hat = m / (1 - beta1 ** (step + 1))
        v_hat = v / (1 - beta2 ** (step + 1))
        delta = delta + lr * m_hat / (np.sqrt(v_hat) + eps)
        p = phi_at(delta)
        if p > best_phi:
            best_phi = p
            best_U = to_U(delta)
            rho_U = best_U.conj().T @ rho0 @ best_U
            Y_U = best_U.conj().T @ Y0 @ best_U
            best_obj = objective(rho_U, Y_U, N, n_a, tau, LAMBDA)
    if best_obj is None:
        rho_U = best_U.conj().T @ rho0 @ best_U
        Y_U = best_U.conj().T @ Y0 @ best_U
        best_obj = objective(rho_U, Y_U, N, n_a, tau, LAMBDA)
    return best_U, best_phi, best_obj


def schmidt_pmax(psi, N, n_a):
    d_a, d_b = 2 ** n_a, 2 ** (N - n_a)
    _, s, _ = np.linalg.svd(psi.reshape((d_a, d_b)))
    return float((s ** 2).max())


def dominant_eigenvector(rho):
    evals, evecs = np.linalg.eigh(rho)
    return evecs[:, np.argmax(evals)]


def horizontal_basis_for(N, d_A, d_B):
    """Full horizontal basis, or a random subspace for large N."""
    H_basis = horizontal_basis(d_A, d_B)
    n_h = H_basis.shape[2]
    sub = SUBSPACE_DIM.get(N)
    if sub is None or sub >= n_h:
        return H_basis, n_h, False
    rng = np.random.default_rng(SUBSPACE_SEED)
    W = rng.standard_normal((n_h, sub))
    Q, _ = np.linalg.qr(W)  # (n_h, sub) with orthonormal columns
    # g_k = sum_m Q[m, k] G_m  ->  H_sub[:, :, k] = einsum('mk,ijm->ijk', Q, H_basis)
    H_sub = np.einsum("mk,ijm->ijk", Q, H_basis)
    return H_sub, sub, True


def run_d4a(N, state_name, cut_key):
    n_a, d_A, d_B = CUTS[N][cut_key]
    H, gvec, L = build_system(N)
    states = build_states(H, N)
    rho0 = states[state_name]
    Y0 = (L @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
    H_basis, n_h, is_sub = horizontal_basis_for(N, d_A, d_B)
    steps, seeds_n = STEPS[N], SEEDS[N]

    print(f"=== D4-A N={N} cut={cut_key} state={state_name} "
          f"(dim {n_h}{' [subspace]' if is_sub else ''}, {steps} steps, {seeds_n} seeds) ===")
    cont = objective(rho0, Y0, N, n_a, "j0", LAMBDA)
    print(f"  contiguous: Phi={cont[0]:.4f} I={cont[1]:.4f} J={cont[2]:.4f}")

    seeds = []
    t0 = time.time()
    for s in range(seeds_n):
        seed = SEED0 + s
        rng = np.random.default_rng(seed)
        theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
        U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))
        U_best, phi_best, obj = adam_optimize(U0, rho0, Y0, Y0, N, n_a, "j0",
                                              H_basis, steps)
        psi0 = dominant_eigenvector(rho0)
        p_max = schmidt_pmax(U_best.conj().T @ psi0, N, n_a)
        seeds.append(dict(seed=seed, Phi=phi_best, I=obj[1], J_dyn0=obj[2],
                          C_F_sq=obj[3], p_max=p_max,
                          U_real=U_best.real.tolist(), U_imag=U_best.imag.tolist()))
    phis = np.array([d["Phi"] for d in seeds])
    success = float(np.mean([d["Phi"] > cont[0] for d in seeds]))
    result = dict(
        N=N, cut=cut_key, state=state_name, protocol=dict(
            functional="J_dyn^(0) (Spec v2.2 canonical)",
            quotient_dim=n_h, subspace=is_sub, steps=steps, seeds=seeds_n,
            gamma=gamma_vec(N), couplings=COUPLINGS[N]),
        contiguous=dict(Phi=cont[0], I=cont[1], J=cont[2]),
        seeds=seeds,
        best_Phi=float(phis.max()), median_Phi=float(np.median(phis)),
        std_Phi=float(phis.std()), success_rate=success,
        min_C_F_sq=float(min(d["C_F_sq"] for d in seeds)),
        max_p_max=float(max(d["p_max"] for d in seeds)),
        elapsed=time.time() - t0,
    )
    out = OUTDIR / f"d4_a_n{N}_cut{cut_key}_{state_name}.json"
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"  best Phi={result['best_Phi']:.4f} std={result['std_Phi']:.4f} "
          f"success={success:.2f} minC0={result['min_C_F_sq']:.4f} "
          f"pmax={result['max_p_max']:.4f} ({time.time()-t0:.0f}s)")
    print(f"Saved to {out}")


def run_d4b(N, cut_key, tau_b=0.1, tau_a=1e-3, reps_only="AB"):
    n_a, d_A, d_B = CUTS[N][cut_key]
    H, gvec, L = build_system(N)
    states = build_states(H, N)
    rho0 = states["haar"]
    H_basis, n_h, is_sub = horizontal_basis_for(N, d_A, d_B)
    steps, seeds_n = STEPS[N], SEEDS[N]

    print(f"=== D4-B N={N} cut={cut_key} Haar "
          f"(dim {n_h}{' [subspace]' if is_sub else ''}, {steps} steps, {seeds_n} seeds) ===")

    def prep(tau):
        rho_tau = (sla.expm(tau * L) @ rho0.reshape(-1, order="F")).reshape(rho0.shape, order="F")
        return rho_tau

    reps = {}
    for label, tau in (("A", tau_a), ("B", tau_b)):
        if label not in reps_only:
            continue
        rho_tau = prep(tau)
        seeds = []
        t0 = time.time()
        for s in range(seeds_n):
            seed = SEED0 + s
            rng = np.random.default_rng(seed)
            theta0 = rng.standard_normal(H_basis.shape[2]) * 1.0
            U0 = sla.expm(np.einsum("ijk,k->ij", H_basis, theta0))
            U_best, phi_best, obj = adam_optimize(U0, rho0, rho_tau, rho_tau, N, n_a,
                                                  tau, H_basis, steps)
            psi0 = dominant_eigenvector(rho0)
            p_max = schmidt_pmax(U_best.conj().T @ psi0, N, n_a)
            seeds.append(dict(seed=seed, Phi=phi_best, I=obj[1], J_tau=obj[2],
                              C_F_sq_0=obj[3], p_max=p_max,
                              U_real=U_best.real.tolist(), U_imag=U_best.imag.tolist()))
        phis = np.array([d["Phi"] for d in seeds])
        b = max(seeds, key=lambda d: d["Phi"])
        reps[label] = dict(
            tau=tau, seeds=seeds, best_Phi=float(phis.max()),
            median_Phi=float(np.median(phis)), std_Phi=float(phis.std()),
            best_U_real=b["U_real"], best_U_imag=b["U_imag"],
            elapsed=time.time() - t0,
        )
        print(f"  F_{label} (tau={tau}): best {reps[label]['best_Phi']:.4f} "
              f"std {reps[label]['std_Phi']:.4f} ({time.time()-t0:.0f}s)")

    if "A" not in reps or "B" not in reps:
        # partial run: save what we have and exit (caller runs the other rep)
        out = OUTDIR / f"d4_b_n{N}_cut{cut_key}_haar.json"
        if out.exists():
            prev = json.loads(out.read_text(encoding="utf-8"))
            prev.update({k: v for k, v in reps.items()})
            out.write_text(json.dumps(prev, indent=2, ensure_ascii=False), encoding="utf-8")
            print(f"Updated {out}")
        return

    F_A = np.array(reps["A"]["best_U_real"]) + 1j * np.array(reps["A"]["best_U_imag"])
    F_B = np.array(reps["B"]["best_U_real"]) + 1j * np.array(reps["B"]["best_U_imag"])
    dAB = tps_distance(F_A, F_B, d_A, d_B)
    print(f"  d_F(F_A, F_B) = {dAB:.4f}")

    # DeltaPhi on a dense tau grid (objective-level, fixed reps)
    taus = [1e-3, 2e-3, 3e-3, 5e-3, 8e-3, 1e-2, 1.5e-2, 2e-2, 3e-2, 4e-2,
            5e-2, 7e-2, 1e-1, 1.5e-1, 2e-1, 3e-1, 5e-1, 7e-1, 1.0]
    dphi = {}
    for tau in taus:
        rho_tau = prep(tau)
        oA = objective(F_A.conj().T @ rho0 @ F_A, F_A.conj().T @ rho_tau @ F_A,
                       N, n_a, tau, LAMBDA)
        oB = objective(F_B.conj().T @ rho0 @ F_B, F_B.conj().T @ rho_tau @ F_B,
                       N, n_a, tau, LAMBDA)
        dphi[tau] = oB[0] - oA[0]

    tau_c = None
    for i in range(len(taus) - 1):
        t0, t1 = taus[i], taus[i + 1]
        d0, d1 = dphi[t0], dphi[t1]
        if d0 * d1 <= 0 and d1 != d0:
            tau_c = t0 + (t1 - t0) * (-d0) / (d1 - d0)
            break
    print("  DeltaPhi(tau) = Phi_tau(F_B) - Phi_tau(F_A):")
    for t in taus:
        print(f"    tau={t}: {dphi[t]:+.4f}")
    print(f"  tau_c estimate = {tau_c}")

    result = dict(
        N=N, cut=cut_key, protocol=dict(
            functional="J_dyn^(tau) = -(C_F^2(tau)-C_F^2(0))/(2tau)",
            quotient_dim=n_h, subspace=is_sub, steps=steps, seeds=seeds_n,
            gamma=gamma_vec(N), couplings=COUPLINGS[N],
            reps="F_A @ tau_a=%.0e, F_B @ tau_b=%.0e" % (tau_a, tau_b)),
        F_A=reps["A"], F_B=reps["B"], d_F_FA_FB=dAB,
        DeltaPhi=dphi, tau_c=tau_c,
    )
    out = OUTDIR / f"d4_b_n{N}_cut{cut_key}_haar.json"
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Saved to {out}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)

    pa = sub.add_parser("d4a")
    pa.add_argument("--n", type=int, required=True)
    pa.add_argument("--state", choices=["ground", "haar"], required=True)
    pa.add_argument("--cut", default="21")

    pb = sub.add_parser("d4b")
    pb.add_argument("--n", type=int, required=True)
    pb.add_argument("--cut", default="22")
    pb.add_argument("--tau-a", type=float, default=1e-3)
    pb.add_argument("--tau-b", type=float, default=0.1)
    pb.add_argument("--reps", default="AB", help="which reps to run: A, B, or AB")

    args = ap.parse_args()
    OUTDIR.mkdir(exist_ok=True)
    if args.mode == "d4a":
        run_d4a(args.n, args.state, args.cut)
    else:
        run_d4b(args.n, args.cut, args.tau_b, args.tau_a, args.reps)


if __name__ == "__main__":
    main()
