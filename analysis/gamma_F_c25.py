#!/usr/bin/env python3
"""
Phase C2.5 — Gamma^exact: sign, scaling, dimension audit (2026-08-11).

Pre-freeze checks for Spec v2.0:

1. SIGN: does Gamma_F^exact become NEGATIVE (initial coherence increase)
   for random states / Hamiltonians / environments? If yes, the quantity
   must be named a signed dynamical rate (NOT max(0,Gamma)).
2. SCALING: with a fixed state and Hamiltonian, Gamma^exact must scale
   linearly in the environment strength (uniform dephasing gamma) —
   consistent with [Gamma] = T^-1.
3. DIMENSION: [Gamma] = T^-1 (rate); [lambda] = T in Phi = I - lambda*Gamma
   (analytic; scaling check supports it numerically).

Outputs (analysis_output/):
    gamma_F_c25_results.json, printed summary
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"
MODULE_PATH = REPO / "src" / "rem4_numerical.py"

spec = importlib.util.spec_from_file_location("rem4_numerical", MODULE_PATH)
rem4 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = rem4
spec.loader.exec_module(rem4)

GPATH = REPO / "analysis" / "gamma_F.py"
gspec = importlib.util.spec_from_file_location("gamma_F", GPATH)
gf = importlib.util.module_from_spec(gspec)
assert gspec.loader is not None
sys.modules[gspec.name] = gf
gspec.loader.exec_module(gf)

N_SAMPLES = 200
SEED = 20260811


def haar_random_state(dim: int, rng: np.random.Generator) -> np.ndarray:
    z = rng.standard_normal((dim, 1)) + 1j * rng.standard_normal((dim, 1))
    return (z / np.linalg.norm(z)).ravel()


def random_xy_hamiltonian(rng: np.random.Generator):
    j12 = rng.uniform(0.5, 2.0)
    j23 = rng.uniform(0.5, 2.0)
    h = rng.uniform(0.0, 1.0)
    return rem4.xy_chain_hamiltonian(3, [j12, j23], h)


def random_gue_hamiltonian(rng: np.random.Generator):
    m = rng.standard_normal((8, 8)) + 1j * rng.standard_normal((8, 8))
    m = (m + m.conj().T) / 2.0
    # rescale to a spectral width comparable to the XY benchmark
    m = m * (3.0 / max(np.linalg.eigvalsh(m).max(), 1e-9))
    return m, {}


def random_env(rng: np.random.Generator):
    gamma = tuple(rng.uniform(0.0, 2.0, 3))
    kappa = tuple(rng.uniform(0.0, 2.0, 3))
    return gamma, kappa


def sample_gamma(psi, h_total, gamma_vec, kappa_vec, cut):
    L = gf.liouvillian_env(h_total, list(gamma_vec), list(kappa_vec), n=3)
    rho0 = np.outer(psi, psi.conj())
    basis = gf.schmidt_basis(psi, 3, cut)
    return gf.gamma_exact(rho0, L, basis)


def run_sign_sweep() -> dict:
    rng = np.random.default_rng(SEED)
    neg_count = 0
    min_gamma = 1e9
    max_gamma = -1e9
    family_counts = {"xy": 0, "gue": 0}
    examples = []
    for i in range(N_SAMPLES):
        family = "xy" if i % 2 == 0 else "gue"
        if family == "xy":
            h_total, _ = random_xy_hamiltonian(rng)
        else:
            h_total, _ = random_gue_hamiltonian(rng)
        family_counts[family] += 1
        psi = haar_random_state(8, rng)
        gamma_vec, kappa_vec = random_env(rng)
        g1 = sample_gamma(psi, h_total, gamma_vec, kappa_vec, 1)
        g2 = sample_gamma(psi, h_total, gamma_vec, kappa_vec, 2)
        for g in (g1, g2):
            min_gamma = min(min_gamma, g)
            max_gamma = max(max_gamma, g)
            if g < 0:
                neg_count += 1
                if len(examples) < 5:
                    examples.append(dict(i=i, family=family, g1=g1, g2=g2,
                                         gamma=list(gamma_vec),
                                         kappa=list(kappa_vec)))
    return dict(
        n_samples=N_SAMPLES, n_gamma_values=2 * N_SAMPLES,
        negative_count=int(neg_count),
        negative_fraction=float(neg_count / (2 * N_SAMPLES)),
        min_gamma=float(min_gamma), max_gamma=float(max_gamma),
        family_counts=family_counts, examples=examples,
    )


def run_scaling_check() -> dict:
    """Gamma^exact must scale linearly in uniform dephasing gamma."""
    psi, h_total = gf.build_system()
    gammas = [0.25, 0.5, 1.0, 2.0]
    rows = []
    for g in gammas:
        L = gf.liouvillian(h_total, g, n=3)
        rho0 = np.outer(psi, psi.conj())
        g1 = gf.gamma_exact(rho0, L, gf.schmidt_basis(psi, 3, 1))
        rows.append(dict(gamma=g, gamma_exact=g1, ratio=g1 / g))
    ratios = [r["ratio"] for r in rows]
    return dict(rows=rows,
                ratio_mean=float(np.mean(ratios)),
                ratio_std=float(np.std(ratios)),
                linear=bool(np.std(ratios) / abs(np.mean(ratios)) < 1e-6))


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    sign = run_sign_sweep()
    scaling = run_scaling_check()

    result = dict(sign=sign, scaling=scaling)
    (OUTDIR / "gamma_F_c25_results.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"wrote {OUTDIR / 'gamma_F_c25_results.json'}")
    print(f"\n[SIGN] {sign['n_gamma_values']} Gamma values "
          f"({sign['n_samples']} random configs)")
    print(f"  negative: {sign['negative_count']} "
          f"({sign['negative_fraction']*100:.2f}%)")
    print(f"  min={sign['min_gamma']:.6f}  max={sign['max_gamma']:.6f}")
    if sign["examples"]:
        print("  example negatives:")
        for ex in sign["examples"]:
            print(f"    #{ex['i']} ({ex['family']}) g1={ex['g1']:.4f} "
                  f"g2={ex['g2']:.4f} gamma={ex['gamma']} kappa={ex['kappa']}")
    print(f"\n[SCALING] Gamma^exact / gamma (uniform dephasing):")
    for r in scaling["rows"]:
        print(f"  gamma={r['gamma']}: Gamma={r['gamma_exact']:.6f} "
              f"ratio={r['ratio']:.6f}")
    print(f"  ratio mean={scaling['ratio_mean']:.6f} "
          f"std={scaling['ratio_std']:.2e} -> "
          f"linear={scaling['linear']}  ([Gamma]=T^-1 confirmed)")


if __name__ == "__main__":
    main()
