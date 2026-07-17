#!/usr/bin/env python3
"""
REM4 numerical scaffold.

Purpose
-------
Implements a concrete exact-diagonalization workflow for REM4 based on an
asymmetric XY/Heisenberg-like chain,

    H = J12 (X1 X2 + Y1 Y2) + J23 (X2 X3 + Y2 Y3) + h sum_i Zi

with optional extension to N = 4, 6, 8 sites.  The script evaluates the REM
variational functional

    Phi(F) = Phi_S(F) + lambda * Phi_H(F)

for candidate contiguous cuts F, constructs the lambda* crossover, generates a
2D phase diagram over (J12/J23, lambda), and estimates a pointer-basis mismatch
angle between the REM-selected and dynamics-selected structures.

This file is deliberately explicit and auditable rather than minimal.
"""
from __future__ import annotations

import argparse
import dataclasses
import math
import os
from typing import Dict, List, Sequence, Tuple

import numpy as np

try:
    import scipy.linalg as sla
except Exception as exc:  # pragma: no cover
    raise RuntimeError("scipy is required for rem4_numerical.py") from exc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


Array = np.ndarray

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_n(ops: Sequence[Array]) -> Array:
    out = np.array([[1.0 + 0.0j]])
    for op in ops:
        out = np.kron(out, op)
    return out


def local_op(n: int, site: int, op: Array) -> Array:
    ops = [I2] * n
    ops[site] = op
    return kron_n(ops)


def two_site_op(n: int, i: int, op_i: Array, j: int, op_j: Array) -> Array:
    ops = [I2] * n
    ops[i] = op_i
    ops[j] = op_j
    return kron_n(ops)


def xy_chain_hamiltonian(n: int, couplings: Sequence[float], h: float) -> Tuple[Array, Dict[Tuple[int, int], Array]]:
    """Build H = sum_i J_i (X_i X_{i+1} + Y_i Y_{i+1}) + h sum_i Z_i."""
    if len(couplings) != n - 1:
        raise ValueError(f"Expected {n-1} couplings, got {len(couplings)}")

    dim = 2 ** n
    h_total = np.zeros((dim, dim), dtype=complex)
    bond_terms: Dict[Tuple[int, int], Array] = {}

    for i, coupling in enumerate(couplings):
        term = coupling * (
            two_site_op(n, i, X, i + 1, X)
            + two_site_op(n, i, Y, i + 1, Y)
        )
        h_total += term
        bond_terms[(i, i + 1)] = term

    for i in range(n):
        h_total += h * local_op(n, i, Z)

    return h_total, bond_terms


def computational_basis_state(bitstring: str) -> Array:
    idx = int(bitstring, 2)
    vec = np.zeros(2 ** len(bitstring), dtype=complex)
    vec[idx] = 1.0
    return vec


def normalize(psi: Array) -> Array:
    norm = np.linalg.norm(psi)
    if norm == 0:
        raise ValueError("Zero vector cannot be normalized.")
    return psi / norm


def initial_state_family(n: int, alpha: float, beta: float = 0.0) -> Array:
    if n < 3:
        raise ValueError("n must be >= 3")

    mid = n // 2
    basis_a = computational_basis_state("1" + "0" * (n - 1))
    bits_mid = ["0"] * n
    bits_mid[mid] = "1"
    basis_b = computational_basis_state("".join(bits_mid))
    basis_c = computational_basis_state("0" * (n - 1) + "1")

    psi = (
        math.cos(alpha) * basis_a
        + math.sin(alpha) * math.cos(beta) * basis_b
        + math.sin(alpha) * math.sin(beta) * basis_c
    )
    return normalize(psi)


def ground_state(h_total: Array) -> Tuple[float, Array]:
    evals, evecs = sla.eigh(h_total)
    idx = np.argmin(evals.real)
    return float(evals[idx].real), evecs[:, idx]


def time_evolved_state(h_total: Array, psi0: Array, t: float) -> Array:
    u = sla.expm(-1j * h_total * t)
    return normalize(u @ psi0)


def partial_trace(rho: Array, keep: Sequence[int], n: int) -> Array:
    keep = sorted(keep)
    trace_out = [i for i in range(n) if i not in keep]
    reshaped = rho.reshape([2] * n + [2] * n)
    current_n = n
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    dim_keep = 2 ** len(keep)
    return reshaped.reshape((dim_keep, dim_keep))


def von_neumann_entropy(rho: Array, base: float = 2.0) -> float:
    evals = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    evals = np.maximum(evals.real, 0.0)
    nz = evals[evals > 1e-12]
    if len(nz) == 0:
        return 0.0
    log = np.log2 if base == 2.0 else np.log
    return float(-np.sum(nz * log(nz)))


def mutual_information_for_cut(psi: Array, n: int, cut: int, base: float = 2.0) -> float:
    rho = np.outer(psi, psi.conj())
    rho_left = partial_trace(rho, list(range(cut)), n)
    rho_right = partial_trace(rho, list(range(cut, n)), n)
    return (
        von_neumann_entropy(rho_left, base=base)
        + von_neumann_entropy(rho_right, base=base)
        - von_neumann_entropy(rho, base=base)
    )


def boundary_energy(psi: Array, bond_term: Array) -> float:
    return float(np.real(np.vdot(psi, bond_term @ psi)))


def schmidt_data(psi: Array, n: int, cut: int) -> Tuple[Array, Array, Array]:
    mat = psi.reshape((2 ** cut, 2 ** (n - cut)))
    u, s, vh = np.linalg.svd(mat, full_matrices=False)
    return s, u, vh.conj().T


def dominant_schmidt_angle(psi: Array, n: int, cut_a: int, cut_b: int) -> float:
    _, u1, v1 = schmidt_data(psi, n, cut_a)
    _, u2, v2 = schmidt_data(psi, n, cut_b)
    vec1 = normalize(np.kron(u1[:, 0], v1[:, 0]))
    vec2 = normalize(np.kron(u2[:, 0], v2[:, 0]))
    overlap = np.clip(np.abs(np.vdot(vec1, vec2)), 0.0, 1.0)
    return float(np.arccos(overlap))


@dataclasses.dataclass
class CutMetrics:
    cut: int
    phi_s: float
    phi_h: float
    phi: float
    boundary_energy: float


def evaluate_cuts(
    psi: Array,
    n: int,
    bond_terms: Dict[Tuple[int, int], Array],
    lambda_value: float,
    base: float = 2.0,
) -> List[CutMetrics]:
    out: List[CutMetrics] = []
    for cut in range(1, n):
        bond = bond_terms[(cut - 1, cut)]
        mi = mutual_information_for_cut(psi, n, cut, base=base)
        ebd = boundary_energy(psi, bond)
        phi_h = -ebd
        phi = mi + lambda_value * phi_h
        out.append(CutMetrics(cut, mi, phi_h, phi, ebd))
    return out


def rem_selected_cut(metrics: Sequence[CutMetrics]) -> int:
    return max(metrics, key=lambda x: x.phi).cut


def std_selected_cut(metrics: Sequence[CutMetrics]) -> int:
    return max(metrics, key=lambda x: x.phi_h).cut


def info_selected_cut(metrics: Sequence[CutMetrics]) -> int:
    return max(metrics, key=lambda x: x.phi_s).cut


def lambda_star_for_two_cuts(c1: CutMetrics, c2: CutMetrics) -> float | None:
    dphi_s = c1.phi_s - c2.phi_s
    dphi_h = c2.phi_h - c1.phi_h
    if abs(dphi_h) < 1e-12:
        return None
    return dphi_s / dphi_h


def first_crossover_lambda(psi: Array, n: int, bond_terms: Dict[Tuple[int, int], Array]) -> float | None:
    metrics0 = evaluate_cuts(psi, n, bond_terms, lambda_value=0.0)
    i_cut = info_selected_cut(metrics0)
    i_metric = next(m for m in metrics0 if m.cut == i_cut)
    lambdas = []
    for metric in metrics0:
        if metric.cut == i_cut:
            continue
        lam = lambda_star_for_two_cuts(i_metric, metric)
        if lam is not None and lam > 0:
            lambdas.append(lam)
    return float(min(lambdas)) if lambdas else None


def softmax(x: Array) -> Array:
    z = x - np.max(x)
    ez = np.exp(z)
    return ez / np.sum(ez)


def gradient_flow_on_cut_simplex(metrics: Sequence[CutMetrics], steps: int = 400, lr: float = 0.08) -> Tuple[Array, Array]:
    phi = np.array([m.phi for m in metrics], dtype=float)
    eta = np.zeros_like(phi)
    history = []
    for _ in range(steps):
        p = softmax(eta)
        log_p = np.log(np.clip(p, 1e-12, None))
        avg_phi = np.sum(p * phi)
        avg_log_p = np.sum(p * log_p)
        grad = p * (phi - avg_phi - log_p + avg_log_p)
        eta += lr * grad
        history.append(p.copy())
    return softmax(eta), np.array(history)


def save_phase_diagram(outpath: str, ratios: Array, lambdas: Array, region: Array, boundary: Array, title: str) -> None:
    plt.figure(figsize=(6.0, 4.6))
    plt.imshow(
        region,
        origin="lower",
        aspect="auto",
        extent=[ratios.min(), ratios.max(), lambdas.min(), lambdas.max()],
    )
    mask = ~np.isnan(boundary)
    if np.any(mask):
        plt.plot(ratios[mask], boundary[mask], linewidth=1.5)
    plt.xlabel(r"$J_{12}/J_{23}$")
    plt.ylabel(r"$\lambda$")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def save_scaling_plot(outpath: str, ns: Array, lambdas_star: Array) -> None:
    plt.figure(figsize=(5.5, 4.0))
    plt.plot(ns, lambdas_star, marker="o")
    plt.xlabel(r"$N$")
    plt.ylabel(r"$\lambda^\ast(N)$")
    plt.title("Finite-size crossover estimate")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def make_three_qubit_phase_scan(
    outdir: str,
    j23: float,
    h: float,
    lambda_max: float,
    grid_r: int,
    grid_l: int,
    alpha: float,
    beta: float,
    state_mode: str,
    evolve_time: float,
) -> Dict[str, float]:
    ratios = np.linspace(0.25, 4.0, grid_r)
    lambdas = np.linspace(0.0, lambda_max, grid_l)
    region = np.zeros((grid_l, grid_r), dtype=float)
    boundary = np.full(grid_r, np.nan, dtype=float)
    sample_lambda_star = None
    sample_delta_theta = None

    for ix, ratio in enumerate(ratios):
        h_total, bond_terms = xy_chain_hamiltonian(3, [ratio * j23, j23], h)
        if state_mode == "ground":
            _, psi = ground_state(h_total)
        elif state_mode == "initial":
            psi = initial_state_family(3, alpha=alpha, beta=beta)
        elif state_mode == "evolved":
            psi0 = initial_state_family(3, alpha=alpha, beta=beta)
            psi = time_evolved_state(h_total, psi0, t=evolve_time)
        else:
            raise ValueError(f"Unknown state_mode={state_mode}")

        lam_star = first_crossover_lambda(psi, 3, bond_terms)
        if lam_star is not None:
            boundary[ix] = lam_star

        for iy, lam in enumerate(lambdas):
            metrics = evaluate_cuts(psi, 3, bond_terms, lambda_value=lam)
            rem_cut = rem_selected_cut(metrics)
            std_cut = std_selected_cut(metrics)
            region[iy, ix] = 1.0 if rem_cut == std_cut else 0.0
            if abs(ratio - 1.5) < 0.015 and abs(lam - 0.2) < (lambda_max / grid_l) and rem_cut != std_cut:
                sample_delta_theta = dominant_schmidt_angle(psi, 3, std_cut, rem_cut)

        if abs(ratio - 1.5) < 0.015:
            sample_lambda_star = lam_star

    phase_png = os.path.join(outdir, "REM4_phase_exact.png")
    save_phase_diagram(
        phase_png,
        ratios,
        lambdas,
        region,
        boundary,
        "REM4 phase diagram (exact diagonalization scaffold)",
    )
    return {
        "phase_png": phase_png,
        "sample_lambda_star_at_ratio_1p5": sample_lambda_star if sample_lambda_star is not None else float("nan"),
        "sample_delta_theta_rad": sample_delta_theta if sample_delta_theta is not None else float("nan"),
    }


def coupling_profile(n: int, left: float, right: float, profile: str = "linear") -> Array:
    if n < 2:
        raise ValueError("n must be >= 2")
    m = n - 1
    if profile == "linear":
        return np.linspace(left, right, m)
    if profile == "step":
        k = m // 2
        return np.array([left] * k + [right] * (m - k), dtype=float)
    if profile == "weak_center":
        arr = np.full(m, left, dtype=float)
        arr[m // 2] = right
        return arr
    raise ValueError(f"Unknown coupling profile: {profile}")


def finite_size_scan(
    outdir: str,
    ns: Sequence[int],
    h: float,
    j_scale_left: float,
    j_scale_right: float,
    alpha: float,
    beta: float,
    state_mode: str,
    evolve_time: float,
    coupling_profile_name: str = "weak_center",
) -> Dict[int, float]:
    del alpha, beta, state_mode, evolve_time
    lambdas_star = {}
    for n in ns:
        couplings = coupling_profile(n, j_scale_left, j_scale_right, profile=coupling_profile_name)
        h_total, bond_terms = xy_chain_hamiltonian(n, couplings, h)
        _, psi = ground_state(h_total)
        lam_star = first_crossover_lambda(psi, n, bond_terms)
        lambdas_star[n] = np.nan if lam_star is None else lam_star

    ns_arr = np.array(list(lambdas_star.keys()), dtype=float)
    lam_arr = np.array(list(lambdas_star.values()), dtype=float)
    mask = ~np.isnan(lam_arr)
    if np.any(mask):
        save_scaling_plot(os.path.join(outdir, "REM4_scaling.png"), ns_arr[mask], lam_arr[mask])
    return lambdas_star


def print_metrics_table(metrics: Sequence[CutMetrics]) -> None:
    print("cut\tPhi_S(bits)\tPhi_H\tPhi\t<E_boundary>")
    for metric in metrics:
        print(
            f"{metric.cut}\t{metric.phi_s:.6f}\t{metric.phi_h:.6f}\t"
            f"{metric.phi:.6f}\t{metric.boundary_energy:.6f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="REM4 exact diagonalization scaffold")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--mode", choices=["three", "scaling", "all"], default="three")
    parser.add_argument("--state-mode", choices=["ground", "initial", "evolved"], default="ground")
    parser.add_argument("--alpha", type=float, default=0.85)
    parser.add_argument("--beta", type=float, default=0.35)
    parser.add_argument("--time", type=float, default=0.8)
    parser.add_argument("--j12", type=float, default=1.2)
    parser.add_argument("--j23", type=float, default=1.0)
    parser.add_argument("--h", type=float, default=0.2)
    parser.add_argument("--lambda", dest="lambda_value", type=float, default=0.2)
    parser.add_argument("--lambda-max", type=float, default=1.0)
    parser.add_argument("--grid-r", type=int, default=80)
    parser.add_argument("--grid-l", type=int, default=80)
    parser.add_argument("--coupling-profile", choices=["linear", "step", "weak_center"], default="weak_center")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    h_total, bond_terms = xy_chain_hamiltonian(3, [args.j12, args.j23], args.h)
    if args.state_mode == "ground":
        e0, psi = ground_state(h_total)
        print(f"Ground-state energy: {e0:.8f}")
    elif args.state_mode == "initial":
        psi = initial_state_family(3, alpha=args.alpha, beta=args.beta)
    else:
        psi0 = initial_state_family(3, alpha=args.alpha, beta=args.beta)
        psi = time_evolved_state(h_total, psi0, args.time)

    metrics = evaluate_cuts(psi, 3, bond_terms, lambda_value=args.lambda_value)
    print_metrics_table(metrics)
    info_cut = info_selected_cut(metrics)
    std_cut = std_selected_cut(metrics)
    rem_cut = rem_selected_cut(metrics)
    print(f"info-selected cut: {info_cut}")
    print(f"std-selected cut:  {std_cut}")
    print(f"rem-selected cut:  {rem_cut}")

    if rem_cut != std_cut:
        angle = dominant_schmidt_angle(psi, 3, std_cut, rem_cut)
        print(f"pointer mismatch angle proxy Δθ = {angle:.6f} rad ({math.degrees(angle):.3f} deg)")
    else:
        print("pointer mismatch angle proxy Δθ = 0.000000 rad (cuts coincide)")

    p_final, _ = gradient_flow_on_cut_simplex(metrics)
    print(f"gradient-flow relaxed cut weights = {p_final}")

    if args.mode in ("three", "all"):
        result = make_three_qubit_phase_scan(
            args.outdir,
            args.j23,
            args.h,
            args.lambda_max,
            args.grid_r,
            args.grid_l,
            args.alpha,
            args.beta,
            args.state_mode,
            args.time,
        )
        print("three-qubit scan:", result)

    if args.mode in ("scaling", "all"):
        scaling = finite_size_scan(
            args.outdir,
            [3, 4, 6, 8],
            args.h,
            1.5,
            0.6,
            args.alpha,
            args.beta,
            args.state_mode,
            args.time,
            args.coupling_profile,
        )
        print("finite-size lambda* estimates:", scaling)


if __name__ == "__main__":
    main()
