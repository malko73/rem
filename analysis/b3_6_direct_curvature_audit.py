#!/usr/bin/env python3
"""
B3.6: Direct Curvature Audit — reconstruct Hessian from directional second derivatives.

Purpose: Verify Gate 1 condition 3 by reconstructing the Hessian from direct
directional measurements, bypassing the problematic finite-difference Hessian
implementation that produced spurious positive eigenvalues in B3.4.

Method:
  - For each of 45 horizontal basis directions e_i, compute D²_{e_i}Φ via central difference
  - For each of 990 pair directions v_{ij} = (e_i + e_j)/√2, compute D²_{v_{ij}}Φ
  - Reconstruct H_ii = D²_{e_i}Φ
  - Reconstruct H_ij = D²_{v_{ij}}Φ - (H_ii + H_jj)/2
  - Compute eigenvalues of H_direct
  - If λ_max(H_direct) ≤ ε, then local/flat maximum confirmed

This avoids the finite-difference Hessian entirely.
"""
import numpy as np
import json
import sys
from pathlib import Path
import scipy.linalg as sla

# Load modules
REPO = Path(__file__).parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "analysis"))

from rem4_numerical import xy_chain_hamiltonian, ground_state
from quotient_geometry import horizontal_basis
from b2_canonical_open import canonical_open_objective, apply_tps, liouvillian_dephasing

# Setup
N, N_A = 3, 2
D_A, D_B = 4, 2
DIM = 8
LAMBDA = 0.2

h_total, _ = xy_chain_hamiltonian(N, [1.5, 0.6], 0.2)
_, psi0 = ground_state(h_total)
H_basis = horizontal_basis(D_A, D_B)

# Load A1 converged U* from B3.4
conv_path = REPO / "analysis_output" / "b3_4_gate1_convergence.json"
with open(conv_path) as f:
    conv4 = json.load(f)["results"]

U = np.array(conv4["A1_deph_0512"]["U_real"]) + 1j * np.array(conv4["A1_deph_0512"]["U_imag"])
L = liouvillian_dephasing(h_total, [0.5, 1.0, 2.0])

def phi(U):
    rho_U, Y_U = apply_tps(psi0, L, U)
    obj = canonical_open_objective(rho_U, Y_U, N_A)
    return obj["Phi"]

phi_star = phi(U)
print(f"A1: Φ* = {phi_star:.8f}")

# Direct curvature measurement
EPS = 1e-3
n_dim = 45

print(f"\nMeasuring D²_{{e_i}}Φ for {n_dim} diagonal directions (ε={EPS})...")
D2_diag = np.zeros(n_dim)
for i in range(n_dim):
    # Direction e_i
    H_i = H_basis[:, :, i]
    U_p = U @ sla.expm(EPS * H_i)
    U_m = U @ sla.expm(-EPS * H_i)
    phi_p = phi(U_p)
    phi_m = phi(U_m)
    D2_diag[i] = (phi_p - 2*phi_star + phi_m) / EPS**2

print(f"  Diagonal curvatures: min={D2_diag.min():.4e}, max={D2_diag.max():.4e}")

print(f"\nMeasuring D²_{{v_ij}}Φ for {n_dim*(n_dim-1)//2} pair directions...")
D2_pair = np.zeros((n_dim, n_dim))
count = 0
for i in range(n_dim):
    for j in range(i+1, n_dim):
        # Direction v_{ij} = (e_i + e_j)/√2
        H_i = H_basis[:, :, i]
        H_j = H_basis[:, :, j]
        H_ij = (H_i + H_j) / np.sqrt(2)
        U_p = U @ sla.expm(EPS * H_ij)
        U_m = U @ sla.expm(-EPS * H_ij)
        phi_p = phi(U_p)
        phi_m = phi(U_m)
        D2_pair[i, j] = (phi_p - 2*phi_star + phi_m) / EPS**2
        D2_pair[j, i] = D2_pair[i, j]
        count += 1
        if count % 100 == 0:
            print(f"  Progress: {count}/{n_dim*(n_dim-1)//2}")

print(f"  Pair curvatures: min={D2_pair.min():.4e}, max={D2_pair.max():.4e}")

# Reconstruct Hessian
print("\nReconstructing H_direct...")
H_direct = np.zeros((n_dim, n_dim))

# Diagonal: H_ii = D²_{e_i}Φ
for i in range(n_dim):
    H_direct[i, i] = D2_diag[i]

# Off-diagonal: H_ij = D²_{v_{ij}}Φ - (H_ii + H_jj)/2
for i in range(n_dim):
    for j in range(i+1, n_dim):
        H_direct[i, j] = D2_pair[i, j] - (D2_diag[i] + D2_diag[j]) / 2
        H_direct[j, i] = H_direct[i, j]

# Compute eigenvalues
eigs = np.linalg.eigvalsh(H_direct)
lambda_max = eigs.max()
lambda_min = eigs.min()

print(f"\nH_direct eigenvalue spectrum:")
print(f"  λ_max = {lambda_max:.6e}")
print(f"  λ_min = {lambda_min:.6e}")
print(f"  n_positive (>1e-3) = {np.sum(eigs > 1e-3)}")
print(f"  n_zero (|λ|≤1e-3)  = {np.sum(np.abs(eigs) <= 1e-3)}")
print(f"  n_negative (<-1e-3) = {np.sum(eigs < -1e-3)}")

# Verdict
if lambda_max <= 1e-3:
    verdict = "PASS: local/flat maximum confirmed"
else:
    verdict = f"FAIL: λ_max = {lambda_max:.4e} > 1e-3"

print(f"\nVerdict: {verdict}")

# Save results
out_path = REPO / "analysis_output" / "b3_6_direct_curvature_audit.json"
results = {
    "environment": "A1_deph_0512",
    "phi_star": phi_star,
    "epsilon": EPS,
    "H_direct": H_direct.tolist(),
    "eigenvalues": eigs.tolist(),
    "lambda_max": float(lambda_max),
    "lambda_min": float(lambda_min),
    "n_positive": int(np.sum(eigs > 1e-3)),
    "n_zero": int(np.sum(np.abs(eigs) <= 1e-3)),
    "n_negative": int(np.sum(eigs < -1e-3)),
    "verdict": verdict
}

with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(f"\nSaved to {out_path}")
