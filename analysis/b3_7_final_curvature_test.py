#!/usr/bin/env python3
"""
B3.7: Final curvature test — direct measurement along v_+ from H_direct.

Purpose: Verify whether the positive eigenvalue from B3.6's reconstructed Hessian
is real or an artifact by directly measuring the curvature along its eigenvector.

Method:
  - Load H_direct from B3.6
  - Compute eigenvalues/eigenvectors
  - Extract v_+ (eigenvector for λ_max = 2.70)
  - For ε = 3e-2, 1e-2, 3e-3, 1e-3, 3e-4:
      U±(ε) = exp(±ε K_{v_+}) U*
      D²_{v_+}Φ = (Φ(U+) - 2Φ(U*) + Φ(U-)) / ε²
  - If D²_{v_+}Φ → +2.7: true saddle (Gate 1 FAIL)
  - If D²_{v_+}Φ ≤ 0: Hessian reconstruction error (need alternative method)
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

# Load H_direct from B3.6
b36_path = REPO / "analysis_output" / "b3_6_direct_curvature_audit.json"
with open(b36_path) as f:
    b36 = json.load(f)

H_direct = np.array(b36["H_direct"])
lambda_max_b36 = b36["lambda_max"]
print(f"\nB3.6 H_direct: λ_max = {lambda_max_b36:.6f}")

# Eigendecomposition
eigs, vecs = np.linalg.eigh(H_direct)
idx_max = np.argmax(eigs)
v_plus = vecs[:, idx_max]
lambda_max = eigs[idx_max]

print(f"Recomputed: λ_max = {lambda_max:.6f} (eigenvector index {idx_max})")
print(f"|v_+| = {np.linalg.norm(v_plus):.6f} (should be 1.0)")

# Construct K_{v_+} = Σ_k v_+[k] H_k
K_vplus = np.einsum('ijk,k->ij', H_basis, v_plus)
print(f"|K_{{v_+}}|_F = {np.linalg.norm(K_vplus, 'fro'):.6f}")

# Direct curvature measurement
print(f"\nDirect measurement: D²_{{v_+}}Φ(ε)")
print(f"ε        | Φ(U+)    | Φ(U-)    | ΔΦ+      | ΔΦ-      | D²_{{v_+}}Φ | 予測 ½λε²")

results = []
for eps in [3e-2, 1e-2, 3e-3, 1e-3, 3e-4]:
    U_p = U @ sla.expm(eps * K_vplus)
    U_m = U @ sla.expm(-eps * K_vplus)
    phi_p = phi(U_p)
    phi_m = phi(U_m)
    dp = phi_p - phi_star
    dm = phi_m - phi_star
    d2 = (dp + dm) / eps**2
    pred = 0.5 * lambda_max * eps**2
    print(f"{eps:.0e} | {phi_p:.8f} | {phi_m:.8f} | {dp:+.4e} | {dm:+.4e} | {d2:+.6e} | {pred:+.6e}")
    results.append({
        "epsilon": eps,
        "phi_plus": float(phi_p),
        "phi_minus": float(phi_m),
        "delta_phi_plus": float(dp),
        "delta_phi_minus": float(dm),
        "d2_curvature": float(d2),
        "prediction": float(pred)
    })

# Verdict
d2_small_eps = [r["d2_curvature"] for r in results if r["epsilon"] <= 1e-3]
d2_mean = np.mean(d2_small_eps)

print(f"\nMean D²_{{v_+}}Φ for ε≤1e-3: {d2_mean:+.6e}")
print(f"Expected if true saddle: +{lambda_max:.6f}")

if abs(d2_mean - lambda_max) / lambda_max < 0.2:  # within 20%
    verdict = "CASE A: True saddle confirmed (Gate 1 FAIL)"
else:
    verdict = "CASE B: Hessian reconstruction error (need alternative method)"

print(f"\nVerdict: {verdict}")

# Save results
out_path = REPO / "analysis_output" / "b3_7_final_curvature_test.json"
output = {
    "environment": "A1_deph_0512",
    "phi_star": float(phi_star),
    "lambda_max_b36": float(lambda_max_b36),
    "lambda_max_recomputed": float(lambda_max),
    "v_plus_norm": float(np.linalg.norm(v_plus)),
    "K_vplus_norm": float(np.linalg.norm(K_vplus, 'fro')),
    "measurements": results,
    "d2_mean_small_eps": float(d2_mean),
    "verdict": verdict
}

with open(out_path, "w") as f:
    json.dump(output, f, indent=2)

print(f"\nSaved to {out_path}")
