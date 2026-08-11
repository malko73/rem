#!/usr/bin/env python3
"""
B3.5c: Directional second derivative verification for uniform and A2.

Purpose: Verify whether the Hessian's positive eigenvalue is real or a numerical artifact.
Method: Compute D_v²Φ directly along v_max direction and compare with Hessian prediction.
"""
import numpy as np
import json
import sys
import importlib.util
from pathlib import Path
import scipy.linalg as sla

# Load modules
REPO = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location("rem4_numerical", REPO / "src" / "rem4_numerical.py")
rem4 = importlib.util.module_from_spec(spec)
sys.modules["rem4_numerical"] = rem4
spec.loader.exec_module(rem4)

spec2 = importlib.util.spec_from_file_location("quotient_geometry", REPO / "src" / "quotient_geometry.py")
qgeom = importlib.util.module_from_spec(spec2)
sys.modules["quotient_geometry"] = qgeom
spec2.loader.exec_module(qgeom)

spec3 = importlib.util.spec_from_file_location("gamma_F", REPO / "analysis" / "gamma_F.py")
gamma_F = importlib.util.module_from_spec(spec3)
sys.modules["gamma_F"] = gamma_F
spec3.loader.exec_module(gamma_F)

spec4 = importlib.util.spec_from_file_location("b2_canonical_open", REPO / "analysis" / "b2_canonical_open.py")
b2 = importlib.util.module_from_spec(spec4)
sys.modules["b2_canonical_open"] = b2
spec4.loader.exec_module(b2)

# Setup
h_total, _ = rem4.xy_chain_hamiltonian(3, [1.5, 0.6], 0.2)
_, psi0 = rem4.ground_state(h_total)
H_basis = qgeom.horizontal_basis(4, 2)

def phi(U, L):
    rho_U, Y_U = b2.apply_tps(psi0, L, U)
    return b2.canonical_open_objective(rho_U, Y_U, 2)['Phi']

# Load converged U* from B3.4
conv4 = json.loads((REPO / "analysis_output" / "b3_4_gate1_convergence.json").read_text())['results']

envs = {
    'uniform_dephasing': b2.liouvillian_dephasing(h_total, [1.0, 1.0, 1.0]),
    'A2_deph_2105':      b2.liouvillian_dephasing(h_total, [2.0, 1.0, 0.5]),
}

results = {}

for env_name, L in envs.items():
    print(f'\n===== {env_name} =====')
    U = np.array(conv4[env_name]['U_real']) + 1j*np.array(conv4[env_name]['U_imag'])
    phi_star = phi(U, L)
    print(f'Φ* = {phi_star:.8f}')

    # Hessian (ε=3e-2)
    n = 45
    H_EPS = 3e-2
    H = np.zeros((n, n))
    f0 = phi(U, L)
    for i in range(n):
        Hi = H_basis[:, :, i]
        H[i,i] = (phi(U @ sla.expm(H_EPS*Hi), L) - 2*f0 + phi(U @ sla.expm(-H_EPS*Hi), L)) / H_EPS**2
    for i in range(n):
        Hi = H_basis[:, :, i]
        for j in range(i+1, n):
            Hj = H_basis[:, :, j]
            H[i,j] = (phi(U@sla.expm(H_EPS*Hi)@sla.expm(H_EPS*Hj), L) - phi(U@sla.expm(H_EPS*Hi)@sla.expm(-H_EPS*Hj), L)
                      - phi(U@sla.expm(-H_EPS*Hi)@sla.expm(H_EPS*Hj), L) + phi(U@sla.expm(-H_EPS*Hi)@sla.expm(-H_EPS*Hj), L)) / (4*H_EPS**2)
            H[j,i] = H[i,j]

    eigs, vecs = np.linalg.eigh(H)
    idx = np.argmax(eigs)
    v_max = vecs[:, idx]
    lam_max = eigs[idx]
    print(f'λ_max(H) = {lam_max:.6f}')

    # Directional second derivative
    K_v = np.einsum('ijk,k->ij', H_basis, v_max)
    print(f'ε        | ΔΦ+      | ΔΦ-      | D_v²Φ     | Hessian予測')
    d2_results = []
    for eps in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4]:
        U_p = U @ sla.expm(eps * K_v)
        U_m = U @ sla.expm(-eps * K_v)
        dp = phi(U_p, L) - phi(U, L)
        dm = phi(U_m, L) - phi(U, L)
        d2 = (dp + dm) / eps**2
        pred = 0.5 * lam_max * eps**2
        print(f'{eps:.0e} | {dp:+.4e} | {dm:+.4e} | {d2:+.4e} | {pred:+.4e}')
        d2_results.append({'eps': eps, 'd2': d2, 'pred': pred})

    results[env_name] = {
        'phi_star': phi_star,
        'lambda_max_hessian': lam_max,
        'directional_d2': d2_results
    }

# Save results
out_path = REPO / "analysis_output" / "b3_5c_directional_verification.json"
out_path.write_text(json.dumps(results, indent=2))
print(f'\nSaved to {out_path}')
