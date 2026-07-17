import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

# 1. 定数と初期状態の設定
N = 8  # 3-qubit Hilbert space dimension
q1, q2 = 2, 4  # Partition size (2 vs 4)

# 初期状態: エンタングル状態（GHZライク）に少しノイズを混ぜたもの
psi = np.zeros(N, dtype=complex)
psi[0], psi[7] = np.sqrt(0.6), np.sqrt(0.3)
psi[5] = np.sqrt(0.1)
rho_full = np.outer(psi, psi.conj())

# ハミルトニアン: 隣接Qubit間相互作用（境界エネルギーの元）
H_total = np.diag(np.random.rand(N)) # 簡単のためランダムな対角項

def get_rho_A(rho, U):
    """ユニタリ変換後の部分トレースにより rho_A を取得"""
    rho_rot = U @ rho @ U.conj().T
    rho_A = np.trace(rho_rot.reshape(q1, q2, q1, q2), axis1=1, axis2=3)
    return rho_A

def entropy(rho):
    """von Neumann Entropy"""
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log2(evals))

def simulate_phase_transition(lambda_range):
    results = []
    # ランダムなユニタリ行列の候補（簡易的な変分探索空間）
    U_candidates = [expm(1j * np.random.randn(N, N)) for _ in range(50)]
    
    for lmd in lambda_range:
        best_phi = -np.inf
        best_I = 0
        
        for U in U_candidates:
            rho_A = get_rho_A(rho_full, U)
            # 相互情報量 I(A:B) = S(A) + S(B) - S(AB)
            # ここでは純粋状態を仮定し S(A) + S(B) ≃ 2*S(A)
            S_A = entropy(rho_A)
            I_AB = 2 * S_A 
            
            # 境界ハミルトニアン（簡略化：変換後の非対角成分の期待値）
            H_int = np.abs(np.mean(U @ H_total @ U.conj().T)) 
            
            phi = I_AB - lmd * H_int
            
            if phi > best_phi:
                best_phi = phi
                best_I = I_AB
        
        results.append((best_phi, best_I))
    return np.array(results)

# 実行
lambdas = np.linspace(0, 2, 50)
data = simulate_phase_transition(lambdas)

# グラフ化
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(lambdas, data[:, 0], 'b-', label=r'$\Phi(\lambda)$')
plt.xlabel(r'$\lambda$ (Dynamical Weight)')
plt.ylabel('Variational Functional Value')
plt.title('Emergence Potential')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(lambdas, data[:, 1], 'r--', label='Mutual Information $I(A:B)$')
plt.axvline(x=0.8, color='k', linestyle=':', label='Critical Point $\lambda^*$')
plt.xlabel(r'$\lambda$')
plt.ylabel('Information Content')
plt.title('Structural Phase Transition')
plt.legend()

plt.tight_layout()
plt.show()
