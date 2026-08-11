"""
Phase B0: Quotient geometry for the factorization manifold.

Provides the geometric decomposition of the TPS quotient space
  F = U(D) / Im(U(d_A) × U(d_B))
via a real SVD of the coefficient matrix of raw local generators.

Key construction (per master's decision 2026-08-11):
  1. Build an orthonormal REAL basis of u(D) (dimension D²) using the
     real inner product  <X, Y> = Re Tr(X† Y).
  2. Build the RAW local/gauge generators (d_A² + d_B² of them), WITHOUT
     manually removing the U(1) kernel:
         A_a ⊗ I_{d_B}  (a = 1..d_A²),   I_{d_A} ⊗ B_b  (b = 1..d_B²)
  3. Map each raw generator to a D²-dimensional real coefficient vector.
  4. SVD: M = U Σ V^T,  M ∈ R^{D² × (d_A² + d_B²)}.
     The rank of M is automatically d_A² + d_B² - 1 (the U(1) kernel is
     detected as a linear dependence, NOT removed by hand).
  5. vertical basis   = U[:, 0:rank]
     horizontal basis = U[:, rank:D²]
     which guarantees V^T H = 0 and dim V = rank, dim H = D² - rank.

Everything here is geometry-only: no optimizer, no Hamiltonian evaluation.
"""
from __future__ import annotations

import numpy as np


def u_basis(d: int) -> np.ndarray:
    """
    Orthonormal REAL basis of u(d) (anti-Hermitian matrices).

    Returns
    -------
    basis : ndarray, shape (d, d, d²), complex128
        basis[:, :, k] is the k-th generator, orthonormal w.r.t.
        <X, Y> = Re Tr(X† Y).

    Construction
    ------------
    For each ordered pair of indices (p, q) we take two real directions:
      * symmetric real:  (E_{pq} + E_{qp}) / sqrt(2)   for p <= q
      * antisymmetric imaginary: i (E_{pq} - E_{qp}) / sqrt(2)  for p < q
      * diagonal imaginary: i E_{pp}
    This yields d² real directions, all anti-Hermitian.
    """
    basis = np.zeros((d, d, d * d), dtype=complex)
    k = 0
    # Diagonal imaginary: i E_{pp}
    for p in range(d):
        basis[p, p, k] = 1j
        k += 1
    # Off-diagonal real ANTIsymmetric (anti-Hermitian): (E_{pq} - E_{qp})/sqrt2
    for p in range(d):
        for q in range(p + 1, d):
            basis[p, q, k] = 1.0 / np.sqrt(2.0)
            basis[q, p, k] = -1.0 / np.sqrt(2.0)
            k += 1
    # Off-diagonal imaginary SYMMETRIC (anti-Hermitian): i (E_{pq} + E_{qp})/sqrt2
    for p in range(d):
        for q in range(p + 1, d):
            basis[p, q, k] = 1j / np.sqrt(2.0)
            basis[q, p, k] = 1j / np.sqrt(2.0)
            k += 1
    assert k == d * d, f"u_basis produced {k} generators, expected {d*d}"
    return basis


def _coeff_vector(X: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """
    Expand an anti-Hermitian matrix X in the orthonormal basis.

    Returns the real coefficient vector c with X ≈ Σ_k c_k basis[:,:,k].
    c_k = <basis_k, X> = Re Tr(basis_k† X) = -Re Tr(basis_k X)
    (since basis_k is anti-Hermitian, basis_k† = -basis_k).
    """
    d = basis.shape[0]
    c = np.zeros(basis.shape[2])
    Xc = X.astype(complex)
    for k in range(basis.shape[2]):
        c[k] = np.real(np.trace(basis[:, :, k].conj().T @ Xc))
    return c


def raw_vertical_generators(d_a: int, d_b: int) -> np.ndarray:
    """
    RAW local/gauge generators: {A_a ⊗ I_{d_B}} ∪ {I_{d_A} ⊗ B_b}.

    The U(1) kernel (i·I_{d_A·d_B}) is intentionally NOT removed here;
    it appears as a linear dependence among the d_A² + d_B² generators.

    Returns
    -------
    gens : ndarray, shape (D, D, d_A² + d_B²), complex128
    """
    D = d_a * d_b
    u_a = u_basis(d_a)  # (d_A, d_A, d_A²)
    u_b = u_basis(d_b)  # (d_B, d_B, d_B²)
    gens = np.zeros((D, D, d_a * d_a + d_b * d_b), dtype=complex)
    k = 0
    for a in range(d_a * d_a):
        gens[:, :, k] = np.kron(u_a[:, :, a], np.eye(d_b, dtype=complex))
        k += 1
    for b in range(d_b * d_b):
        gens[:, :, k] = np.kron(np.eye(d_a, dtype=complex), u_b[:, :, b])
        k += 1
    return gens


def vertical_basis(d_a: int, d_b: int) -> np.ndarray:
    """
    Orthonormal vertical basis of the gauge subspace.

    Returns
    -------
    basis : ndarray, shape (D, D, rank), complex128
        rank = d_A² + d_B² - 1  (19 for d_A=4, d_B=2).
        Derived from SVD column space of the raw generators' coefficient matrix.
    """
    D = d_a * d_b
    u_full = u_basis(D)
    raw = raw_vertical_generators(d_a, d_b)
    n_raw = raw.shape[2]

    M = np.zeros((D * D, n_raw))
    for r in range(n_raw):
        M[:, r] = _coeff_vector(raw[:, :, r], u_full)

    U, S, _ = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(S > 1e-10))
    # Vertical basis = left singular vectors spanning M's column space,
    # expressed in the u(D) basis: V_col = u_full_flat @ U[:, :rank]
    V = (u_full.reshape(D * D, D * D) @ U[:, :rank]).reshape(D, D, rank)
    return V


def horizontal_basis(d_a: int, d_b: int) -> np.ndarray:
    """
    Orthonormal horizontal basis (orthogonal complement of vertical).

    Returns
    -------
    basis : ndarray, shape (D, D, D² - rank), complex128
        For d_A=4, d_B=2: 64 - 19 = 45.
    """
    D = d_a * d_b
    u_full = u_basis(D)
    raw = raw_vertical_generators(d_a, d_b)
    n_raw = raw.shape[2]

    M = np.zeros((D * D, n_raw))
    for r in range(n_raw):
        M[:, r] = _coeff_vector(raw[:, :, r], u_full)

    U, S, _ = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(S > 1e-10))
    # Horizontal basis = remaining left singular vectors,
    # expressed in the u(D) basis: H_col = u_full_flat @ U[:, rank:]
    H = (u_full.reshape(D * D, D * D) @ U[:, rank:]).reshape(D, D, D * D - rank)
    return H


def project_vertical(K: np.ndarray, d_a: int, d_b: int) -> np.ndarray:
    """Project anti-Hermitian K onto the vertical (gauge) subspace."""
    V = vertical_basis(d_a, d_b)
    D = d_a * d_b
    out = np.zeros((D, D), dtype=complex)
    Kc = K.astype(complex)
    for k in range(V.shape[2]):
        c = np.real(np.trace(V[:, :, k].conj().T @ Kc))
        out += c * V[:, :, k]
    return out


def project_horizontal(K: np.ndarray, d_a: int, d_b: int) -> np.ndarray:
    """Project anti-Hermitian K onto the horizontal subspace."""
    H = horizontal_basis(d_a, d_b)
    D = d_a * d_b
    out = np.zeros((D, D), dtype=complex)
    Kc = K.astype(complex)
    for k in range(H.shape[2]):
        c = np.real(np.trace(H[:, :, k].conj().T @ Kc))
        out += c * H[:, :, k]
    return out


def geometry_report(d_a: int = 4, d_b: int = 2) -> dict:
    """Return the geometry invariants as a dict (for tests and diagnostics)."""
    D = d_a * d_b
    u_full = u_basis(D)
    raw = raw_vertical_generators(d_a, d_b)
    V = vertical_basis(d_a, d_b)
    H = horizontal_basis(d_a, d_b)

    # Projectors
    P_V = np.einsum("abk,cdk->abcd", V, V.conj()).reshape(D * D, D * D)
    P_H = np.einsum("abk,cdk->abcd", H, H.conj()).reshape(D * D, D * D)
    Id = np.eye(D * D, dtype=complex)

    report = {
        "dim_u": u_full.shape[2],                 # D²
        "raw_vertical": raw.shape[2],             # d_A² + d_B²
        "rank_vertical": V.shape[2],              # d_A² + d_B² - 1
        "dim_horizontal": H.shape[2],             # D² - rank
        "VdagH_max": _max_abs_VdagH(V, H),
        "PV2_minus_PV": float(np.max(np.abs(P_V @ P_V - P_V))),
        "PH2_minus_PH": float(np.max(np.abs(P_H @ P_H - P_H))),
        "PV_plus_PH_minus_I": float(np.max(np.abs(P_V + P_H - Id))),
        "V_orthonorm": _orthonormality_error(V),
        "H_orthonorm": _orthonormality_error(H),
    }
    return report


def _max_abs_VdagH(V: np.ndarray, H: np.ndarray) -> float:
    """max_{i,j} |<V_i, H_j>| where <A,B> = Re Tr(A† B)."""
    maxv = 0.0
    for i in range(V.shape[2]):
        for j in range(H.shape[2]):
            val = abs(np.real(np.trace(V[:, :, i].conj().T @ H[:, :, j])))
            if val > maxv:
                maxv = val
    return float(maxv)


def _orthonormality_error(B: np.ndarray) -> float:
    """max |<B_i, B_j> - δ_ij|."""
    maxv = 0.0
    for i in range(B.shape[2]):
        for j in range(B.shape[2]):
            val = abs(np.real(np.trace(B[:, :, i].conj().T @ B[:, :, j])) - (1.0 if i == j else 0.0))
            if val > maxv:
                maxv = val
    return float(maxv)
