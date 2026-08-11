"""
Phase B0 — Geometry-layer unit tests (master's gate list).

Locks the SVD-based quotient decomposition:
  dim u(8)       = 64
  raw vertical   = 20
  rank vertical  = 19
  dim horizontal = 45
  V†H            ≈ 0
  P_V²           = P_V
  P_H²           = P_H
  P_V + P_H      = I
"""
import numpy as np
import pytest

from quotient_geometry import (
    u_basis,
    raw_vertical_generators,
    vertical_basis,
    horizontal_basis,
    project_vertical,
    project_horizontal,
    geometry_report,
)


def test_dim_u8_is_64():
    basis = u_basis(8)
    assert basis.shape == (8, 8, 64)


def test_raw_vertical_is_20():
    raw = raw_vertical_generators(4, 2)
    assert raw.shape == (8, 8, 20)


def test_rank_vertical_is_19():
    V = vertical_basis(4, 2)
    assert V.shape == (8, 8, 19)


def test_dim_horizontal_is_45():
    H = horizontal_basis(4, 2)
    assert H.shape == (8, 8, 45)


def test_VdagH_zero():
    V = vertical_basis(4, 2)
    H = horizontal_basis(4, 2)
    maxv = 0.0
    for i in range(V.shape[2]):
        for j in range(H.shape[2]):
            val = abs(np.real(np.trace(V[:, :, i].conj().T @ H[:, :, j])))
            maxv = max(maxv, val)
    assert maxv < 1e-10, f"V†H max = {maxv}"


def test_projectors_idempotent_and_complete():
    rep = geometry_report(4, 2)
    assert rep["PV2_minus_PV"] < 1e-10, f"P_V²≠P_V: {rep['PV2_minus_PV']}"
    assert rep["PH2_minus_PH"] < 1e-10, f"P_H²≠P_H: {rep['PH2_minus_PH']}"
    assert rep["PV_plus_PH_minus_I"] < 1e-10, f"P_V+P_H≠I: {rep['PV_plus_PH_minus_I']}"


def test_bases_orthonormal():
    rep = geometry_report(4, 2)
    assert rep["V_orthonorm"] < 1e-10, f"V not orthonormal: {rep['V_orthonorm']}"
    assert rep["H_orthonorm"] < 1e-10, f"H not orthonormal: {rep['H_orthonorm']}"


def test_projections_agree_with_basis_expansion():
    rng = np.random.default_rng(7)
    K = 1j * rng.standard_normal((8, 8))
    K = K - K.conj().T  # anti-Hermitian

    Pv = project_vertical(K, 4, 2)
    Ph = project_horizontal(K, 4, 2)
    # Decomposition is exact: K = Pv + Ph
    assert np.max(np.abs((Pv + Ph) - K)) < 1e-10
