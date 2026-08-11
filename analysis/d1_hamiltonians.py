#!/usr/bin/env python3
"""
Phase D1: Hamiltonian Families Implementation

Implement the 5 Hamiltonian families defined in D0 protocol:
1. Asymmetric XY Chain (existing benchmark)
2. Transverse-Field Ising Model
3. Heisenberg XXZ Model
4. XYZ Model
5. Random Local Hamiltonian
"""

import numpy as np
from typing import Dict, List, Tuple

def pauli_matrices() -> Dict[str, np.ndarray]:
    """Return Pauli matrices."""
    return {
        'I': np.eye(2, dtype=complex),
        'X': np.array([[0, 1], [1, 0]], dtype=complex),
        'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
        'Z': np.array([[1, 0], [0, -1]], dtype=complex)
    }

def tensor_product(operators: List[np.ndarray]) -> np.ndarray:
    """Compute tensor product of a list of operators."""
    result = operators[0]
    for op in operators[1:]:
        result = np.kron(result, op)
    return result

def build_operator_term(pauli_string: str, coeff: complex, n_sites: int) -> np.ndarray:
    """
    Build a term in the Hamiltonian from a Pauli string.
    
    Args:
        pauli_string: String like 'XIZ' or 'ZZI'
        coeff: Coefficient
        n_sites: Total number of sites
    
    Returns:
        Hamiltonian term as 2^n_sites x 2^n_sites matrix
    """
    paulis = pauli_matrices()
    operators = []
    
    for char in pauli_string:
        if char in paulis:
            operators.append(paulis[char])
        else:
            raise ValueError(f"Invalid Pauli character: {char}")
    
    # Pad with identity if needed
    while len(operators) < n_sites:
        operators.append(paulis['I'])
    
    return coeff * tensor_product(operators)

def asymmetric_xy_chain(n_sites: int = 3, J_12: float = 1.5, J_23: float = 0.6, h: float = 0.2) -> np.ndarray:
    """
    Asymmetric XY Chain: H = J₁₂(X₁X₂ + Y₁Y₂) + J₂₃(X₂X₃ + Y₂Y₃) + h(Z₁ + Z₂ + Z₃)
    
    This is the existing benchmark from Phase B/C.
    """
    H = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    
    # J₁₂(X₁X₂ + Y₁Y₂)
    H += J_12 * build_operator_term('XXI', 1.0, n_sites)
    H += J_12 * build_operator_term('YYI', 1.0, n_sites)
    
    # J₂₃(X₂X₃ + Y₂Y₃)
    H += J_23 * build_operator_term('IXX', 1.0, n_sites)
    H += J_23 * build_operator_term('IYY', 1.0, n_sites)
    
    # h(Z₁ + Z₂ + Z₃)
    H += h * build_operator_term('ZII', 1.0, n_sites)
    H += h * build_operator_term('IZI', 1.0, n_sites)
    H += h * build_operator_term('IIZ', 1.0, n_sites)
    
    return H

def transverse_field_ising(n_sites: int = 3, J: float = 1.0, h: float = 0.5) -> np.ndarray:
    """
    Transverse-Field Ising Model: H = -J Σ Z_i Z_{i+1} - h Σ X_i
    """
    H = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    
    # -J Σ Z_i Z_{i+1}
    for i in range(n_sites - 1):
        pauli_string = 'I' * i + 'ZZ' + 'I' * (n_sites - i - 2)
        H += -J * build_operator_term(pauli_string, 1.0, n_sites)
    
    # -h Σ X_i
    for i in range(n_sites):
        pauli_string = 'I' * i + 'X' + 'I' * (n_sites - i - 1)
        H += -h * build_operator_term(pauli_string, 1.0, n_sites)
    
    return H

def heisenberg_xxz(n_sites: int = 3, J: float = 1.0, Delta: float = 0.5) -> np.ndarray:
    """
    Heisenberg XXZ Model: H = J Σ (X_i X_{i+1} + Y_i Y_{i+1} + Δ Z_i Z_{i+1})
    """
    H = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    
    for i in range(n_sites - 1):
        # X_i X_{i+1}
        pauli_string = 'I' * i + 'XX' + 'I' * (n_sites - i - 2)
        H += J * build_operator_term(pauli_string, 1.0, n_sites)
        
        # Y_i Y_{i+1}
        pauli_string = 'I' * i + 'YY' + 'I' * (n_sites - i - 2)
        H += J * build_operator_term(pauli_string, 1.0, n_sites)
        
        # Δ Z_i Z_{i+1}
        pauli_string = 'I' * i + 'ZZ' + 'I' * (n_sites - i - 2)
        H += J * Delta * build_operator_term(pauli_string, 1.0, n_sites)
    
    return H

def xyz_model(n_sites: int = 3, J_x: float = 1.0, J_y: float = 0.8, J_z: float = 0.6) -> np.ndarray:
    """
    XYZ Model: H = Σ (J_x X_i X_{i+1} + J_y Y_i Y_{i+1} + J_z Z_i Z_{i+1})
    """
    H = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    
    for i in range(n_sites - 1):
        # J_x X_i X_{i+1}
        pauli_string = 'I' * i + 'XX' + 'I' * (n_sites - i - 2)
        H += J_x * build_operator_term(pauli_string, 1.0, n_sites)
        
        # J_y Y_i Y_{i+1}
        pauli_string = 'I' * i + 'YY' + 'I' * (n_sites - i - 2)
        H += J_y * build_operator_term(pauli_string, 1.0, n_sites)
        
        # J_z Z_i Z_{i+1}
        pauli_string = 'I' * i + 'ZZ' + 'I' * (n_sites - i - 2)
        H += J_z * build_operator_term(pauli_string, 1.0, n_sites)
    
    return H

def random_local_hamiltonian(n_sites: int = 3, seed: int = 42, 
                             h_range: Tuple[float, float] = (-0.5, 0.5),
                             J_range: Tuple[float, float] = (0.3, 0.7)) -> np.ndarray:
    """
    Random Local Hamiltonian: H = Σ h_i · σ_i + Σ J_ij σ_i · σ_j
    
    Uses fixed seed for reproducibility.
    """
    rng = np.random.default_rng(seed)
    H = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    
    paulis = ['X', 'Y', 'Z']
    
    # Single-site terms: h_i · σ_i
    for i in range(n_sites):
        for pauli in paulis:
            h_i = rng.uniform(h_range[0], h_range[1])
            pauli_string = 'I' * i + pauli + 'I' * (n_sites - i - 1)
            H += h_i * build_operator_term(pauli_string, 1.0, n_sites)
    
    # Two-site terms: J_ij σ_i · σ_j
    for i in range(n_sites - 1):
        for pauli in paulis:
            J_ij = rng.uniform(J_range[0], J_range[1])
            pauli_string = 'I' * i + pauli + pauli + 'I' * (n_sites - i - 2)
            H += J_ij * build_operator_term(pauli_string, 1.0, n_sites)
    
    return H

def get_hamiltonian(name: str, n_sites: int = 3, **kwargs) -> np.ndarray:
    """
    Get Hamiltonian by name.
    
    Args:
        name: One of 'asymmetric_XY', 'transverse_field_ising', 'heisenberg_xxz', 'xyz', 'random_local'
        n_sites: Number of sites (default 3)
        **kwargs: Additional parameters for the Hamiltonian
    
    Returns:
        Hamiltonian matrix
    """
    hamiltonians = {
        'asymmetric_XY': asymmetric_xy_chain,
        'transverse_field_ising': transverse_field_ising,
        'heisenberg_xxz': heisenberg_xxz,
        'xyz': xyz_model,
        'random_local': random_local_hamiltonian
    }
    
    if name not in hamiltonians:
        raise ValueError(f"Unknown Hamiltonian: {name}. Choose from {list(hamiltonians.keys())}")
    
    return hamiltonians[name](n_sites=n_sites, **kwargs)

def get_ground_state(H: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Compute ground state energy and state.
    
    Returns:
        (E_ground, psi_ground)
    """
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    E_ground = eigenvalues[0]
    psi_ground = eigenvectors[:, 0]
    return E_ground, psi_ground

if __name__ == "__main__":
    # Test all Hamiltonians
    print("Testing Hamiltonian implementations...")
    print()
    
    for name in ['asymmetric_XY', 'transverse_field_ising', 'heisenberg_xxz', 'xyz', 'random_local']:
        H = get_hamiltonian(name, n_sites=3)
        E_ground, psi_ground = get_ground_state(H)
        
        print(f"{name}:")
        print(f"  Dimension: {H.shape}")
        print(f"  Hermitian: {np.allclose(H, H.conj().T)}")
        print(f"  Ground state energy: {E_ground:.6f}")
        print()
