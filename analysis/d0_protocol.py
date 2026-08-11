#!/usr/bin/env python3
"""
Phase D0: Protocol Freeze

Define and freeze the evaluation protocol for Phase D (generality validation).
This ensures reproducibility and prevents post-hoc adjustments.
"""

import json
from pathlib import Path

# ============================================================================
# Evaluation Metrics
# ============================================================================

METRICS = {
    "Phi_star": {
        "description": "Optimized objective function value",
        "formula": "max_F [I_ρ(F) - λ C_dyn^(τ)(F)]",
        "unit": "dimensionless"
    },
    "d_F": {
        "description": "Gauge-invariant factorization distance",
        "formula": "||P_{A_F1} - P_{A_F2}||_F",
        "unit": "dimensionless",
        "range": "[0, 1]"
    },
    "Gamma_F_exact": {
        "description": "Instantaneous decoherence rate at t=0",
        "formula": "-2 Re<Q_F ρ, Q_F L(ρ)> / ||Q_F ρ||²",
        "unit": "T⁻¹"
    },
    "C_dyn_tau": {
        "description": "Finite-time structure selection functional",
        "formula": "-(1/τ) log(C_F(τ)/C_F(0))",
        "unit": "T⁻¹"
    }
}

# ============================================================================
# Fixed Parameters
# ============================================================================

PARAMETERS = {
    "lambda": {
        "value": 0.2,
        "justification": "Benchmark value from Phase B/C",
        "unit": "T"
    },
    "tau": {
        "value": 0.1,
        "justification": "Used in C3.3b finite-time prediction",
        "unit": "T"
    },
    "optimizer": {
        "value": "Adam",
        "justification": "Consistent with Phase B/C",
        "hyperparameters": {
            "learning_rate": 0.01,
            "beta1": 0.9,
            "beta2": 0.999,
            "epsilon": 1e-8
        }
    },
    "steps": {
        "value": 200,
        "justification": "Consistent with Phase B/C"
    },
    "seeds": {
        "value": 6,
        "justification": "Used in C3.3a/C3.3b"
    },
    "system_size": {
        "value": 3,
        "justification": "N=3 for D1-D3, N=4,5 for D4"
    }
}

# ============================================================================
# Hamiltonian Families (D1)
# ============================================================================

HAMILTONIAN_FAMILIES = {
    "asymmetric_XY": {
        "name": "Asymmetric XY Chain",
        "formula": "H = J₁₂(X₁X₂ + Y₁Y₂) + J₂₃(X₂X₃ + Y₂Y₃) + h(Z₁ + Z₂ + Z₃)",
        "parameters": {
            "J_12": 1.5,
            "J_23": 0.6,
            "h": 0.2
        },
        "status": "existing benchmark (Phase B/C)"
    },
    "transverse_field_ising": {
        "name": "Transverse-Field Ising Model",
        "formula": "H = -J Σ Z_i Z_{i+1} - h Σ X_i",
        "parameters": {
            "J": 1.0,
            "h": 0.5
        },
        "status": "new"
    },
    "heisenberg_xxz": {
        "name": "Heisenberg XXZ Model",
        "formula": "H = J Σ (X_i X_{i+1} + Y_i Y_{i+1} + Δ Z_i Z_{i+1})",
        "parameters": {
            "J": 1.0,
            "Delta": 0.5
        },
        "status": "new"
    },
    "xyz": {
        "name": "XYZ Model",
        "formula": "H = Σ (J_x X_i X_{i+1} + J_y Y_i Y_{i+1} + J_z Z_i Z_{i+1})",
        "parameters": {
            "J_x": 1.0,
            "J_y": 0.8,
            "J_z": 0.6
        },
        "status": "new"
    },
    "random_local": {
        "name": "Random Local Hamiltonian",
        "formula": "H = Σ h_i · σ_i + Σ J_ij σ_i · σ_j",
        "parameters": {
            "seed": 42,
            "h_range": [-0.5, 0.5],
            "J_range": [0.3, 0.7]
        },
        "status": "new"
    }
}

# ============================================================================
# Environment Set
# ============================================================================

ENVIRONMENT = {
    "type": "pure_dephasing",
    "parameters": {
        "gamma": [0.5, 1.0, 2.0]
    },
    "justification": "A1 environment from Phase C, non-uniform dephasing"
}

# ============================================================================
# Success Criteria
# ============================================================================

SUCCESS_CRITERIA = {
    "D1_hamiltonian_generality": {
        "description": "Non-trivial F* exists for each Hamiltonian family",
        "criterion": "Phi_star > contiguous_best for all 5 families",
        "metric": "Phi_star"
    },
    "D2_state_generality": {
        "description": "F* is stable across state types",
        "criterion": "d_F(F*_ground, F*_thermal) < 0.1",
        "metric": "d_F"
    },
    "D3_timescale_structure": {
        "description": "Characterize F*(τ) behavior",
        "categories": [
            "constant: d_F(F*(τ₁), F*(τ₂)) < 0.05 for all τ",
            "continuous_deformation: 0.05 < d_F < 0.3, smooth variation",
            "crossover: sharp change in d_F at specific τ",
            "discontinuous_transition: d_F > 0.5 at transition point"
        ],
        "metric": "d_F vs τ"
    },
    "D4_N_scaling": {
        "description": "Similar structure emerges for N=4,5",
        "criterion": "Non-trivial F* exists, comparable Phi_star",
        "metric": "Phi_star, d_F"
    }
}

# ============================================================================
# Protocol Summary
# ============================================================================

PROTOCOL = {
    "version": "D0_v1",
    "date": "2026-08-12",
    "metrics": METRICS,
    "parameters": PARAMETERS,
    "hamiltonian_families": HAMILTONIAN_FAMILIES,
    "environment": ENVIRONMENT,
    "success_criteria": SUCCESS_CRITERIA
}

def save_protocol():
    """Save protocol to JSON file."""
    output_dir = Path(__file__).parent.parent / "analysis_output"
    output_dir.mkdir(exist_ok=True)
    
    output_path = output_dir / "d0_protocol.json"
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(PROTOCOL, f, indent=2, ensure_ascii=False)
    
    print(f"Protocol saved to: {output_path}")
    return PROTOCOL

def print_protocol_summary():
    """Print human-readable protocol summary."""
    print("=" * 70)
    print("Phase D0: Protocol Freeze")
    print("=" * 70)
    print()
    
    print("Evaluation Metrics:")
    for name, info in METRICS.items():
        print(f"  - {name}: {info['description']}")
    print()
    
    print("Fixed Parameters:")
    for name, info in PARAMETERS.items():
        if isinstance(info.get('value'), dict):
            print(f"  - {name}: {info['value']}")
        else:
            print(f"  - {name} = {info['value']} ({info['unit'] if 'unit' in info else 'dimensionless'})")
    print()
    
    print("Hamiltonian Families (D1):")
    for name, info in HAMILTONIAN_FAMILIES.items():
        print(f"  {info['name']}:")
        print(f"    {info['formula']}")
        print(f"    Parameters: {info['parameters']}")
        print(f"    Status: {info['status']}")
    print()
    
    print("Environment:")
    print(f"  Type: {ENVIRONMENT['type']}")
    print(f"  Parameters: {ENVIRONMENT['parameters']}")
    print()
    
    print("Success Criteria:")
    for phase, criteria in SUCCESS_CRITERIA.items():
        print(f"  {phase}:")
        print(f"    {criteria['description']}")
        if 'criterion' in criteria:
            print(f"    Criterion: {criteria['criterion']}")
        elif 'categories' in criteria:
            print(f"    Categories:")
            for cat in criteria['categories']:
                print(f"      - {cat}")
    print()
    
    print("=" * 70)

if __name__ == "__main__":
    print_protocol_summary()
    save_protocol()
