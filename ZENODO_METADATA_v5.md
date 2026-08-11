# Zenodo Metadata — REM_lambda v5

## Title
Competing Informational and Dynamical Criteria for Subsystem Structure:
Relational Emergence Model — Revised Version 5

## Version
v5

## Publication type
Preprint

## Creator
Yoshifumi Maruko

## Description
REM_lambda v5 is a corrective revision of v4 of the Relational Emergence
Model (REM).

This revision addresses two independent issues identified during
post-publication validation.

First, the continuous tensor-product-structure optimization used in v4
contained a non-unitary parameterization. The generator was anti-Hermitian
but was exponentiated as exp(iG), rather than exp(G). The corrected
parameterization is explicitly unitary, and all continuous-optimization
results have been recomputed. The previously reported continuous optimum
Phi ≈ 1.05 and the associated 1.95x improvement are therefore withdrawn.
The dimension of the relevant U(8)/(U(4) x U(2)) factorization space is
also corrected from 44 to 45 after accounting for the one-dimensional
U(1) kernel.

Second, open-system falsification tests show that the closed-system
quantity

    C_H^closed(F) = <H_boundary(F)^2>

does not generally predict factorization-dependent decoherence stability.
It is therefore retained only as a closed-system structural surrogate.

The canonical REM variational form is generalized to

    Phi(F; lambda, rho, L, tau)
      = I_rho(F) - lambda C_dyn(F; rho, L, tau),

where C_dyn is an environment-dependent signed dynamical functional.

The first operational open-system realization introduced in v5 is the
signed structural decoherence rate

    C_Gamma^(0)(F; rho, L)
      = Gamma_F^exact(0).

Numerical Lindblad tests show that structural stability depends on the
environment: symmetric noise can produce factorization-independent initial
rates, while non-uniform dephasing can reverse the preferred ordering
between factorizations.

Accordingly, the previously reported lambda* ≈ 0.165 is retained only as
a structural crossover of the closed-system second-moment surrogate model;
it is no longer interpreted as a demonstrated open-system physical
crossover.

For C_dyn with dimensions of inverse time, lambda has dimensions of time
and is interpreted as a characteristic timescale.

This release is synchronized with REM Spec v2.0 and with the corrected,
reproducible numerical implementation.

## Major corrections from v4

1. Corrected the continuous TPS parameterization to preserve unitarity.
2. Recomputed all affected continuous-optimization results.
3. Withdrawn the previous Phi ≈ 1.05 / 1.95x continuous-improvement claim.
4. Corrected factorization-space dimension 44 -> 45.
5. Reclassified <H_boundary^2> as a closed-system structural surrogate.
6. Added open-system falsification results.
7. Introduced the environment-dependent signed dynamical functional C_dyn.
8. Introduced Gamma_F^exact(0) as the first operational open-system realization.
9. Reinterpreted lambda* ≈ 0.165 as a closed-system surrogate crossover.
10. Synchronized the paper with REM Spec v2.0.

## Reproducibility

Source repository:
https://github.com/malko73/rem

Reference code commit:
a9fb272

Publication-package commit:
56b2c61

Tests:
50 passed

LaTeX build:
GitHub Actions — green
- REM_lambda v5: 6 pages
- REM Spec v2.0: 7 pages
- fatal errors: 0
- undefined references: 0
- missing files: 0

## Files

- REM_lambda_v5.pdf
- REM_lambda_v5_source.tex
- REM_spec_v2_0.pdf
- REM_spec_v2_0.tex
- CHANGELOG.md

## Keywords

Relational Emergence Model
REM
quantum foundations
tensor product structure
subsystem structure
open quantum systems
Lindblad dynamics
decoherence
relational quantum mechanics
variational principle
quantum information
emergent structure

## License

CC BY 4.0 (same as REM_lambda v4)

## Related identifiers

Previous version:
REM_lambda v4
DOI: 10.5281/zenodo.21427776

Repository:
https://github.com/malko73/rem

## Notes

This version supersedes the numerical continuous-TPS results and the
open-system interpretation of the dynamical proxy reported in v4.
The earlier version is retained as part of the publication history.
