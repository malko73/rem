# Analysis index

The `analysis/` directory contains several generations of numerical work. **Not every script in this directory represents the current canonical REM formulation.** Use this index together with [CURRENT_STATUS.md](../CURRENT_STATUS.md), [REM Spec v2.2](../REM_spec_v2_2.md), and the [Phase D validation report](PHASE_D_FINAL_REPORT.md).

## Current primary validation chain

The present paper-level evidence is organized around the falsification-driven D-series and its Phase E mechanism analysis.

| Phase | Purpose | Primary scripts | Primary outputs |
|---|---|---|---|
| D2.1 | Falsify the normalized global TPS objective | `d2_1_haar_singularity_audit.py` | `../analysis_output/d2_1_haar_singularity_audit.json` |
| D2.2-A | Test denominator-free $J_{\mathrm{dyn}}^{(0)}$ | `d2_2a_unnormalized_instantaneous.py` | `../analysis_output/d2_2a_unnormalized_instantaneous.json` |
| D2.2-B | Check finite-time extension and $\tau\to0$ consistency | `d2_2b_finite_time.py`, `d2_2b_gates.py` | `../analysis_output/d2_2b_tau_*.json`, `../analysis_output/d2_2b_gates.json` |
| D1-R | Revalidate Hamiltonian response under the corrected canonical functional | `d1r_canonical_revalidation.py` | `../analysis_output/d1r_canonical_revalidation.json` |
| D2-R | Revalidate state response under the corrected canonical functional | `d2r_state_revalidation.py` | `../analysis_output/d2r_state_revalidation.json` |
| D3 | Identify the fixed-representative finite-time crossover | `d3_timescale.py` | `../analysis_output/d3_timescale.json` |
| D4 | Test persistence for the reported $N=3$–$5$ finite systems | `d4_nscaling.py`, `d4_gates.py` | `../analysis_output/d4_*.json`, `../analysis_output/d4_finite_size.json` |
| E1 | Explain the crossover scale through inter-basin expansion | `e1_tauc_mechanism.py`, `e1_sensitivity.py` | `../analysis_output/e1_tauc_mechanism.json`, `../analysis_output/e1_sensitivity.json` |

The exact frozen protocol, gates, commits, caveats, and reported values are consolidated in [PHASE_D_FINAL_REPORT.md](PHASE_D_FINAL_REPORT.md).

## Supporting D-series scripts

The following files support protocol construction, generality checks, or intermediate validation:

- `d0_protocol.py`
- `d1_generality.py`
- `d1_hamiltonians.py`
- `d2_state_generality.py`

Where a corrected revalidation script exists (`d1r_*`, `d2r_*`), the **revalidation result is the relevant current evidence** for Spec v2.2.

## Earlier / exploratory analysis

Scripts with prefixes such as `b*`, `c*`, `gamma_F*`, `m2_vs_var*`, `compare_*`, and `robustness*` document earlier stages of the research program, optimizer comparisons, proxy tests, and falsification work that led to the current formulation.

They are retained for traceability, but **do not infer that an older script's objective or terminology is canonical simply because the file remains in the repository**. Check [CHANGELOG.md](../CHANGELOG.md) and [CURRENT_STATUS.md](../CURRENT_STATUS.md) before using an earlier output as a current REM claim.

## Machine-readable outputs

`../analysis_output/` contains committed JSON, text tables, and figures. These files serve two purposes:

1. preserve the numerical record that underlies the validation report;
2. allow third parties to inspect reported values without immediately rerunning expensive optimizations.

For current claim-to-output mapping, use [Phase D §9](PHASE_D_FINAL_REPORT.md#9-reproducibility) or [REPRODUCING.md](../REPRODUCING.md).

## Reproduction rule for D3/D4

The finite-time crossover is defined from fixed basin representatives:

$$
\Delta\Phi(\tau)=\Phi_\tau(F_B)-\Phi_\tau(F_A),
\qquad
\Delta\Phi(\tau_c)=0.
$$

An optimizer jumping from one basin to another is **not by itself** the crossover definition. This distinction is required when independently checking D3 or D4.

## Interpretation hierarchy

When two files appear to conflict, use this order:

1. `../REM_spec_v2_2.md` — canonical definitions;
2. `PHASE_D_FINAL_REPORT.md` — current numerical evidence and limitations;
3. `../CURRENT_STATUS.md` — repository-level interpretation map;
4. D-series / E-series scripts and committed outputs;
5. earlier B/C/REM4-era analyses — historical or supporting context.
