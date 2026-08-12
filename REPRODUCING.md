# Reproducing the current REM results

This guide separates three different tasks that are easy to confuse in this repository:

1. **Integrity check** — does the current codebase still pass its automated tests?
2. **Evidence inspection** — can a reader trace each paper claim to the committed script and raw output without recomputing it?
3. **Numerical reproduction** — can the relevant Phase D / Phase E calculation be rerun from the frozen protocol?

The canonical definition is [REM Spec v2.2](REM_spec_v2_2.md). The primary evidence record is [Phase D](analysis/PHASE_D_FINAL_REPORT.md), including its Phase E addendum.

## 1. Environment

Python 3.11+ is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## 2. Fast integrity check

```bash
pytest -q
```

The Phase D report records **74 passed + 1 xfailed** for the validation state closed on 2026-08-12. Treat that number as a historical checkpoint, not as a permanent promise that the test count will never change.

A passing test suite checks code invariants and regressions. It is **not** equivalent to rerunning every expensive Phase D / D4 / Phase E experiment.

## 3. Inspect the committed evidence first

Before launching a numerical rerun, locate the claim in [CURRENT_STATUS.md](CURRENT_STATUS.md), then follow it to the corresponding section of the Phase D report.

The primary committed outputs are under `analysis_output/`. The report's [§9 reproducibility table](analysis/PHASE_D_FINAL_REPORT.md#9-reproducibility) records the script, output, and commit associated with each validation phase.

> **Do not rerun current analysis scripts directly on the authoritative evidence archive.** Several scripts write to tracked files under `analysis_output/`; perform reruns on a dedicated working branch or disposable clone/worktree, keep the starting tree clean, and review any `analysis_output/` diff before preserving new results.

This allows a third party to audit the evidence chain without assuming that a newly rerun optimizer must reproduce every floating-point value bit-for-bit.

## 4. Current claim-to-script map

Run these commands from the repository root.

| Research question | Script | Main committed output |
|---|---|---|
| Why the normalized objective fails | `python analysis/d2_1_haar_singularity_audit.py` | `analysis_output/d2_1_haar_singularity_audit.json` |
| Does the denominator-free instantaneous functional remove the failure mode? | `python analysis/d2_2a_unnormalized_instantaneous.py` | `analysis_output/d2_2a_unnormalized_instantaneous.json` |
| Does the finite-time extension approach the instantaneous limit without recurrence of the singularity? | `python analysis/d2_2b_finite_time.py` and `python analysis/d2_2b_gates.py` | `analysis_output/d2_2b_tau_*.json`, `analysis_output/d2_2b_gates.json` |
| Does the selected TPS respond to Hamiltonian choice? | `python analysis/d1r_canonical_revalidation.py` | `analysis_output/d1r_canonical_revalidation.json` |
| Does the selected TPS respond to state choice? | `python analysis/d2r_state_revalidation.py` | `analysis_output/d2r_state_revalidation.json` |
| Is there a fixed-representative finite-time crossover? | `python analysis/d3_timescale.py` | `analysis_output/d3_timescale.json` |
| Does the crossover persist at the tested finite sizes? | `python analysis/d4_nscaling.py` and `python analysis/d4_gates.py` | `analysis_output/d4_*.json`, `analysis_output/d4_finite_size.json` |
| What sets the crossover scale? | `python analysis/e1_tauc_mechanism.py` | `analysis_output/e1_tauc_mechanism.json` |

The D4 full-quotient $N=5$ calculation uses the documented reduced protocol and is materially more expensive than the fast test suite. See Phase D §7 before rerunning it.

## 5. What must be recorded in an independent reproduction

For a result intended to count as an independent check, record at minimum:

- repository commit;
- Python and dependency versions;
- random seeds;
- $N$ and bipartition dimensions;
- full-quotient vs restricted-subspace coverage;
- optimizer steps, learning rate, and stopping rule;
- finite-time grid when estimating $\tau_c$;
- convergence diagnostics;
- any deviation from the frozen protocol.

For D3/D4, estimate a crossover from the **fixed-representative objective difference**

$$
\Delta\Phi(\tau)=\Phi_\tau(F_B)-\Phi_\tau(F_A),
$$

and its sign change. Optimizer basin hopping alone is not the crossover definition.

## 6. Current open replication targets

- [Issue #5 — D3 independent finite-time crossover reproduction](https://github.com/malko73/rem/issues/5)
- [Issue #4 — D4 finite-size persistence stress test](https://github.com/malko73/rem/issues/4)

Negative or non-replication results are valid outcomes and should be preserved with the protocol details above.

## 7. Legacy reproduction path

The root script

```bash
./reproduce.sh
```

runs the **earlier REM4 asymmetric-XY benchmark** implemented in `src/rem4_numerical.py`. It does **not** reproduce the current Spec v2.2 Phase D validation chain.

That script remains in the repository so earlier REM4 numerical work can still be reproduced. Do not use its output as evidence for the current denominator-free open-system claims unless the relevant current analysis explicitly references it.

## 8. CI interpretation

The GitHub Actions CI currently runs the automated test suite and a reduced-grid smoke test of the earlier REM4 scaffold. Therefore:

- **green CI** means the tested code paths are passing;
- it does **not** mean all D-series and Phase E optimization runs were recomputed in CI;
- the committed `analysis_output/` files plus the Phase D report remain the evidence map for the full numerical record.

## 9. Paper and specification regeneration

The repository also contains publication tooling:

- `tools/gen_spec.py` — specification generation support;
- `tools/gen_paper.py` — paper generation support;
- `tools/gen_figures.py` — figure generation support;
- `tools/qc_paper.py` — manuscript quality-control checks.

These tools manage research artifacts; they are separate from the numerical validation scripts listed above.
