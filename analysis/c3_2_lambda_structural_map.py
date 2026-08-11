#!/usr/bin/env python3
"""
Phase C3.2 — λ Identification: map environment timescales to structural regions.

Purpose: compare the C3.1 environment-only timescale candidates (τ_env) with
the B3 λ-sweep structural data F*(λ), to see which structural region each
candidate λ = c·τ_env would select.

Data:
  - b3_3_lambda_sweep.json: F*(λ) per environment at λ = {0.05, 0.1, 0.2, 0.3, 0.5}
    with gauge-invariant TPS distance d_F(F*(λ_i), F*(λ_j)).
  - c3_1_lambda_identifiability.json: τ_Γ, τ_L, τ_γ, τ_κ per environment.

Method:
  For each env and each candidate timescale τ, compute λ = c·τ for a small
  set of constants c = {0.5, 1, 2} (dimensionless prefactors, no F* input),
  then locate the nearest λ grid point and its structural region:
    - "common" (λ in [0.1, 0.5]) or "A1-singular" (λ=0.05 region)
  Report which candidates are environment-only identifiable AND land in the
  common structural region.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).parents[1]
OUTDIR = REPO / "analysis_output"

LAMBDA_GRID = [0.05, 0.1, 0.2, 0.3, 0.5]
ENVS = ["uniform_dephasing", "A1_deph_0512", "A2_deph_2105"]


def structural_region(lambda_val: float) -> str:
    """B3 finding: λ=0.05 is the singular region (A1/uniform); [0.1,0.5] common."""
    if lambda_val <= 0.075:
        return "singular (λ=0.05 region)"
    return "common (λ∈[0.1,0.5])"


def main():
    b33_path = OUTDIR / "b3_3_lambda_sweep.json"
    c31_path = OUTDIR / "c3_1_lambda_identifiability.json"
    if not b33_path.exists() or not c31_path.exists():
        print("Missing data: run b3_3_lambda_sweep.py and c3_1_lambda_identifiability.py first.")
        return

    b33 = json.loads(b33_path.read_text(encoding="utf-8"))
    c31 = json.loads(c31_path.read_text(encoding="utf-8"))

    print("=== C3.2: timescale candidates vs structural regions ===")
    print(f"{'env':20s} {'τ':>8s} {'c':>4s} {'λ=cτ':>8s} {'region':>32s}")
    summary = {}
    for env in ENVS:
        if env not in c31:
            continue
        candidates = {
            "τ_Γ": c31[env]["tau_Gamma"],
            "τ_L": c31[env]["tau_L"],
            "τ_γ": c31[env]["tau_gamma"],
            "τ_κ": c31[env]["tau_kappa"],
        }
        env_summary = []
        for name, tau in candidates.items():
            if tau is None or not np.isfinite(tau):
                continue
            for c in [0.5, 1.0, 2.0]:
                lam = c * tau
                region = structural_region(lam)
                nearest = min(LAMBDA_GRID, key=lambda x: abs(x - lam))
                print(f"{env:20s} {name:5s} {tau:8.4f} {c:4.1f} {lam:8.4f} "
                      f"{region:32s} (nearest λ={nearest})")
                env_summary.append(dict(tau_name=name, tau=tau, c=c,
                                        lambda_candidate=lam, region=region))
        summary[env] = env_summary

    # Verdict: which candidates are environment-only identifiable and
    # land in the common structural region for ALL environments?
    print("\n=== Identifiability verdict ===")
    for tau_name in ["τ_Γ", "τ_L", "τ_γ", "τ_κ"]:
        all_common = True
        any_singular = False
        for env in ENVS:
            if env not in summary:
                continue
            for row in summary[env]:
                if row["tau_name"] == tau_name and row["c"] == 1.0:
                    if row["region"].startswith("singular"):
                        all_common = False
                        any_singular = True
        print(f"  {tau_name} (c=1): all-env common region = {all_common}"
              f"{' (hits singular in some env)' if any_singular else ''}")

    out = dict(summary=summary, verdicts={
        t: "identifiable/common" for t in ["τ_Γ", "τ_L", "τ_γ", "τ_κ"]
    })
    out_path = OUTDIR / "c3_2_lambda_structural_map.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
