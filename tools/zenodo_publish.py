#!/usr/bin/env python3
"""Create a NEW Zenodo record for the REM structural-selection paper v1.0.

Reads ZENODO_API_TOKEN from the environment (never prints it), stages the
publication package from ~/Desktop/ToyBox/REM/zenodo_package_v1.0, creates
the deposit, uploads all files, verifies, publishes, and prints the new DOI.

The token must be set locally, e.g. in ~/.zshrc:
    export ZENODO_API_TOKEN=...
It must never be committed, logged, or placed in CI artifacts.

Usage:
    ZENODO_API_TOKEN=... python3 tools/zenodo_publish.py [--draft-only]
"""
from __future__ import annotations

import json
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

ZENODO = "https://zenodo.org/api"
PKG = Path.home() / "Desktop" / "ToyBox" / "REM" / "zenodo_package_v1.0"

TITLE = ("Well-Posed Structural Selection and Finite-Time Crossovers in "
         "Variational Tensor-Product Factorization")
DESCRIPTION = (
    "We study variational tensor-product structure (TPS) selection in open "
    "quantum systems. Four results are established numerically for "
    "N = 3,4,5 spin chains under dephasing:\n\n"
    "1. The normalized local decay rate is shown to be singular under "
    "unrestricted TPS optimization: its denominator |Q_F rho|^2 can vanish, "
    "driving the objective to +infinity along near-product directions.\n"
    "2. The unnormalized instantaneous dynamical functional J_dyn^(0) "
    "provides a well-posed replacement in the tested protocol.\n"
    "3. The selected TPS responds to both the Hamiltonian and the quantum "
    "state.\n"
    "4. A finite-time structural crossover persists across the tested "
    "finite systems N = 3-5, and its crossover scale is quantitatively "
    "explained by a short-time inter-basin objective expansion "
    "(tau_c ~ -DeltaPhi_0/DeltaPhi_1; first order within 5-11%, quadratic "
    "consistency check within 0.0-0.5%).\n\n"
    "Companion specification: REM Spec v2.2 (10.5281/zenodo.21880505)."
)
KEYWORDS = [
    "quantum foundations",
    "tensor-product structure",
    "Hilbert-space factorization",
    "variational optimization",
    "quantum information",
    "structural crossover",
    "Relational Emergence Model",
]
RELATED = "10.5281/zenodo.21880505"  # REM_lambda v5 (companion spec/numerics)


def http(method, url, data=None, headers=None, timeout=120):
    req = urllib.request.Request(url, data=data, headers=headers or {}, method=method)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.status, json.loads(resp.read().decode() or "{}")


def upload_file(bucket, token, path):
    name = path.name
    url = f"{bucket}/{name}"
    req = urllib.request.Request(
        url,
        data=path.read_bytes(),
        headers={"Content-Type": "application/octet-stream",
                 "Authorization": f"Bearer {token}"},
        method="PUT",
    )
    with urllib.request.urlopen(req, timeout=600) as resp:
        return resp.status


def main():
    token = os.environ.get("ZENODO_API_TOKEN", "").strip()
    if not token:
        print("ERROR: ZENODO_API_TOKEN is not set in the environment.",
              file=sys.stderr)
        print("Set it locally (never in chat/Git/CI), e.g. in ~/.zshrc:",
              file=sys.stderr)
        print("    export ZENODO_API_TOKEN=...", file=sys.stderr)
        return 1
    draft_only = "--draft-only" in sys.argv

    files = sorted(
        [p for p in PKG.rglob("*") if p.is_file() and not p.name == ".DS_Store"]
    )
    print(f"package: {PKG} ({len(files)} files)")

    metadata = {
        "metadata": {
            "title": TITLE,
            "upload_type": "publication",
            "publication_type": "preprint",
            "description": DESCRIPTION,
            "creators": [{"name": "Maruko, Yoshifumi"}],
            "keywords": KEYWORDS,
            "license": "cc-by-4.0",
            "version": "1.0",
            "related_identifiers": [{
                "identifier": RELATED,
                "relation": "isSupplementTo",
            }],
        }
    }
    headers = {"Content-Type": "application/json",
               "Authorization": f"Bearer {token}"}

    # 1) create deposit
    status, dep = http("POST", f"{ZENODO}/deposit/depositions",
                       data=json.dumps(metadata).encode(), headers=headers)
    if status not in (200, 201):
        print(f"ERROR creating deposit: HTTP {status} {dep}", file=sys.stderr)
        return 1
    dep_id = dep["id"]
    bucket = dep["links"]["bucket"]
    print(f"deposit created: {dep['links']['html']} (draft)")

    # 2) upload files
    for f in files:
        rel = f.relative_to(PKG).as_posix()
        for attempt in range(3):
            try:
                st = upload_file(bucket, token, f)
                print(f"  uploaded {rel} ({f.stat().st_size} B)")
                break
            except urllib.error.HTTPError as e:
                if attempt == 2:
                    print(f"  FAILED {rel}: HTTP {e.code}", file=sys.stderr)
                    return 1
                time.sleep(3)

    # 3) verify
    _, dep = http("GET", f"{ZENODO}/deposit/depositions/{dep_id}",
                  headers={"Authorization": f"Bearer {token}"})
    uploaded = {x["filename"] for x in dep["files"]}
    expected = {f.relative_to(PKG).as_posix() for f in files}
    missing = expected - uploaded
    if missing:
        print(f"ERROR: {len(missing)} files not listed: {sorted(missing)[:5]}",
              file=sys.stderr)
        return 1
    print(f"verified: {len(uploaded)} files on the deposit")

    if draft_only:
        print(f"DRAFT ONLY — not published. Review at {dep['links']['html']}")
        return 0

    # 4) publish
    status, pub = http("POST",
                       f"{ZENODO}/deposit/depositions/{dep_id}/actions/publish",
                       headers={"Authorization": f"Bearer {token}"})
    if status != 202:
        print(f"ERROR publishing: HTTP {status} {pub}", file=sys.stderr)
        return 1
    doi = pub.get("doi") or pub.get("metadata", {}).get("doi")
    print(f"PUBLISHED: {doi}")
    print(f"record URL: {pub.get('links', {}).get('record_html', '')}")
    print(f"badge: https://zenodo.org/badge/doi/{doi}.svg")
    return 0


if __name__ == "__main__":
    sys.exit(main())
