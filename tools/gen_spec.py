#!/usr/bin/env python3
"""Regenerate papers/REM_spec_v2_2.tex from REM_spec_v2_2.md.

Reproducible generation: takes the pandoc-generated body of the current spec
Markdown and assembles it with the v2.0 pandoc preamble plus the pdflatex
fixes required by this document:
  - newunicodechar mappings for text-mode Unicode (Φ λ ρ τ − ⊗ ≈ — and the
    v2.2 additions × Γ →)
  - longtable/booktabs/array + calc (the D-series evidence table)

Usage: python3 tools/gen_spec.py [--check]
  --check: exit non-zero if the .tex on disk differs from the generated one.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MD = REPO / "REM_spec_v2_2.md"
V20 = REPO / "papers" / "REM_spec_v2_0.tex"
OUT = REPO / "papers" / "REM_spec_v2_2.tex"

NEWUNICODE_V22 = (
    "\\newunicodechar{⊗}{\\ensuremath{\\otimes}}\n"
    "\\newunicodechar{×}{\\ensuremath{\\times}}\n"
    "\\newunicodechar{Γ}{\\ensuremath{\\Gamma}}\n"
    "\\newunicodechar{→}{\\ensuremath{\\to}}\n"
)
TABLE_PKGS = (
    "% Table support (pandoc longtable output)\n"
    "\\usepackage{longtable,booktabs,array}\n"
    "\\usepackage{calc}\n"
)


def generate() -> str:
    v20 = V20.read_text(encoding="utf-8")
    marker = "\\begin{document}"
    preamble = v20[: v20.index(marker) + len(marker)]

    body = subprocess.run(
        ["pandoc", str(MD), "--from", "markdown", "--to", "latex"],
        capture_output=True, text=True, check=True,
    ).stdout

    # insert v2.2-specific preamble additions after the v2.0 newunicodechar block
    anchor = "\\newunicodechar{⊗}{\\ensuremath{\\otimes}}\n"
    preamble = preamble.replace(anchor, anchor + NEWUNICODE_V22 + TABLE_PKGS, 1)
    return preamble + "\n\n" + body + "\n\\end{document}\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    tex = generate()
    if args.check:
        if OUT.read_text(encoding="utf-8") != tex:
            print("DRIFT: papers/REM_spec_v2_2.tex is out of date; run tools/gen_spec.py")
            return 1
        print("OK: papers/REM_spec_v2_2.tex matches REM_spec_v2_2.md")
        return 0

    OUT.write_text(tex, encoding="utf-8")
    print(f"wrote {OUT} ({len(tex.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
