#!/usr/bin/env python3
"""Regenerate papers/REM_paper_v0_1.tex from papers/REM_paper_v0_1.md.

Reproducible generation: pandoc body + the Spec v2.2 preamble (which carries
the newunicodechar mappings and the longtable/booktabs/array/calc table
support). Usage: python3 tools/gen_paper.py [--check]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MD = REPO / "papers" / "REM_paper_v0_1.md"
SPEC_TEX = REPO / "papers" / "REM_spec_v2_2.tex"
OUT = REPO / "papers" / "REM_paper_v0_1.tex"

PAPER_NEWUNICODE = (
    # text-mode unicode used by the paper but not mapped in the spec preamble
    "\\newunicodechar{Δ}{\\ensuremath{\\Delta}}\n"
    "\\newunicodechar{Ĩ}{\\ensuremath{\\tilde{I}}}\n"
    "\\newunicodechar{§}{\\S}\n"
    "\\newunicodechar{–}{--}\n"
)



def _strip_figure_alt(body: str) -> str:
    """Remove pandoc's alt={...} attribute from includegraphics.

    Pandoc 3.x emits [keepaspectratio,alt={...}] where the alt text is a
    text-mode linearization of the caption containing backslash commands.
    That breaks pdflatex's keyval parsing (undefined control sequences /
    brace imbalance). The visible \caption carries the same content in
    proper math mode, so the alt attribute is dropped.
    """
    out = []
    i = 0
    while True:
        j = body.find("alt={", i)
        if j < 0:
            out.append(body[i:])
            break
        out.append(body[i:j])
        k = j + 5
        depth = 1
        while k < len(body) and depth > 0:
            ch = body[k]
            if ch == "\\":
                k += 2
                continue
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
            k += 1
        i = k
    return "".join(out)


def generate() -> str:
    spec = SPEC_TEX.read_text(encoding="utf-8")
    marker = "\\begin{document}"
    preamble = spec[: spec.index(marker) + len(marker)]

    body = subprocess.run(
        ["pandoc", str(MD), "--from", "markdown", "--to", "latex"],
        capture_output=True, text=True, check=True,
    ).stdout
    body = _strip_figure_alt(body)

    anchor = "\\newunicodechar{→}{\\ensuremath{\\to}}\n"
    if anchor in preamble:
        preamble = preamble.replace(anchor, anchor + PAPER_NEWUNICODE, 1)
    # The spec preamble has no figures, so it neither loads graphicx nor
    # defines \pandocbounded (pandoc >= 3.1 wraps includegraphics in it).
    if "\\usepackage{graphicx}" not in preamble:
        preamble = preamble.replace(
            "\\begin{document}",
            "\\usepackage{graphicx}\n"
            "\\providecommand{\\pandocbounded}[1]{#1}\n"
            "\\begin{document}",
            1,
        )
    return preamble + "\n\n" + body + "\n\\end{document}\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    tex = generate()
    if args.check:
        if OUT.read_text(encoding="utf-8") != tex:
            print("DRIFT: papers/REM_paper_v0_1.tex is out of date; run tools/gen_paper.py")
            return 1
        print("OK: papers/REM_paper_v0_1.tex matches papers/REM_paper_v0_1.md")
        return 0

    OUT.write_text(tex, encoding="utf-8")
    print(f"wrote {OUT} ({len(tex.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
