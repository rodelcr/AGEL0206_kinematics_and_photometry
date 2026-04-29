#!/usr/bin/env python3
"""
Build TESTS_AND_DIAGNOSTICS.pdf from the markdown source.

Pipeline: pandoc (md -> tex) -> patch preamble for unicode fallbacks +
table-friendly typography -> xelatex (tex -> pdf, twice for refs).

Leaves TESTS_AND_DIAGNOSTICS.md untouched.
"""
from __future__ import annotations
import os
import re
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(REPO, "TESTS_AND_DIAGNOSTICS.md")
TEX = os.path.join(REPO, "TESTS_AND_DIAGNOSTICS.tex")
PDF = os.path.join(REPO, "TESTS_AND_DIAGNOSTICS.pdf")


PREAMBLE_PATCH = r"""
% --- patch added by scripts/build_tests_pdf.py ---
\usepackage{newunicodechar}
\usepackage[htt]{hyphenat}    % allow hyphenation inside \texttt{}
\usepackage{pifont}           % \ding{51} = check mark
\providecommand{\cmark}{{\color{green!55!black}\ding{51}}}
\newunicodechar{✅}{\cmark}
\newunicodechar{→}{\ensuremath{\rightarrow}}
\newunicodechar{⊕}{\ensuremath{\oplus}}
\newunicodechar{∈}{\ensuremath{\in}}
\newunicodechar{≤}{\ensuremath{\leq}}
\newunicodechar{≥}{\ensuremath{\geq}}
\newunicodechar{≈}{\ensuremath{\approx}}
\newunicodechar{²}{\ensuremath{{}^{2}}}
\newunicodechar{Σ}{\ensuremath{\Sigma}}
\newunicodechar{σ}{\ensuremath{\sigma}}
\newunicodechar{α}{\ensuremath{\alpha}}
\newunicodechar{β}{\ensuremath{\beta}}
\newunicodechar{Δ}{\ensuremath{\Delta}}
\newunicodechar{×}{\ensuremath{\times}}
\newunicodechar{±}{\ensuremath{\pm}}
\newunicodechar{−}{\ensuremath{-}}    % unicode minus
\newunicodechar{★}{\ensuremath{\star}}
\newunicodechar{☉}{\ensuremath{\odot}}
\newunicodechar{”}{\textquotedblright}
\newunicodechar{“}{\textquotedblleft}
% Make tabular cells more forgiving of long monospace strings.
\renewcommand{\arraystretch}{1.05}
% --- end patch ---
"""


def run(cmd, **kw):
    print("$", " ".join(cmd))
    return subprocess.run(cmd, check=True, **kw)


def main():
    if shutil.which("pandoc") is None:
        sys.exit("pandoc not on PATH")
    if shutil.which("xelatex") is None:
        sys.exit("xelatex not on PATH")

    # 1) MD -> TeX (standalone). Tight margins + small monospace make the
    # wide test matrix fit.
    run([
        "pandoc", MD, "-o", TEX,
        "--standalone",
        "--pdf-engine=xelatex",
        "--columns=200",
        "-V", "geometry:margin=0.55in",
        "-V", "documentclass=article",
        "-V", "fontsize=10pt",
        "-V", "colorlinks=true",
        "-V", "linkcolor=blue",
        "-V", "urlcolor=blue",
        "-V", "mainfont=STIX Two Text",
        "-V", "monofont=Menlo",
        "-V", "monofontoptions=Scale=0.78",
    ])

    # 2) Patch the preamble: unicode fallbacks + tt-hyphenation + sloppy.
    with open(TEX, "r", encoding="utf-8") as f:
        tex = f.read()

    if r"% --- patch added by scripts/build_tests_pdf.py ---" not in tex:
        tex = tex.replace(
            r"\begin{document}",
            PREAMBLE_PATCH + "\n" + r"\begin{document}" + "\n" + r"\sloppy",
            1,
        )

    # Pandoc emits longtable with `>{\raggedright\arraybackslash}p{...}`
    # for table cells. Add `\hangafter=0\hangindent=0pt` so multi-line
    # cells in narrow columns don't acquire phantom indents. Belt & braces.
    tex = re.sub(
        r"\\begin\{longtable\}",
        r"\\begingroup\\small\\setlength{\\tabcolsep}{4pt}\\begin{longtable}",
        tex,
    )
    tex = re.sub(
        r"\\end\{longtable\}",
        r"\\end{longtable}\\endgroup",
        tex,
    )

    with open(TEX, "w", encoding="utf-8") as f:
        f.write(tex)

    # 3) xelatex twice for cross-references / longtable widths.
    for _ in range(2):
        subprocess.run(
            ["xelatex", "-interaction=nonstopmode", "-halt-on-error",
             os.path.basename(TEX)],
            cwd=REPO,
            check=False,  # don't abort on the "may have changed" warning
        )

    # 4) Cleanup aux files.
    for ext in ("aux", "log", "out", "toc"):
        p = os.path.join(REPO, f"TESTS_AND_DIAGNOSTICS.{ext}")
        if os.path.exists(p):
            os.remove(p)

    print(f"OK -> {PDF}")


if __name__ == "__main__":
    main()
