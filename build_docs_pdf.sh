#!/usr/bin/env bash
# Canonical md -> PDF build for the project docs (DRAFTING_FACTS, TESTS_AND_DIAGNOSTICS,
# METHODS_AND_SYSTEMATICS, ...). ALWAYS use this, not a bare `pandoc` call:
# omitting the header include drops the STIX Two Math glyph fallback and math/Greek
# symbols (→ ≤ ⊕ σ Γ …) render as missing/?. (Learned the hard way 2026-06-11.)
#
#   .pandoc_header.tex       — STIX Two Math fallback for glyphs absent from STIX Two
#                              Text + table padding/font + \BreakablePath (required by
#                              the lua filter).
#   .pandoc_break_paths.lua  — rewrites inline `code` paths with break opportunities so
#                              long paths don't overflow table cells.
#
# Usage:  ./build_docs_pdf.sh                       # rebuild the default doc set
#         ./build_docs_pdf.sh FILE.md [FILE2.md ...] # rebuild specific files
set -euo pipefail
cd "$(dirname "$0")"

DEFAULT_DOCS=(DRAFTING_FACTS_paper_2026-05-29.md TESTS_AND_DIAGNOSTICS.md)
DOCS=("$@"); [ ${#DOCS[@]} -eq 0 ] && DOCS=("${DEFAULT_DOCS[@]}")

for md in "${DOCS[@]}"; do
  [ -f "$md" ] || { echo "skip (missing): $md"; continue; }
  pdf="${md%.md}.pdf"
  pandoc "$md" -o "$pdf" \
    --pdf-engine=xelatex \
    -V mainfont="STIX Two Text" \
    -V geometry:margin=1in \
    --include-in-header=.pandoc_header.tex \
    --lua-filter=.pandoc_break_paths.lua
  # report any glyphs that still failed to render
  miss=$(pandoc "$md" -o /dev/null --pdf-engine=xelatex -V mainfont="STIX Two Text" \
         -V geometry:margin=1in --include-in-header=.pandoc_header.tex \
         --lua-filter=.pandoc_break_paths.lua 2>&1 | grep -ci "missing character" || true)
  echo "built $pdf  (missing glyphs: $miss)"
done
