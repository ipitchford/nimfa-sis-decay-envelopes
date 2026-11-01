#!/usr/bin/env bash
# Rebuild nimfa_envelopes.pdf from manuscript.md (Pandoc -> LaTeX -> Tectonic/XeTeX).
set -euo pipefail
cd -- "$(dirname -- "$0")"
pandoc manuscript.md \
  --from=markdown+tex_math_dollars+pipe_tables \
  --standalone --include-in-header=paper_header.tex \
  --variable=geometry:margin=22mm --variable=papersize:a4 --variable=fontsize:11pt \
  --variable=colorlinks:true \
  --resource-path=. --output=nimfa_envelopes.tex
tectonic --chatter minimal --keep-logs nimfa_envelopes.tex > pdf_build.log 2>&1
