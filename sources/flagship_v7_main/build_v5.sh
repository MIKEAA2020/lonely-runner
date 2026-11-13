#!/bin/bash
# Build the v5 flagship paper: compile (2 passes), regenerate cover,
# validate, merge, QA. Run from scripts/paper2_v5.
set -e
cd /home/z/my-project/scripts/paper2_v5

echo "== [1/6] tectonic pass 1+2"
tectonic main.tex > /dev/null 2>&1 || tectonic main.tex
tectonic main.tex 2>&1 | tail -2
echo "pages: $(python3 -c "from pypdf import PdfReader; print(len(PdfReader('main.pdf').pages))")"

echo "== [2/6] check-html on cover"
python3 /home/z/my-project/skills/pdf/scripts/poster_validate.py check-html cover.html | tail -5

echo "== [3/6] cover_validate.js"
node /home/z/my-project/skills/pdf/scripts/cover_validate.js cover.html 2>&1 | tail -3 || true

echo "== [4/6] render cover (html2poster)"
node /home/z/my-project/skills/pdf/scripts/html2poster.js cover.html --output cover.pdf --width 794px 2>&1 | tail -2

echo "== [5/6] merge"
python3 merge_cover.py

echo "== [6/6] pdf_qa (skip cover, formulas)"
python3 /home/z/my-project/skills/pdf/scripts/pdf_qa.py --skip-cover --formulas final.pdf 2>&1 | tail -12
