#!/bin/bash
# lit_search_lrc.sh — Step 1 (literature check) for the Reduction Theorem / (TAU) novelty question.
# Runs a battery of web searches and saves raw JSON per query.
set -u
OUT=/home/z/my-project/scripts/lit
mkdir -p "$OUT"

run() {
  local file="$1"; local query="$2"; local num="${3:-8}"
  echo "== $file: $query"
  z-ai function -n web_search -a "{\"query\": \"$query\", \"num\": $num}" -o "$OUT/$file.json" \
    || echo "FAILED: $query"
}

run q01_complementary  "lonely runner conjecture complementary speeds pair"
run q02_mirror         "lonely runner conjecture pair of speeds sum mirror reduction proof"
run q03_cyclic         "lonely runner conjecture finite cyclic group formulation residues"
run q04_cusick         "Cusick view obstruction lonely runner conjecture"
run q05_bhk            "Bohman Holzman Kleitman six lonely runners"
run q06_barajas        "Barajas Serra lonely runner conjecture seven runners proof gap"
run q07_tight          "tight sets lonely runner conjecture classification"
run q08_smallcases     "lonely runner conjecture proof four runners elementary"
run q09_induction      "lonely runner conjecture reduce number of runners induction"
run q10_survey         "lonely runner conjecture survey history Wills"
run q11_spectrum       "lonely runner conjecture loneliness spectrum"
run q12_recent         "lonely runner conjecture arxiv 2024 2025 new results"
run q13_three          "lonely runner conjecture three runners proof Wills Cusick 1973"
run q14_chromatic      "lonely runner conjecture chromatic number distance graph Chen"

echo "ALL DONE"
ls -la "$OUT"
