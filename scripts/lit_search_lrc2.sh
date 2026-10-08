#!/bin/bash
# lit_search_lrc2.sh — targeted follow-up batch on the novelty question.
set -u
OUT=/home/z/my-project/scripts/lit
mkdir -p "$OUT"

run() {
  local file="$1"; local query="$2"; local num="${3:-8}"
  echo "== $file: $query"
  z-ai function -n web_search -a "{\"query\": \"$query\", \"num\": $num}" -o "$OUT/$file.json" \
    || echo "FAILED: $query"
}

run q15_pairsum    "lonely runner conjecture sum of two speeds time rational"
run q16_betke      "Betke Wills lonely runner four speeds proof diophantine"
run q17_kravitz    "Kravitz lonely runner maximum loneliness tight sets"
run q18_tao        "Tao some remarks lonely runner conjecture 2018 theorem"
run q19_sungk      "Sungkawichai lonely runner finite checking computer"
run q20_product    "lonely runner conjecture speed product primitive set bound"
run q21_rosenfeld  "Rosenfeld lonely runner nine ten runners proof technique"
run q22_jtimes     "lonely runner conjecture equally spaced polygon time proof"

echo "ALL DONE"
