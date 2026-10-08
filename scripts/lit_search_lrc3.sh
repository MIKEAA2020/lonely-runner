#!/bin/bash
# lit_search_lrc3.sh — Reviewer-directed round 2: new phrasings aimed at settling
# whether the EQUALITY M(V) = max_pair tau / pair-sum-lattice attainment is known.
# Distinct from rounds 1-2 (22 queries): targets denominator structure, finite
# cyclic reformulations, exact-value formulas, and the 2024-26 paper ring.
OUT=/home/z/my-project/scripts/lit
mkdir -p $OUT

declare -A Q
Q[r01_denominator_sum]='lonely runner conjecture optimal time denominator "sum of two speeds"'
Q[r02_finite_cyclic]='lonely runner conjecture "finite cyclic" reformulation residue classes equivalent'
Q[r03_dilate_cover]='covering cyclic group dilates arithmetic progression "lonely runner" arc'
Q[r04_exact_formula]='"lonely runner" "maximum loneliness" exact value formula attained'
Q[r05_zhang]='Zhang "lonely runner" tight set "one entry" modification regular'
Q[r06_perarnau_serra]='Perarnau Serra lonely runner survey "lonely runner" problem survey'
Q[r07_rosenfeld_10]='Rosenfeld Trakulthongchai lonely runner ten runners proof 2025'
Q[r08_kravitz_tight]='Kravitz Giri "tight sets" lonely runner classification maximum loneliness'
Q[r09_spectrum_fan]='Fan "lonely runner" spectrum conjecture amending 2026'
Q[r10_lattice_attained]='lonely runner "supremum attained" lattice rational time "speed sum" structure optimal'
Q[r11_arxiv_2026]='arxiv lonely runner conjecture 2026 computer assisted fifteen runners'
Q[r12_bohus_wong]='Bohus Wong lonely runner four speeds proof'
Q[r13_mirror_pair]='"lonely runner" mirror "v_i + v_j" time "1/(v_i+v_j)" lonely'
Q[r14_gap_formulation]='lonely runner "maximum gap" circle formulation equivalent speeds'
Q[r15_cusick_covering]='Cusick view obstruction "covering" discrete sphere "multiplicative"'
Q[r16_regular_time]='lonely runner regular speeds tight "1/(n+1)" equality characterization'

for key in $(echo "${!Q[@]}" | tr ' ' '\n' | sort); do
  echo "### $key : ${Q[$key]}"
  z-ai function -n web_search -a "{\"query\": \"${Q[$key]}\", \"num\": 8}" -o "$OUT/${key}.json" 2>/dev/null
  sleep 1
done
echo "DONE"
