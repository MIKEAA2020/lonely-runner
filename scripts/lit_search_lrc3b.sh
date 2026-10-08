#!/bin/bash
# lit_search_lrc3b.sh — retry round-2 queries without embedded double quotes.
OUT=/home/z/my-project/scripts/lit
mkdir -p $OUT

run() {  # run <file-key> <query>
  local key="$1"; local query="$2"
  if [ -s "$OUT/${key}.json" ]; then echo "skip $key (exists)"; return; fi
  echo "### $key : $query"
  z-ai function -n web_search -a "{\"query\": \"$query\", \"num\": 8}" -o "$OUT/${key}.json" >/dev/null 2>&1
  sleep 1
}

run r01_denominator_sum 'lonely runner conjecture optimal time denominator sum of two speeds'
run r03_dilate_cover 'covering cyclic group dilated arithmetic progression lonely runner arc'
run r04_exact_formula 'lonely runner maximum loneliness exact value formula attained'
run r05_zhang 'Zhang lonely runner tight set one entry modification regular speeds'
run r06_perarnau_serra 'Perarnau Serra lonely runner survey problem'
run r08_kravitz_tight 'Kravitz tight instances lonely runner maximum loneliness classification'
run r09_spectrum_fan 'Fan lonely runner spectrum conjecture amending'
run r10_lattice_attained 'lonely runner supremum attained rational time speed sum structure optimal'
run r13_mirror_pair 'lonely runner mirror speeds sum time lattice pole'
run r14_gap_formulation 'lonely runner maximum gap circle formulation equivalent speeds'
run r15_cusick_covering 'Cusick view obstruction covering multiplicative discrete'
run r16_regular_time 'lonely runner regular speeds tight bound equality characterization'
# extra reviewer-directed phrasings
run r17_scholar_denom 'optimal time denominator equals speed sum lonely runner'
run r18_cyclic_reduction 'lonely runner finite cyclic group reduction Z_N'
run r19_pair_structure 'lonely runner pair of speeds denominator optimum proof'
run r20_four_exact 'lonely runner four runners exact lonely time formula Cusick'
echo "DONE"
