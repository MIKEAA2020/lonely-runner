#!/bin/bash
# Tight-zoo hunts: n=2..6 @ V=50 with theta = 1/(n+1) (tight value).
# v4 probe filter: skips sets with gap > theta; tight sets never exceed theta.
cd /home/z/my-project/scripts
mkdir -p out
for n in 2 3 4 5 6; do
  V=50
  theta_den=$((n+1))
  echo "== n=$n V=$V theta=1/$theta_den =="
  ./lrc_ilp_check_v4 hunt $n $V out/tz${n} 1 0 1 $theta_den nodump 2>&1 | tail -2
  echo "tight zoo lines: $(wc -l < out/tz${n}_tight.txt)"
done
echo DONE
