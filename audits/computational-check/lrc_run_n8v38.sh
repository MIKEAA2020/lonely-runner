#!/bin/bash
# lrc_run_n8v38.sh -- run all 16 shards of the n=8, V=38 sweep, two at a
# time (2 cores), fully detached. C(38,8) = 48,903,492 vectors.
# Progress marker: out/n8_V38_progress.txt
cd /home/z/my-project/scripts
mkdir -p out
: > out/n8_V38_progress.txt
for k in 0 1 2 3 4 5 6 7; do
  for i in $((2*k)) $((2*k+1)); do
    ./lrc_ilp_check_v3 shard 8 38 out/n8_V38_s$i 16 $i nodump \
        > out/n8_V38_s$i.out 2> out/n8_V38_s$i.log &
  done
  wait
  echo "pair $k done $(date +%H:%M:%S)" >> out/n8_V38_progress.txt
done
echo "ALL DONE $(date +%H:%M:%S)" >> out/n8_V38_progress.txt
