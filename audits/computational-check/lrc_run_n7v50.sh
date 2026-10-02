#!/bin/bash
# lrc_run_n7v50.sh <k> -- run shard pair (2k, 2k+1) of the n=7, V=50 sweep.
# 16 shards total, C(50,7) = 99,884,400 vectors, 6,242,775 per shard.
# Two shards run in parallel (2 cores); each call stays well under the
# 10-minute sandbox cap.
set -e
cd /home/z/my-project/scripts
k=$1
for i in $((2*k)) $((2*k+1)); do
  ./lrc_ilp_check shard 7 50 out/n7_V50_s$i 16 $i \
      > out/n7_V50_s$i.out 2> out/n7_V50_s$i.log &
done
wait
cat out/n7_V50_s$((2*k))_summary.txt out/n7_V50_s$((2*k+1))_summary.txt
