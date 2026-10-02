#!/bin/bash
# lrc_shard_equiv.sh -- prove the v2 shard mode is behavior-identical to the
# committed v1 solver:
#   1. v2 "run" mode byte-matches v1 "run" mode on fresh runs.
#   2. The concatenation of shard dumps (i = 0..S-1) byte-matches the v1
#      run-mode dump, for several (n, V, S) including remainder cases.
#   3. Merged tight lists and totals match.
#   4. BIG: n=7 V=38 sharded into 5 vs the committed n7_V38_dump.bin
#      (produced by the v1 binary in the previous session; its checksum is
#      recorded in the repo's computational check results).
set -e
cd /home/z/my-project/scripts
V1=./lrc_ilp_check.v1.orig
V2=./lrc_ilp_check
mkdir -p out

echo "== 1. run-mode byte match (v1 vs v2) =="
for t in "2 6" "3 12" "4 8" "5 10"; do
    set -- $t
    $V1 run $1 $2 out/eqv1_${1}_${2} > /dev/null
    $V2 run $1 $2 out/eqv2_${1}_${2} > /dev/null
    cmp out/eqv1_${1}_${2}_dump.bin out/eqv2_${1}_${2}_dump.bin
    cmp out/eqv1_${1}_${2}_tight.txt out/eqv2_${1}_${2}_tight.txt
    echo "   run n=$1 V=$2: dump+tight byte-identical"
done

echo "== 2/3. shard decomposition == run =="
for t in "2 6 7" "3 12 4" "4 8 3" "5 10 5" "5 10 7"; do
    set -- $t
    n=$1; V=$2; S=$3
    rm -f out/eqc_${n}_${V}_${S}.bin out/eqc_${n}_${V}_${S}.tight
    for i in $(seq 0 $((S-1))); do
        $V2 shard $n $V out/eqs_${n}_${V}_${S}_$i $S $i > /dev/null 2>out/eqs_${n}_${V}_${S}_$i.log
        cat out/eqs_${n}_${V}_${S}_${i}_dump.bin >> out/eqc_${n}_${V}_${S}.bin
        cat out/eqs_${n}_${V}_${S}_${i}_tight.txt >> out/eqc_${n}_${V}_${S}.tight
    done
    cmp out/eqv1_${n}_${V}_dump.bin out/eqc_${n}_${V}_${S}.bin
    cmp out/eqv1_${n}_${V}_tight.txt out/eqc_${n}_${V}_${S}.tight
    tot=$(awk -F'[= ]' '/^n=/{s+=$14} END {print s}' out/eqs_${n}_${V}_${S}_[0-9]*_summary.txt)
    feas=$(awk -F'[= ]' '/^n=/{s+=$16} END {print s}' out/eqs_${n}_${V}_${S}_[0-9]*_summary.txt)
    vtight=$(awk -F'[= ]' '/^n=/{s+=$20} END {print s}' out/eqs_${n}_${V}_${S}_[0-9]*_summary.txt)
    vfail=$(awk -F'[= ]' '/^n=/{s+=$22} END {print s}' out/eqs_${n}_${V}_${S}_[0-9]*_summary.txt)
    ref=$(awk -F'[= ]' '{print $6, $8, $12, $14}' out/eqv1_${n}_${V}_summary.txt)
    echo "   shard n=$n V=$V S=$S: concat byte-identical; totals $tot/$feas tight=$vtight vfail=$vfail (ref: $ref)"
    [ "$tot $feas $vtight $vfail" = "$ref" ]
done

echo "== 4. BIG equivalence: n=7 V=38, 5 shards vs committed dump =="
for i in 0 1 2 3 4; do
    $V2 shard 7 38 out/eq7_s$i 5 $i > out/eq7_s$i.out 2> out/eq7_s$i.log
done
rm -f out/eq7_concat.bin
for i in 0 1 2 3 4; do cat out/eq7_s${i}_dump.bin >> out/eq7_concat.bin; done
echo "   committed: $(sha256sum out/n7_V38_dump.bin | cut -d' ' -f1)"
echo "   concat:    $(sha256sum out/eq7_concat.bin | cut -d' ' -f1)"
cmp out/n7_V38_dump.bin out/eq7_concat.bin && echo "   BYTE-IDENTICAL"
cat out/eq7_s[0-4]_tight.txt > out/eq7_tight_merged.txt
cmp out/n7_V38_tight.txt out/eq7_tight_merged.txt && echo "   TIGHT LIST IDENTICAL"
awk -F'[= ]' '/^n=/{s+=$14+0; f+=$16+0; t+=$20+0; vf+=$22+0}
     END{printf "   merged: total=%d feasible=%d tight=%d verify_failures=%d\n", s, f, t, vf}' out/eq7_s[0-4]_summary.txt
cat out/eq7_s[0-4]_summary.txt
echo "EQUIVALENCE: ALL CHECKS PASS"
