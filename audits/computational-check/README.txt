Computational check of the Lonely Runner polytope integrality (n <= 5, v_n <= 50)
Extended to n = 6 (v_n <= 50) and n = 7 (v_n <= 38)
Extended to n = 7 (v_n <= 50, sharded) and Theorem 6 machine-checked
================================================================================

Contents
  lrc_ilp_check.c   exhaustive solver (exact integer arithmetic), three modes:
                      run   n V prefix    full sweep, writes dump + tight list
                                        + summary
                      shard n V prefix S I [nodump]
                                        the vectors whose lexicographic rank
                                        lies in the I-th of S contiguous
                                        intervals partitioning C(V,n); same
                                        records, same order as run (the
                                        concatenation of shard dumps is
                                        byte-identical to the run dump --
                                        proven on n = 7, V = 38 against the
                                        committed output). Also writes a
                                        reduced gap histogram
                                        <prefix>_spectrum.txt and every
                                        99733rd record to
                                        <prefix>_sample.txt.
                      brute n            literal k-box enumeration, stdin
                    Build: gcc -O2 -o lrc_ilp_check lrc_ilp_check.c
                    (v2 adds shard mode; run/brute are behavior-identical to
                    the v1 binary -- verified byte-for-byte.)
  lrc_verify.py     validation driver for the n <= 5 check (see "computational
                    check results.txt" sections 3-8).
  lrc_extend67.py   extension driver for n = 6 / n = 7 (V = 38): witness
                    re-verification, reference solver, brute, MILP, censuses
                    -> stats67.json.
  lrc_theory_experiments.py
                    experiments behind "integrality analysis.txt" (E1-E4)
                    -> theory_experiments.json.
  lrc_gap_family_check.py
                    machine verification of Theorem 6
                    (gap(1,2,...,n-1,2n) = 2/(2n+1) for every n >= 2) and of
                    its two lemmas: L1 exact family gap on the full candidate
                    set (n <= 200 exhaustive, samples to 1000), L2 residue
                    kill-lemma (n <= 250 exhaustive, samples to 800), L3
                    committed-solver cross-check, L4 pure-Python reference,
                    L5 dense grid -> thm6_family_check.jsonl.
  lrc_shard_equiv.sh
                    byte-equivalence proofs: v2 run mode vs the v1 binary,
                    and shard decomposition vs run (incl. n = 7, V = 38 vs
                    the committed dump) -> equiv_console.txt.
  lrc_run_n7v50.sh  runs one pair of the 16 shards of the n = 7, V = 50
                    sweep (two processes, 2 cores, well under any 10-minute
                    per-process limit).
  lrc_n7v50_validate.py
                    validation layers A-F for the V = 50 sweep: A full numpy
                    re-verification of every dump record against all polytope
                    constraints + rank geometry; B aggregate/spectrum/census
                    checks; C independent reference on all sampled and census
                    vectors; D HiGHS MILP; E Theorem 6 anchors; F sample/dump
                    consistency -> n7v50_validation.jsonl.
  stats.json        machine-readable results, n <= 5.
  stats67.json      machine-readable results, n = 6 / n = 7 (V = 38).
  stats_n7v50.json  machine-readable results, n = 7, V = 50 (per-shard times,
                    spectrum bottom, V38->V50 reconciliation, tight /
                    inter-rung / rung-2 censuses).
  theory_experiments.json   machine-readable E1-E4 data.
  thm6_family_check.jsonl   Theorem 6 layer results (L1-L5).
  n7v50_validation.jsonl    V = 50 validation layer results (B-F).
  equiv_console.txt         shard-equivalence console output.
  tight_n*.txt      complete tight-set census (gap = 1/(n+1)) with witnesses:
                    n = 1..5 (v_n <= 50), n = 6 (v_n <= 50), n = 7 (v_n <= 38)
                    and n = 7 (v_n <= 50: tight_n7_V50.txt, 14 vectors).
  tight_n7_V50_rung2.txt    all 29 vectors with gap exactly 2/15 at n = 7,
                    v_n <= 50 (the rung-2 population, Theorem 6 family incl.).
  inter_rung_n7_V50.txt     all 4 vectors with gap strictly in (1/8, 2/15)
                    at n = 7, v_n <= 50 (the known 3/23 counterexamples).
  gap_spectrum_n7_V50.txt   complete reduced gap-value histogram of the
                    99,884,400 vectors (930 distinct values).

Reproducing the full sweep (about 30 seconds, single core, n <= 5):
  gcc -O2 -o lrc_ilp_check lrc_ilp_check.c
  mkdir -p out
  for n in 1 2 3 4 5; do ./lrc_ilp_check run $n 50 out/n${n}_V50; done
  python3 lrc_verify.py

Reproducing the n = 6 / n = 7 (V = 38) extension (about 11 minutes):
  ./lrc_ilp_check run 6 50 out/n6_V50
  ./lrc_ilp_check run 7 38 out/n7_V38
  python3 lrc_extend67.py
  python3 lrc_theory_experiments.py

Reproducing the n = 7, V = 50 sharded sweep (~56 min CPU, ~30 min wall on
2 cores; each call < 10 min) and its validation:
  ./lrc_shard_equiv.sh                      # byte-equivalence proofs
  for k in 0 1 2 3 4 5 6 7; do ./lrc_run_n7v50.sh $k; done
  python3 lrc_gap_family_check.py           # Theorem 6 layers L1-L5
  python3 lrc_n7v50_validate.py A 0 8
  python3 lrc_n7v50_validate.py A 8 16
  python3 lrc_n7v50_validate.py B C D E F

Dumps (74 MB + 508 MB + 404 MB + 16 x 200 MB) are not committed; SHA-256 of
every reconstructed full dump is recorded in "computational check results.txt"
(sections 9 and 10). Seeds: 20261002.

Note: the Python drivers as committed read/write absolute paths under
/home/z/my-project/scripts/out; adjust OUT/BIN at the top when running
elsewhere.
