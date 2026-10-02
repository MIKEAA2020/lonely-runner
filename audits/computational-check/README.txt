Computational check of the Lonely Runner polytope integrality (n <= 5, v_n <= 50)
Extended to n = 6 (v_n <= 50) and n = 7 (v_n <= 38)
Extended to n = 7 (v_n <= 50, sharded) and Theorem 6 machine-checked
Extended to the rung-2 ladder: R1 proved and machine-checked, zoos at
n = 5/6 (V = 50) and n = 8 (V = 30 and V = 38, sharded) -- see
"../rung-2 ladder and R1.txt"
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
  lrc_ilp_check_v3.c
                    v3: identical enumeration, argmax and exact gap
                    (solve() untouched; dumps and tight lists byte-
                    identical to v2), plus two census outputs in run and
                    shard modes: <prefix>_rung2.txt (gap exactly
                    2/(2n+1), with recorded argmax and binding speeds)
                    and <prefix>_interrung.txt (1/(n+1) < gap <
                    2/(2n+1)). Regression: n = 7, V = 20 and V = 38
                    reproduce the committed censuses exactly (tight
                    4/10, rung2 8/21, interrung 2/4).
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
  lrc_gap_lib.py            exact gap machinery, two independent
                    implementations (integer arithmetic mirroring the C
                    solver + Fraction reference), all-argmax extraction.
  lrc_thm_r1_check.py       Theorem 1 / Theorem 2 (R1) machine check,
                    n = 2..40: gap value, unique argmax pair, binders
                    {1, 2n}, slack >= 3/(2n+1), equality configuration.
  lrc_rung2_classify.py     n = 7 rung-2 classification: Lemma A
                    verification on all argmaxes of the 29 vectors,
                    Schur census, minimal cores ->
                    rung2_n7_classification.json.
  lrc_core_extension.py     hub-core decomposition and extension lists
                    X(S) for every 6-subset of every zoo member.
  lrc_n56_check.py          Lemma A verification for zoo_5 / zoo_6.
  lrc_rung_family_check.py  Theorem 1' (rung family gap(1..n-1,kn) =
                    k/(kn+1)) check, n = 2..12, k = 2..6.
  lrc_n8_analyze.py         n = 8 zoo analysis layers V1-V5 (shard
                    reconciliation, reference re-verification, structure,
                    ladder/hub decomposition, spectrum) ->
                    rung2_n8_classification.json.
  lrc_n8_validate.py        n = 8, V = 38 sweep validation: S1 spectrum
                    sums to C(38,8); S2 496/496 sampled records match the
                    Python reference; S3 census sizes and values.
  lrc_ladder_table.py       ladder tables at V = 38, cross-checks vs the
                    committed censuses, single-speed lifts 5->6->7->8.
  lrc_run_n8v38.sh          sharded n = 8, V = 38 driver (8 shards, two at
                    a time; adapt to a foreground chunk runner if the
                    sandbox reaps background processes).
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
  n5_V50_rung2.txt  zoo_5: all 32 vectors with gap exactly 2/11 at n = 5,
                    v_n <= 50 (6 primitive classes x scalars).
  n6_V50_rung2.txt  zoo_6: all 32 vectors with gap exactly 2/13 at n = 6,
                    v_n <= 50 (8 primitive classes x scalars).
  n8_V30_rung2.txt  n = 8, V = 30: all 10 vectors with gap exactly 2/17.
  n8_V38_rung2.txt  n = 8, V = 38: all 14 vectors with gap exactly 2/17
                    (7 primitive classes x scalars; binding pairs mod 17).
  n8_V38_tight.txt  n = 8, V = 38: all 4 tight vectors (gap 1/9) --
                    exactly {1,...,8} x {1,2,3,4} (tight zoo closed).
  rung2_n7_classification.json / rung2_n8_classification.json
                    machine-readable classifications (argmaxes, binders,
                    residues, Schur triples, cores).

Reproducing the rung-2 ladder extension (~25 min on 2 cores):
  gcc -O2 -o lrc_ilp_check_v3 lrc_ilp_check_v3.c
  ./lrc_ilp_check_v3 run 5 50 out/n5_V50
  for i in 0 1 2 3; do ./lrc_ilp_check_v3 shard 6 50 out/n6_V50_s$i 4 $i nodump & done; wait
  ./lrc_ilp_check_v3 run 8 30 out/n8_V30
  ./lrc_run_n8v38.sh
  python3 lrc_thm_r1_check.py
  python3 lrc_rung2_classify.py
  python3 lrc_core_extension.py
  python3 lrc_n56_check.py
  python3 lrc_rung_family_check.py
  python3 lrc_n8_analyze.py
  python3 lrc_n8_validate.py
  python3 lrc_ladder_table.py

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

=== rung family note additions (Oct 2026) ===

lrc_ilp_check_v4.c   hunt mode: probe-lemma filter (gap >= value at any
                     candidate time; vectors certified gap > theta = 2/15
                     are skipped -- exact census below theta preserved),
                     low-gap census file, rung-k zoo emission.
                     solve()/enumeration/records unchanged from v3.
                     Validation: v3 vs v4 censuses byte-identical at
                     V = 38 (full range + shard level) and at V = 50
                     (sub-range, 105x speedup).
lrc_xs_batch.c       batch filler brute force for lrc_xs_theory.py
                     (solve() copied verbatim).
lrc_xs_theory.py     THE closed-form filler characterization X(S):
                     kill (window condition per escape interval) AND
                     touch (clean at a level time), radius bound
                     R = floor(2g/L_max).  Verified: 612 zoo cores
                     (n=5..8) + 480 random cores + degenerate regimes,
                     X_theory == X_brute everywhere (lrc_xs_batch.c,
                     gap_int, gap_frac triple-checked).
lrc_n8v50_analyze.py analysis of the n=8 V=50 hunt: H1-H8 (accounting,
                     5384/5384 sample validation, zoos, Lemma A,
                     inter-rung emptiness, X(core) ladder, spectrum).
lrc_rung_family_check.py  extended: n=2..16, k=1..12 -- value, EXACT
                     argmax set (unique pair for k>=2, units grid for
                     k=1), binders {1,kn}, slack >= (k+1)/(kn+1).
note_rung_family.tex source of audits/"rung family note.pdf".
n8_V50_*.txt         census files of the n=8, V=50 hunt
                     (C(50,8) = 536,878,650 vectors): tight 6
                     (consecutive x 1..6), rung2 21 (same 7 primitive
                     classes as V=38, scalars to 4), interrung 0 (EMPTY),
                     rung3/4/5/6 = 11/10/5/3, lowgap 2806 (gap <= 2/15).
n8_V50_analysis.txt  full hunt report; all checks pass.
lrc_xs_theory_output.txt  filler theorem verification output.

Reproduction (2 cores):
  gcc -O2 -o lrc_ilp_check_v4 lrc_ilp_check_v4.c
  for i in 0 1 2 3 4 5 6 7; do ./lrc_ilp_check_v4 hunt 8 50 out/n8_V50_s$i 8 $i 2 15 nodump; done
  python3 lrc_n8v50_analyze.py
  gcc -O2 -o lrc_xs_batch lrc_xs_batch.c && python3 lrc_xs_theory.py
  python3 lrc_rung_family_check.py
