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

=== 2026-10-02: tight-set decomposition via X(S) (the rung-2 spine) ===
audit: audits/"tight-set decomposition via X(S).txt"
lrc_tight_hunt.sh     zoo hunts: tight n=2..6 @ theta=1/(n+1); tight+rung2+
                      rung3 at m=3..6 @ theta=3/(3m+1) (v4, seconds each).
tight_n{2..6}_V50.txt complete tight zoos n<=6 (match committed censuses).
rz{3..6}_{tight,rung2,rung3}.txt, rz6_interrung.txt  level-m zoo layers.
lrc_tight_decomp.py   T1 theorem-vs-brute on all 463 (T,x) pairs of all
                      tight sets n=2..8 (0 failures); T4 X-closure at every
                      core (zoo extensions == X(S,1/(n+1)) cut to [1,50]);
                      T3 naive-ladder test (fails exactly on the mod
                      families); T5 core spectrum -> tight_decomp.json.
lrc_tight_atlas.py    canonical min-gap X-chains (every node a zoo member),
                      seed taxonomy, productivity census, 13/13 predictions
                      beyond V=50, 324/324 scaling covariance ->
                      tight_atlas.json.
lrc_tight_spines.py   spine verification + forward chain presentation.
tight_decomp_output.txt / tight_atlas_output.txt / tight_spines_output.txt
Headline: the mod families of the tight zoo are seeded by rung-2 ladders
({1}->{1,4}->{1,3,4}->{1,3,4,5}->{1,3,4,5,9}; [5]+12=+rung2@6->(1,2,3,4,5,7,12);
{1}->{1,6}->{1,5,6}->{1,4,5,6}->{1,4,5,6,7}->{1,4,5,6,7,11}->tight@7);
n=6 and n=8 have no mods because no rung-2@(n-1) member (V<=50) has a
nonempty X(S,1/(n+1)).

Reproduction (2 cores):
  gcc -O2 -o lrc_ilp_check_v4 lrc_ilp_check_v4.c
  bash lrc_tight_hunt.sh
  gcc -O2 -o lrc_xs_batch lrc_xs_batch.c
  python3 lrc_tight_decomp.py && python3 lrc_tight_atlas.py && python3 lrc_tight_spines.py

=== 2026-10-02: silent killers and the induction step (the NSK census) ===
audit: audits/"silent killers and the induction step.txt"
nsk_census/            self-contained subdirectory:
  lrc_nsk_census.c     the no-silent-kill census.  Per (m-1)-core S of
                       [1..B] at level g=1/(m+1): probe pass on the
                       crossing-grid midpoints of s_max with the
                       killer-intersection certificate (if no x in
                       [1,R_found]\S kills all FOUND escapes then
                       D(S,g) is empty -- sound because D <= R <=
                       R_found and D must kill every escape); exact
                       event-scan fallback for near-tight cores; D/X
                       by the General Filler Theorem trichotomy for
                       x in [1,min(R,cap)].  check mode prints R/D/X
                       for arbitrary rational g (the Python port's
                       twin).  Flags: ON (gap=g), SUB (gap<g), RFLAG
                       (R>cap), SILENT.
  lrc_nsk_theory.py    NT1 trichotomy on 2865 (S,g,x) triples; NT2
                       every-order decomposition of all 443 (T,x)
                       pairs of every tight set n=2..8; NT3 reduction
                       equivalence at m=2..5 + negative control
                       (S={1}, g=2/5, x=2 silent); NT4 scaling
                       containment + level-2 strictness errata.
  lrc_nsk_validate.py  V-A 480/480 C-vs-Python; V-B 68 trichotomy
                       pairs vs gap_int; V-C m=4 B=14 census ==
                       Python enumeration.
  lrc_nsk_analyze.py   tier tallies (accounting C(B,m-1) exact),
                       productive-core classification (ladder /
                       mod-seeds / every-order others), ALL 2137
                       X-fillers re-verified exact, 120 negative
                       samples > g -> nsk_census_report.json.
  cores_m{2..10}_*.txt merged per-tier productive cores (R, D, X).
  summary_m{2..10}_*.txt per-shard summaries.
  nsk_census_report.json  machine-readable totals and maps.

Headline: 1,709,221,699 cores across m=2..10 (B=500/300/200/150/120/
80/60/50/40), SILENT=0, D=X at every core, no flags.  By Theorem R:
LRC certified for every m-set whose second-largest speed is <= B_m --
including m=9 (B=50: independent computational check of the announced
nine-speed case) and m=10 (B=40: the open frontier; 273,438,880
nine-speed cores, self-contained -- any 9-core at/below 1/11 would
have been flagged).  Productive maps: pure ladder family at
m=6,8,9,10 (rung-2 barrenness extended to B=120/60/50/40; tight zoos
at n=9,10 predicted pure scalings within these universes); ladder +
mod seeds + their every-order cores at m=4,5,7 (B=200/150/80,
reconciling the Task-3 productivity census exactly).

Reproduction (2 cores, ~70 min):
  gcc -O2 -o lrc_nsk_census nsk_census/lrc_nsk_census.c
  python3 nsk_census/lrc_nsk_theory.py && python3 nsk_census/lrc_nsk_validate.py
  for spec in "2 500 1" "3 300 1" "4 200 1" "5 150 1" "6 120 2" "7 80 4" \
              "8 60 6" "9 50 8" "10 40 6"; do set -- $spec; \
    for i in $(seq 0 $(($3-1))); do \
      ./lrc_nsk_census $1 $2 out/nsk_m$1_s$i $3 $i 3000 & done; wait; done
  python3 nsk_census/lrc_nsk_analyze.py

=== DROP EXPERIMENT (the drop question) ===
Files: ../"the drop question.txt" (audit), drop_experiment/lrc_drop_experiment.py,
drop_experiment/lrc_drop_analyze.py, drop_experiment/drop_records.jsonl (17,784
pairs), drop_experiment/drop_analysis.txt.

The advisor's drop question (kill-not-clean at level g = gap(S) = no-touch
droppers; complementary to the NSK fixed-level census).  Theorem D (drop
formula): gap(S+x) = max(Q*, C*) -- Q* = S's own critical ladder cleared by
x (q-corner landings, machine-equal), C* = best crossing value (mixed
arithmetic); dichotomy crossing/q-corner proved and machine-verified on
every pair; residue law (binders u,w at t=a/b: u = +-w mod b) checked
28,390/28,390.  Findings: margin >= 0 on all 8,922 drops (SILENT=0, the
induction step of LRC certified in-universe); 105 floor drops = tight
extensions; drops QUANTIZED -- min positive margin = rung-2 margin
1/((2m+1)(m+1)) at aug m in {2,3,4,5,6,8}, inter-rung window EMPTY there;
inhabited at aug=7 by {1,2,3,4,5,7}+18 -> 3/23 (the committed census
vector, now with mechanism: tight core + overshooting filler).  LRC residue
isolated: the C-half (crossing values >= 1/(k+2)).  Plus the second-language
(kernel/covering) dictionary as the bridge prospectus.

Reproduction (~6 min, pure Python):
  python3 drop_experiment/lrc_drop_experiment.py --stage count
  for i in 0 1 2 3; do python3 drop_experiment/lrc_drop_experiment.py \
      --stage complete --shard $i/4; done
  for i in 0 1; do python3 drop_experiment/lrc_drop_experiment.py \
      --stage zoo --shard $i/2; done
  python3 drop_experiment/lrc_drop_analyze.py
