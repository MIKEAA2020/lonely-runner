# N9 — The V=56 data point and the independent equality verification

Reviewer's pre-submission checklist, executed. Two items: (1) the optional
V=56 run ("the single experiment that decides between 'stable coverage
with growing classification burden' and 'fundamental decay'"), and (5) an
independent verification of the equality theorem ("a from-scratch
implementation, different code path, on a random subset of sets").

## 1. The V=56 run (`scripts/N9_v56_extension.py`, log `scripts/n9_run.log`)

Battery byte-identical to the committed N6/N7/N8 (verified: the function
block from `RIGID_PRIMES` through `report()` is character-identical to
N8; only ranges extended: tables/capability/sigcap to N ≤ 111, corpus to
V = 56). This environment reaps background processes between tool calls,
so the run is staged (six foreground stages, per-corpus state persisted
to `scripts/n9_state.json`, hard controls re-asserted at every stage).

HARD CONTROLS — all PASS exactly:
- n=4: 1651/1661/1744 of 1745, OPEN 1, open set (3,5,8,13).
- same battery: V=16 4179/4193/4311; V=24 35778/36362/41240/41656
  (unflagged 62,089 = 19,852 P4-target + 42,237 kernel); V=32
  156461/161596/196751; V=40 480454/500280/641166 (unflagged 534,098);
  V=48 1200777/1253258/1665356 (unflagged 1,041,235 = 358,745 P4-target
  + 532,700 + 137,484 + 12,306).
- augmented: V=24 38611/38924/41656; V=32 173337/176376/196751; V=40
  526607/540415/641166 (unflagged 338,126); V=48 1295478/1335857/1665356
  (unflagged 682,490 = 532,700 + 137,484 + 12,306).
- capability: zones {7,11,13} (n=4) and {7,13,17,19,37} with counts
  20/68/16/64/32; 48..95 empty (committed N8); sigcap 48..79 (N6) and
  80..95 (N8) dictionaries exact.
- V=56 corpus size: 3,712,576 primitive 5-sets (Mobius count; the same
  formula reproduces the committed V=48 total exactly).

RIGID ZONE: no capable modulus anywhere in 96..111; primes 97, 101, 103,
107, 109 safe (inside the committed (47,150] scan); new odd composites
99 (19 non-unit classes, 18 reachable), 105 (28/25), 111 (19/18) all
PROVED safe by the reachability mechanism — the micro-case family stays
{49, 77, 91}. Verified composite range: 48..111.

MEASUREMENT (n=5, v ≤ 56, 3,712,576 sets):
- same battery: ports 2,640,683 = 71.13%; +D3 2,757,665 = 74.28%;
  +V 100%; OPEN 0; ground truth (TAU-5) 3,712,576/3,712,576.
- augmented battery: ports 2,805,505 = 75.57%; +D3 2,901,615 = 78.16%;
  +V 100%; OPEN 0.
- unflagged decomposition (augmented): 1,205,429 kernel-world pairs =
  861,988 @ N≤47 + 267,331 @ 48..79 + 65,626 @ 80..95 + 10,484 @
  96..111; self-check unflagged allinv @19/37 == 0 (residual purely
  kernel-world, as at every earlier augmented range).
- soundness: P1 witness fails 0; m4 mismatches 0; classification
  exactness 0; flag=>cert sweeps: 96..111 random 0/64,000; N=99
  exhaustive 0/4,249,575; N=111 exhaustive 0/6,672,876.

THE HONEST READING (reviewer's framing, enforced):
- augmented +D3 series: 93.44 → 89.64 → 84.29 → 80.21 → 78.16;
  residual 6.56 → 10.36 → 15.71 → 19.79 → 21.84; increments
  +3.80, +5.36, +4.07, +2.06 (mean 3.82 pts per 8 V) — roughly linear
  growth with fluctuation; five points cannot distinguish deceleration
  from a noisy linear trend.
- Level fact: V=48 was a marginal pass of the 80% continuation bar
  (+0.21 pts); V=56 falls BELOW it (−1.84 pts). Proved coverage did not
  stabilize in the tested range.
- Within the tested range the reformulation remains productive (OPEN 0,
  +V 100%, ground truth exact, rigid zone constant); the residual is
  entirely kernel-world; closing it requires the kernel-world structure
  theorem — the single open mathematical problem in the reformulation.
- Per the reviewer's dichotomy: the reading is "productive in the tested
  range with fundamental decay / unresolved asymptotics", NOT "stable
  coverage with growing classification burden". The 60% pivot level is
  not near (mean-increment extrapolation crosses it around V ≈ 96); no
  pivot to R2 is indicated by the level rule (60 ≤ 78.16 < 80).

## 2. Independent equality verification (`scripts/EQ_independent_verify.py`,
log `scripts/eq_independent_run.log`)

From-scratch reimplementation, different code path (time domain, exact
integer arithmetic, written against the theorem statement without
reference to the pipeline or to V1):
- Route A (unrestricted folklore candidates): max of g over breakpoints
  j/(2v_p) ∪ ALL pairwise crossings c/(v_p ± v_q) — differences included.
- Route B: max over pair-sum lattices.
- C1 equality: 2,870/2,870 (100.000%) across n=2..6 (300–1,200 random
  primitive sets per rung; 1,200 in the paper's domain v ≤ 56; a
  beyond-corpus probe at v ≤ 200).
- C2 LRC bound: 2,870/2,870.
- Anchors exact: (1,2)→1/3, (1,2,3)→1/4, (1,2,3,4)→1/5, (1,…,5)→1/6,
  (1,…,6)→1/7, (1,2,3,5)→1/4, (1,3,4,7)→1/5, (1,3,4,5,9)→1/6; and the
  n=4 open set (3,5,8,13) has exact M = 5/21 (new committed datum);
  (2,6,8,10,11) has M = 2/9.
- Honest nuance: on 962/2,870 sets a difference time ALSO attains the
  optimum (never needed) — consistent with Lemma 2's sufficiency form;
  recorded in the paper (§2.5) to keep the sharpening's phrasing exact.

## 3. Paper updates delivered

- Title: "The Lonely Runner **Problem** on Pair-Sum Lattices: ..."
  (cover, metadata; the word "Conjecture" no longer leads).
- Abstract: reformulation still first; V=56 numbers; decline framing
  ("does not establish stabilization; asymptotics unresolved").
- §5.4 rewritten as the honest reading (residual series, increments,
  below-threshold level fact, kernel-world-only residual).
- Rifford paragraph (§1.1): explicit shared/different statement.
- §2.5: the independent verification (2,870 sets) + the difference-time
  nuance.
- §4: composite range ≤ 111 (table +3 rows), successor-conjecture range.
- §6: decay open problem updated (below-bar fact); new open problem 9
  (bridges: additive combinatorics on multiplicative arcs; zonotope-
  matroid / log-concave polynomial encoding).
- §7: soundness rows extended (2,870-set independent check; exhaustive
  99/111; ground truth v ≤ 56; N9 + EQ artifacts).
- Figures regenerated with the V=56 point and the honest annotations.

Build record: tectonic 2 runs, ZERO TeX warnings; cover validated
(poster_validate + cover_validate PASS) and rendered via html2poster.js;
merged with page-size normalization; metadata set; pdf_qa residual
warning = the documented TOC false positive (31 link annotations
verified present via pypdf); font.check findings identical to the prior
build's visually-refuted set (8 notdef + 24 control chars, inherited).
