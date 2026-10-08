# The n=5 Leverage Experiment: The Decisive Measurement

**Program:** LRC residue-language reformulation (pair-sum finitization).
**Charge (reviewer recommendation 1, verbatim intent):** at n = 5, enumerate
all sets, run the same cyclic-lemma battery, and measure the fraction closed
by proved lemmas alone. Comparable to the committed n = 4 figure (94.6%) ⟹
the reformulation is genuinely productive; a sharp drop (say to 50%) ⟹ it
closes n = 4 by accident of dimension and the plain regime is where the
difficulty lives. **Decision rule (recommendation 4):** ≥ 80% closure by
proved lemmas ⟹ continue the program; < 60% ⟹ pivot to R2 with the
reformulation as a bridge.
**Committed scripts:** `scripts/N5_leverage_experiment.py` (the experiment),
`scripts/N5_prime_scan_crosscheck.py` (independent verification of the
rigid-zone extension).

---

## 1. Design

**Corpus.** All primitive 5-subsets of [1, 24] — 41,656 sets, the E2b base.
Ground truth re-verified independently this turn: every set has ≥ 1
certifying pair (41,656/41,656; bitmask method, E2b used grid arithmetic) —
(TAU-5) holds on the full base.

**Covering framework (n generic, T = n + 1, k = n − 1).** For a pair (p, q)
with N = v_p + v_q, the k effective residues rs = (v_p, others) mod N
(v_q ≡ −v_p shares v_p's bad set). The pair certifies iff some k has
min_w ||wk||_N ≥ N/T, iff the k bad sets

  B_w = { k ∈ Z_N : T·||wk||_N < N } = { k : ||wk||_N ≤ h },  h = ⌊(N−1)/T⌋

do NOT cover Z_N. N ≤ 47 < 64 bits, so every covering question is a bitmask
OR — the entire corpus is decidable exactly and fast, with the battery
running on top of it.

**Battery (status discipline as at n = 4: PROVED lemmas count toward the
headline; per-modulus exhaustive VERIFICATIONS are reported separately).**

| flag | statement | status |
|---|---|---|
| P1 j-gon | j \| N, 2 ≤ j ≤ T, all rs ≢ 0 mod j ⟹ k = N/j certifies | PROVED (Lemma J port; j ≤ 6 at n = 5) |
| P2 counting | Σ exact sizes − (k−1) < N ⟹ no covering | PROVED (B0/B5 port + Lemma S5 below) |
| P4 classifications | at N ∈ {7, 13, 17} (n = 5): covering iff the class condition holds | PROVED (K7/K13/K17 below) |
| P5 D3 | ±-coincidence among residues + invertible + counting ⟹ certifies | PROVED (new this turn; Lemma D3 below) |
| V1 safe modulus | all residues invertible and N not unit-covering-capable | VERIFIED (exhaustive per modulus) |
| V2 safe signature | the pair's order-signature is never-covering at N | VERIFIED (exhaustive per (N, sig)) |

**Positive control.** The whole pipeline was first run at n = 4 (v ≤ 16,
1,745 sets, T = 5, 3 bad sets) with the n = 4 ports (P4 at {7, 11, 13}).
It reproduces the committed record **exactly**: proved-only 1,651/1,745
(94.61%), full battery leaving exactly one open set, and that set is
V = (3, 5, 8, 13) — the committed open set, identified, not just counted.

## 2. The n = 5 rigid zone (new phenomenon)

**Unit covering capability** (four unit bad sets, T = 6 arcs, exhaustive
scale-normalized enumeration):

| range | capable moduli |
|---|---|
| n = 4, 3 bad sets, N ≤ 31 (control) | {7, 11, 13} — committed result reproduced |
| n = 5, 4 bad sets, N ≤ 47 | **{7, 13, 17, 19, 37}** |
| n = 5, primes 53–150 (targeted extension) | **none** |

Covering-configuration counts: 20 at 7, 68 at 13, 16 at 17, 64 at 19, 32 at
37. The extension scan (mask-based) was cross-checked by an independent
set-based implementation at N ∈ {17, 19, 37, 53, 59, 61, 67, 71, 73, 79}:
counts agree exactly (16/64/32 and 0s). Hand-verified anchor: the residues
(1, 4, 6, 7) on Z_17 cover (B_1 ∪ B_4 ∪ B_6 ∪ B_7 = {0,±1,±2} ∪ {0,±4,±9} ∪
{0,±3,±6} ∪ {0,±5,±10} = Z_17), realizable as V = (1, 4, 6, 7, 16), pair (1, 16).

Three structural facts:

1. **The n = 4 no-covering conjecture does not port.** Primes 17, 19, 37 are
   covering-capable at n = 5. The landscape is not "covering confined to
   {7, 11, 13}"; it is *scattered* — and apparently **terminal**: within the
   verified range the last capable prime is 37. The natural n = 5 conjecture:
   *for prime N ≥ 41, four unit bad sets never cover Z_N* (verified to 150).
   In class language (prime N): four multiplicative translates of
   S_h = {classes of 1, …, h} never cover C_{(N−1)/2} for h ≥ 7 except at
   N = 37 (h = 6). Note the size bound never obstructs (4(2h+1) − 3 ≥ N for
   all these N), so the obstruction is genuinely multiplicative.
2. **No composite modulus ≤ 47 is unit-capable.** Even composites are
   provably safe: all residues odd ⟹ k = N/2 certifies (j-gon, j = 2).
   Odd composites ≤ 47 are safe by a non-unit counting argument (unit bad
   sets reach non-unit elements only through j·w^{-1} with gcd(j, N) > 1,
   j ≤ h, and the non-unit part of Z_N is too large to cover) — hand-checked
   at 9, 15, 21, 25, 35. Beyond the range: open.
3. **Three-set coverings exist at n = 5** — base moduli {7, 8, 13} plus
   exactly their lifts {14, 16, 21, 24, 26, 28, 32, 35, 39, 40, 42}
   (e.g. Z_7 residues (1, 2, 3): B_1 ∪ B_2 ∪ B_3 = {0,±1} ∪ {0,±4} ∪
   {0,±5} = Z_7). The lift law (verified, 2,727 random samples, 0
   violations) explains the composite list completely.

**Consequence — the n = 4 B2 coincidence port is DEAD.** At n = 4,
"two effective residues coincide or are antipodal ⟹ certifies" held
(verified N ≤ 80). At n = 5 it is false: a pair with eff ≡ (1, 1, 2, 3) on
Z_7 fails (its three distinct bad sets cover). The language self-corrected:
Lemma D3 (below) is the sound replacement, and it is *stronger where it
applies* — it contributes 584 set-closures at n = 5 versus 10 at n = 4.

## 3. New proved lemmas

**Lemma S5 (unified exact size formula).** For a residue w ∈ Z_N of additive
order j = N/gcd(w, N),

  |B_w| = (N/j) · (2·⌊(j−1)/T⌋ + 1).

*Proof.* Multiplication by w maps Z_N onto 〈w⌋ = {t·N/j : 0 ≤ t < j}, each
image value having exactly N/j preimages, so |B_w| = (N/j)·|A ∩ 〈w⌋| where
A = {x : ||x||_N ≤ h}. Now ||tN/j||_N = (N/j)·min(t, j−t) ≤ h·(N/j)... more
precisely (N/j)·min(t, j−t) < N/T ⟺ min(t, j−t) < j/T, which holds for
exactly 2⌊(j−1)/T⌋ values of t ∈ {1, …, j−1} (t and j−t paired). With t = 0
this gives |A ∩ 〈w⌋| = 2⌊(j−1)/T⌋ + 1. ∎*

Verified: 0 violations over all N ≤ 47 (T = 6) and N ≤ 31 (T = 5), all
residues. This subsumes the n = 4 Lemma S (orders ≤ 5 give |B_w| = N/j;
invertible gives 2h+1) and makes every counting bound hand-computable from
the orders alone.

**Lemma C5 (class reduction, prime N).** For prime N and invertible w, B_w
is antipodal-saturated, and its nonzero part is the union of the h antipodal
classes {±j·w^{-1} : 1 ≤ j ≤ h}. Hence k unit bad sets cover Z_N iff the
corresponding k class-sets X_w·S_h (X_w = class(w^{-1})) cover all (N−1)/2
nonzero classes. *Proof:* k ∈ B_w ⟺ wk ≡ ±j (j ≤ h) ⟺ k ≡ ±j·w^{-1};
0 ∈ B_w always. ∎ (At composite N the reduction is invalid — unit bad sets
need not reach non-unit elements — and is not used there.)

**Theorem K7.** At N = 7 (h = 1): four unit bad sets cover Z_7 iff their
antipodal classes {±w^{-1}} cover all three nonzero classes of Z_7^*.
*Proof:* Lemma C5; S_1 is a singleton. ∎ [20 covering 4-multisets.]

**Theorem K13 (domino tiling).** At N = 13 (h = 2; class(2) generates
C_6 since 2^6 ≡ −1 mod 13): each B_w covers the class domino {X_w, X_w·g},
g = class(2). Four dominos {t, t+1} cover C_6 iff T ∪ (T+1) = Z_6 (T the
start-exponents) iff the complement of T contains two cyclically adjacent
elements. *Proof:* Lemma C5 plus: p ∉ T ∪ (T+1) ⟺ p, p−1 ∉ T. ∎
[68 covering 4-multisets.]

**Theorem K17 (perfect matching).** At N = 17 (h = 2; c = class(2) has
order 4 since 2^4 ≡ −1 mod 17): each B_w covers the domino {X_w, X_w·c},
and dominos live inside the two cosets of the order-4 subgroup 〈c〉 (each a
4-cycle under ·c). Four dominos cover C_8 iff each coset contains exactly
two of the four starts and the two are non-adjacent in that 4-cycle (i.e.
differ by c²). *Proof:* Lemma C5; within a coset the edges {X, Xc}, {Y, Yc}
cover the 4-cycle iff disjoint iff Y ∉ {X, Xc, Xc^{-1}}. ∎ [16 covering
4-multisets — the sparsity is the matching constraint.]

All three classifications were verified against exhaustive enumeration over
all unit 4-multisets at those moduli (predicate ⟺ mask-covering):
0 mismatches.

**Lemma D3 (coincidence counting; sound replacement of the dead B2).** Let
the k = n−1 effective residues at modulus N span pm antipodal classes with
pm ≤ k − 1 (some two residues equal or antipodal), all invertible. Then the
pair certifies whenever pm(2h+1) − (pm−1) < N. In particular at n = 5:
pm = 2 always certifies (4h + 1 < N), and pm = 3 certifies iff N ≢ 1
(mod 6). *Proof.* B_r = B_{r'} when r' ≡ ±r, so the union has at most pm ≤
k − 1 distinct members, each of size exactly 2h+1 (Lemma S5, invertible).
k′ sets all containing 0 have union of size at most the sum of sizes minus
(k′−1) (each set after the first adds at most |A| − 1 elements, as 0 is
already present); the bound is increasing in k′, so |∪| ≤ pm(2h+1) −
(pm−1). For pm = 2 this is 4h+1 ≤ (2/3)(N−1) + 1 < N; for pm = 3 at T = 6
it is 6h+1, which is < N exactly when N ≢ 1 (mod 6). ∎*

Hand-checked anchor: N = 25, rs = (5, 7, 18, 20): classes {5}, {7, 18},
{5, 20} give pm = 2; B_5 = {0, 5, 10, 15, 20}, B_7 = B_18 = 7-arc of 9
elements; union ≤ 13 < 25; indeed k = 1 gives min(||5||, ||7||) = 5 ≥
⌈25/6⌉. The N ≡ 1 (mod 6) exclusion is sharp: Z_7's (1, 1, 2, 3) is the
falsifier of the naive port.

## 4. The measurement

**Primary corpus (n = 5, v ≤ 24, 41,656 sets):**

| battery | sets closed | fraction |
|---|---|---|
| PROVED ports only (P1 + P2 + P4) | 35,778 | **85.89%** |
| all PROVED (+ P5 = Lemma D3) | 36,362 | **87.29%** |
| + VERIFIED (V1 safe moduli) | 41,240 | 99.00% |
| + VERIFIED (V2 safe signatures) | 41,656 | **100.00%** |
| OPEN | 0 | 0% |

**Matched-range comparison (v ≤ 16 — same modulus range N ≤ 31):**

| battery | n = 4 (1,745 sets) | n = 5 (4,311 sets) |
|---|---|---|
| proved ports only | 94.61% | **96.94%** |
| all proved (+ D3) | 95.19% | **97.26%** |
| + verified | 99.94% (1 open: (3,5,8,13)) | 100.00% |

**Per-flag pair counts (primary corpus):** j-gon j = 2/3/4/5/6:
7,920 / 22,656 / 22,880 / 35,549 / 12,771; P2 counting 35,850;
classifications: K7 1,290, K13 5,522, K17 10,262; D3 21,714; V1 85,648;
V2 288,909.

**The plain regime, measured correctly.** Two different notions must not be
conflated. E2b's "plain 73%" was a *best-pair* census (the single
highest-τ pair lacks j-gon structure). The battery-relevant notion is the
*any-pair* census: no pair of the set has j-gon structure. On that basis:

| | n = 4 (v ≤ 16) | n = 5 (v ≤ 24) |
|---|---|---|
| plain sets (any-pair basis) | 669 (38.4%) | 14,307 (34.3%) |
| plain sets closed by proved lemmas | — | 9,013 (63.0% of plain) |
| plain sets closed only by verification | — | 5,294 (37.0% of plain) |

The plain share does **not** grow from n = 4 to n = 5 on the battery basis —
it shrinks. And the plain regime is not monolithic: proved lemmas beyond
j-gon (counting, classifications, D3) close 63% of it. Every one of the
5,294 sets that need the verified tier is plain — the difficulty at n = 5
lives exactly there, but "there" is 12.7% of the corpus, not 73%.

## 5. Where the n = 5 difficulty lives (precisely)

1. **The kernel world.** Mixed order-signatures at composite moduli — e.g.
   (5, 25, 25, 25) at N = 25, (13, 13, 26, 26) at 26, (7, 21, 21, 21) at 21.
   Counting has slack (an order-d residue has |B| = N/d, large), so uniform
   proofs need structure, not size. This is the same "kernel-arc covering"
   residual already present at n = 4 (the B6 class), now bigger. The
   capable signatures at composite N are lifts and mixed-order
   configurations, fully mapped per modulus by the signature scan.
2. **The h ≥ 3 capable primes: 19 and 37.** All-invertible certifying pairs
   there lack classifications (the segment-covering problem: 4 translates
   of S_h in C_{(N−1)/2}, h = 3 and 6). K7/K13/K17 are the h ≤ 2 cases; 19
   and 37 are the natural next classifications.
3. **The pair-selection law is untouched** — why a certifying pair always
   exists is the (TAU) content. The experiment adds an empirical law on
   top: *every one of the 41,656 sets has a certifying pair sitting at a
   never-covering (N, order-signature) configuration* (that is exactly why
   V2 closes 100%). The per-pair hard core — certifying pairs at capable
   configurations — exists (62,089 pairs corpus-wide) but is never
   load-bearing.

## 6. Verdict

**Favorable. Both proved readings clear the 80% threshold:**

- proved ports only: **85.89%** (n = 4 comparator: 94.61%);
- all proved lemmas: **87.29%**;
- at matched range (v ≤ 16): **96.94%, above n = 4's 94.61%** — the
  strongest single datum against "n = 4 closed by accident of dimension".

Per the reviewer's decision rule: **continue the reformulation program; do
not pivot to R2.** The language did at n = 5 exactly what it did at n = 4:
when a ported tool died (B2), the framework generated its sound replacement
(D3, proved); when the rigid zone changed shape, it produced the new map
(scattered capable primes, terminal at 37 in range) and three new
classifications (K7/K13/K17). The paper (reviewer recommendation 3) stands
as framed: *a lossless finitization of LRC via pair-sum lattices, with a
cyclic covering obstruction that scales with n* — with the folklore
citation (Tao 2015 blog comment) in the introduction, and the frontier
noted honestly (LRC proved through 10 runners; the value of this work is
structural, not new small cases).

**Honest caveats.** (i) The proved-port decay 94.6% → 85.9% over the full
range is real and range-driven: larger moduli bring the h ≥ 3 capable
primes and the kernel world into play; at matched range there is no decay.
(ii) 12.7% of sets close only by per-modulus verification — that gap is the
concrete to-do list (classifications at 19/37; kernel-arc structure).
(iii) The trajectory continues: at n = 6 (T = 7, 5 bad sets) the rigid zone
will be denser and the kernel world larger; the program's method — generate
the classifications the new level demands — is exactly what this turn did
for n = 5.

**Do not attack the no-covering conjecture yet** (reviewer recommendation
2, still in force): the n = 5 successor conjecture ("no unit covering at
primes ≥ 41", verified ≤ 150) is now *stated*, and it is a side result, not
the main line. The main line is the paper (recommendation 3) and then the
n = 5 residual classifications as method validation.

## 7. Soundness record

| check | result |
|---|---|
| Lemma S5 size formula, all N ≤ 47 (T=6), N ≤ 31 (T=5) | 0 violations |
| K7/K13/K17 (and n=4 7/11/13) predicate ⟺ mask-covering, exhaustive over unit multisets | 0 mismatches |
| flag ⟹ certification sweep, exhaustive N ≤ 31 | 0 / 46,371 (n=4), 0 / 324,626 (n=5) |
| flag ⟹ certification sweep, random N ∈ 32–47 | 0 / 64,000 |
| P1 constructive witness (bit N/j clear in the union) | 0 failures |
| direct grid-max cross-check (independent m4* code) | 0 mismatches |
| lift law, random lifts | 0 / 2,727 |
| ground truth (every set has a certifying pair) | 41,656 / 41,656 |
| n = 4 positive control vs committed record | exact (1651; open set = (3,5,8,13)) |
| rigid-zone extension, independent set-based scan | agrees exactly at 17/19/37/53–79 |
| Z_17 covering example, D3 example | hand-verified (in §2, §3) |

**Four bugs caught this turn** (the discipline record continues):
(1) `capable_units` enumerated one bad set too many in both branches —
caught because the n = 4 control's capable list contradicted the committed
{7, 11, 13}; (2) a flag-string slice off-by-one (crash); (3) the V1 flag
leaked into the "proved + P5" tier via a loose `if fl` — caught because
"100%" appeared where 87% was expected; (4) the open-pair diagnostic
accumulated corpus-wide instead of open-tier-only — caught because 180,101
exceeded the 4,160 arithmetic maximum. All four were fixed before any
number in this note was recorded.

## 8. Committed artifacts and next steps

- `scripts/N5_leverage_experiment.py` — the full experiment (tables, rigid
  zone, battery, control, sweeps, corpora, diagnostics, prime extension).
- `scripts/N5_prime_scan_crosscheck.py` — independent verification of the
  rigid-zone extension.
- This note: `download/n5_leverage_experiment.md`.

**Next steps (in order):** (1) write the paper per recommendation 3
(framing fixed this turn; folklore citation; open problems: the
capable-prime characterization, the pair-selection law, kernel-arc
coverings); (2) as method validation, the K19/K37 classifications (the
h ≥ 3 segment problem) and the kernel-world structure — these are also the
paper's §open-problem content; (3) the plain regime and n ≥ 6 remain
untouched by decision, not by neglect.
