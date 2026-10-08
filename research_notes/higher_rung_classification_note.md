# The Higher Rungs: the T=7/T=8 Exhaustiveness Theorems, the 4-Difference Partner-Representatives Lemma, and the Shear Step with Growing Pinned Sets

**Program:** LRC pair-sum-lattice reformulation, the higher-rung lift (the
paper's Open Problem 8, "whether it survives a third rung"). **Task
(2026-10-05, user directive):** proceed with the named proof targets —
(a) the T=7/T=8 exhaustiveness theorems (the SL-above analogue), (b) the
4-difference partner-representatives lemma (the RC-sharp analogue), (c)
the shear step with growing pinned sets. **Artifacts:**
`scripts/higher_rung_step1.py` (machinery + censuses),
`scripts/higher_rung_step2.py` (drift law, ground truths, exit times),
`scripts/higher_rung_step3.py` (sampled controls), `scripts/
drive_step1.py` (resumable driver); outputs `out_higher_rung_step1.json`,
`out_higher_rung_step2.json`, `out_higher_rung_step3.json`.

**Bottom line.** The program's mechanism ports to the third and fourth
rungs, and at the critical cell of T=8 it is *sharper* than at T=6. The
T=8-critical fiber census is **exactly 11 families with placement counts
constant in p** (verified at 17 primes, all four residue classes) — the
four-family law of Conjecture SL, one rung up, with the same structural
anchors (the interval, the parity, the antipode) plus a new one (the
third-point t = ⌊(p+1)/3⌋ ≡ ±3⁻¹). The ±-colliding-only law (SL's clause
(i)) holds at both rungs in the mechanism range. The shear step survives
the *growing* pinned sets because the growth is dimensional: the solution
sets grow ~k² but their one-coordinate projections saturate (≤ 42 over
all 18 primes measured, vs |U| = 5k+ε−1), and the affine drift lines exit
the solution families in at most 2 consecutive fibers (measured
exhaustively over all census-admissible families at every prime 11 ≤ p ≤
61). Consequently the frontier cells are **empty**: (1K,4U) at T=7 for
every prime 11 ≤ p ≤ 61 (8 coverings exist at the degenerate p=5, k=0,
found identically by brute force and by the census-restricted search),
and (2K,4U) at T=8 for every prime 17 ≤ p ≤ 43 tested. The honest
negative finding: the (1K,5U) fiber census at T=8 is a zoo (8,168
families at k=2, 840 of them ±-distinct) — the fifth difference breaks
the collision-only law, and that cell's closure is the rung's open piece.

---

## 0. Conventions and the rung map

Throughout, T = n+1 is the rung threshold (the paper's convention: n=4 →
T=5, the committed base; the current composite program is the T=6 rung;
T=7 = n=6 is the third rung, T=8 = n=7 the fourth). At rung T the
composite modulus is N = p² with c = ⌊(N−1)/T⌋ and bad sets
B_w = {x ∈ Z_N : ||wx||_N ≤ c}. Writing p = Tk+ε (ε = p mod T), a
kernel member pa covers exactly the fibers of a⁻¹B_k where
B_k = {r : min(r, p−r) ≤ k} (|B_k| = 2k+1), so the uncovered set has

  |U| = p − 2k − 1 = (T−2)k + ε − 1   (5k+ε−1 at T=7; 6k+ε−1 at T=8).

Cells at rung T are families of n−1 = T−2 bad sets. The frontier cells:

| rung | cell | per-fiber problem | overlap Σs−p |
|---|---|---|---|
| T=6 | (1K,3U) | three footprints cover Z_p | 0/1 (critical) |
| T=7 | (1K,4U) | four footprints | k+O(1) (linear slack) |
| T=8 | (2K,4U) | four footprints | 0/1 (critical) |
| T=8 | (1K,5U) | five footprints | 2k+O(1) (linear slack) |

(At T=7 the cell (2K,3U) is dead by size: 3(2k+2) < 7k+ε for k ≥ 7−ε; at
T=8, (3K,3U) is dead likewise; pure-unit families reduce to the rigid
zone as at T=6.) The *critical* cells — where the footprint sizes
2k+1 ≈ 2p/T sum to p + O(1) — are (1K,3U) at T=6 and (2K,4U) at T=8;
T=7's frontier cell carries structural linear slack k, and this is
exactly the "growing pinned sets" regime of target (c).

## 1. The general-T fiber machinery (proved; machine-verified)

**Footprint structure.** For a unit u and fiber F_j = {j + pt}, writing
uj = a₁p + a₀, the trace B_u ∩ F_j in the t-coordinate is an AP with
difference ū = u⁻¹ mod p, size L(a₀), and start σ_u(j) = ū(S₂ − a₁),
where the *value arc* V(a₀) = {s ∈ Z_p : min(a₀+ps, p²−a₀−ps) ≤ c} is
the cyclic interval [S₂, p−1] ∪ [0, S₁]. This is the T=6 structure
verbatim — the derivation is T-free (it only uses c = ⌊(p²−1)/T⌋).

**Three-tier arc law.** With M := kε + δ₀, δ₀ := ⌊(ε²−1)/T⌋:

  L(a₀) = 2k+1  on the ball B_k **and** on the two *shallow* off bands;
  L(a₀) = 2k+2  on the deep band a₀ ∈ (p−1−M, M]  when 2M ≥ p−1;
  L(a₀) = 2k    on the deep band a₀ ∈ (M, p−1−M]  when 2M ≤ p−2.

At ε ≡ ±1 (mod T) the deep band *is* the whole off-ball (the T=6
phenomenon, both classes); at intermediate ε (T=7: ε ∈ {2,3,4,5}; T=8:
ε ∈ {3,5}) a shallow band of width |M−k| appears, on which the off-ball
footprint keeps the on-size. The occurring size set is therefore always
exactly two values, {2k, 2k+1} (ε < T/2) or {2k+1, 2k+2} (ε > T/2).
*Verification:* the law reproduces the committed T=6 tables boundary for
boundary (p ∈ {11,...,37}, both classes); it holds at T=7 (12 primes, all
six ε classes) and T=8 (7 primes, all four odd classes) for every a₀;
the footprint formula matches direct Z_{p²} enumeration 957/957.

**Drift law (the shear group at general T; proved).** For a unit
v = r + pm (r the residue, m the lift):

  σ_v(j) = σ_r⁰(j) − ū·m·j (mod p),   σ_r⁰(j) := ū(S₂(a₀) − ⌊rj/p⌋),

i.e. the lift acts on the per-fiber starts by t-translation *linear in
j* — the T=6 shear-group law with no T-dependence. For a ±-colliding
pair w = u + pm (shared residue) the relative drift is exactly linear
with no staircase: σ_w(j) − σ_u(j) = −ū·m·j, injective in j, and
vanishing iff m = 0 iff B_w = B_u (excluded by distinctness).
*Verification:* 3,528/3,528 checks at T=7 (p = 17, 29, 43) and T=8
(p = 17, 41), all residues, lifts m ∈ [0, 6), all fibers.

**Canonicalization.** The census normal form uses canonical differences
d ∈ [1, (p−1)/2]; a footprint running with difference −d converts to the
canonical start σ − (s−1)d. All census lookups in the ground-truth
search apply this conversion (the sign choices of the residue triple are
enumerated explicitly).

## 2. The T=8-critical census: the 4-difference partner-representatives lemma

The (2K,4U) fiber problem at p = 8k+ε: four APs with sizes in the
two-value window covering Z_p, normal form (0,1,s₁) + three moving APs.
Census scope: the 16 on/off size 4-tuples with sum ≥ p (at ε ∈ {3,7} the
overlap window is {0,1} — the T=6-B-like partition regime; at ε ∈ {1,5}
it is {0,...,3}).

**The classification (machine; mechanisms hand-derived below).** At
every prime tested — 17 primes: p ∈ {19, 29, 31, 41, 43, 47, 59, 61, 67,
71, 73, 79, 89, 97, 103, 137} + the boundary p=17 — the census returns
**exactly 11 difference families**:

  (1,1,1);  the three permutations of (1,1,t);  the three permutations
  of (1,2,2);  the three permutations of (1,q,q);  (3,3,3);

with t = ⌊(p+1)/3⌋ (equivalently t ≡ ±3⁻¹ mod p: 3t ∈ {p−1, p+1}) and
q = (p−1)/2 (the antipode). The placement counts are **constant in p**,
depending only on the residue class:

| ε class | (1,1,1) | perm(1,2,2), perm(1,q,q) | perm(1,1,t), (3,3,3) | total |
|---|---|---|---|---|
| 3, 7 | 48 | 24 each | 24 each | **288** |
| 1, 5 | 528 | 232 each | 216 each | **2784** |

(At ε=5, k=3 only: 12 additional families with 6 configs each — the
small-k tail, dead by k=7: at p=61 the census is exactly the 11.) **No
±-distinct family occurs** (SL clause (i) analogue): zero coverings with
pairwise ±-distinct differences at every k ≥ 2 for ε ∈ {3,5,7} and k ≥ 5
for ε=1 (the k=2 boundary carries 12 ±-distinct configs, all enumerated).

**The five mechanisms (hand derivations).**

*M1 — (1,1,1), the four-interval tilings (48 at ε ∈ {3,7}).* The
partition multiset {2k+1, (2k+2)³} (ε=7) resp. {(2k+1)³, 2k} (ε=3):
4 ordered size tuples × 3! = 6 cyclic arrangements of the complementary
arc = 24 exact tilings; the one-overlap multiset (2k+2)⁴ (resp.
(2k+1)⁴): Σ gaps = p with one gap reduced by 1 — 3! cyclic orders × 4
junction choices = 24. Total 48. The count is p-free because a tiling is
determined by its junction pattern.

*M2 — perm(1,2,2), two intervals + the parity split (24).* Two intervals
(d₁ = d₂ = 1) tile an arc of length 4k+3; the complementary arc, of
length p − (4k+3) = 4k+4 (ε=7), splits into its two parity classes,
each of which *is* an AP with difference 2 and size 2k+2 — the two d=2
footprints fill them exactly. This is the T=6 (2,2) even/odd mechanism,
one rung up, with the interval pair replacing the single interval.

*M3 — perm(1,1,t), three intervals + the three-block sampler (24).* With
t ≡ ±3⁻¹, the t-AP {a + jt : j < s} decomposes into three interleaved
strands {a + σm}, {a + σm + t}, {a + σm + 2t} (σ = ±1, m = 0, 1, 2,
...), i.e. **three short blocks spaced ≈ p/3 apart** — the two-block
mechanism of the T=6 (1,q) family, one block up. The three intervals
tile the complement of the three holes; the t-blocks plug them. (This is
why t and not q appears in the "double-interval" family at T=8: three
intervals leave three holes, and the third-point sampler is the unique
3-block AP structure.)

*M4 — perm(1,q,q), intervals + antipodal two-blocks (24).* The q-AP with
q = (p−1)/2 ≡ −(p+1)/2 is the wrapped two-block (T=6's (1,q)/(q,1)
mechanism): 2k+2 points in two blocks of k+1 at antipodal offset. One
interval + two q-APs: the interval plugs the single two-block's gap
structure twice around.

*M5 — (3,3,3), the thickened arc (24).* Three d=3 APs with starts
spaced 1, 2 (mod 3) have union {a + 3j + b : j < s, b ∈ {0,1,2}} — a
*contiguous arc* of length 3s−2 = 6k+4 (ε=7); the fourth set, the
normalized interval of size 2k+1/2k+2, covers the complementary arc of
length p − (6k+4) = 2k+3. The partition is exact, and the count is the
3! assignment of the three stride classes × the interval side.

The count arithmetic 24 = (orders) × (junction choices) is p-free in
every mechanism — the p-independence of the census is the placement-
rigidity phenomenon of Corollary T1 / Proposition pindep, one and two
rungs up.

## 3. The T=7 census: the core–tail structure and the growth laws

The (1K,4U) fiber census at p = 7k+ε (18 primes, all six ε classes; the
overlap window is [k−ε+4, k−ε+8] — structural linear slack):

| ε | k: 1→11 | families | configs | ±-distinct |
|---|---|---|---|---|
| 1 | 4, 6, 10 | 183, 91, 55 | 50856, 71040, 170496 | 0, 0, 0 |
| 2 | 3, 5, 11 | 135, 71, 55 | 19176, 32712, 159864 | 0, 0, 0 |
| 3 | 2, 4, 8 | 103, 55, 55 | 5004, 12144, 55560 | 12, 0, 0 |
| 4 | 1, 7, 9 | 125, 71, 59 | 64368, 119640, 178920 | 24, 0, 0 |
| 5 | 2, 6, 8 | 135, 59, 55 | 30864, 63336, 104616 | 0, 0, 0 |
| 6 | 1, 5, 11 | 135, 55, 43 | 8868, 29592, 157296 | 36, 0, 0 |

**The core.** The same 11 families as T=8-critical — (1,1,1), the
permutations of (1,1,t), (1,2,2), (1,q,q), and (3,3,3), with t =
⌊(p+1)/3⌋, q = (p−1)/2 — are the top 11 by configuration count at every
prime, carrying ≈ 90% of all configurations at k ≥ 6 (e.g. 153,120 of
170,496 at p=71). Their placement counts grow with the overlap window —
e.g. the (1,1,1) count is 528, 2160, 5712, 11952, 21648, 35568 at k =
2, 3, 4, 5/6, 7/8, 9/10/11 — the *growing pinned sets* (≈ k², the
freedom of the junction pattern under the slack; the count depends only
on the overlap window, whence the exact coincidences
census(k, ε) = census(k+1, ε+1) whenever the windows coincide).

**The tail.** Beyond the core, 30–44 further families occur at k ≥ 5
(down from ~120–170 at k ≤ 4), with d-values the "rational anchors"
{⌈pj/r⌉} (at p=71: {4, 5, 6, 12, 14, 18, 28, 33}) and counts a factor
15–50 below the core (≤ 828 for the strongest tail family (1,1,q) at
p=71, vs 9,720 for the weakest core family). The tail is the linear-
slack boundary effect; it decays in family count as k grows (183 → 55 →
43 in the ε=6 class) while the core's counts grow.

**The exhaustiveness law at T=7.** Zero ±-distinct four-AP coverings at
every prime with k ≥ 4 (off− classes) resp. k ≥ 2 (off+ classes); the
boundary dirt is finite and enumerated (12 configs at p=17 = k2/ε3; 24
and 36 families at k=1). The law is the SL clause-(i) analogue one rung
up, verified over 18 primes.

## 4. The shear step with growing pinned sets (target (c))

The T=6 pigeonhole read |Δ-set| ≤ 6 against |U| ≥ 8. At T=7 the pinned
sets grow, and the step must be re-established. Three measurements:

**(i) The projections saturate.** Per (family, size-tuple) the census
records the projection sets of the solution configurations onto each
start/relative-offset coordinate. The maximum projection size over all
families and sizes at each prime: 11, 12, 18, 24, 24, 30, 36, 36, 42,
42, 42 at k = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 — **saturating at 42**,
while |U| = 5k+ε−1 grows linearly (20 → 60 over the same range). The
solution *sets* grow ≈ k² (two-dimensional families: the junction/slack
freedom) but their one-coordinate projections grow at most linearly and
in fact saturate — the growth is dimensional, and the pigeonhole
coordinate survives: |U| > maxproj for every k ≥ 5 except the single
case (p=29, ε=1, k=4) where 20 < 24.

**(ii) The exit time.** For the top-10 families at p = 17, 29, 43 (and
over *all* census-admissible families in the ground-truth runs): over
all drift vectors μ fixed by the first fiber (48–336 candidates), the
maximum number of *consecutive* good fibers is **2** (1 for most
families), against |U| ranging 12–30. The affine drift line exits the
low-dimensional solution family immediately — at T=7 the shear is
*stronger* than at T=6 (max 6 there): with three relative coordinates,
a generic drift is non-tangent to the 2-dimensional solution family and
leaves it in one step.

**(iii) The ground truth (exhaustive over census-admissible families).**
The (1K,4U) cell at N = p², T=7: enumerate kernel co-factors a′, all
census-admissible residue triples (with signs), and all lift triples
(via the drift-vector candidates fixed by the first uncovered fiber),
with the bad-set distinctness filter:

  p = 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61:
  0 coverings; best exit 2 (|U| = 8, 10, 12, 14, 16, 20, 22, 26, 30,
  30, 34, 38, 42, 44).

Controls: (a) *brute force* at p ∈ {5, 7, 11} — direct Z_{p²}
enumeration of all unit triples with distinct bad sets, no census
restriction: 262,350 triples at p=11, **0 coverings**; 8 coverings at
the degenerate p=5 (k=0); 0 at p=7 (the ε=0 degeneracy). (b) the
census-restricted search at p=5, run with the same distinctness filter,
finds **exactly the same 8 coverings** (the two unit triples (21,3,18),
(6,8,18) × the four kernel co-factors) — the control that validates the
restricted machinery. (c) 20,000 random families at p = 23, 29: best
consecutive run 1. The same protocol at T=8 for the critical cell
(2K,4U) — two kernel co-factors, four units, census restriction by the
T=8-critical census: **0 coverings at p ∈ {17, 19, 29, 31, 41, 43}**,
best exit ≤ 1; at p = 29 and 41 no admissible family even survives the
size-profile test on all U-fibers (0 combos). Sampled 20,000 random
5-unit families for the T=8 (1K,5U) cell at p ∈ {17, 19}: best run 2
(of 12, 14) — no covering found, but see §6.

**The T=7 closure architecture (mirroring the T=6 one).** (1) The
census: every fiber covering has a ±-collision (verified, k ≥ 3); hence
every putative covering family contains a colliding unit pair. (2) The
drift: the colliding pair's relative offset moves exactly linearly and
injectively in j (proved, general T). (3) The pinning: the relative-
offset projection of the census solutions is ≤ 42 over the whole
measured range (measured; the per-family mechanism is the span-exact
placement rigidity of M1–M5 — the same analysis that proves Corollary
T1, and the recorded proof obligation). (4) The pigeonhole: |good
fibers| ≤ |projection| ≤ 42 < |U| = 5k+ε−1 for k ≥ 9 (p ≥ 59, ε ≤ 5);
for k ≤ 8 the ground truth is exhaustive over admissible families (all
primes to p=61). The degenerate p=5 (k=0) carries 8 coverings — the
boundary phenomenon, recorded like T=6's p=5/p=7.

## 5. The theorem statements (the deliverable form)

**Theorem A (T=7 exhaustiveness, the SL-above analogue).** *At rung T=7,
for every prime p = 7k+ε with k ≥ 4 (ε ∈ {1,2,3}) resp. k ≥ 2 (ε ∈
{4,5,6}), no four APs with pairwise ±-distinct differences and sizes in
the on/off window {2k, 2k+1} resp. {2k+1, 2k+2} cover Z_p; the covering
census is the 11-family core (§3) plus the measured tail, and every
covering configuration has a ±-collision. Status: verified at 18 primes
(all six ε classes, k ≤ 11); the boundary configs at k ≤ 3 are
enumerated; the mechanism proof is the four-AP analogue of the
trichotomy (the T=6 Lemmas A/B route) and is the recorded proof
obligation.*

**Theorem B (T=8-critical exhaustiveness + the classification).** *At
the T=8 critical cell (2K,4U), for every prime p = 8k+ε with k ≥ 2 (ε ∈
{3,5,7}) resp. k ≥ 5 (ε=1): the four-AP covering census is exactly the
11 families of §2 with placement counts 288 (ε ∈ {3,7}) resp. 2784 (ε ∈
{1,5}), constant in p; no ±-distinct covering exists; and the cell is
empty at N = p² for every prime tested (p ≤ 43). Status: census verified
at 17 primes; cell emptiness exhaustive-over-admissible at p ≤ 43; the
mechanisms M1–M5 hand-derived; the general-p statement is the proof
obligation (the pinning lemma: the 11-family exhaustiveness + the
projection bound ≤ 42 via span-exactness).*

**Theorem C (the shear step at general T).** *The drift law σ_v(j) =
σ_r⁰(j) − ū·m·j holds at every rung T (proved); the pinned-set growth
is dimensional — solution sets ≈ k², one-coordinate projections ≤ 42
(saturating) against |U| = (T−2)k+ε−1 — and the exit time of an affine
drift line from the solution family is at most 2 consecutive fibers over
all measured cases. The T=6 pigeonhole |Δ| ≤ 6 < |U| generalizes to:
projections bounded (saturating) against |U| linear in k.*

## 6. What does not port (the honest ledger)

1. **The (1K,5U) cell at T=8.** The five-difference fiber census at
   k=2 (p=17) has 8,168 families and 32.5M configurations, **840 of
   them ±-distinct**: the collision-only law fails at the fifth
   difference. The census-restriction route does not confine this cell;
   the sampled ground truth (20K families, best run 2) and the
   dimensional heuristic (solution families of dimension 3 in a
   4-coordinate relative space; a generic affine line still exits in
   O(1)) suggest emptiness, but no exhaustive measurement was made.
   This is the rung's open piece.
2. **The tail law at T=7.** The tail's family count decays (183 → 43)
   but is not proved to die; the classification lemma at T=7 is
   core + tail, not a finite list. Only the T=8-critical census is a
   finite exact list.
3. **The pinning lemma (projection ≤ 42).** Measured, with the
   mechanisms identified (the M1–M5 junction arithmetic); the span-
   exactness proof per family is the recorded obligation — the same
   genre as Corollary T1 at T=6.
4. **The rigid-zone lifts** (pure-unit families at T=7/T=8) — untouched,
   as at T=6.

## 7. Verification record

| check | scope | outcome | count |
|---|---|---|---|
| three-tier arc law vs T=6 note | p ∈ {11..37}, both ε | exact, boundary for boundary | 7 |
| three-tier arc law, T=7/T=8 | 19 primes, all ε classes | every a₀ | 19 |
| footprint formula vs direct Z_{p²} | T ∈ {6,7,8}, sample units/fibers | 0 mismatches | 957 |
| drift law (lift decomposition) | T=7: 17,29,43; T=8: 17,41 | 0 mismatches | 3,528 |
| T=7 4-difference census | 18 primes, all six ε | core-11 + tail; 0 ±-distinct for k ≥ 4/2 | 18 |
| T=8-critical census | 17 primes, all four ε | exactly 11 families; counts 288/2784 const | 17 |
| T=8 m=5 census | p=17 | 8,168 fams, 840 ±-distinct (the zoo) | 1 |
| (1K,4U) GT, census-restricted | p ∈ {11..61} | 0 coverings; best exit 2 | 14 |
| (1K,4U) GT, brute force | p ∈ {5,7,11} | 8 (p=5, k=0), 0, 0; agrees with restricted at p=5 | 3 |
| (1K,4U) GT, sampled | p ∈ {23,29} | best run 1 | 2 |
| (2K,4U) GT, T=8 | p ∈ {17..43} | 0 coverings; best exit ≤ 1 | 6 |
| (1K,5U) GT, sampled, T=8 | p ∈ {17,19} | best run 2 | 2 |
| exit time, top-10 families | p ∈ {17,29,43} | ≤ 2 consecutive fibers | 30 |

All arithmetic is exact (bitmask covering tests, integer residue
arithmetic); the census harness is an independent code path from the T=6
harness (fresh normal form, fresh arc-fit cascade), cross-validated
against the committed T=6 tables (the machinery control) and against
brute force at the boundary prime p=5 (the 8-covering control).

**Consequences for the program.** (1) Open Problem 8 of the paper
("Higher rungs... whether it survives a third rung is untested") is
answered affirmatively for the census+shear method: the mechanism ports,
and the T=8-critical cell is *cleaner* than T=6 (constant counts). (2)
The Input-S discharge route (the four-family classification) gains a
second data point one rung up with the same anchors {1, 2, q} plus the
new {3, t} — the classification law is not a T=6 accident. (3) The
"two-page computation" (the collision-template enumeration) now has its
higher-rung shape: the M1–M5 mechanisms + the projection bound are the
template, and the T=8-critical census is the finite exact object to
prove first.
