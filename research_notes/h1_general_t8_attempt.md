# H1-General Attacked at T=8: the Volume Census, the Renormalization Lemma, and the Reduction to the Strand Census

**Program:** LRC pair-sum-lattice reformulation. **Directive (reviewer,
2026-10-06, post-T-19):** "Attempt H1-general at T=8 only. If H1-general
can be proved at the smallest rung where the pattern is visible, the proof
structure will suggest whether it generalizes. If it can't be proved at
T=8, the general case is very unlikely to go through, and the
reformulation is a bounded reduction." No T=11, no paper edit. **Artifacts
(this session):** `scripts/h1_t8_census.py` (+ state), `scripts/
h1_t8_structure.py`, `scripts/h1_t8_prime.py` (+ state), `scripts/
h1_t8_path.py`, outputs `out_h1_t8_census.json`,
`out_h1_t8_structure.json`, `out_h1_t8_prime.json`.

**Bottom line.** H1-general is **not proved at T=8 this session**, but the
attack changed its shape in three permanent ways. (1) A new structural
lemma is proved — **Lemma V (renormalization invariance)**: the covering
volume of a difference family, summed over size assignments, is invariant
under re-choosing which arc is the interval; the volume is a function of
the dilation class of the difference multiset, not of the family key.
Verified exactly at p=17 (35 sets collapse to 7 orbits) and p=19 (70 sets,
14 orbits in 9 volume classes). (2) The volume law is measured across the
k-range at T=8: the worst distinct family reaches **0.21 / 0.083 / 0.013**
of ov⁴ at p = 17 / 19 / 29 — the law |Sol| ≤ c₁·ov^{m−1} holds with
c₁ = 1 and a margin that **grows** with k (the distinct class sits 10× /
43× / 160× below the chain class at equal ov). No family approaches the
chain-class constant; the empirical obstruction the reviewer asked about
is absent, in the direction that matters. (3) The proof obligation is
reduced to a named bounded object — the **strand census**: after the
interval-spine reduction (cover the complement interval I = [s₁, p) of
length ≈ (T−2)k by the I-restrictions of the m−1 moving APs), every
covering is a cyclic arrangement of strand handoffs, and the volume is
(arrangement count) × (phase freedom per arrangement). What does not yet
exist is the uniform count over resonance classes — the same genre as the
committed Lemma-E machinery, and a real obligation, not a measurement.

---

## 0. Scope and compliance

The cell: T=8 (1K,5U), m=5 arcs, sizes {2k, 2k+1} (off−) or {2k+1, 2k+2}
(off+), census normal form (normalizer interval d₁=1 at the origin, ±
distinct differences in [2, D], D=(p−1)/2). Primes: p=17 (the committed
zoo point — full census), p=19 (k=2, ε=3 — full census; a committed
pool-law prime from T-19), p=29 (k=3, ε=5 — targeted: all binary-cascade
families, the doubling-path quadruples, anchor families, and a 120-set
random sample; a committed pool-law prime from T-19). The covering
convention is the cascade's (bystanders included — the m=4 GT convention
of T-18/T-19). Every decision procedure was validated against an
independent direct brute force (plain p⁴ enumeration, no cascade, no
pruning) before any number is cited below. No T=11, no paper edit, no new
cascade-depth measurement.

## 1. The census

**p=17 (k=2, ε=1, off−; sizes {4,5}; 32 size tuples; ov ∈ [3,8]).** All
840 ±-distinct families censused across all 32 size tuples.

- Admissible: **840/840** — the committed pool identity reproduced by an
  independent path (the T-19 saturation checker vs this session's
  counting cascade).
- Per-family totals (over 32 size tuples): median 2,981, p90 9,633,
  max 9,633. Per-(family,size) maxima: median 410, p90 876, max 876.
- The top size tuple is (5,5,5,5,5) (ov=8) for every family; the volume
  profile in ov is unimodal (top family: 48 / 507 / 1,954 / 3,434 / 2,814
  / 876 at ov = 3…8).

**p=19 (k=2, ε=3, off−; sizes {4,5}; ov ∈ [1,6]).** All 1,680 distinct
families censused.

- Admissible: **960/1,680 = 0.571** — exactly the committed T-19 pool-law
  saturation figure at p=19 (an independent reproduction).
- Totals: median 102, p90 645, max 801; per-size max 108 at (5,5,5,5,5)
  (ov=6).

**p=29 (k=3, ε=5, off+; sizes {7,8}; ov ∈ [6,11]).** Targeted census:
every 4-set containing {2,4,8} (all orders), the doubling-path
quadruples, the anchor sets, and a 120-set random sample — 138 sets /
710 families measured (the full space is C(13,4)·24 = 17,160 families;
not claimed).

- Admissible among measured: 20/138 sets (14.5%) — consistent with the
  committed full-scope saturation 0.133 at p=29.
- Top: total 1,311, per-size max **186 at (8,8,8,8,8) (ov=11)**.

**GT record (independent brute force, direct covering test over all p⁴
start tuples):** p=17: 876/876, 410/410, 0/0, 8,880/8,880 (the chain
control) — 4/4 keys; p=19: 108/108 (distinct top), 4,680/4,680 (chain) —
2/2; p=29: 186/186 ({2,4,7,8} top), 186/186 ({4,8,13,3}), 29,760/29,760
(chain) — 3/3. **9/9 MATCH, zero mismatches.**

## 2. Lemma V — renormalization invariance (proved)

> **Lemma V.** Let D = {1} ∪ S be the difference multiset of a census
> family (S a ±-distinct (m−1)-set), and for δ ∈ D let
> renorm_δ(S) = {canon(d·δ⁻¹ mod p) : d ∈ D, d ≠ δ}, where
> canon(x) = min(x, p−x). Then the total covering count summed over all
> size assignments is equal for S and renorm_δ(S):
> Σ_sz |Sol(S, sz)| = Σ_sz |Sol(renorm_δ(S), sz)|.

*Proof.* The covering condition is invariant under the global dilation
φ(x) = δ⁻¹·x of Z_p (a bijection; it carries an arc a + d·[0,s) onto the
arc δ⁻¹a + (d·δ⁻¹)·[0,s)). Apply φ to a covering placement of the family
S with size assignment (s_d)_{d∈D}: the image is a covering placement of
the family with difference multiset δ⁻¹·D, in which the arc that had
difference δ now has difference 1 — i.e. it is the normalizer of the
transformed system. Re-anchor the transformed system at that arc's start
(covering is translation-invariant) and canonicalize each difference
d ↦ canon(d) (an arc with difference −d′ is the reflection, x ↦ −x
composed with a translation, of an arc with difference d′; reflections
preserve coverings). The result is a census-normal-form covering placement
of the family renorm_δ(S), with the sizes carried along by the arcs.
Summing over all size tuples removes the bookkeeping dependence on which
arc received which size (the admissible size-tuple set is
permutation-closed), and the map is a bijection at each size assignment.
Hence the totals agree. ∎

**Verification.** p=17: the 35 difference-sets partition into exactly 7
renormalization orbits (each of size 5, one per choice of δ ∈ D), and the
7 volume classes are precisely these orbits (all 24 orderings of a set
share the total — the trivial relabeling invariance — and the 5 sets of
an orbit share it by Lemma V). p=19: 70 sets = 8 positive orbits + 6
zero-orbits; the 9 volume classes are exactly the unions of orbits
(8 singleton classes + the 30-set zero class = 6 orbits). Machine-checked
set-by-set in both primes. **The volume is a function of the dilation
class of D alone.** This cuts the family space 5-to-1 in general and
identifies the correct equivalence for any future classification: at T=8
the census object is not C(D−2, m−1) sets but their renormalization
orbits.

*A consequence for the cascade.* The two-fiber cut at dilation
λ = j/j₀ compares Sol₀ with a dilate of Sol_j — and Lemma V says the
volumes involved are dilation-class invariants. The cut law's input
volumes are therefore orbit functions; H1g needs constants uniform over
orbits, not over keys.

## 3. The volume law at T=8 (the H1g measurement)

Benchmark: the chain family (1,1,1,1) — the proved class (B1′/B2) —
measured at the same top size tuple, GT'd, with the ov⁴-normalized
constant stable at 2.0–3.6 across primes (the chain law
|Sol| ≈ c·ov^{m−1} with c ≈ 2–3.6 is the right scale).

| prime | k, ε | top distinct orbit | distinct max/size @ ov | ov⁴-normalized | chain @ same tuple | ratio distinct/chain |
|---|---|---|---|---|---|---|
| 17 | 2, 1 | orbit of {1,2,3,4,8} | 876 @ ov=8 | 0.214 | 8,880 | 0.099 |
| 19 | 2, 3 | orbit of {1,2,4,5,8} | 108 @ ov=6 | 0.083 | 4,680 | 0.023 |
| 29 | 3, 5 | orbit of {1,2,4,7,8} | 186 @ ov=11 | 0.013 | 29,760 | 0.0063 |

**Reading.** (i) The uniform bound |Sol| ≤ ov⁴ (c₁ = 1) holds at every
measured prime, with margin 4.7× / 12× / 79× — **the margin grows with
k** (and tightens with the ε-window at fixed k). (ii) The distinct class
is uniformly *thinner* than the chain class at equal ov, by a factor that
grows from 10× to 160×. (iii) The direction of risk for H1g is therefore
not a fat counterexample family but proof coverage — the opposite of an
obstruction. (iv) Honest caveats: three primes, two k-values, both ε
classes at k=2 but only one at k=3; the top orbit at each prime was found
by census at p=17/19 but only within the targeted scope at p=29 (the
binary-cascade + anchor + doubling-path + random sets); a family outside
that scope beating 186 at p=29 is possible in principle and would not
change the law's direction unless it beat ov⁴ outright, which the pool
collapse makes unlikely (only 20/138 measured sets are admissible at all).

## 4. The mechanism of the top class (the structure probes)

Probes on the top families at the top size tuple (p=17: (2,4,6,8), 876
solutions; plus the min-class (2,3,5,6) and a mid-class for contrast).

- **All arcs essential.** Every one of the 876 solutions needs all five
  arcs (essential-arc pattern (1,1,1,1,1) in 876/876). No spectator
  redundancy: the fourth difference genuinely participates. The
  p=17 orbit accident (all {2,4,8,x} share the total 9,633) is Lemma V —
  those sets are renormalizations of each other, not spectators.
- **The overlap budget is exact and distributed.** The identity
  Σ_{i<j} O_ij = ov + Σ_x C(o(x),2) holds with 0 violations in all
  solutions probed. Pairwise overlaps are small: mean 0.5–1.5, max 2–3,
  and every pair respects the spacing cap O_ij ≤ s_j/b_min(λ_ij, s_i)+1
  with room. No pair carries the covering; the arrangement is collective.
- **Projections are pinned.** Coordinate projections of the solution set:
  [16, 17, 14, 12] (top family), [9, 9, 7, 8] (min class) — the
  T=7-style pinning law persists at the zoo cell.
- **The block geometry.** Sample solutions show the mechanism: the
  interval covers [0, s₁); the d=2 arc covers every other point of a
  window of span 2s−1; the d=4 arc covers every fourth point of a span
  4(s−1); the d=8 arc is the antipodal two/three-block sampler; the
  circle is closed by interleaving these against the interval and each
  other. The top class is the **binary cascade** {1, 2, 4, 8} plus a
  fifth difference — at p=17 the orbit swallows every fifth value
  {3,5,6,7}; at p=19 the fifth is 5; at p=29 it is 7 (observed pattern:
  fifth = ε+2 at the three primes measured — three points, flagged as
  such, not a law).
- The doubling path at p=29 (2 is a primitive root: the canonicalized
  doubling orbit 2→4→8→13→3→6→12→5→10→9→11→7→14 visits every difference)
  shows that only cascade-containing quadruples survive among path
  quadruples — resonance structure is necessary for admissibility at
  k=3, consistent with the pool-law collapse.

## 5. The reduction to the strand census (where H1g now lives)

**Step 1 — the interval-spine reduction (elementary).** The complement of
the normalizer interval is the interval I = [s₁, p) of length
ℓ = p − s₁ ≈ (T−2)k. Every covering of Z_p by [interval + (m−1) APs]
restricts to a covering of I by the I-restrictions of the moving APs
(plus the two boundary points at the wrap). The moving arcs' total size
is (m−1)·2k ≈ (T−3)·2k against ℓ ≈ (T−2)k — a near-exact cover with
slack ov. So H1g at T=8 is: *count the ways to cover an interval of
length ≈ 6k by the I-restrictions of 4 APs (steps d₂..d₅, sizes ≈ 2k),
times the boundary/bystander corrections.*

**Step 2 — the strand structure.** The I-restriction of an AP with step
d is a union of *strands*: maximal runs of consecutive t whose points
a+t·d stay inside I; within a strand the points sit d apart. Covering
consecutive points of I requires strands of *different* steps to
interleave — a pair (d, d′) interleaves only at compatible phases, and
the phase-compatibility condition per pair is a three-distance-type
statement about λ = d′/d (exactly the Lemma-E genre, on the pair lines
t·d − u·d′ ≡ δ that already govern the pairwise overlaps O_ij).

**Step 3 — the arrangement form of the volume.** Every covering induces a
cyclic *arrangement*: the maximal runs of I and, for each run, the subset
of strands carrying it; the handoffs between runs are pair-phase
conditions. The volume decomposes as

  |Sol| = Σ_arrangements (phase freedom of the arrangement),

and the data says the answer is ≤ ov⁴·(small, orbit-dependent constant)
with the constant *decaying* in k. What is missing for a proof is the
uniform count: bound the number of arrangements (the handoff patterns)
and the phase freedom per arrangement, uniformly over the resonance
classes of the difference set — i.e., over renormalization orbits
(Lemma V reduces exactly to this). The chain class is the degenerate
resonance class where all steps coincide and the strand picture collapses
to the B1′/B2 telescoping; the binary cascade is the pure ratio-2 class;
the general orbit mixes ratio classes.

**Why this is the honest reduction.** The three necessary-condition
families already available (the overlap budget Σ_{i<j}O_ij ≥ ov; the
pair spacing caps; the last-arc fit count ≤ s₅+1) do NOT discriminate —
random placements typically satisfy the budget (the expected total
pairwise overlap ≈ 12 exceeds ov at p=17). The discriminating structure
is the pointwise no-holes condition along I, which is what the
arrangement/strand language encodes. A proof of H1g at T=8 must count
arrangements; that is a bounded, concrete combinatorial obligation — but
it is class-by-class work (per resonance type), which is precisely the
fractal pattern the reviewer named.

## 6. The verdict

- **Is H1g true at T=8?** Every measurement says yes with a growing
  margin: c₁ = 1 works at all three primes; the worst distinct family
  decays relative to both ov⁴ and the chain benchmark as k grows.
- **Is H1g proved at T=8?** No. The proved pieces are: Lemma V
  (dilation-class reduction — new this session), the exact overlap-budget
  identity and pair caps (verified on the full solution sets), and the
  chain class (prior). The open core is the uniform strand census.
- **Does the T=8 proof structure suggest generalization?** The reduction
  chain (renormalization → interval-spine → strand census) is T-free in
  form and ports verbatim to any rung; the strand census itself is where
  class-by-class work lives, and its difficulty grows with the number of
  resonance classes per orbit, which grows with D = (p−1)/2. The T=8
  evidence (margin growing in k) suggests the general law is comfortably
  true; the T=8 proof attempt suggests the *technique* generalizes only
  as fast as the pair-line/three-distance machinery covers resonance
  classes uniformly — the same bottleneck as H2m's sharp constants.
- **What would change the verdict.** A family at some (T=8, k) whose
  per-size volume exceeds ov⁴ — the measured pool collapse makes this
  increasingly unlikely (admissibility itself dies at 13% by p=29) — or a
  strand-census proof at T=8 with uniform constants, which would make
  H1g a theorem at the rung and supply the template for general T.

## 7. Verification record

| check | scope | outcome |
|---|---|---|
| counting cascade vs direct brute force | 9 keys across p=17/19/29 (incl. the chain controls 8,880 / 4,680 / 29,760 and zero keys) | 9/9 exact MATCH |
| admissibility reproduction | p=17: 840/840; p=19: 960/1,680 = 0.571 | matches the committed T-19 pool census exactly (independent code path) |
| Lemma V orbit law | p=17 (35 sets), p=19 (70 sets) | 7 resp. 14 orbits; every volume class is a union of size-5 orbits; all 24 orders share totals |
| budget identity Σ_{i<j}O_ij = ov + ΣC(o(x),2) | 876 + 30 + ~400 solutions probed | 0 violations |
| pair spacing caps O_ij ≤ s_j/b_min(λ_ij,s_i)+1 | all 10 pairs × probed families | holds with room (max 2–3 vs caps 3–6) |
| essential-arc patterns | top/mid/min families | (1,1,1,1,1) in all solutions — no spectators |
| order-independence of volume | all 840 (p=17), 1,680 (p=19) families | totals constant across the 24 orderings |
| solutions_m5 list rebuild vs counting cascade | 4 families × 32 size tuples | 0 mismatches |

**Session hygiene note.** The p=29 targeted census is partial by design
(138 of 715 sets); the full p=29 census was out of scope-time and is not
claimed anywhere. The top-orbit identification at p=29 rests on the
targeted scope (all cascade-containing sets, the doubling path, anchors,
random 120); a missed top outside that scope would have to exceed ov⁴ to
affect the law's direction, and the admissibility collapse argues against
it. This is recorded as a caveat, not a claim.
