# The Multi-Character Rigidity Developed: the Arrangement–Phase Decomposition of the Joint Survivor Set, the Phase Freedom Closed by a Packing Lemma, and the Residual Compressed to the Arrangement Count

**Program:** LRC pair-sum-lattice reformulation. **Directive (reviewer,
2026-10-07, post-T-25):** "proceed with Multi-character rigidity
(calibrated at T-22, partially developed)." **Reading adopted:** the
partially developed object is the JOINT structure the T-22 calibration
named and no route has attacked — the marginals are full, the sparsity
is joint — so the development is the decomposition of the joint
survivor set along the character system X: the arrangement space (the
6 pairwise-difference characters, dimension m−2 = 3 at T=8) versus the
phase direction (the diagonal). **Artifacts (this session):**
`scripts/t26_multichar_rigidity_session.py`, output
`scripts/out_t26_session.json`, log `scripts/t26_run.log`.

**Bottom line.** The session ran under the pre-committed criterion (a
uniform-in-k margin for |Sol| against the trivial box, or the
decomposition with the margin localized to one named factor). The
criterion was not attained; the decomposition landed, and the margin is
now localized to exactly one factor. Three permanent results. (1)
**Lemma D (the phase-packing bound), proved**: the phase set of every
arrangement class is contained in a single cyclic interval of length
≤ s1 − span(H_u) + 1 ≤ ov + 1 (H_u the hole set of the moving union) —
machine-verified class-by-class at all seven committed top keys, zero
violations — so |Sol| ≤ (ov+1)·N_arr and the phase factor of the
margin is provably ≥ p/(ov+1) → T/(T−6), k-free. (2) **The measured
factorization**: the committed margins 95.3× / 1,206.7× / 3,802.6×
factor exactly as (arrangement compression) × (phase factor) =
15.75×6.05 / 90.25×13.37 / 283.6×13.41 (262.2×14.50), while the chain
benchmark sits at 1.95/3.89/3.34 × 4.82/7.16/7.11 — **the
distinct-class thinness IS the arrangement compression; the chain's
arrangement space is full.** In overlap units: N_arr/ov³ =
0.609 / 0.352 / 0.065 / 0.070 (distinct, decaying) versus
4.92 / 8.17 / 5.48 (chain, bounded); the chain volume constants
2.17 / 3.61 / 2.03 reproduce the committed 2.0–3.6 record exactly.
(3) **The residual compressed**: the strand-census obligation "bound
arrangements × phase freedom uniformly over renormalization orbits" is
discharged on the phase side by Lemma D and reduced on the other side
to a single named count — **the arrangement-count law**: N_arr =
o(ov^{m−2}) uniformly in k and over orbits. The paper is unchanged
this turn (the T-22 precedent for attack turns); the edit queue is at
the end of this note.

---

## 0. Scope and compliance

Charter fixed before the work: no transfer operators and no spectral
quantities (Lemma M's boundary respected — the session decomposes the
STATIC per-fiber solution sets); not a fourth marginal/projection
variation (no single-character law is sought; marginals computed only
as GT anchors against the T-22 record); pre-committed criterion and
three-way outcome tree; committed primes and families only (p = 17, 19,
29; the seven committed keys: the zoo chain (1,1,1,1) and top-distinct
(2,3,4,8) at p=17; chain and (2,4,5,8) at p=19; (2,4,7,8), (3,4,8,13),
chain at p=29); no T=11; no new cells; every number below is an exact
census over the stated scope (no sampling, no intervals). GT: all
committed counts reproduced — totals 9,633 / 801 / 1,311 / 1,311, top
keys 876 / 108 / 186 / 186, chains 8,880 / 4,680 / 29,760, the three
T-22 ten-functional marginal profiles exact, and covers()
cross-validation on samples of every top key (7/7 families, 0
mismatches).

## 1. Lemma D (the phase-packing bound)

**Setup.** Fix a committed key (family, size tuple). Solutions are
placements a = (a2..a5) of the four moving arcs such that the five arcs
(anchor [0, s1) + four moving APs) cover Z_p with pairwise overlaps
≤ ov = Σs_i − p. The arrangement map α(a) = (a2−a3, a2−a4, a2−a5) is
the coset of the diagonal (1,1,1,1) — the phase direction; classes
C_u = Sol ∩ (u + diagonal); |Sol| = Σ_u |C_u| exactly.

> **Lemma D.** For every arrangement u with hole set H_u = Z_p \
> ⋃(moving arcs) (phase-invariant), the phase set of C_u is contained
> in a single cyclic interval of length at most s1 − span(H_u) + 1,
> hence at most s1 − h_u + 1 ≤ ov + 1 (for h_u ≥ 1; if h_u = 0 the
> covering bound is vacuous and the class is bounded only by p).

*Proof.* Covering at phase t is exactly H_u + t ⊆ [0, s1) (the holes of
the moving union must land inside the anchor arc). The set of such t is
⋂_{x∈H_u} {t : (t+x) mod p < s1} — an intersection of cyclic intervals
of common length s1 ≤ p/2 (in-range: s1/p ≤ (2k+2)/Tk < 1/2 at T ≥ 5),
which is a single cyclic interval; its length is s1 minus the circular
span of H_u, and span(H_u) ≥ h_u − 1 ≥ s1 − ov − 1. ∎

**Machine record (class-by-class, all seven committed top keys):**
zero violations of any clause; zero cases of equality with the proved
cap (the measured max class sizes sit 3–6× below it: 5 / 2 / 3 / 2
distinct, versus caps 9 / 7 / 12 / 12 — the program's known
measured-to-proved constant gap, here on the phase side); h-ranges
[1,4] / [2,4] / [3,5] / [3,4] at the distinct keys — **the covering
binds at every distinct arrangement** (h_u = 0 never occurs), while the
chains carry 120 / 24 / 120 full-phase classes (h_u = 0, the moving
intervals alone tile the circle at every phase).

**The interval law (measured sharpening).** Not only is each phase set
contained in an interval — it IS an interval: every multi-member class
at every committed key has all phase gaps equal to 1 (564 / 32 / 100 /
93 / 6,360 / 22,464 / … consecutive-phase classes, 100%). The further
constraints (arc-5's AP containment; anchor overlaps) trim interval
ends and never punch interior holes. Proved: containment. Measured:
equality. This is the static-phase analog of the T-25 drift-orbit run
structure ("maximal lattice runs ≤ 2"), now with the interval shape
proved rather than observed.

## 2. The factorization, measured

| key (top size tuple) | \|Sol\| | margin p⁴/\|Sol\| | N_arr | arr-factor p³/N_arr | phase-factor p/avg | avg / max class | c1 = \|Sol\|/ov⁴ | N_arr/ov³ |
|---|---|---|---|---|---|---|---|---|
| 17 (2,3,4,8) | 876 | 95.3× | 312 | **15.75×** | 6.05 | 2.81 / 5 | 0.214 | **0.609** |
| 19 (2,4,5,8) | 108 | 1,206.7× | 76 | **90.25×** | 13.37 | 1.42 / 2 | 0.083 | **0.352** |
| 29 (2,4,7,8) | 186 | 3,802.6× | 86 | **283.59×** | 13.41 | 2.16 / 3 | 0.0127 | **0.065** |
| 29 (3,4,8,13) | 186 | 3,802.6× | 93 | **262.25×** | 14.50 | 2.00 / 2 | 0.0127 | **0.070** |
| 17 chain | 8,880 | 9.4× | 2,520 | 1.95× | 4.82 | 3.52 / 17 | 2.168 | 4.92 |
| 19 chain | 4,680 | 27.8× | 1,764 | 3.89× | 7.16 | 2.65 / 19 | 3.611 | 8.17 |
| 29 chain | 29,760 | 23.8× | 7,296 | 3.34× | 7.11 | 4.08 / 29 | 2.033 | 5.48 |

Readings, graded. (i) *Measured:* the margin factorizes exactly as
(arrangement compression) × (phase factor) — an algebraic identity
once the classes are defined; the content is the values. (ii)
*Measured:* **the arrangement compression carries the distinct-class
growth** — 15.75 → 90.25 → 283.6 (262.2) across (k,eps) =
(2,1)/(2,3)/(3,5) — while the chain sits flat at 1.95/3.89/3.34. The
distinct/chain arrangement ratio grows 8.1 → 23.2 → ~85. (iii)
*Measured:* the phase factor grows like p/avg with avg ∈ [1.4, 2.8]
(no trend, three primes), for both classes alike — the chain's own
margin growth (9.4 → 27.8 → 23.8) is its phase factor, its
arrangements being full. (iv) *Proved (Lemma D):* the phase factor is
≥ p/(ov+1) = 1.89 / 2.71 / 2.42 in-range, → T/(T−6) = 4 asymptotically
— k-free. **Therefore the provable growth of the margin is carried by
the arrangement compression alone**, and H1g (c1 → 0) is implied by
N_arr/ov³ → 0: c1 = |Sol|/ov⁴ ≤ (1 + 1/ov)·N_arr/ov³. Conversely
N_arr ≤ |Sol|, so the law and H1g differ by at most the class-size
factor — measured ≤ 4.1, proved ≤ ov+1. (v) *The k-trend of
N_arr/ov³ (0.609 → 0.352 → 0.065) rests on three primes at two
k-values with mixed eps classes — the same caveat the volume law
carried at T-20 and T-22, restated here rather than inherited
silently. The p=17 point (0.609) is not small; the decay claim is the
trend, not any single point.

## 3. The arrangement space: the pairwise ladder is empty

|Arr_adm| — the gauge-fixed arrangement triples satisfying all six
moving-moving overlap constraints (the t-invariant part of the census
pruning) — equals **p³ identically at every committed top key**
(4,913 / 6,859 / 24,389). At the operating point the all-on size
tuples carry ov ≥ every arc size, so no pairwise overlap constraint
can bind. Consequences: (a) the three-level ladder p³ → |Arr_adm| →
N_arr collapses — the entire arrangement compression (15.75× / 90.25× /
283.6×) is the covering-realizability bite, i.e., exactly the pointwise
no-holes condition the T-20 reduction named; (b) the arrangement count
is a pure joint-realizability count — no pairwise/marginal reduction of
it exists at the operating point, which is the T-22 "joint, not
marginal" calibration confirmed at the arrangement level.

## 4. The orbit test (the C10 anchor-mixing, made concrete)

p=19, S=(2,4,5,8) → S2=(3,4,6,8) (δ=5, the committed Lemma V
bijection): totals 801 = 801 and top keys 108 = 108 (the bijection
preserves everything global), but **N_arr = 76 → 32** and the average
class 1.42 → 3.375 — the arrangement/phase split is NOT
orbit-invariant. The renormalization maps the source diagonal to the
direction (0, 15, 0, 0) — a single-slot direction in the target (the
new normalizer's slot), not the target diagonal — and only 44/76
source classes map into single target classes. This is the C10
mechanism (slot singletons ↔ anchor pairs) operating on the split: the
orbit mixes arrangement into phase and back. The invariant content is
|Sol| (hence the margin, hence c1), the size-tuple multiset, and C10's
ten-marginal multiset; the FACTORIZATION is gauge-dependent, and the
uniform-over-orbits form of the residual must be stated as: the
arrangement-count law holds in every normalization along the orbit
(equivalently, for the worst normalization — measured spread across
one renormalization step: 76 vs 32).

## 5. What this closes, and what it does not

**Closed (proved, this session):**
- The phase freedom is sublinear: |Sol| ≤ (ov+1)·N_arr with the phase
  sets single intervals (Lemma D, + the class-by-class machine record).
  The phase half of the named obligation "arrangements × phase
  freedom" is discharged.
- The distinct/chain contrast is now structural: the thinness of the
  distinct class is exactly its arrangement compression (measured,
  exhaustively, on the committed keys); the chain — the proved
  benchmark — has a full arrangement space and its margin is its phase
  factor alone.
- The pairwise ladder is empty at the operating point: the arrangement
  count is pure joint realizability (the T-20 no-holes object, now the
  ONLY object on the arrangement side).

**Not closed (the residual, at its finest grain):**
- **The arrangement-count law**: N_arr = o(ov^{m−2}), uniformly in k
  and over renormalization orbits — equivalently, the fraction of
  arrangement classes admitting a covering phase vanishes. Measured
  0.609 → 0.352 → 0.065 (three primes, the standing caveat). This is
  a pure counting statement — census-shaped — and it is exactly the
  type the program's toolkit has never produced (the type-mismatch
  reading, third verification, now by construction: this session's one
  new proved lemma is a packing closure).
- The interval equality (end-trimming only) is measured, not proved.
- The sharp class-size constant (measured ≤ ~4.1, proved cap ov+1) is
  open — the phase-side analog of the program's measured-to-proved
  constant gaps.

**Pre-committed verdict: outcome (b).** The criterion (a uniform-in-k
margin) was not attained; the decomposition landed with the margin
localized to one named factor. The session does not close H1g; it
compresses the strand-census obligation from a product of two freedoms
to a single count, with the other freedom closed by a proved packing
lemma. Nothing here touches H2m (the correlation side); Theorem SC
remains conditional on H1g ∧ H2m with H1g now reduced, up to a
provably-sublinear factor, to the arrangement count.

## 6. Edit queue (paper; NOT executed this turn, per the T-22 precedent for attack turns — awaiting the reviewer's read)

1. `secscale.tex` — a new subsection after the transfer-matrix session
   paragraph: "The arrangement–phase decomposition (T-26)": Lemma D
   with proof; the class-by-class record; the factorization table (the
   seven committed keys); the Arr_adm = p³ datum; the orbit-mixing
   record; the residual restated as the arrangement-count law; the
   type-audit line.
2. `secsc.tex` — the H1g ledger row: the remaining-distance text
   extended (the phase freedom discharged by Lemma D; the obligation
   reduced to the arrangement count; the fourth route-closure framing
   NOT used — this route was developed, not closed).
3. `sec1.tex` — one sentence in the abstract and one in intro item 5:
   the decomposition result and the residual at its finest grain.

## 7. Verification record

| check | scope | outcome |
|---|---|---|
| committed totals / top counts / chain anchors | 7 families, 3 primes | 9,633 / 801 / 1,311 / 1,311; 876 / 108 / 186 / 186; 8,880 / 4,680 / 29,760 — all exact, 0 mismatches |
| T-22 ten-functional marginal profiles | 3 distinct top keys | exact match (singletons and pairwise differences) |
| covers() independent covering path | samples of every top key | 7/7 families, 0 mismatches |
| Lemma D, class-by-class | all classes of all 7 top keys | 0 violations; 0 equalities with the cap; h-ranges as recorded; zero-h classes only in the chains |
| interval law (contiguity) | every multi-member class, all 7 keys | 100% contiguous (all phase gaps = 1) |
| chain volume constants (GT) | 3 chain top keys | 2.17 / 3.61 / 2.03 — the committed 2.0–3.6 record reproduced |
| orbit bijection totals | p=19, S → S2 | 801 = 801; top 108 = 108; split 76 → 32 (the mixing, as predicted by C10's mechanism) |

**Scope honesty.** The decomposition is stated and verified at the
committed top keys (the margin convention of the committed record);
the class-size and N_arr statistics of non-top size tuples are not
part of this session's claims. The k-trend statements rest on three
primes at two k-values with mixed eps classes (the standing caveat).
The orbit test covers one committed renormalization step (the only one
with a committed explicit map). No new cells, no T=11, no cascade
runs, no spectral quantities.
