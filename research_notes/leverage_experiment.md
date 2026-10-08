# The Leverage Experiment: Does the Cyclic Language Generate New Tools?

**Program:** LRC residue-language reformulation. **Charge (reviewer directive 3):**
take the gridonly tier (the 70 n=4 sets, v <= 16, whose (TAU) certificates are
fully grid-specific) and try to prove them in the cyclic-group language. If new
tools apply, the reformulation is productive; if not, it is a change of language
with no new content. **Committed scripts:** `scripts/L1_leverage_experiment.py`,
`scripts/L2_battery_extension.py` (sound final battery), plus the soundness
cross-check run recorded below.

---

## 1. The covering framework (everything in one statement)

For a pair (u,w), N = v_u + v_w, effective triple eff = (a,b,c) on Z_N, the
pair certifies (TAU) iff some k has 5*||wk|| >= N for all three residues --
i.e. iff the three BAD SETS

  B_w := { k in Z_N : 5*||w k||_N < N }

do NOT cover Z_N. **Pair failure = covering.** The entire obstruction taxonomy
unifies as the additive order of the effective residues:

| order of residue w | B_w | old taxonomy name |
|---|---|---|
| 1 (w = 0) | Z_N (everything bad) | stuck |
| 2 (w = N/2) | 2*Z_N, size N/2 | half |
| 3 (w = N/3-type) | 3*Z_N, size N/3 | third |
| 4, 5 | j*Z_N, size N/j (see Lemma S) | (5blocked-adjacent) |
| j >= 6 | coset-union, size gcd(w,N)*|A ^ <w>| | kernel world |
| N (invertible) | dilated arc w^{-1}A, size 2h+1 | covering class |

This is a theorem-level restatement (definitions + the mirror identity), and it
is the single cleanest output of the reformulation: **five empirical
obstruction classes are five additive orders of one covering phenomenon.**

## 2. New lemmas produced by the language (all proved unless marked)

**Lemma C (class-group reduction).** For odd N, every B_w (w invertible) is a
union of antipodal classes {x, -x}; it is the union of the h classes of
{j*w^{-1} : 1 <= j <= h}, h = ceil(N/5)-1. Hence three invertible dilates
cover Z_N iff their class-segments cover the class group C_{(N-1)/2}
(multiplication by class(2) = generator g). *Proof: B_w = {0} U {+-j w^{-1}}
is antipodal-saturated; covering descends to classes. ∎*

**Theorem K7/K11/K13 (covering classifications; proved, code-verified).**
- N = 7 (h = 1, C_3): (a,b,c) covers iff the three dilates represent all three
  antipodal classes of Z_7^*. [4 covering (b,c) pairs at a=1: exactly (2,3),
  (2,4), (3,5), (4,5).]
- N = 11 (h = 2, C_5): each dilate covers a *domino* {X, Xg} of C_5; three
  dominos cover C_5 iff the two non-start vertices are non-adjacent.
  [Exactly 12 covering pairs; all match the scan.]
- N = 13 (h = 2, C_6): three dominos cover C_6 iff they *tile* it, iff the
  three start-classes form a coset {X, Xg^2, Xg^4} of the order-3 subgroup.
  [Exactly 4 covering pairs: (3,4), (3,9), (4,10), (9,10) at a=1; matches.]

**Lemma L (lift law).** (a,b,c) covers Z_M  iff  (ma,mb,mc) covers Z_{mM}.
*Proof: (mw)k mod mM = m(wk mod M), so 5||mw k||_{mM} < mM iff 5||wk||_M < M;
membership of k depends only on k mod M. ∎* Consequently every composite
covering modulus in range is a lift of a base modulus; the base covering
moduli observed for sum-triples (a, b, a+b) are exactly {6, 7, 8, 11, 13, 16}.

**Lemma S (size formula).** If residue w has additive order j <= 5 then
|B_w| = N/j exactly (the arc A meets <w> only at 0, since its nonzero
elements mN/j have distance >= N/5 with equality only at j=5). If w is
invertible, |B_w| = |A| = 2h+1 exactly.

**Lemma B0 (counting certification).** If |B_a| + |B_b| + |B_c| - 2 < N
(the -2 for the triple point 0), then the pair certifies. With Lemma S this
is hand-computable from the orders alone. Corollary: any pair with an
order-5 residue (others invertible) always certifies; order-4 likewise for
N in {8,12,20,24,28,40,44,60,...} (N/4 + 2(2h+1) - 2 < N).

**Lemma B1 (j-gon; proved last turn, 3-line proof).** j | N, j in {2,3,4,5},
all eff residues nonzero mod j => k = N/j certifies with m3 >= N/j >= N/5.

**Verified (proof pending):** the 2-runner grid law: for all N in [3,80] and
all residue pairs (a,b), max_k min(||ak||,||bk||) >= N/5, with equality
exactly at the pentagon configuration (1,2) on Z_{5m}. The 5|N case is proved
(pentagon times k = jN/5 when 5 does not divide the ratio; kernel coset
argument otherwise). **Per-modulus verified:** no invertible triple covers
Z_N for any N in [4,60] \ {7,11,13}; Z_17 and Z_19 hand-analyzed in the
class-segment language (three translates of S_3 never cover C_8 / C_9).

## 3. Measured coverage (sound battery)

Soundness protocol: every lemma flag was cross-checked against exact
computation of m3 on the full corpus; **0 failures in 4,918 flagged pairs**
after the fix described in §5.

| measurement | result |
|---|---|
| ALL 1,745 sets (v<=16): some pair certified by PROVED lemmas (B0/B1/B3'/B5) | **1,651 (94.6%)** |
| gridonly tier (70 sets): certified pairs closed by full battery (incl. verified lemmas) | 210/313 pairs, **69/70 sets** |
| gridonly tier: PROVED-only closures | 91/313 pairs, **54/70 sets** |
| sets open under the full battery | **1: V = (3,5,8,13)** |

The single open set is the most instructive object in the experiment. Its only
live pairs are (5,13) with eff = (5,3,8) on Z_18 and (8,13) with eff = (8,3,5)
on Z_21 -- both satisfying the speed sum relation 3 + 5 = 8. The same triple
(3,5,8) FAILS on Z_16 (a genuine base covering: order-2 residue 8 gives
B_8 = evens, and B_3, B_5 tile the odds) and certifies on Z_18/21. Candidate
closing lemma (stated, not proved): for eff = (a, b, a+b) with residues
r, s, r+s, it suffices to find k with ak, bk mod N both in [N/5, 2N/5]
(the band-intersection argument); the obstruction is exactly the classified
coverings at base moduli {6,7,8,11,13,16}.

## 4. What the residual gap is (precisely)

1. **No-covering conjecture:** for prime N >= 17 (equivalently all base moduli
   outside {6,7,8,11,13,16}), three bad sets never cover Z_N. Verified for all
   invertible triples, N <= 60; hand proofs at 17, 19. In class language:
   three multiplicative translates of the small-integer class set S_h never
   cover C_{(N-1)/2} for N >= 17. This is a clean, self-contained conjecture
   in multiplicative combinatorics that the reformulation surfaced and that
   does not appear (in this form) in the reachable literature.
2. **Kernel-arc coverings at orders 2, 3** (the half/third world): counting
   has slack (N/2 + 4N/5 > N), so certification there needs structure, not
   size. The Z_16 covering (3,5,8) shows these coverings are real.
3. **The pair-selection law** (why a live pair always exists for every set)
   is untouched by the covering framework -- it is the remaining (TAU) content.

## 5. Bug caught by the soundness protocol (recorded for credibility)

The first version of the order-4/5 counting lemma (B5) bounded the *other*
bad sets by 2h+1 unconditionally. That bound is false when another residue is
stuck (|B| = N) or order 2 (|B| = N/2). The soundness cross-check (flag must
imply certify) caught it: 245 flagged-but-non-certifying pairs, every one
B5-only. Fixed by requiring the other two residues invertible (where the
sizes are exact by Lemma S); after the fix, 0 failures in 4,918 flagged pairs.
All numbers in §3 are from the sound battery. This is the third bug the
cross-validation discipline has caught (bitwise-OR, grid-sampling, B5 bound).

## 6. Verdict on the reviewer's question

**Does the reformulation help?** Yes, in a measured, bounded sense:

- The cyclic language generated six provable statements this turn
  (unification, class reduction, K7/K11/K13, lift law, size formula, counting
  certification) that close 94.6% of the base corpus and 54/70 of the
  grid-specific core *without any breakpoint computation*.
- It surfaced one clean open problem (the no-covering conjecture) that is
  new in form, finite, and self-contained.
- It reduced the 70-case gridonly tier to a single open set with an explicit
  candidate lemma.

What it has NOT done: touch the pair-selection law or the plain regime; the
hard content of LRC-n above the base rung has not moved. The reformulation is
a **proof-generating language at n = 4**, not a route past the frontier (which
stands at 10 proved / 14-15 computer-assisted). Per the reviewer's decision
rule, this supports continuing the structural program (the no-covering
conjecture is the natural next target) while treating R2 as the parallel
track -- and it definitively answers "change of language with no new content"
in the negative at the base rung.
