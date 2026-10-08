# The Strand Census Attacked with the Multi-Character Rigidity: the Character System Proved, the Marginal Route Measured Dead, and the Threshold Form of H1g

**Program:** LRC pair-sum-lattice reformulation. **Directive (reviewer,
2026-10-06, post-T-21):** "Attack the strand census. ... The next attack
needs a different tool: the multi-character rigidity that was calibrated
but not developed, or an entirely new angle." **Reading adopted:** the
multi-character rigidity is the projection-pinning / p-independent-counts
phenomenon, calibrated at four committed sites and never developed into a
lemma: (i) the T-7 zoo 42-saturation of one-coordinate projections; (ii)
the T-8c census counts 24/312 constant in p over seventeen primes; (iii)
the pinning-size lemma input (relative-offset sets <= 6, counts 28/68);
(iv) the T-8 zoo structure probes ([16,17,14,12] top vs [9,9,7,8] min at
p=17). If the reviewer meant a different object, the correction is one
line and this note's structure survives. **Artifacts (this session):**
`scripts/strand_rigidity_probe.py`,
`scripts/strand_rigidity_debug.py`,
`scripts/strand_rigidity_probe2.py`,
`scripts/strand_vprime_settle.py`,
`scripts/strand_vprime2_check.py`,
`scripts/strand_char10_check.py`, outputs
`out_strand_rigidity.json`, `out_strand_rigidity_debug.json`,
`out_strand_rigidity2.json`, `out_strand_gamma_p{17,19,29}.json`.

**Bottom line.** The attack was executed with the different tool, and the
tool **does not close the census**. Three permanent results came out of
the attempt. (1) **The character system is now a theorem** (Lemma C10):
the correct multi-character object on the moving-start space is the
extended system X = {4 start singletons} + {6 pairwise differences} (10
functionals at T=8); it is closed under the renormalization action up to
affine shifts, and its marginal-multiset is an orbit invariant —
machine-verified on the full top orbits at p=17 and p=19. The naive
readings (singleton profiles, pairwise profiles) are both **false**,
caught by machine check, with the mixing mechanism identified. (2) **The
uniform marginal rigidity law is measured dead at the zoo cell**: across
1,678 admissible keys (exhaustive at p=17/p=19, targeted at p=29, all
GT'd 1,654/1,654 against the committed T-20 census), the marginals are
essentially full — envelope = p at both k=2 primes, 27/29 at k=3;
fullness fractions 1.000 / 1.000 / 0.931 — and the sparsity that carries
the volume law is **joint** (median density in the marginal box 1.5% /
7.4% / 6.9%), which is exactly the arrangement-level structure the strand
census names. (3) **The threshold form of H1g** (Proposition T): the
trivial bound |Sol| <= p^{m-1} sits *exactly at* the cascade threshold —
p/ov -> T/(T-6), so c1_triv = (T/(T-6)+o(1))^{m-1} = 256 / 243 / 244.1
at T = 8/9/10, the value at which Theorem SC's depth denominator
vanishes. H1g's content is precisely that the distinct class beats the
trivial bound by a factor growing in k — measured 95x / 1,207x / 3,803x
at (k,eps) = (2,1)/(2,3)/(3,5). The paper is unchanged in scope (the
reviewer's "if it doesn't close" branch); the record gains two
route-closures and one sharpened target.

---

## 0. Scope and compliance

Cells: the committed zoo cell T=8 (1K,5U), census normal form, sizes
{2k,2k+1} (off-) at p=17/19 and {2k+1,2k+2} (off+) at p=29. Primes:
p=17 (full, 35 sets = 7 orbits), p=19 (full, 70 sets = 14 orbits), p=29
(targeted: every set containing the binary-cascade core {2,4,8} plus the
doubling-path quadruples — 19 of the committed 138 targeted sets;
disclosed). Enumeration is solution-set post-processing of committed
cells via `solutions_m5` (itself GT'd 9/9 against direct brute force at
T-20); no new cascade-depth measurements, no T=11, no sampling (hence no
intervals: every number below is an exhaustive census over the stated
scope). GT against the committed T-20 census states: **1,654/1,654
matched keys, 0 mismatches** (p=17: 1,010; p=19: 620; p=29: 24, with 24
skips = ordered-family keys outside the committed 710). No new fitted
parameters: every statistic below is a measured envelope, not a constant
entering any formula.

## 1. The tool developed: the character system (Lemma C10)

**The false first readings, and why they fail.** The natural first
developments of the pinning phenomenon — (V') the multiset of the 4
singleton marginal sizes is an orbit invariant; (V'') the same for the 6
pairwise-difference marginals — are both **false**, each caught by
machine check (V': profile multisets differ across orbit members at both
full primes; V'': 0/32 sigma-matched size tuples agree at p=19). The
mechanism, pinned by an explicit construction of the renormalization map
(p=19, S=(2,4,5,8), delta=5 -> S'=(3,4,6,8)): the bijection acts on
start coordinates as

  b_slot   = inv*(a_rho - a_delta) + const(sizes)   [moving slots]
  b_normal = inv*a_delta                             [the new normalizer]

so every moving slot is coupled to the anchor coordinate a_delta. The
image's *slot singletons* are the source's *anchor pairs*; the image's
*pairwise differences* (anchor cancels) are the source's non-anchor
singletons and pairs. Neither subsystem is closed; **the union is**.

> **Lemma C10 (the extended character system).** Let
> X = {e_i} U {e_i - e_j : i < j} on the moving-start space Z_p^{m-1}
> (|X| = (m-1) + C(m-1,2); 10 functionals at T=8). Under the
> renormalization bijection of Lemma V, each functional of X maps to a
> nonzero scalar multiple of a functional of X plus a size-dependent
> constant. Consequently the multiset of marginal sizes
> {|pi_chi(Sol(S,szt))| : chi in X}, aggregated over the
> permutation-closed set of size tuples, is a renormalization-orbit
> invariant.

*Verification.* The construction above is machine-verified end-to-end:
**801/801** images of Sol((2,4,5,8), szt) land in Sol((3,4,6,8),
szt∘sigma), injectively, at **every** size tuple; per-key counts match
under the explicit sigma with **0 mismatches**. The C10 invariance then
holds on the full top orbits: p=17 orbit of (2,3,4,8) (5 members, 9,633
solutions per set) and p=19 orbit of (2,4,5,8) (5 members, 801 per set) —
aggregate 10-functional multisets **identical across members at both
primes**. Independent Lemma V confirmations: the two admissible p=29
probe sets, (2,4,7,8) and (3,4,8,13), are members of one orbit (delta=7)
with equal totals 1,311.

**Consequences.** (a) The committed Lemma V proof *stands and is now
verified at a stronger level than committed* (per-assignment bijection
with explicit sigma; the prior record verified only the summed totals).
(b) Any future class-level (orbit-level) law for the census must be
stated on X — not on absolute starts, not on pairwise phases alone. (c)
The census object's symmetry analysis is complete: the orbit reduction
(Lemma V) plus the character closure (C10) together fix the coordinate
system in which "uniform over resonance classes" has meaning.

## 2. The threshold form of H1g (Proposition T)

> **Proposition T.** At every cell in the measured family, |Sol(f,sz)|
> <= p^{m-1} trivially, and p/ov = (Tk+eps)/((T-6)k+O(1)) -> T/(T-6).
> Hence the trivial bound realizes c1 = (T/(T-6)+o(1))^{m-1} — exactly
> the value at which Theorem SC's depth denominator
> (m-1)log(T/(T-6)) - log c1 vanishes.

Numerically the threshold is **256 / 243 / 244.1 at T = 8/9/10** —
near-constant across rungs: the cascade does not soften at higher T, and
the volume law's job does not get easier. Three consequences:

- **H1g, sharpened:** its content is *beat the trivial bound by a factor
  growing in k* — a joint-sparsity statement, not a marginal one.
  Measured sparsity below the trivial box: **95x / 1,207x / 3,803x** at
  p = 17/19/29 (top family, top size tuple: 17^4/876, 19^4/108,
  29^4/186). The committed margins over ov^4 (4.7x / 12x / 79x) are the
  same fact in volume-law units.
- **Unification with the sharp-H2m question:** the chain class's *proved*
  constants (Theorem B3: (m!)^2 * 4m * s_max, loose by orders of
  magnitude) face the same ~250 threshold; the chain's *measured*
  constant (2.0-3.6) clears it by ~70-100x. The sharp-constants open
  problem and H1g are one problem: every volume constant in the program
  must beat ~(T/(T-6))^{T-4}, and every measured one does.
- **The marginal route's precise burden:** a sufficient marginal law
  would be |pi_chi| <= gamma0 * ov on X with gamma0^{m-1} < threshold.
  This requires the fullness fraction |pi_chi|/p to fall below
  gamma0*(T-6)/T asymptotically (e.g. below 75% for gamma0 = 3). The
  measured fractions are 1.000 / 1.000 / 0.931 — no downward trend in
  evidence across the three committed primes.

## 3. The calibration (the rigidity law, measured on the correct system)

Per prime: the pinning envelope (max marginal size over all admissible
keys and all 10 characters), the fullness fraction, the gamma-envelope
(max |Pi|/ov over keys with ov >= 2) against the threshold 4, and the
top family's 10-functional profile at its top size tuple:

| prime | k, eps | envelope | fullness | gamma-env (thr 4) | top family @ top szt: singletons / pairs | n |
|---|---|---|---|---|---|---|
| 17 | 2, 1 | **17 (= p)** | 1.000 | 4.25 (ov=4 key) | (2,3,4,8)@(5,5,5,5,5): [12,17,17,16] / [17,10,17,15,14,12] | 876 |
| 19 | 2, 3 | **19 (= p)** | 1.000 | 5.00 (ov=2 key) | (2,4,5,8)@(5,5,5,5,5): [4,4,9,8] / [4,8,11,8,12,16] | 108 |
| 29 | 3, 5 | **27 / 29** | 0.931 | 2.70 (ov=10 key) | (2,4,7,8)@(8,8,8,8,8): [9,7,10,27] / [7,6,17,6,15,14] | 186 |

Chain-contrast control (the proved full-dimensional class), same probe:
singletons [17,17,17,17] / [19,19,19,19] / [29,29,29,29]; volumes
8,880 / 4,680 / 29,760 (the committed GT anchors reproduced exactly).

**Readings, graded.** (i) *Measured:* the marginals are essentially full
at every prime — the constant-envelope pinning law (the direct analog of
the T-7 42-saturation) does not exist at the zoo cell. (ii) *Measured:*
the top-family profiles are non-monotone in (k, eps) — [12,17,17,16] ->
[4,4,9,8] -> [9,7,10,27] — family- and window-dependent, not a law; the
p=29 profile (three pinned characters 9/7/10 and one free 27, the free
one the d=8 antipodal-sampler arc of the committed mechanism picture) is
recorded as an observation, not a pattern. (iii) *Measured:* the
gamma-envelopes at k=2 exceed the threshold through the additive
small-ov slack (keys with ov = 2..4 and near-full marginals); at k=3 the
envelope sits below threshold (2.70 < 4) — but with the fullness
fraction still at 0.93, the asymptotic law has no supporting data point
(the burden of §2 stands). (iv) *Measured:* the only invariant sparsity
is joint — median density in the marginal box 1.5% / 7.4% / 6.9%, max
1.0 only at singleton-solution keys. (v) *Flagged as a reading, not a
law (three committed sites, no fit):* cell-level marginal pinning is a
property of the other committed cell types (the T-7 zoo-type
42-saturation; the T-8c 24/312 p-independent counts) and of per-fiber
cascade objects (the 28/68 pinning-size input) — it does not port to the
T=8 (1K,5U) zoo cell. Whether the boundary is the cell type or the
number of moving arcs is not determined by the committed record.

## 4. The verdict

- **Does the multi-character rigidity close the strand census?** No. As
  a uniform marginal law it is false at the zoo cell (measured,
  exhaustively, on the correct character system). As a structural tool
  it delivered the symmetry analysis (C10) and nothing that bounds the
  arrangement count. The volume law's sparsity is joint — the pointwise
  no-holes condition along I, exactly the object the T-20 reduction
  named. The obligation after T-22 is unchanged in form: bound
  arrangements x phase freedom uniformly over renormalization orbits.
- **What the attack closed.** Two routes, with data: (a) marginal
  pinning (this note); (b) the necessary-condition families (T-20,
  committed). What remains live: the arrangement count itself. The
  sharpened target: constants must beat (T/(T-6))^{T-4} ~ 250
  (Proposition T) — the measured margins over that bar grow 95x ->
  3,803x with k, and the obstruction direction remains empirically
  absent.
- **The paper.** Unchanged in scope, per the reviewer's framing of the
  risk. Two items queued for the next revision (not executed now, per
  scope discipline): the H1g ledger row restated in the threshold form
  of Proposition T; a two-paragraph character-system subsection (C10)
  replacing the informal "projections pinned" phrasing in the
  mechanism-probes part of section 9.
- **No new lemma of the fractal kind was created.** C10 is a closure
  theorem for the symmetry analysis (it retires questions rather than
  opening them: singleton and pairwise profiles are provably the wrong
  objects; X is the right one). The obligation count is unchanged: the
  strand census stands, now with its coordinate system fixed, two
  route-closures recorded, and its success criterion numerical (beat
  ~250 with margin growing in k).

## 5. Verification record

| check | scope | outcome |
|---|---|---|
| solution counts vs committed T-20 census states | 1,654 keys across p=17/19/29 | 1,654/1,654 exact, 0 mismatches (24 p=29 skips: ordered-family keys outside the committed 710) |
| Lemma V per-assignment bijection (explicit construction) | p=19, S=(2,4,5,8) -> (3,4,6,8), delta=5, every size tuple | 801/801 images in target, injective, 0 per-key count mismatches under the explicit sigma |
| V' singleton-profile orbit invariance | full top orbits, p=17 and p=19 | FALSE (profile multisets differ across members) — mechanism: anchor mixing |
| V'' pairwise-profile orbit invariance | p=19, 32 sigma-matched size tuples | FALSE (0/32) — same mechanism |
| **Lemma C10 aggregate invariance** | full top orbits (5 members each), p=17 and p=19 | TRUE at both primes; totals 9,633 / 801 per set |
| Lemma V orbit law at k=3 | p=29 probe sets (2,4,7,8) and (3,4,8,13) | one orbit (delta=7), equal totals 1,311 |
| chain control | (1,1,1,1) at top size tuples, three primes | 8,880 / 4,680 / 29,760 — the committed GT anchors reproduced |

**Scope honesty.** The p=29 calibration covers the committed targeted
classes {2,4,8,x} + doubling path (19 sets, 2 admissible, both in one
orbit) — not the full 17,160-family space; the p=29 envelope (27) and
fullness (0.931) are statistics of that scope, disclosed as such. The
k-trend statements rest on three primes at two k-values with mixed eps
classes — the same caveat the volume law carried at T-20, restated here
rather than inherited silently.
