# The Kernel World Is a Fibration — Step 1 of the Kneser Bridge

**Program:** LRC pair-sum-lattice reformulation, bridge track (strategy memo of
2026-10-04). **Task:** the memo's first concrete step — "pick N = 49, compute
the stabilizer structure of B_u + B_v + B_w for all triples; if a Kneser-type
bound suffices to show non-coverage, generalize to N = p², then to products."
**Artifact:** `scripts/kneser_bridge_step1.py`, `scripts/out_kneser_step1.json`
(runtime 21.7 s). **Result:** the proposed kernel-world theorem is FALSE as
stated; the correct theorem is a fibration/transfer structure, proved and
verified; Kneser itself carries no slack anywhere covering-relevant; the true
residual is the *unit-assisted* cells, which are empty throughout the verified
range.

Conventions (byte-compatible with the committed battery N5..N9): n = 5 speeds,
T = n+1 = 6; pair (p,q) with N = v_p + v_q; bad set
B_w = {k ∈ Z_N : 6·min(wk mod N, N − wk mod N) < N} (strict); the pair
certifies τ(N) ≥ N/6 iff the FOUR bad sets (pair residue + three other speeds)
do not cover Z_N. "Kernel speed": gcd(w, N) > 1. The 3-set maps below (the
memo's literal "three bad sets") are reported alongside the 4-set τ-maps.

---

## 1. The proposed theorem is refuted at its first test case

> **Proposed (memo):** "For every composite N and every triple (u,v,w) of
> residues with non-trivial gcd structure, the three bad sets B_u, B_v, B_w
> do not cover Z_N."

**False at N = 49.** B_7 = {t : t ≡ 0,1,6 (mod 7)}, B_14 = {t : t ≡ 0,3,4},
B_21 = {t : t ≡ 0,2,5} (each = three full 7-fibers, 21 elements), and

    {0,1,6} ∪ {0,3,4} ∪ {0,2,5} = Z_7  ⟹  B_7 ∪ B_14 ∪ B_21 = Z_49.

Hand-verified twice (residue arithmetic and direct products, e.g. 14·17 =
238 ≡ 42 ∈ arc); asserted as a hard anchor in the script. The 4-set form
fails too: any pair with v_p + v_q = 49 and other speeds ≡ 7·{1,2,3}, e.g.
V = (20,29,7,14,21), has τ(49) < 49/6 — the pair does not certify (the set V
itself satisfies LRC via another pair, e.g. (7,21) on Z_28; losslessness is
unaffected). **Covering kernel configurations exist at 49 — but they are
exactly the transferred Z_7 ones.** That is the corrected statement, and it
is a theorem, not a census:

## 2. The fibration theorem (four lemmas, all proved, all verified)

**L-FIB (per-speed pullback — the heterogeneous lift law).** Let d = gcd(u,N)
> 1 and M = N/d. Then for every threshold T:

    B_u^{(N)} = π_M^{-1}( B_{u/d}^{(M)} ),   π_M : Z_N → Z_M reduction mod M.

*Proof.* u = dk with gcd(k,M) = 1. ut mod dM = d·(kt mod M); min(dc, dM−dc)
= d·min(c, M−c); so T·dist(ut, N) < N ⟺ T·dist(kt, M) < M ⟺ (kt mod M) ∈
arc(M) ⟺ (t mod M) ∈ B_k^{(M)}. ∎

Every kernel bad set is the pullback of a **unit** bad set at the reduced
modulus: the kernel world is a fibration tree over the divisor lattice of N.
Verified: 62 kernel speeds at {49,77,91,121,169} and **280 kernel speeds at
N = 1001 = 7·11·13** (three-prime tree), 0 failures.

**L-HOM (homogeneous transfer = the committed lift law, k-set version).**
If all k speeds share divisor d, covering at N ⟺ covering at N/d (preimages
distribute over union; π_M surjective). Verified: 56/126 (49), 276/841 (77),
420/1491 (91), 220/715 (121), 364/1365 (169) — 0 mismatches, k = 3 and 4.

**L-CRT (transverse covering at N = pq).** A kernel-only family of p-speeds
and q-speeds covers Z_{pq} ⟺ the p-speeds' transferred Z_q-images cover Z_q
OR the q-speeds' transferred Z_p-images cover Z_p. *Proof.* (⇐) immediate.
(⇒) contrapositive: if neither side covers, pick r ∉ ∪S_i (Z_q) and
s ∉ ∪T_j (Z_p); the CRT gives t ≡ r (mod q), t ≡ s (mod p); t is uncovered. ∎
Verified: 816/3876 (77), 1140/5985 (91) — 0 mismatches.

**L-FC (fiber counting at N = p², p ≥ 7).** Let c_p = max number of elements
of a unit bad set in any p-fiber (= max multiplicity of the arc's residue
distribution mod p; c_7 = 3, c_11 = 4, c_13 = 5). If a family has ≥ 2 kernel
speeds, then it covers ⟺ its kernel speeds alone cover (⟺ their transferred
sets cover Z_p): a non-covering kernel family leaves ≥ 2 full p-fibers gap,
and (k − #kernel) ≤ 2 unit sets contribute ≤ 2c_p < p elements per gap fiber.
*2c_p < p holds for every p ≥ 7* (c_p ≈ p/3 + 1). For 3-sets, 1 kernel + 2
units also never covers (2c_p < p). Verified with **0 prediction mismatches
over all 3.16 M combinations** at {49,77,91,121,169}, both arities.

## 3. The verified map (what actually covers at the composite moduli)

Distinct bad sets at N (up to the identity B_u = B_{−u}), all 3- and 4-set
combinations, classified by speed-type cell. "Covering" = union = Z_N.
**unit-assisted** = a covering whose kernel subfamily alone does NOT cover
(units needed). **unit-assisted = 0 everywhere.**

| N | structure | 3-set coverings | 4-set coverings | what they are |
|---|---|---|---|---|
| 49 = 7² | 24 sets (21 U, 3 K7) | 1 | 24 | kernel speeds 7k with k spanning the three antipodal classes of Z_7 — the K7 classification, transferred. Example: (7,14,21); 4-set e.g. (1,7,14,21). |
| 77 = 7·11 | 38 sets (30 U, 5 K7, 3 K11) | 1 | 38 | 11-speeds 11k spanning the Z_7 classes, e.g. (11,22,33). The 7-side never covers: Z_11 is not capable at n=5 (L-HOM). |
| 91 = 7·13 | 45 sets (36 U, 6 K7, 3 K13) | 3 | 138 | EITHER 13-speeds spanning Z_7 classes, e.g. (13,26,39), OR 7-speeds forming a transferred Z_13 covering: (7,21,28) and (14,35,42) — the two C₆ parity tilings of K13, pulled back. |
| 121 = 11² | 60 sets (55 U, 5 K11) | **0** | **0** | no covering configuration at all (Z_11 not capable ⇒ transferred world empty; unit world safe). |
| 169 = 13² | 84 sets (78 U, 6 K13) | 2 | 171 | 13-speeds forming a transferred Z_13 covering: (13,39,52), (26,65,78) — again the two parity tilings. |

**The composite kernel world is not a new world — it is the rigid-zone world
seen through the fibration Z_N → Z_{N/gcd}.** Covering at composite N ⟺ the
kernel subfamily covers ⟺ (L-FIB/L-HOM/L-CRT) a classified prime-level (or
product-side) unit configuration covers.

Re-verified against the committed record, exactly: rigid-zone normalized
counts {7:20, 13:68, 17:16, 19:64, 37:32}; primes 11,23,29,31,41,43,47 safe;
composite unit-safety at 49/77/91 (the committed micro-cases). **New safety
data beyond the committed scan range (≤ 95): N = 121 and N = 169 have no
covering configuration of any kind** (unit, kernel, or mixed), including the
cap-failing cells (K,U,U,U) — 146,300 combos at 121 and 492,960 at 169, all
empty.

## 4. The Kneser/stabilizer computation — and its honest verdict

The memo's literal request, computed exhaustively at N = 49 (all 300
distinct-set pairs; all 2,600 distinct-set triples) and all pairs at
77/91/121/169. H = translation stabilizer of the sumset; Kneser bound
|X+Y| ≥ |X+H| + |Y+H| − |H|.

| pair type (N=49) | outcome |
|---|---|
| (K7,K7), same class | H = 7Z_49 (order 7), \|X+Y\| = 35 = 21+21−7 — **Kneser-tight** (the transferred Z_7 arcs are APs: Cauchy–Davenport equality) |
| (K7,K7), distinct classes | X+Y = Z_49, H = Z_49 — trivially saturated |
| (K7,U) — all 63 pairs | X+Y = Z_49 (kernel set contains the subgroup 7Z_49; unit set meets every fiber) — saturated |
| (U,U) — ±-pairs | H trivial, \|X+Y\| = 33 = 2·17−1 — **Kneser-tight** (arc AP-structure) |
| (U,U) — distinct | X+Y = Z_49 — saturated (e.g. A + 2A = Z_49 by direct computation) |
| **all 2,600 triples** | B_u + B_v + B_w = Z_49, H = Z_49 — **everything saturates** |

At 77/91/121/169 the pattern repeats (e.g. (K11,K11) at 77: tight (11,55,55)
or saturated; (K13,K13) at 169: tight (13,117,117) or saturated). The only
pairs with nonzero Kneser slack live strictly inside already-classified
non-covering cells (transferred Z_11 arc sums at 77; a few unit pairs with
|A+λA| = 73 < 77).

**Verdict: Kneser's inequality is tight or trivially saturated on 100% of
covering-relevant configurations. There is no additive slack to convert into
non-coverage. The sumset route cannot be the bridge — the covering
obstruction is multiplicative, not additive.** (This is itself a structural
finding: the kernel-world sets sit at exact Kneser equality — they are as
"additively structured" as possible — which is precisely why their covering
behavior is controlled by the fibration instead.)

## 5. The multiplicative spectrum — where the bridge actually lives

|A ∩ λA| over all units λ (this equals |B_u ∩ B_v| for v/u = λ — all pairwise
intersections of unit bad sets; spectrum is even, λ ↔ −λ, verified 0
violations):

| N | \|A\| | min \|A ∩ λA\| | argmin (±-pairs) |
|---|---|---|---|
| 7 | 3 | 1 | {±2, ±3} — and covering EXISTS at Z_7 (three classes disjoint outside 0) |
| 49 | 17 | 3 | {±5, ±10}, 5·10 ≡ 1 — an inverse pair |
| 77 | 25 | 5 | {±5, ±31}, 5·31 ≡ 1 |
| 91 | 31 | 7 | {±18, ±23} ∪ {4,5} |
| 121 | 41 | 9 | {±5, ±24} |
| 169 | 57 | 11 | {±5, ±34} |

Covering by four unit sets needs total overlap exactly 4|A| − N (a perfect
near-tiling: 68−49 = 19 at N=49). The spectrum shows individual pairs can be
near-minimal, but min > 1 at every composite (vs min = 1 at Z_7, where
covering exists): **the multiplicative arithmetic at composites never
permits the simultaneous six-way near-disjointness a tiling requires.** That
is the quantitative content behind the committed micro-cases 49/77/91 — and
now behind their extension to 121/169 and to every mixed cell measured.

## 6. The corrected target for the 3-month program

> **Conjecture (no unit-assisted covering at composite moduli).** For
> composite N, no covering family exists in which the kernel subfamily does
> not already cover. Equivalently: units never complete a partial kernel
> covering, and the unit world at composite N never covers.
> [Verified: 0 occurrences in 3.16 M combinations at N ∈ {49,77,91,121,169},
> 3-set and 4-set; consistent with the committed 49/77/91 verification.]

Combined with the proved lemmas, this conjecture yields the uniform composite
theorem: covering at composite N ⟺ the fibration tree's prime-level worlds
cover — reducing the entire composite covering problem to the classified
rigid zone plus the unit worlds. Program, in priority order:

1. **Prove the no-unit-assisted conjecture.** The tools are multiplicative,
   not additive: lower bounds on |A ∩ λA| forced by composite structure, or
   direct transversal arguments (the witness formulation: ∃t with ut, vt, wt
   in the good arc G for all unit triples). This is the memo's direction 2
   (Lev/Konyagin/Ruzsa), now with a precise cell structure and first data
   (§5). The cap-failing cells (1 kernel + 3 units at p²; any kernel + units
   at pq) are the exact open targets.
2. **Extend the fibration tree to three-prime moduli.** L-FIB verified at
   1001; the recursive CRT lemma (mixed 7/11/13-speed families) is the next
   proof target, same method as L-CRT.
3. **Battery augmentation (not run, per the directive of no further
   V-extensions).** L-FC/L-CRT prove certification for every battery pair at
   p²/pq moduli whose four speeds include ≥ 2 kernel speeds not spanning the
   transferred classes — and prove failure for the spanning ones. This
   converts a named slice of the kernel-world residual from verified to
   proved; the 1-kernel and 0-kernel slices await item 1.

## 7. Soundness record

- Hand anchors asserted in-code: (7,14,21) covers Z_49; (7,14,35) does not;
  B_7/B_14/B_21 fiber descriptions exact.
- Committed-record anchors re-verified exactly: rigid-zone counts
  20/68/16/64/32 (pair-normalized 4-tuples, conventions reconciled); primes
  11,23,29,31,41,43,47 safe; 49/77/91 unit-safe.
- Lemma verification: L-FIB 0 failures (incl. 280 checks at 1001); L-HOM,
  L-CRT 0 mismatches; **0 prediction mismatches over all 3,160,535 mapped
  combinations**; unit-assisted coverings: 0.
- Spectrum evenness λ ↔ −λ: 0 violations (an argmin print truncation was
  traced and is cosmetic only).
- Full per-cell tables, examples, Kneser details, and spectra:
  `scripts/out_kneser_step1.json`.
