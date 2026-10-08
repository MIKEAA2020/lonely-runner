# No Unit Assistance — The Corrected Theorem, Proved in Its Provable Part

**Program:** LRC pair-sum-lattice reformulation, kernel-world track. **Task
(T-11, directive of 2026-10-04):** Priority 1 — prove the no-unit-assistance
statement ("units never complete a partial kernel covering") that replaces the
refuted "no covering triple at composites" conjecture; fill in the counting
constants; handle the small primes; extend to N = pq. Priority 2 — state the
kernel-world structure theorem. No censuses: proofs, plus targeted
verification of the exact statements the proofs rely on (the established
discipline). **Artifacts:** `scripts/nua_proof_verify.py`,
`scripts/out_nua_verify.json` (runtime ≈ 55 s, all assertions pass).
**Result:** the sketch of 2026-10-04 is right in spirit, wrong in two
constants and one argument; with the corrections, the theorem is provable for
**at most two unit sets, for three pure unit sets, and for pure-unit
families reduced to the rigid zone** — at every prime p ≥ 5 (with p = 3
genuinely excepted: the blanket statement is FALSE at N = 9). The residual is
one named cell family plus rigid-zone lifts, now with sharp per-fiber budgets.
This closes Priority 1 in its provable part and delivers Priority 2 as a
stated structure theorem with an explicit proved/conditional split.

Conventions (byte-compatible with T-10 / the committed battery): n = 5, T = 6,
N composite, c = ⌊(N−1)/6⌋, bad set
B_w = {k ∈ Z_N : 6·min(wk mod N, N − wk mod N) < N}, arc A = B_1 = [−c, c].
**Kernel speed:** gcd(w, N) > 1; **unit speed:** gcd(w, N) = 1. Covering
families are multisets of ≤ 4 distinct bad sets (the pair convention: the
pair's two speeds contribute one bad set with multiplicity two). For a prime
M | N, the **M-fibers** are the cosets of ker(Z_N → Z_M), each of size N/M;
ρ_p := ⌊p/6⌋ is the transferred ball radius at modulus p (this is also the
Z_p bad-set radius at T = 6, since ⌈p/6⌉ − 1 = ⌊p/6⌋ for prime p ≥ 5 — the
fibration is threshold-compatible, which is L-FIB).

---

## 0. Executive summary

| statement | status |
|---|---|
| Transferred kernel ball at p² has "3 residues mod p" (the sketch) | **wrong for p ≥ 11**: it has 2ρ_p + 1 residues (3 only at p = 7, 11) |
| Case 0+3 "reduces via the p-adic fibration to three scaled balls on Z_p" (the sketch) | **invalid as stated**: every unit bad set surjects onto Z_p under the mod-p projection (§1.2); the correct reduction for pure-unit families runs through the *p-multiples* (Theorem 4) |
| Blanket "no unit assistance at composite N" | **false at N = 9** (§2); true and proved for p ≥ 5 in the cells with ≤ 2 unit sets (Theorems 2, 3) |
| Three unit bad sets never cover at N ≥ 31, N ≢ 0 mod 6 | **proved** (Theorem 1; new tool: a Dirichlet-approximation lower bound on arc intersections, Lemma 2) |
| Pure-unit 4-covering at p² exists only if p is in the rigid zone | **proved** (Theorem 4); hence no pure-unit covering at N = 77 — a formerly verified-only fact, now proved |
| (1 kernel + 3 units) and rigid-zone lifts | **open** — the unique residual cell family; occupied at p = 3 (N = 9), verified empty at p = 5 and p ∈ {7, 11, 13} |
| Kernel-world structure theorem (Priority 2) | **stated** with the proved/conditional split (§7) |

The honest one-line verdict on the directive's success criterion: *the
counting constants close for ≤ 2 units and for 3 pure units; they do not
close for ≥ 3 units together with ≤ 1 kernel set.* That residue is no longer
a fog — it is one cell family with a sharp per-fiber near-tiling budget
(§6), an occupied boundary instance at p = 3, and an identified attack
route.

## 1. Corrections to the sketch

### 1.1 The transferred radius is ρ_p = ⌊p/6⌋, not 1

For a kernel speed w = pa (p ∤ a) at N = p²: ||pak||_{p²} = p·||ak||_p, so

    B_{pa}^{(p²)} = {k : ||ak||_p ≤ ⌊c/p⌋} = π^{-1}(B̃_a^{(p)}),
    B̃_a^{(p)} = a^{-1}·[−ρ_p, ρ_p],   ρ_p = ⌊p/6⌋,

using ⌊(p²−1)/(6p)⌋ = ⌊p/6 − 1/(6p)⌋ = ⌊p/6⌋. The transferred ball has
2ρ_p + 1 residues: **3 at p = 7 and 11, 5 at p = 13, 7 at p = 19, 11 at
p = 31.** The sketch's "in general |B_a mod p| = 3" holds only at p ≤ 11;
all counting below uses 2ρ_p + 1. (This is why 169 = 13² transfers the two
C₆ parity tilings as 5-element balls — as T-10 measured — and why a
3-kernel covering at p² needs 3(2ρ_p+1) ≥ p, which fails for every
p ≡ 5 mod 6 and holds for p ≡ 1 mod 6, restricting 3-kernel coverings to
p ≡ 1 (mod 6) primes: within the committed classification, p ∈ {7, 13}.)

### 1.2 The mod-p projection of a unit bad set is everything

The sketch's Case 0+3 argues that three unit bad sets at p² "reduce via the
p-adic fibration to a uniform covering of Z_p by three scaled balls." The
reduction step is not available: **π_p(B_u) = Z_p for every unit u and every
p ≥ 5.** Indeed, for a residue class j mod p, the values uk (k ≡ j mod p)
fill the residue class of uj mod p; that class contains an element of
norm-distance ≤ (p−1)/2 from 0, and (p−1)/2 ≤ c holds for every p ≥ 3
((p−1)(p−2) ≥ 0). So each fiber of π_p is met by every unit bad set, and
projection loses all information. The usable reduction is in the other
direction — through the *p-multiples* (the zero fiber of π_p), where the
transfer is exact (Theorem 4).

### 1.3 The blanket statement has a counterexample at p = 3

At N = 9 (c = 1): B_1 = {0,1,8}, B_2 = {0,4,5}, B_3 = {0,3,6} (the unique
kernel-speed bad set), B_4 = {0,2,7}, and

    B_1 ∪ B_2 ∪ B_3 ∪ B_4 = Z_9,   while  B_3 = {0,3,6} alone does not cover.

Units are essential: this is a genuine unit-assisted covering, and it is
pair-realizable (V = (4,5,1,2,3): the pair (4,5) has N = 9 and its four bad
sets are exactly B_4, B_1, B_2, B_3 — the pair fails to certify). Every
theorem below is therefore stated with explicit prime hypotheses, and the
small primes are mapped exactly (§2): p = 2 vacuous (no covering of any kind
at N = 4), p = 3 false (the covering above), p = 5 clean (all cells empty at
N = 25, verified — the open cell included).

## 2. The small primes, settled as boundary data

**N = 4 (p = 2).** c = 0; the only bad sets are {0} (units) and {0,2}
(kernel). No family of ≤ 4 distinct bad sets covers Z_4. The blanket
statement holds vacuously. (Fiber cap cap(2) = 1, but 2·cap(2) = 2 is not <
2, so the ≤ 2-unit theorem does not apply — it does not need to.)

**N = 9 (p = 3).** The four distinct bad sets are listed in §1.3; the
covering map over 3-multisets (20 combos) is empty in every cell; over
4-multisets (35 combos) exactly one covering exists — the family
(1,2,3,4) of cell (K3, U, U, U), unit-assisted. Note the mechanism: the
kernel set is F_0 (ρ_3 = 0), and the three unit balls meet each of the
fibers F_1, F_2 in exactly cap(3) = 1 element each — 3 units × 1 element =
3 = fiber size — an *exact per-fiber tiling by three units at minimum
profile*. The counterexample lives precisely in the cell that §6 leaves
open; any proof of the residual conjecture must therefore use p ≥ 5 (i.e.
ρ_p ≥ 0 with cap ≥ 2, or ρ_p ≥ 1).

**N = 25 (p = 5).** c = 4, |A| = 9, distinct bad sets = 10 unit balls + the
single kernel set F_0 = 5Z (all kernel speeds 5a transfer to the Z_5-ball of
radius ρ_5 = 0, which is {0}). The full maps: 3-multisets (286 combos) —
empty in every cell; 4-multisets (1001 combos) — empty in every cell,
including (K5, U, U, U) (220 combos) and (U, U, U, U) (715 combos). So p = 5
clean like p ≥ 7: no unit-assisted covering exists — but only the parts
below that apply at p = 5 are *proved* ((0,4) is proved via Theorem 4 since
every Z_5 unit bad set is {0}; (1,3) is verified empty, per-fiber budget
3·cap(5) − 5 = 1). The spectrum floor also holds at N = 25 by direct
computation: min_{λ≢±1} |A ∩ λA| = 3 (attained at λ ∈ {3, 4, 6, 8, 9, 11,
14, 16, 17, 19, 21, 22}; the remaining λ give 5 — an 18-line hand check,
recorded in the harness).

## 3. Lemma 1 — the uniform fiber cap (exact, closed form)

> **Lemma 1.** Let N = sM with s, M primes (s = M allowed; s ≥ 2), and let
> F_j = {k ≡ j (mod M)} be an M-fiber (of size s). For every unit u:
>
>     B_u ∩ F_j = u^{-1}·(A ∩ F_{uj mod M})   and
>     |B_u ∩ F_j| ≤ cap(s) := 2⌊s/6⌋ + 1 + [s ≡ 5 (mod 6)].
>
> The bound is sharp (equality is attained for some unit and fiber), and
> cap(s) < s for every prime s ≥ 3, with 2·cap(s) < s for every prime s ≥ 5
> and 2·cap(3) = 2 < 3.

*Proof.* The identity: k ∈ B_u ∩ F_j ⟺ uk ∈ A and k ≡ j (mod M) ⟺ (u is a
unit mod M since gcd(u, N) = 1) uk ≡ uj (mod M) and uk ∈ A ⟺ k = u^{-1}r
with r ∈ A ∩ F_{uj mod M}. For the count, A ∩ F_a (a ∈ [0, M)) = {a + Mt :
t ∈ [0, s)} ∩ ([0, c] ∪ [N−c, N)), so

    |A ∩ F_a| = ⌊(c−a)/M⌋ + 1 + ⌊(c+a)/M⌋ = 2γ + 1 − [a > δ] + [a ≥ M−δ],

with γ = ⌊c/M⌋ and δ = c mod M. The identity holds for every a ∈ [0, M):
the first term is 0 exactly when a > c, which is the case ⌊(δ−a)/M⌋ = −1
(since δ − a ∈ (−M, M)); the last term uses ⌊(δ+a)/M⌋ = [δ + a ≥ M].
Now γ = ⌊(sM−1)/(6M)⌋ = ⌊s/6 − 1/(6M)⌋ = ⌊s/6⌋ = ρ_s (subtracting 1/(6M)
never crosses an integer for prime s ≥ 5; s ∈ {2, 3} check directly). The
maximum over a is 2γ + 2 iff some a satisfies a ≤ δ and a ≥ M − δ, i.e. iff
2δ ≥ M; otherwise 2γ + 1. It remains to evaluate 2δ ≥ M. Write r = N mod 6
∈ {1, 5} (products of primes ≥ 5), so 6c = N − r exactly, hence
6c ≡ −r (mod M) and c ≡ −r·μ_M (mod M) with μ_M = 6^{-1} mod M. Case
analysis on (s mod 6, M mod 6) ∈ {1, 5}²:

| (s, M) mod 6 | r | μ_M | c mod M | 2δ ≥ M? |
|---|---|---|---|---|
| (1, 1) | 1 | (5M+1)/6 | (M−1)/6 = ρ_M | no |
| (1, 5) | 5 | (M+1)/6 | (M−5)/6 = ρ_M | no |
| (5, 1) | 5 | (5M+1)/6 | (5M−5)/6 = 5ρ_M | yes |
| (5, 5) | 1 | (M+1)/6 | (5M−1)/6 | yes |

The "yes" cases are exactly s ≡ 5 (mod 6) — and for s = M = p this
reproduces the direct p² computation (δ = ρ_p at p ≡ 1; δ = (5p−1)/6 at
p ≡ 5). So the cap is 2ρ_s + 1 + [s ≡ 5 mod 6], depending only on the fiber
size s, not on the co-prime. Sharpness: the maximum is attained (a = 0 always
gives 2ρ_s + 1; when s ≡ 5 mod 6 any a in [M−δ, δ] gives 2ρ_s + 2).
Finally 2·cap(s) < s: for s ≡ 1 (mod 6), s ≥ 7: 4ρ_s + 2 = (2s+4)/3 < s;
for s ≡ 5 (mod 6), s ≥ 5: 4ρ_s + 4 = (2s+2)/3 < s; s = 3: 2 < 3. ∎

*Verification.* The harness asserts equality of the computed maximum with
the closed form for **every** prime-indexed fibration of every modulus
tested — {4, 9, 25, 49, 77, 91, 121, 169} — and the profile identity
B_u ∩ F_j = u^{-1}(A ∩ F_{uj}) over **all** units and fibers (e.g. 2028
checks at 169, 1210 at 121, both fibrations at 77 and 91): 0 failures. The
committed T-10 caps c_7 = 3, c_11 = 4, c_13 = 5 are the closed form at
s = 7, 11, 13; at N = 77 the two fibrations give caps 3 (fibers of size 7)
and 4 (fibers of size 11) — the co-prime independence, confirmed.

Lemma 1 is the corrected, closed form of T-10's L-FC cap (there computed
numerically per modulus). Two consequences are used repeatedly: a single
unit set never covers a full s-fiber (cap(s) < s), and two unit sets never
cover a full s-fiber for s ≥ 5 (2cap(s) < s).

## 4. Lemma 2 and Theorem 1 — three units never cover

> **Lemma 2 (Dirichlet arc-intersection floor).** Let N ≥ 31 with
> N ≢ 0 (mod 6). For every unit λ mod N: |A ∩ λA| ≥ 2.

*Proof.* c = ⌊(N−1)/6⌋ ≥ 5. By Dirichlet's approximation theorem with
Q = c there exist s ∈ [1, c] and k ∈ Z with |sλ/N − k| ≤ 1/(c+1), so
|sλ − kN| ≤ N/(c+1). Writing N = 6m + r (1 ≤ r ≤ 5, c = m):
N/(c+1) = (6m+r)/(m+1) = 6 − (6−r)/(m+1) < 6, so the integer |sλ − kN| ≤ 5
≤ c. Set x := sλ mod N. Then ||x||_N ≤ 5 ≤ c, so x ∈ A; x = λ·s with
||s||_N = s ≤ c, so x ∈ λA; and x ≠ 0 because λ is a unit and
0 < s ≤ c < N. Together with 0 ∈ A ∩ λA this gives |A ∩ λA| ≥ 2. ∎

(The committed spectra are far stronger — min |A ∩ λA| = 3, 5, 7, 9, 11 at
49, 77, 91, 121, 169, all re-verified exactly by the harness — but those
minima are computed facts; the floor of 2 is what a three-line proof gives,
and it is exactly what Theorem 1 needs. The lemma's honest range is
N ≥ 31 with N ≢ 0 mod 6, plus N ≥ 42 when 6 | N; it fails at N ≤ 30 and
N = 36 — the harness verifies the boundary, and at N = 25 the floor still
holds by the direct check of §2, while at N = 9 it fails with min = 1,
which is precisely why Z_9 is the exceptional world.)

> **Theorem 1 (no three unit balls cover).** Let N ≥ 31, N ≢ 0 (mod 6) — in
> particular N = p² with p ≥ 7, and N = pq with 5 ≤ p < q. No three unit
> bad sets cover Z_N. (Distinct or not: if two coincide the union has at
> most 2(2c+1) < N elements.)

*Proof.* Suppose B_{u_1} ∪ B_{u_2} ∪ B_{u_3} = Z_N, the three distinct.
Every B_u is u^{-1}A with |A| = 2c+1, and |B_i ∩ B_j| = |A ∩ λ_{ij}A| with
λ_{ij} = u_i u_j^{-1} a unit. Inclusion–exclusion gives

    Σ_{i<j} |B_i ∩ B_j| − |B_1 ∩ B_2 ∩ B_3| = 3(2c+1) − N =: Δ.

Since 0 lies in every bad set, the triple intersection t ≥ 1 and t ≤ each
pairwise intersection. By N mod 6:

* **N ≡ 4, 5, 0:** 3(2c+1) ≤ (N−4) + 3 < N — the union is too small.
* **N ≡ 2:** Δ = 1. If t ≥ 2 then Σ pairwise ≥ 3t > 1 + t = Σ — impossible;
  if t = 1 then Σ pairwise = 1 < 3 ≤ Σ (each pairwise ≥ 1) — impossible.
* **N ≡ 3:** Δ = 0, and the same two steps give Σ pairwise = t with
  t ≥ 3t or t = 1 < 3 — impossible.
* **N ≡ 1:** Δ = 2. If t ≥ 2 then Σ pairwise = 2 + t < 3t ≤ Σ — impossible;
  if t = 1 then Σ pairwise = 3, but each pairwise intersection is ≥ 2 by
  Lemma 2, so Σ ≥ 6 — impossible. ∎

Note p² ≡ 1 (mod 6) always — the square case is exactly the Lemma-2 case —
and pq ≡ 1 or 5. At N = 9 (≡ 3 mod 6) the escalation still applies without
Lemma 2, which is why three units do not cover Z_9 (verified); the Z_9
counterexample of §1.3 needs all four sets. Theorem 1 also subsumes the
3-set pure-unit cells of T-10's maps at 49–169 (harness-asserted) and at 25
(by the direct spectrum check).

## 5. The no-unit-assistance theorems

> **Theorem 2 (N = p², at most two units never essential).** Let p be an
> odd prime and let F be a family of at most 4 distinct bad sets at
> N = p² whose unit members number υ ≤ 2. Then F covers Z_{p²} only if its
> kernel subfamily covers.

*Proof.* If υ = 0 there is nothing to prove, so assume υ ≥ 1, and let
κ ≤ 4 − υ ≤ 3 be the number of kernel members. Each kernel
member is π^{-1} of a Z_p-ball of radius ρ_p (§1.1) — a union of exactly
2ρ_p + 1 complete p-fibers. If the kernel balls cover Z_p, the kernel
subfamily covers and we are done. Otherwise some complete fiber F_j lies in
no kernel member. The fiber has p elements; by Lemma 1 each unit member of F
meets F_j in at most cap(p) elements, so the υ ≤ 2 units cover at most
2·cap(p) < p of them (p ≥ 3). Some element of F_j is covered by nothing —
F does not cover. ∎

This extends T-10's L-FC in two ways: the (3 kernel, 1 unit) case is handled
(the three kernel balls may already cover all p fibers — then the unit is
not needed — or leave a full fiber — then it is helpless), and the cap is
the exact closed form. At p = 5 the theorem applies verbatim (2cap(5) =
4 < 5); at p = 3 it also applies, and indeed the Z_9 counterexample has
υ = 3.

> **Theorem 3 (N = pq, at most two units never essential).** Let N = pq with
> 5 ≤ p < q primes, and let F be a family of at most 4 distinct bad sets
> with υ ≤ 2 unit members, κ_p kernel-p members (multiples of p) and κ_q
> kernel-q members. Then F covers Z_{pq} only if its kernel subfamily
> covers — **with exactly one exception**: the cell
> (κ_p, κ_q, υ) = (1, 1, 2) survives the cap counting.

*Proof.* Work in the CRT picture Z_{pq} ≅ Z_p × Z_q. A kernel-p speed pa
satisfies ||pak||_{pq} = p·||ak||_q ≤ c ⟺ ||ak||_q ≤ ρ_q, so its bad set is
π_q^{-1}(Z_q-ball) — a union of 2ρ_q + 1 complete rows (q-coordinate
classes, each of p elements). Symmetrically kernel-q members are unions of
2ρ_p + 1 complete columns. Let Y ⊆ Z_q be the union of the κ_p row-balls and
X ⊆ Z_p the union of the κ_q column-balls. The kernel subfamily covers
exactly when X = Z_p or Y = Z_q (the complement of its union is the
rectangle X^c × Y^c). Suppose neither, so X^c, Y^c ≠ ∅ and the units must
cover the rectangle. Fix a row y ∈ Y^c: its |X^c| cells must lie in unit
bad sets, and a unit meets any row (a π_q-fiber of size p) in ≤ cap(p)
elements (Lemma 1), so |X^c| ≤ υ·cap(p); similarly |Y^c| ≤ υ·cap(q). On the
other hand |X^c| ≥ p − κ_q(2ρ_p + 1) and |Y^c| ≥ q − κ_p(2ρ_q + 1). Now:

* **κ_q = 0** (all kernel on the p-side): |X^c| = p > 2cap(p) ≥ υ·cap(p) —
  contradiction. (κ_p = 0 symmetric, using 2cap(q) < q for q ≥ 7.)
* **κ_q ≥ 1, υ = 1:** |X^c| ≥ p − (2ρ_p + 1) > cap(p) — the inequality
  p − (2ρ_p+1) > cap(p) reads (2p−2)/3 > (p+2)/3 for p ≡ 1 mod 6 (true for
  p ≥ 7) and (2p+2)/3 > (p+1)/3 for p ≡ 5 mod 6 (always) — contradiction.
  (If κ_q ≥ 2 then κ_p ≤ 1, and the same inequality on the q-side,
  q − (2ρ_q+1) > cap(q) for q ≥ 7, applies with κ_p = 1; κ_p = 0 is the
  previous case.)
* **κ_q ≥ 1, κ_p ≥ 1, υ = 2** (forcing κ_p = κ_q = 1): the necessary
  conditions become p − (2ρ_p+1) ≤ 2cap(p) and q − (2ρ_q+1) ≤ 2cap(q) —
  both consistent (for p ≡ 5 mod 6 the first is an *equality*:
  p − (2ρ_p+1) = (2p+2)/3 = 2cap(p)). No contradiction: this is the
  (1, 1, 2) cell. ∎

*Verification.* The harness reproduces the committed T-10 maps at
77 and 91 exactly (unit-assisted = 0 in every cell, both arities), which
covers the pq instances of all killed cells; the (1,1,2) cell is empty
there too (verified, not proved — consistent with Theorem 3's statement,
which does not claim it).

> **Theorem 4 (pure-unit reduction to the rigid zone).** Let N = sM with
> s, M primes ≥ 3, and suppose unit bad sets B_{u_1}, …, B_{u_m} cover Z_N.
> Then for each prime p | N, the reduced speeds cover the reduced world:
>
>     ∪_i B̃^{(N/p)}_{u_i mod (N/p)} = Z_{N/p},
>
> where B̃^{(m)} is the bad set at modulus m, T = 6. In particular, at
> N = p² a pure-unit covering can exist only if Z_p admits a covering by
> unit bad sets at T = 6 — i.e. only if p lies in the n = 5 rigid zone
> {7, 13, 17, 19, 37} within the verified range (primes ≤ 150).

*Proof.* The p-multiples {pj : j ∈ Z_{N/p}} must be covered. Now
u pj mod N = p·((u j) mod (N/p)) and ||p·y||_N = p·||y||_{N/p} for
y ∈ [0, N/p), so ||upj||_N ≤ c ⟺ ||uj||_{N/p} ≤ ⌊c/p⌋. And
⌊c/p⌋ = ⌊(N−1)/(6p)⌋ = ⌊(N/p)/6 − 1/(6p)⌋ = ρ_{N/p}, which is exactly the
Z_{N/p} bad-set radius at T = 6. So pj ∈ B_u ⟺ j ∈ B̃_{u mod (N/p)}, and
covering the p-multiples is covering Z_{N/p} by the reduced bad sets. ∎

*Verification.* The harness asserts the set identity
{k ∈ B_u : p | k} = p·B̃_{u mod (N/p)} for **every** unit u and **both**
prime divisors at N ∈ {25, 49, 77, 91, 121, 169} (e.g. 156 checks at 169):
0 failures. Base cases: every unit bad set at Z_3 and Z_5 is {0}.

**Corollaries.** (i) No pure-unit 4-covering exists at N = 77: the
reduction demands a covering of Z_11, and 11 is not rigid-zone — a
committed verified-only fact, now *proved*. (ii) The same holds at 121
(reduction to Z_11 — consistent with T-10's "no covering configuration of
any kind" there, which the corollary re-derives for the pure-unit cells).
(iii) At N = 25 no pure-unit covering exists (Z_5 unit bad sets are {0}).
(iv) At N = 49, 91, 169 (and 289, 361, 1369 if ever needed) the pure-unit
question is *reduced* to lifting the classified rigid-zone configurations
through the fibration — a finite, explicit question per prime (§6).

**What is now proved, assembled.** At N = p² (p ≥ 5 odd) and N = pq
(5 ≤ p < q): a covering family of at most 4 distinct bad sets with at most
two unit members has its kernel subfamily covering (Theorems 2, 3); three
unit balls never cover at all (Theorem 1); a pure-unit 4-family covers only
above a rigid-zone prime (Theorem 4). The remaining possibilities for unit
assistance are exactly: **(1 kernel + 3 units) at p²; (1 kernel-p, 1
kernel-q, 2 units), (1 kernel + 3 units), and (4 units above rigid-zone
primes) at pq; and rigid-zone lifts of (0, 4) at p².**

## 6. The residual, precisely: why it resists, and the attack route

The surviving cells share one mechanism, and it is the same one that made
the N = 9 counterexample possible: **when ≥ 3 units act on a fiber that
kernel sets leave fully uncovered, the counting budget is positive but
O(1):**

| configuration | per-uncovered-fiber budget Σ cap − p | status |
|---|---|---|
| (1K, 3U) at p², p ≡ 1 mod 6 | 3cap(p) − p = (p−1) + 3 − p = **2** | open (empty at p = 5, 7, 11, 13: verified) |
| (1K, 3U) at p², p ≡ 5 mod 6 | 3(2ρ+2) − (6ρ+5) = **1** | open; occupied at p = 3 (budget 3·1−3 = 0!) |
| (1, 1, 2) at pq, p ≡ 5 mod 6 | row-side slack p − (2ρ_p+1) − 2cap(p) = **0** (exact) | open |
| (0, 4) at p², p rigid-zone | (reduced, not counted) | open as lifts |

So in every surviving cell the units must *near-tile* each uncovered fiber
with total overlap ≤ 2 (usually ≤ 1, at N = 9 exactly 0), for every
uncovered fiber simultaneously. The within-fiber structure (obtained by
writing k = j + pt and u = u_0 + pu_1): a unit's footprint in fiber j is

    T_u(j) = α_u · S_{a_u(j)} + μ_u · j,    α_u = (u mod p)^{-1},
    μ_u = −α_u·u_1,   a_u(j) = (u mod p)·j,

an arithmetic progression in Z_p with step α_u and size n_{a} ∈ {2ρ, 2ρ+1,
2ρ+2} (the profile of Lemma 1), whose position moves *affinely in j*. The
covering condition per uncovered fiber j is: three APs, sizes ~p/3, cover
Z_p with overlap ≤ 2; and this must hold for all j in the far arc
|J| = p − 2ρ − 1 ≈ 2p/3, with the three positions moving at generically
different rates μ_u.

This is a rigidity problem — near-tilings of Z_p by three short APs that
must persist under a moving linear reparameterization — and it is the
identified attack route. A worked indication of the mechanism (partial,
labeled as such): at p ≡ 5 (mod 6) with budget 1, the three within-fiber
sizes must be (2ρ+2, 2ρ+2, 2ρ+1) or (2ρ+2, 2ρ+2, 2ρ+2), i.e. per fiber at
most one unit's reduction lands in the thin region; and if the three steps
are α_u = ±1 (intervals), a persistent near-tiling under moving offsets
forces the three slopes μ_u to coincide (a nonzero relative slope sweeps
|J| ≫ 3 admissible offset values), which forces the three speeds to share a
common multiplicative factor u_i = ū_i·w mod p² — collapsing the family to a
common scaling of reduced balls and pushing the question back to the
Z_p-level covering, which the rigid-zone classification controls. The
unproved step is the general-step case: three APs with steps α ∉ {±1}
near-tiling Z_p with overlap ≤ 1 for a whole arc of parameter values. Any
theorem of the form "Z_p admits no near-tiling by three APs of size ~p/3
with distinct steps, stable under affine reparameterization" closes the
(1K, 3U) cell outright, and its pq siblings with it.

Two further reductions sharpen the target list. For (1K, 3U) at p ≡ 1
(mod 6) (budget 2): every uncovered fiber j needs at least one unit
reduction ū_i·j in the kernel ball S = [−ρ, ρ] (the fat region of the
profile) — a covering incidence condition between three scaled copies of S
and the far arc J that already fails to contradict by O(1) in the counting
of §5. For the rigid-zone lifts of (0, 4): by Theorem 4 the reduced speeds
must form one of the *classified* covering configurations at Z_p (K7, K13,
K17, K19, K37 — all proved classifications), so the lift question is finite
per prime: for each classified configuration (ū_1, …, ū_4) and each lift
u_i = ū_i + p·s_i, decide whether the fiber-wise shifted-quasi-ball
conditions hold for all κ̄ ∈ Z_p. No uniform argument is known at present;
the first test-beds are N = 49 and 169 (both measured empty) and, untested
per the no-census directive, N = 289 and 361.

**Summary of the residual conjecture** (the corrected replacement for the
refuted original): *at N = p² and N = pq with p ≥ 5, no covering family of
at most 4 distinct bad sets exists in which the kernel subfamily does not
cover.* Proved for ≤ 2 unit members and for 3 pure unit members and in the
pure-unit 4 case outside the rigid zone; open in the cells above; false at
p = 3. The boundary instances are exact: the conjecture's failure at N = 9
is the unique small-prime counterexample class, and it sits in the (1K, 3U)
cell — the cell every remaining proof must target.

## 7. The kernel-world structure theorem (Priority 2)

Assembling L-FIB, L-HOM, L-CRT (T-10, proved), the rigid-zone
classifications K7/K13/K17/K19/K37 (T-5/T-6, proved), and this note's
Theorems 1–4:

> **Theorem S (kernel-world structure at composite moduli, n = 5).** Let
> N = p² (p ≥ 5 odd) or N = pq (5 ≤ p < q), and let F be a family of at
> most 4 distinct bad sets.
>
> 1. *(Fibration — proved.)* Every kernel bad set is π^{-1} of a unit bad
>    set at the reduced modulus with the *same threshold* (radius ρ); a
>    kernel-only family covers Z_N iff its transferred family covers the
>    reduced world (L-HOM; at pq: iff either side covers — L-CRT).
> 2. *(No unit assistance — proved in part.)* If F covers with at most two
>    unit members, its kernel subfamily covers (Theorems 2, 3; sole
>    exception cell (1,1,2) at pq). Three unit members never cover
>    (Theorem 1). A pure-unit family covers only above rigid-zone primes
>    (Theorem 4).
> 3. *(Consequence, conditional on the residual conjecture.)* Covering at
>    composite N ⟺ the kernel subfamily covers ⟺ the transferred
>    prime-level world covers ⟺ the prime is in the rigid zone with a
>    classified configuration. Within the verified prime range this makes
>    the composite covering problem a *closed* classification: at p²,
>    3-kernel coverings exist iff p ∈ {7, 13} (committed 3-coverable
>    primes) and 4-kernel coverings iff p is rigid-zone; at pq, kernel
>    coverings exist iff one side transfers to a covering (L-CRT).

Item 3 is the "theorem, not census" the program wanted: with the residual
conjecture proved, the composite world contributes *no* covering
configurations beyond the rigid-zone pullbacks. Without it, the statement
holds with four named exceptions (§6) — the honest current form. Note also
what Theorem S does **not** say: it classifies coverings *given* the
rigid-zone classifications, and the successor conjecture (no capable prime
beyond 37) remains a hypothesis outside the verified range. The structure
theorem is a statement about the reformulation's covering problem — a
bounded contribution, in the sense recorded in the program's strategic
honesty notes; it is not a step toward LRC itself.

## 8. Verification record (the harness, `scripts/nua_proof_verify.py`)

* **Anchors.** All committed T-10 covering maps reproduced exactly, per
  cell: 49 — 3-set 1, 4-set 24 (3 + 21); 77 — 1, 38 (3 + 5 + 30);
  91 — 3, 138 (3 + 6 + 36 + 6 + 15 + 72); 121 — 0, 0; 169 — 2, 171
  (15 + 156); unit-assisted = 0 everywhere. Hand anchors: (7,14,21) covers
  Z_49, (7,14,35) does not; the N = 9 counterexample asserted cell-exactly
  (1 covering, unit-assisted, family (1,2,3,4)).
* **Lemma 1.** Closed-form cap equality at every prime fibration of
  {4, 9, 25, 49, 77, 91, 121, 169} (caps: 1, 1, 2, 3, {3, 4}, {3, 5}, 4, 5);
  profile identity over all units × fibers: 6,174 checks, 0 failures.
* **Lemma 2 / spectra.** Committed minima 3/5/7/9/11 at 49/77/91/121/169
  reproduced; Dirichlet case table verified over N ∈ [31, 400) with
  failures exactly at N ≤ 30 and N = 36; N = 25 floor (min 3) by direct
  check; N = 9 min = 1 recorded (the exception).
* **Theorem 1.** All 3-unit cells empty at every modulus tested, including
  the boundary moduli.
* **Theorem 4.** Reduction identity over all units at all six composite
  moduli (both prime sides at 77/91): 592 checks, 0 failures; Z_3, Z_5 base
  cases asserted.
* **Boundary maps.** N = 4: no covering; N = 9: exactly one covering
  (1K + 3U, unit-assisted); N = 25: no covering in any cell (286 + 1001
  combos). Total harness runtime ≈ 55 s; results in
  `scripts/out_nua_verify.json`.

Per the directive, no new moduli were scanned (289/361 are *not* computed),
no V-extension was run, and no capability scan beyond the committed record.

## 9. Consequences for the program

The reformulation's covering problem at composites is now: fibration
(proved) + no-unit-assistance (proved in the ≤ 2-unit, 3-unit, and
pure-unit-reduction parts) + rigid-zone classifications (proved) + **four
named residual cells with sharp budgets and an identified rigidity
attack**. The month-scale picture the reviewer drew is confirmed by the
mathematics: the residual is small, explicit, and genuinely hard — every
counting bound in it sits within O(1) of tight, which is why censuses kept
returning zero without explaining it. The next unit of work on Priority 1
is the AP near-tiling rigidity lemma of §6 (or a transversal/witness
formulation of the same statement); nothing else in the composite world
remains open below the rigid-zone lift question.

Strategically: if the residual cells close, the kernel-world structure
theorem becomes unconditional, and the reformulation is *closed as a
bounded contribution* — a lossless re-encoding of LRC with a complete
classification of its covering obstructions at the tested rungs. That is
the decision point the program's honesty notes anticipated: submit the
structural paper (v2 with the fibration and this note's theorems), and
either end the track there or pivot to R2 with the reformulation as the
bridge between the two languages. What is *not* at the decision point:
more data. The directive stands — proofs, not computations.
