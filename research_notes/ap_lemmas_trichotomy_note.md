# The AP Lemmas after the Wrap-Around Repair — Trichotomy, Case Analyses, Verification

**Program:** LRC pair-sum-lattice reformulation, kernel-world track, p²-cell
closure. **Task (2026-10-04, post-assessment directive):** the human's
AP-lemma sketch has a real gap — the span bound (s−1)|d| ≤ L fails for
wrapping APs (counterexample at p = 23: the step-8 AP of size 8 inside the
arc [5, 22]); the corrected analysis converged on a **trichotomy** for APs
contained in the complement arc. The directive: stop hand-analyzing, run the
exhaustive computation (Lemma A, Lemma B, the trichotomy claim, the
rigid-region data, all primes p ≤ 250), then write the four-section
deliverable — trichotomy lemma / Lemma A case analysis / Lemma B case
analysis / computational verification — and clarify whether the human's
proof structure survives the repair. **Artifacts:**
`scripts/ap_lemmas_verify.py` (+ `scripts/ap_inventory_dump.py`),
`scripts/out_ap_lemmas.json`, `scripts/ap_lemmas_run.log` (runtime 13.3 s;
every anchor passes; the search logic is cross-validated against an
independent brute force at p ∈ {7, 11, 13, 17, 19}). **Result:** both
lemmas verified at scale; the trichotomy is proved and explains the
boundary exceptions exactly; **Lemma B is proved in full for k ≥ 3
(p ≥ 23)**; **Lemma A is proved for k ≥ 7 (p ≥ 43)**, verified for
k ∈ {3, 5, 6} (p = 19, 31, 37), and is **FALSE at p ∈ {7, 13}**, where
±-distinct coverings exist — a correction to the framework's belief that
the boundary solutions are non-±-distinct; those two primes are absorbed
by direct cell checks (re-verified independently here). The human's
architecture survives the repair; Steps 1–2 are subsumed by the
trichotomy (§0.2).

Conventions. p prime, p = 6k+1 (class ε = 1) or p = 6k+5 (ε = 5). An
**AP** of difference d and size s is A = {a + jd mod p : j ∈ [0, s)}; the
**±-rep** of d is min(d mod p, p − d mod p) ∈ [1, (p−1)/2] (A is unchanged
under d → p − d, so ±-reps enumerate AP sets exactly once). Differences
are **±-distinct** when pairwise distinct up to sign. **Normalization:**
scaling by d₁⁻¹ and translating (both preserve sizes, covering, and
±-distinctness) lets us designate one AP — always one of maximum size,
chosen per case — as A₁ = [0, s₁−1] with d₁ = 1; then
I₁ := [s₁, p−1] (the complement arc), q := p − s₁, L := p−1−s₁ (span).
Critical sizes: ε = 5: s ∈ {2k+1, 2k+2}; ε = 1: s ∈ {2k, 2k+1, 2k+2}.
**Lemma A** (ε = 1): no three ±-distinct APs with sizes in the critical
set, Σsᵢ − p ∈ {0, 1, 2} (the five multisets {2k,2k,2k+1},
{2k,2k+1,2k+1}, {2k,2k,2k+2}, {2k,2k+1,2k+2}, {2k+1,2k+1,2k+1}),
covering Z_p. **Lemma B** (ε = 5): same with Σsᵢ − p ∈ {0, 1} (the two
multisets {2k+1,2k+2,2k+2} and {2k+2,2k+2,2k+2}). These are exactly the
t-space statements produced by the fiber reduction (§5): a covering of
Z_{p²} by the kernel bad set plus three unit bad sets forces, in **every**
off-kernel fiber, a covering of the fiber's t-space (≈ Z_p) by three APs
with differences uᵢ⁻¹ mod p (pairwise ±-distinct iff the units are) and
sizes L(cᵢ) ∈ {2k, 2k+1} (ε = 1) resp. {2k+1, 2k+2} (ε = 5) — the
capacity formula, re-verified from scratch in §4.

---

## 0. Executive summary and the subsumption answer

### 0.1 Findings

| statement | status |
|---|---|
| Human's Step 1/Step 2 span bound for APs in an arc | **valid only in the non-wrapping regime** (d ≤ s₁); replaced by Lemma T, which covers both regimes |
| Trichotomy (contained critical-size APs have d = 2 or (p−1)/2) | **proved** for k ≥ 2 (ε = 5) and k ≥ 3 (ε = 1); exceptions exactly at p = 11 (d ∈ {3,4}) and p = 13 (d = 5), both *explained by the proof's own inequalities* (§1.3); p = 7 and p = 17 are clean |
| Two-block structure at d = (p−1)/2 | **proved and verified**: blocks of ⌈s/2⌉, ⌊s/2⌋ consecutive integers, right endpoints offset by exactly (p−1)/2 (0 mismatches at every prime ≤ 250) |
| Lemma B (p ≡ 5 mod 6) | **verified 0 solutions at every prime 11 ≤ p ≤ 250 (25 primes); proved for k ≥ 3 (p ≥ 23)**; p = 11, 17 by exhaustive verification; p = 5 degenerate (§2.4) |
| Lemma A (p ≡ 1 mod 6) | **verified 0 solutions at every prime 19 ≤ p ≤ 250 (23 primes); proved for k ≥ 7 (p ≥ 43)**; p = 19, 31, 37 verified; **FALSE at p = 7 (50 raw configurations) and p = 13 (2)** — ±-distinct coverings exist, exhibited in §3.4 |
| Human's claim "p = 7, 13 exhibit the non-±-distinct solutions" | **corrected**: at p = 7 and p = 13 there exist *±-distinct* coverings (at p = 13 even with reduction sizes {2k+1}³) — the boundary is not merely the rigid families |
| Relaxed (non-±-distinct) landscape | exactly **four families per size multiset at every p ≡ 5** — differences (1,1,1), (1,2,2), (1,2,(p−1)/2), (1,(p−1)/2,(p−1)/2) — *every one contains a ±-collision*; the lemmas live and die by ±-distinctness precisely |
| Fiber capacity + t-space AP structure | **re-verified from scratch** (p ≤ 37, all units × all fibers): sizes match the on/off-ball formula with 0 violations; every fiber bad set is an AP with difference u⁻¹ (0 violations) |
| (1K,3U) cells at p² | **empty, independently re-verified by full brute force over unit ±-classes mod p² for p ≤ 31** (9 primes; the committed T-12c structured search covered p ≤ 47) |
| Rigid zone (n = 5 side: min ±-distinct dilates of B_k covering Z_p) | 3 dilates at {7, 13}; 4 at {17, 19, 37}; ≥ 5 elsewhere — **exactly the committed zone**, re-verified p ≤ 61 |

### 0.2 Does the human's proof structure survive the repair?

The assessment asked precisely this: if the repair requires the trichotomy
anyway, Step 1 is redundant; if it only patches Steps 1–2 locally, the
structure survives. The answer from the completed analysis:

- **Step 1 is subsumed.** Case 1 of Lemma T's proof *is* Step 1, done
  correctly: for d ≤ s₁ the walk cannot wrap (a wrap lands below s₁),
  so the span bound (s−1)d ≤ L is valid there and forces d = 2. The
  human's error was applying the span bound in the wrapping regime, where
  it is false. The trichotomy's Case 2 replaces the invalid half of Step 1
  with the two-lattice descent argument and delivers more: not just
  "which d are possible" but the exact block structure and placement.
- **Step 2 is subsumed.** The gap-sum/parity bookkeeping of A₂ ∩ I₁ is
  replaced by the explicit d = 2 inventory (parity classes and their
  endpoint truncations, with exact placement windows — §1.2), which is
  what the case analyses actually consume.
- **Step 3 survives intact.** The architecture — normalize, pass to the
  complement arc, do overlap accounting (|A₁∩A₂| ≤ ov etc.), split into
  partition / overlap cases, kill each by structure — is exactly the
  backbone of §2 and §3, and it now runs on verified ground.

So the human's framework is reusable at the architecture level; the
implementation of Steps 1–2 is replaced wholesale by Lemma T, and the
lemma statements are re-scoped as above. One load-bearing element becomes
visible in the process: the ±-distinctness hypothesis is essential (the
parity family (1,2,2) exists at every p ≡ 5 and is killed by nothing
else), and its justification sits in the reduction, not the lemma — see
§5 for the precise residual this leaves.

---

## 1. The Trichotomy Lemma

### 1.1 Statement

**Lemma T.** Let p = 6k+5 with k ≥ 2 (resp. p = 6k+1 with k ≥ 3). Let
A ⊆ Z_p be an AP with ±-rep difference d ∈ [2, (p−1)/2] and size s with
2k+1 ≤ s ≤ 2k+2 (resp. 2k ≤ s ≤ 2k+2). Let I₁ = [s₁, p−1] with
2k+1 ≤ s₁ ≤ 2k+2 (resp. 2k ≤ s₁ ≤ 2k+2). If A ⊆ I₁ then:

1. **d = 2**, and A = {a, a+2, …, a+2(s−1)} is a non-wrapping integer AP
   with a ∈ [s₁, p−1−2(s−1)]; or
2. **d = (p−1)/2**, and A is the **two-block set**
   [a−⌈s/2⌉+1, a] ∪ [a+d−⌊s/2⌋+1, a+d] (mod p): two blocks of ⌈s/2⌉ and
   ⌊s/2⌋ consecutive integers whose right endpoints differ by exactly d.

In particular the only differences that ever occur are 2 and (p−1)/2.

### 1.2 Proof

Write L := p−1−s₁ (note ε = 5: L ≤ 4k+3; ε = 1: L ≤ 4k). Walk the AP as
x_j = a + jd mod p, so x₀ = a and consecutive walk points differ by
exactly d (mod p).

**Case 1: d ≤ s₁ (the non-wrapping regime — the corrected Step 1).**
If some non-terminal x_j exceeds p−1−d, the next point is
x_j + d − p ∈ [0, d−1] ⊆ [0, s₁−1], outside I₁ — forbidden. (The terminal
point x_{s−1} may lie in [p−d, p−1]; it has no successor.) Hence the walk
never wraps: x_j = a + jd as integers, all in [s₁, p−1], so

  (s−1)·d = x_{s−1} − a ≤ p−1−s₁ = L.

With s−1 ≥ 2k (ε = 5): 2k·d ≤ 4k+3, so d ≤ 2 + 3/(2k) < 3, i.e. d = 2.
With s−1 ≥ 2k−1 (ε = 1): (2k−1)d ≤ 4k, so d ≤ 2 + 2/(2k−1) < 3 for
k ≥ 3, i.e. d = 2. This proves alternative 1, with the placement window
a ∈ [s₁, p−1−2(s−1)] read off from x₀ ≥ s₁ and x_{s−1} ≤ p−1.

**Case 2: d ≥ s₁+1 (the wrapping regime).** Then d > p/3 (ε = 5:
d ≥ 2k+2 > 2k+5/3; ε = 1: d ≥ 2k+1 > 2k+1/3) and d ≤ (p−1)/2 < p/2, so
2d < p < 3d and

  r := p − 2d  satisfies  1 ≤ r < s₁, r odd.

(The upper bound: r ≤ p−2s₁−2 ≤ 2k−1 < s₁ in both classes.) Since
2d ≡ −r (mod p), the walk satisfies x_{j+2} ≡ x_j − r (mod p): the
even-index and odd-index subsequences each *descend by r per index pair*.

*No-wrap claim.* From a current point y ∈ I₁, the next point of the same
subsequence is y − r (mod p). If y ≥ s₁ + r then y − r ∈ [s₁, p−1] as an
integer (no wrap); if y ∈ [s₁, s₁+r) then y − r ∈ [0, s₁) — outside I₁,
forbidden; and y − r < 0 cannot occur since y ≥ s₁ > r. So both
subsequences descend inside [s₁, p−1] as integers, never wrapping; in
particular every non-terminal element of a subsequence is ≥ s₁ + r.

*Span bound, sub-case (i): the first step does not wrap* (x₁ = x₀+d ≤
p−1). Let M₁ := ⌈s/2⌉−1, M₂ := ⌊s/2⌋−1. Then max A ≥ x₁ = x₀ + d and
min A ≤ x_{2M₁} = x₀ − M₁r (the last even-index point, an integer by the
no-wrap claim), so

  span(A) = max − min ≥ d + M₁r,  hence  d + M₁r ≤ L.

ε = 5: M₁ = k and d ≥ 2k+2 give kr ≤ L − d ≤ 4k+3−2k−2 = 2k+1, so
r ≤ 2 + 1/k < 3, whence r = 1 (r odd). ε = 1: M₁ ≥ k−1 and d ≥ 2k+1
give (k−1)r ≤ 2k−1, so r ≤ 2 + 1/(k−1) < 3 for k ≥ 3, whence r = 1.

*Span bound, sub-case (ii): the first step wraps* (x₀ ∈ [p−d, p−1],
x₁ = x₀+d−p ∈ [s₁, d−1]). Then max A ≥ x₀ and min A ≤ x_{2M₂+1} =
x₁ − M₂r, so

  span(A) ≥ (x₀ − x₁) + M₂r = (p − d) + M₂r.

ε = 5: p−d ≥ 3k+3 and M₂ ≥ k−1 give (k−1)r ≤ k, so r ≤ 1+1/(k−1) < 2,
whence r = 1. ε = 1: p−d ≥ 3k+1 and M₂ ≥ k−1 give (k−1)r ≤ k−1, whence
r = 1.

Both sub-cases force r = 1, i.e. p = 2d+1 and d = (p−1)/2. With r = 1 the
even subsequence is {a − m : m ∈ [0, ⌈s/2⌉)} = [a−⌈s/2⌉+1, a] and the odd
subsequence is [a+d−⌊s/2⌋+1, a+d] (mod p): two blocks of ⌈s/2⌉ and ⌊s/2⌋
consecutive integers whose right endpoints a and a+d differ by exactly d.
This proves alternative 2. ∎

### 1.3 The boundary exceptions, explained by the proof itself

The two inequalities that force d = 2 and r = 1 slacken at small k, and
the slack is exactly where the computational exceptions sit:

- **p = 11 (k = 1, ε = 5):** Case 1 gives d ≤ (4k+3)/(2k) = 3.5, so
  d = 3 is an extra *non-wrapping* possibility; Case 2 at k = 1 allows
  r = 3 (from kr ≤ 2k+1), giving d = (p−3)/2 = 4. Census: contained
  differences {2, 3, 4, 5} — exactly {2, 3, 4, (p−1)/2}. ✔
- **p = 13 (k = 2, ε = 1):** Case 2 at k = 2 allows r = 3 (from
  (k−1)r ≤ 2k−1, i.e. r ≤ 3), giving d = (p−3)/2 = 5; the fit is exact:
  d + M₁r = 5 + 1·3 = 8 = L for s = 2k, s₁ = 2k — the AP
  {4, 7, 9, 12} ⊆ [4, 12] of the census. ✔
- **p = 7 (k = 1, ε = 1):** d ≤ 4k/(2k−1) = 4 and ±-reps stop at 3 =
  (p−1)/2: no exception — the census reads {2, 3} = {2, (p−1)/2}. ✔
- **p = 17 (k = 2, ε = 5):** both inequalities are already tight enough
  (r ≤ 2.5 in sub-case (i), r ≤ 2 in (ii)); the census is exactly
  {2, 8} = {2, (p−1)/2} — the trichotomy *holds* at p = 17, as the
  assessment's hand-check found. ✔

### 1.4 Corollaries consumed by the case analyses

**Corollary T1 (exact placements).** In the situation of Lemma T, at
ε = 5 with s₁ = 2k+2 (L = 4k+2): the size-(2k+2) forms are *unique* in
each class — d = 2: only E := {2k+2, 2k+4, …, 6k+4} (span 4k+2 = L forces
a = 2k+2); d = (p−1)/2: only T := [2k+2, 3k+2] ∪ [5k+4, 6k+4] (the block
windows a ≥ 3k+2 and a+3k+2 ≤ 6k+4 force a = 3k+2). The size-(2k+1) forms:
d = 2: a ∈ {2k+2, 2k+3, 2k+4} (three parity truncations); d = (p−1)/2: the
two block orders, [2k+2, 3k+2] ∪ [5k+5, 6k+4] and
[2k+3, 3k+2] ∪ [5k+4, 6k+4]. (Verified: the p = 239 census reproduces
exactly these 1+1+3+2 forms.)

**Corollary T2 (Run Lemma).** Let A be an AP of size σ with ±-rep
difference d ≢ ±1 (mod p), let 2(σ−1) < p, and suppose A contains ρ ≥ 2
consecutive residues. Then ν(d⁻¹) ≤ (σ−1)/(ρ−1), where ν(·) denotes the
±-rep. In particular ρ = σ−1 forces ν(d⁻¹) = 1, i.e. d ≡ ±1.

*Proof.* The run c, c+1, …, c+ρ−1 lies in A; consecutive-in-value pairs
have AP-index difference δ ≡ d⁻¹ (mod p) with |δ| ≤ σ−1. Since
2(σ−1) < p, at most one of {ι, ι−p} (ι := d⁻¹ mod p) lies in
[−(σ−1), σ−1], so every such δ equals one fixed δ₀ with |δ₀| = ν(d⁻¹) ≥ 1.
The run's indices form an arithmetic progression with step δ₀ inside
[0, σ−1], so (ρ−1)·ν(d⁻¹) ≤ σ−1. ∎

**Corollary T3 (Run Lemma at spacing 2).** With the same hypotheses,
if A contains ρ ≥ 2 residues in arithmetic progression of value-spacing
γ (i.e. c, c+γ, …, c+(ρ−1)γ), then ν(γ·d⁻¹) ≤ (σ−1)/(ρ−1). In
particular, containing a full parity class of ρ = σ−1 points forces
ν(2d⁻¹) = 1, i.e. **d ≡ ±2**.

*Proof.* Identical, with δ ≡ γ·d⁻¹ for each value-consecutive pair of the
spacing-γ run; the "at most one candidate" step needs 2(σ−1) < p, which
holds for all our sizes (4k+2 < 6k+5 and 4k+2 < 6k+1). ∎

T3 is the clean replacement for the previous session's messy parity-walk
attempt at the "one foot in A₁" case: an AP that contains a full parity
class of size σ−1 *must* have difference ±2 — no walk analysis needed.

---

## 2. Lemma B — complete case analysis (p = 6k+5)

**Lemma B.** *Let p = 6k+5 be prime, k ≥ 1. There do not exist three APs
in Z_p with pairwise ±-distinct differences, sizes sᵢ ∈ {2k+1, 2k+2},
Σsᵢ − p ∈ {0, 1}, and A₁ ∪ A₂ ∪ A₃ = Z_p.*

**Status: proved for k ≥ 3 (p ≥ 23); k = 2 and k = 1 (p = 17, 11)
verified by exhaustive search (0 solutions); k = 0 (p = 5) degenerate
(§2.4).** The proof consumes Lemma T and Corollaries T1–T3 only.

### 2.1 Setup and overlap accounting

Both multisets contain a size-(2k+2) AP, so we may normalize one of them:
scale by its difference and translate to get

  A₁ = [0, 2k+1] (d₁ = 1),  I₁ = [2k+2, 6k+4],  q = 4k+3,  L = 4k+2.

Let ov = Σsᵢ − p ∈ {0, 1}. By inclusion–exclusion,
ov = |A₁∩A₂| + |A₁∩A₃| + |A₂∩A₃| − |A₁∩A₂∩A₃|, so **|A₁∩Aᵢ| ≤ ov for
i = 2, 3** (both other pairwise terms dominate the triple term). Write
Bᵢ := Aᵢ ∩ I₁ and tᵢ := |Bᵢ|; then t₂ + t₃ − |B₂∩B₃| = q, tᵢ ≤ sᵢ, and
tᵢ ≥ sᵢ − ov.

### 2.2 The partition case (ov = 0, multiset {2k+1, 2k+2, 2k+2})

Here t₂ = s₂ and t₃ = s₃: both A₂, A₃ ⊆ I₁ and A₂ ⊔ A₃ = I₁, with sizes
{2k+2, 2k+1}. By Lemma T (k ≥ 2 suffices here), d₂, d₃ ∈ {2, 3k+2};
±-distinctness forces {d₂, d₃} = {2, 3k+2}. By Corollary T1 the
size-(2k+2) forms are E (d = 2) and T (d = 3k+2), and the size-(2k+1)
forms are three parity truncations and two two-blocks.

- (2k+2-form = E, 2k+1-form = two-block): the two-block must equal
  I₁ ∖ E = O, the odd class of I₁ — 2k+1 points with *no two consecutive*.
  A two-block has blocks of ⌈(2k+1)/2⌉ = k+1 ≥ 4 consecutive integers.
  Contradiction. ✔
- (2k+2-form = T, 2k+1-form = parity truncation): the parity form must
  equal I₁ ∖ T = [3k+3, 5k+3] — 2k+1 *consecutive* integers. A parity
  truncation has no two consecutive. Contradiction. ✔
- (E, parity) or (T, two-block): both differences equal — ±-collision.
  ✔

All four combinations die; note that the surviving structural pair
(E, O) — the exact relaxed parity family — is killed *only* by
±-distinctness. This is where the hypothesis earns its keep.

### 2.3 The overlap-1 case (multiset {2k+2, 2k+2, 2k+2})

The accounting t₂ + t₃ − |B₂∩B₃| = 4k+3 with tᵢ ∈ {2k+1, 2k+2} leaves
exactly two possibilities:

**(α) t₂ = t₃ = 2k+2, |B₂∩B₃| = 1.** Both contained, so both are E or T;
±-distinctness forces {E, T}. Now |E ∩ T| = #even elements of
T = (evens of [2k+2, 3k+2]) + (evens of [5k+4, 6k+4]) ≥ 2·⌈(k+1)/2⌉ ≥
k+1 ≥ 4 > 1 (k ≥ 3). Contradiction with |B₂∩B₃| = 1. ✔

**(β) one contained, the other with exactly one foot in A₁.** Say A₃ ⊆ I₁
and t₂ = 2k+1, so B₂ ⊔ A₃ = I₁ and A₂ = B₂ ∪ {z} with z ∈ A₁.

- *A₃ = T.* Then B₂ = I₁ ∖ T = [3k+3, 5k+3], so A₂ contains 2k+1 =
  s₂−1 consecutive residues; Corollary T2 forces d₂ ≡ ±1, excluded by
  ±-distinctness from d₁ = 1. ✔
- *A₃ = E.* Then B₂ = I₁ ∖ E = O, the full odd class of I₁ (2k+1 =
  s₂−1 points at spacing 2); Corollary T3 forces d₂ ≡ ±2, a ±-collision
  with d₃ = 2. ✔ (This is the general-k form of the p = 17 rigid family,
  which realizes it with d₂ = d₃ = 2 — excluded here by hypothesis.)

The case split is exhaustive (t₂ = t₃ = 2k+1 would need |B₂∩B₃| = −1).
This completes the proof of Lemma B for k ≥ 3. ∎

### 2.4 The boundary primes

- **k = 2 (p = 17) and k = 1 (p = 11):** exhaustive search (§4), 0
  solutions each. The proof does not close these because at k ≤ 2 the
  *one-foot* census is not yet clean: with |A₁∩A₂| ≤ 1 the viable
  differences at p = 17 are {2, 4, 5, 7, 8} and at p = 11 are
  {2, 3, 4, 5} — the middle values are the small-k slack of Lemma T's
  inequalities applied to APs with one foot (the same slack as §1.3, one
  notch looser because a single point may leave I₁). At k ≥ 3 the
  one-foot census is clean — {2, (p−1)/2} at every prime 23 ≤ p ≤ 250
  (verified; this is the *extended trichotomy*, and it is what lets
  §2.3β run).
- **k = 0 (p = 5):** sizes {1, 2}; a size-1 AP has no defined difference,
  and in any case only two ±-classes exist mod 5, so three units can
  never be pairwise ±-distinct — the hypothesis of the lemma is
  unsatisfiable. The corresponding cell (N = 25) is verified empty
  directly (§4; committed T-11).

---

## 3. Lemma A — case analysis (p = 6k+1)

**Lemma A.** *Let p = 6k+1 be prime. There do not exist three APs in Z_p
with pairwise ±-distinct differences, sizes sᵢ ∈ {2k, 2k+1, 2k+2},
Σsᵢ − p ∈ {0, 1, 2}, and A₁ ∪ A₂ ∪ A₃ = Z_p.*

**Status: FALSE at k ∈ {1, 2} (p = 7, 13) — counterexamples exhibited in
§3.5; verified 0 solutions at k ∈ {3, 5, 6} (p = 19, 31, 37); proved for
k ≥ 7 (p ≥ 43) — the multisets M1–M3 unconditionally (Lemma T +
Corollaries T1–T3 only), M4–M5 modulo Proposition E below (verified).**
Note k = 4, 8, 9 give non-primes, so the primes are covered exactly.

### 3.1 Setup and the containment inventory

Overlap accounting as in §2.1 (|A₁∩Aᵢ| ≤ ov). The five multisets and the
normalization of one maximum-size AP as A₁ = [0, s₁−1]:

| multiset | ov | s₁ | I₁ | q = p−s₁ | L |
|---|---|---|---|---|---|
| M1 {2k,2k,2k+1} | 0 | 2k+1 | [2k+1, 6k] | 4k | 4k−1 |
| M2 {2k,2k+1,2k+1} | 1 | 2k+1 | [2k+1, 6k] | 4k | 4k−1 |
| M3 {2k,2k,2k+2} | 1 | 2k+2 | [2k+2, 6k] | 4k−1 | 4k−2 |
| M4 {2k,2k+1,2k+2} | 2 | 2k+2 | [2k+2, 6k] | 4k−1 | 4k−2 |
| M5 {2k+1,2k+1,2k+1} | 2 | 2k+1 | [2k+1, 6k] | 4k | 4k−1 |

**Containment inventory** (Lemma T plus the exact placement windows;
verified against the p = 241 census, which reproduces it exactly):

- At s₁ = 2k+1 (L = 4k−1): size-(2k) forms: d = 2 with a ∈ {2k+1, 2k+2}
  — the odd class P_o = {2k+1, …, 6k−1} and the even class
  P_e = {2k+2, …, 6k} (span 4k−2); d = 3k with the exact fit
  d + M₁r = 3k + (k−1) = 4k−1 = L forcing a = 3k — the unique two-block
  T_A = [2k+1, 3k] ∪ [5k+1, 6k]. Size-(2k+1) forms: **none** — both
  spans 4k (d = 2) and 4k (d = 3k: 3k + k) exceed L = 4k−1.
- At s₁ = 2k+2 (L = 4k−2): size-(2k) forms: d = 2 with span 4k−2 = L
  exact, a = 2k+2 — the unique even class
  E' = {2k+2, …, 6k}; d = 3k: span 4k−1 > L, **none**. Size-(2k+1)
  forms: **none** (spans 4k).

So at either normalization *no size-(2k+1) AP is ever contained in I₁*,
and at s₁ = 2k+2 the only contained form of any size is E'. This
simplifies everything: in M2 and M4 the size-(2k+1) AP always carries at
least one foot.

### 3.2 M1, M2, M3 — unconditional kills

**M1 (partition, {2k, 2k} ⊔ I₁).** Both contained, so both are in
{P_o, P_e, T_A}; ±-distinctness forces one d = 2 and one d = 3k.
(P, T_A): the two-block must equal the complement of a parity class —
another parity class, which has no two consecutive, while T_A's blocks
are k consecutive. ✘ (T_A, P): the parity form must equal
I₁ ∖ T_A = [3k+1, 5k], 2k consecutive integers. ✘ (P_o, P_e): ±-collision
(both d = 2). ✘ **M1 dead.**

**M2 ({2k, 2k+1}, ov = 1).** The 2k+1-sized AP has t₃ ≤ 2k (never
contained), so the accounting t₂ + t₃ − |B₂∩B₃| = 4k leaves exactly
(t₂, t₃) = (2k, 2k) with B₂ ⊔ B₃ = I₁: A₂ is a contained 2k-form, A₃
has one foot. If A₂ = P (parity): B₃ = the other parity class (2k =
s₃−1 points at spacing 2), so Corollary T3 forces d₃ ≡ ±2 — ±-collision
with d₂ = 2. ✘ If A₂ = T_A: B₃ = [3k+1, 5k] (2k = s₃−1 consecutive), so
Corollary T2 forces d₃ ≡ ±1 — excluded by ±-distinctness from d₁ = 1. ✘
**M2 dead.**

**M3 ({2k, 2k}, ov = 1, s₁ = 2k+2).** Contained 2k-forms: only E'.
The accounting t₂ + t₃ − |B₂∩B₃| = 4k−1 with tᵢ ≤ 2k leaves
(t₂, t₃) = (2k, 2k) with |B₂∩B₃| = 1, or (2k, 2k−1) with ∩ = 0.
- Both contained: both are E' — ±-collision (and A₂ = A₃ = E' would not
  cover I₁ anyway). ✘
- One contained (E'), the other with one foot: the footed AP contains
  B = I₁ ∖ E' = O', the odd class of I₁ — 2k−1 points at spacing 2 with
  σ = 2k, so T3 gives (2k−2)·ν(2d⁻¹) ≤ 2k−1, ν < 2, d ≡ ±2 — collision
  with d₂ = 2. ✘ **M3 dead.**

### 3.3 M4, M5 — kills modulo Proposition E (verified)

**Proposition E (extended trichotomy, one/two feet).** *At p = 6k+1 with
k ≥ 7 (and at p = 6k+5 with k ≥ 3): any AP of a critical size with at
most ov* ≤ 2 points outside I₁ = [s₁, p−1] (s₁ critical) has
d ∈ {2, (p−1)/2}.* — *Verification:* complete census, every prime in
range (§4): the viable-difference set with |A₁∩A| ≤ 2 is exactly
{2, (p−1)/2} for p ≡ 1, p ≥ 43, and with ≤ 1 foot for p ≡ 5, p ≥ 23.
*Proof sketch:* for middle d ∈ [3, s₁], every wrap lands in [0, d−1] ⊆ A₁
and the re-entry climb (steps of d from the landing point up to s₁)
contributes further A₁-points; a three-distance count bounds the number
of A₁-points below by ≈ s·s₁/p − 2, which exceeds 2 precisely in the
range where the census is clean (s·s₁/p ≈ (2k)(2k)/(6k) = 2k/3 > 2 for
k ≥ 4, matching the observed boundary k ≥ 5 for one foot, k ≥ 7 for
two). We flag honestly: the full write-out of this count is not included;
M4/M5 below consume E as a verified input.

**M4 ({2k, 2k+1}, ov = 2, s₁ = 2k+2).** q = 4k−1; tᵢ ≤ 2k (nothing of
size 2k+1 is contained); B₂ ∪ B₃ = I₁ forces t₂ + t₃ ≥ 4k−1, so
(t₂, t₃) ∈ {(2k, 2k), (2k, 2k−1), (2k−1, 2k)}. In every case one of the
two is E' (if t = 2k and the size is 2k, it is contained = E'; if the
size is 2k+1 and t = 2k, it has one foot — and by E, d = 2 with I₁-part
E', the full even class: the d = 2 one-foot 2k+1-forms at s₁ = 2k+2 have
exactly I₁-part E', — the span 4k exceeds L = 4k−2 by
one point, so exactly one endpoint pokes out). Either way the partner AP
must contain O' = I₁ ∖ E' (2k−1 points at spacing 2, σ ≤ 2k+1), and T3
forces its difference ≡ ±2 — a ±-collision with the E'-side difference
d = 2. ✘ **M4 dead** (all three sub-cases at once).

**M5 ({2k+1, 2k+1}, ov = 2, s₁ = 2k+1).** Nothing is contained
(tᵢ ≤ 2k), and t₂ + t₃ − |B₂∩B₃| = 4k forces t₂ = t₃ = 2k with
B₂ ⊔ B₃ = I₁: both one-foot. By E, d₂ ∈ {2, 3k}; by symmetry over the
two sub-cases:

- d₂ = 2: the I₁-part of a d = 2 one-foot 2k+1-form at s₁ = 2k+1 is a
  *full parity class* of I₁ (the span 4k exceeds L = 4k−1 by exactly one
  point — one endpoint pokes out). So B₂ is a parity class, B₃ is the
  other, and T3 applied to A₃ (which contains its full parity class of
  s₃−1 = 2k points) forces d₃ ≡ ±2: ±-collision. ✘
- d₂ = 3k: the two runs of the Case-2 structure are offset by exactly
  d = 3k, so the blocks occupy [2k+1, 3k+1] and its translate
  [5k+1, 6k+1] (clipped to I₁), and the complement B₃ = I₁ ∖ B₂ contains
  the middle interval [3k+2, 5k] — 2k−1 to 2k consecutive integers. With
  σ₃ = 2k+1 this is a run of σ₃−1 or σ₃−2: T2 (and its slack form) forces
  d₃ ≡ ±1 — excluded. ✘ **M5 dead.**

This completes Lemma A for k ≥ 7, with the flagged dependency: M4/M5
consume Proposition E, which is verified computationally over the whole
stated range and whose mechanism is the three-distance foot count. ∎

### 3.4 The near-miss margin (the boundary of truth)

The search records, for every prime, the minimal slack by which a
covering *fails* ("best margin" = min over all configurations of
(minimal containing arc of the residual) − s₃; negative means a
covering exists). The margins grow linearly and confirm that the lemmas
sit far from feasible outside the boundary primes:

  p:      19   31   37   43   61   79   109   151   199   241   (Lemma A)
  margin:  3    8   10   12   18   24    34    48    64    78
  p:      11   17   23   29   41   59   89   131   191   239   (Lemma B)
  margin:  2    5    7    9   13   19   29    43    63    79

At the failure primes the margin is negative (p = 7: −2) or zero
(p = 13): the boundary of truth is exactly {7, 13} for Lemma A and empty
for Lemma B.

### 3.5 The boundary primes p = 7 and p = 13 — Lemma A is FALSE there

The exhaustive search finds genuine ±-distinct coverings — a correction
to the framework's expectation that the boundary solutions are
non-±-distinct:

- **p = 13 (k = 2), sizes (5, 5, 5), ov = 2, differences (1, 3, 4):**
  A₁ = [0, 4], A₂ = {9, 12, 2, 5, 8} (d = 3), A₃ = {7, 11, 2, 6, 10}
  (d = 4). Union = Z₁₃; pairwise differences ±-distinct. All sizes are
  *reduction* sizes ({2k+1}³ with k = 2), so the t-space kill fails at
  p = 13 outright.
- **p = 7 (k = 1):** 50 raw configurations (≤ 2 affine classes per size
  multiset), e.g. sizes (3, 3, 3), differences (1, 2, 3):
  A₁ = [0, 2], A₂ = {1, 3, 6}, A₃ = {1, 4, 5}; and sizes (2, 2, 4),
  differences (1, 2, 3): A₁ = [0, 1], A₂ = {2, 4}, A₃ = {3, 6, 2, 5}.
  The (3, 3, 3) family again uses only reduction sizes ({2k+1}³,
  k = 1).

Both primes are absorbed at the cell level: the (1K, 3U) cells at
p² = 49 and 169 are verified empty — by the committed T-12c structured
search and independently here by full brute force over unit ±-classes
(§4). In the language of the program: p = 7 and p = 13 join the rigid
region {7, 11, 13, 17, 19, 37} as primes where the t-space lemma cannot
close the cell and the direct check must.

---

## 4. Computational verification

Harness: `scripts/ap_lemmas_verify.py` (runtime 13.3 s, single-threaded,
exact integer arithmetic; output `scripts/out_ap_lemmas.json`; log
`scripts/ap_lemmas_run.log`). All primes ≤ 250: 25 of class ε = 1 and 26
of class ε = 5. The lemma search is **exact, not heuristic**:

- *Normalization.* One AP is scaled (d₁ → 1) and translated
  (A₁ = [0, s₁−1]); all size-multiset assignments are enumerated.
- *Exact prune.* |A₁ ∩ A₂| ≤ ov — forced by inclusion–exclusion
  (ov = p₁₂+p₁₃+p₂₃ − triple ≥ p₁₂). Viable (d₂, a₂) pairs are found by
  a difference-array sweep over the forbidden circular arcs.
- *Exact fit test.* R := I₁ ∖ A₂ must sit inside a d₃-AP of size s₃:
  a circular-arc fit check (largest gap of d₃⁻¹R). Covering is then
  automatic (A₁ ∪ A₂ ∪ R = Z_p), and every recorded solution is
  re-verified by direct set construction (0 re-verification failures).
- *Independent cross-validation.* A brute force with the same
  normalization but a₂, a₃ enumerated over all p² placements and direct
  union checks (no prune, no fit test) reproduces the search counts
  exactly at p ∈ {7, 11, 13, 17, 19}.

### 4.1 Results table

| check | range | result |
|---|---|---|
| Lemma B (±-distinct), all size multisets, all differences, all placements | every prime p ≡ 5 mod 6, 11 ≤ p ≤ 250 (25 primes) | **0 solutions** |
| Lemma A (±-distinct) | every prime p ≡ 1 mod 6, 19 ≤ p ≤ 250 (23 primes) | **0 solutions** |
| Lemma A at p = 7 / p = 13 | — | 50 / 2 raw configurations (§3.5) |
| Trichotomy census (contained critical APs) | all 51 primes | d ∈ {2, (p−1)/2} for k ≥ 2 (ε=5) / k ≥ 3 (ε=1); exceptions exactly p = 11: {3,4}, p = 13: {5}; p = 7, 17 clean |
| Two-block structure at d = (p−1)/2 | all 51 primes | 0 mismatches against the exact block form |
| Extended census (≤ 2 feet) | all 51 primes | clean {2, (p−1)/2} for ε=1 at k ≥ 7 (p ≥ 43) and ε=5 at k ≥ 3 (p ≥ 23); dirty exactly at ε=1: {13, 19, 31, 37}, ε=5: {11, 17} — Proposition E's boundary |
| Relaxed search (no ±-distinctness) | all 51 primes | ε=5: exactly 30 raw finds per prime = 4 families × size-multisets (differences (1,1,1), (1,2,2), (1,2,(p−1)/2), (1,q,q) — **every family has a ±-collision**); ε=1: 208 per prime for p ≥ 19 (same four families); richer at 7, 13 |
| Capacity formula L(c) | p ≤ 37, all units mod p × all fibers | sizes match the on/off-ball split with **0 violations**; every fiber bad set is an AP with difference u⁻¹ (**0 violations**) |
| (1K,3U) cells at p², independent brute force over unit ±-classes mod p² | p ≤ 31 (9 primes: 5–31) | **0 coverings** (triples and pairs), confirming the committed T-12c (p ≤ 47) |
| Rigid zone (min ±-distinct dilates of B_k covering Z_p), n = 5 side | 7 ≤ p ≤ 61 | 3 dilates at {7, 13}; 4 at {17, 19, 37}; ≥ 5 elsewhere — the committed zone exactly |

### 4.2 Anchors (hand-verified facts the harness must reproduce)

1. p = 17 rigid family found: A₁ = [0,5], A₂ = {5,7,9,11,13,15},
   A₃ = {6,8,10,12,14,16}, differences (1,2,2), sizes (6,6,6), ov 1. ✔
2. p = 23 wrap-around counterexample: the step-8 AP of size 8 is
   contained in the arc [5, 22] — reproduced as the set
   {5,6,7,13,14,15,21,22}. ✔
3. p = 11 parity family present in the relaxed search; absent in the
   ±-distinct search. ✔
4. p = 13 trichotomy exception d = 5 present (the AP {4,7,9,12} ⊆
   [4,12]). ✔
5. p = 17 census: contained differences exactly {2, 8}. ✔
6. Capacity: p = 17 on/off-ball sizes 5/6; p = 19 sizes 7/6. ✔
7. Cells empty at every tested p (≤ 31). ✔
8. Brute force = pruned search at all five cross-validation primes. ✔
9. Rigid zone: j = 3 exactly at {7, 13}; j = 4 exactly at
   {17, 19, 37}. ✔

All nine anchors pass. The two earlier drafts of this harness had two
bugs, both caught by the anchors before any number was recorded: a
per-column indexing error in the vectorized fit test (caught by the
brute-force mismatch at p = 7), and a wrong "consecutive" test in the
capacity structure check that treated wrap-around arcs as non-arcs
(caught by the p = 31 anchor). This is the cross-validation discipline
doing its job.

---

## 5. Status of the reduction and the named residual

**What the lemmas now buy.** Combining this note with the committed
record, the (1K, 3U) cell at p² — the last open cell of the p² program —
is empty for every prime p ≥ 5, via:

- p ≡ 5 mod 6, p ≥ 23: **Lemma B proved** (this note, §2); p = 11, 17:
  Lemma B verified by exhaustive search; p = 5: degenerate, cell checked.
- p ≡ 1 mod 6, p ≥ 43: **Lemma A proved** (this note, §3 — M1–M3
  unconditional, M4–M5 modulo the verified Proposition E); p = 19, 31,
  37: verified by exhaustive search; p = 7, 13: **lemma false** — closed
  by direct cell verification (committed T-12c to p ≤ 47; independent
  brute force here to p ≤ 31).

**The named residual (stated plainly).** Both lemmas carry the
±-distinctness hypothesis, which the reduction must supply: the t-space
differences are uᵢ⁻¹ mod p, ±-distinct iff the three units are pairwise
±-distinct mod p. The cell definition (distinct speeds mod p²) does
*not* exclude unit pairs with uᵢ ≡ ±uⱼ (mod p), and for such pairs the
t-space APs have ±-colliding differences — the relaxed families exist at
every prime, so the lemma-level kill does not fire. What currently
covers this case: the cell-level verifications (committed T-12c, all
unit triples including mod-p collisions, p ≤ 47; the independent brute
force here, p ≤ 31). What would cover it unconditionally at large p:
either (a) the offset-variation argument — the t-space covering must
hold in *every* off-kernel fiber simultaneously, with per-fiber offsets
that shift with r′, which the static relaxed families cannot tolerate
(natural next proof target; the parity family is destroyed by any
offset change, and the (1, 2, q) and (1, q, q) families deserve the same
check); or (b) at p ≡ 1, the mod-p counting already collapses a
±-colliding triple to a *three-ball* covering problem
(a⁻¹B_k ∪ u₁⁻¹B_k ∪ u₂⁻¹B_k = Z_p), which the rigid-zone data excludes
outside {7, 13} — but that exclusion is itself verified only to
p ≤ 700 (committed T-12b). So: the unconditional-in-p closure of the
p² program now rests on one named obligation, **the ±-colliding unit
case beyond the verified ranges** — a clean, finite-looking target.

**Consequences for the program.**

1. The reviewer's question is answered: the human's proof architecture
   survives; Steps 1–2 are subsumed by the trichotomy; Step 3 is intact.
   The corrected lemma statements and their exact boundaries are now
   pinned: Lemma B true for all k ≥ 1 (proved k ≥ 3); Lemma A true for
   k ≥ 3 at primes (proved k ≥ 7), false at 7 and 13.
2. The paper's §6 open problem on the (1K,3U) cell can be updated: the
   cell is empty for p ≡ 5 unconditionally (Lemma B, proved) and for
   p ≡ 1, p ≥ 43 (Lemma A, proved modulo Proposition E), with the
   residual being the ±-colliding unit case rather than the whole cell.
3. Next targets, in order of value: (i) prove Proposition E's foot count
   (three-distance bookkeeping — removes the last computational input
   from Lemma A's proof); (ii) the offset-variation argument for
   ±-colliding units (removes the residual); (iii) then the pq side of
   the no-unit-assistance program, where the same trichotomy machinery
   should transfer (the pq cells have the same t-space structure per
   CRT component).
