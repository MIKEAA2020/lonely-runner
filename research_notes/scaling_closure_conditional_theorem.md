# The Scaling-Closure Lemma, Stated as a Conditional Theorem

**Program:** LRC pair-sum-lattice reformulation. **Directive (reviewer,
2026-10-06, post-T-19):** consolidate the scattered conditional claims into
ONE named conditional theorem — "Assume H1-general, H2-mid-region, and
H2-resonance. Then for every T ≥ 8 and every p ≥ p₀(T), every (1K, T−3U)
cell at p² is empty, with cascade depth bounded by C(T) as defined" —
listing each hypothesis, its status (proved, measured, conjectured), and
the trace of the derivation. This is that page. It supersedes the
conditional statements scattered through the T-18 and T-19 notes; those
notes remain the evidence records. **No "one more session" commitment is
made or implied anywhere below.**

---

## The theorem

> **Theorem SC (scaling closure; conditional).** *Assume H1-general, H2m,
> H2r, H0, and NF as stated in the ledger below. Let T ≥ 8, p = Tk + ε
> prime with k ≥ 2 (p ≥ p₀(T), the least such prime — p₀(8) = 17), and let
> the cell (1K, T−3U) at N = p² be given, with m = T−3 arcs, size window
> {2k, 2k+1} or {2k+1, 2k+2}, and overlap budget ov = Σs − p ≈ (T−6)k + O(1).
> Then:*
>
> **(i) every census-admissible family dies in the fiber cascade within**
>
> **C(T,k) = log(c₁·ov^{m−1}) / [ (m−1)·log(T/(T−6)) − log c₁ ] + 2Λ₀(T−2)/T + 1 + O(1)**
>
> *fibers — numerically the T-flat plateau ≈ 9 / 9.2 / 10.1 at T = 8 / 9 /
> 10 (k = 2 / 3 / 3) — with the never-exceeded worst case*
> **(T−4) + N_res + O(1) ≤ T + 3, N_res = 1 + 2Λ₀(T−2)/T**;
>
> **(ii) since C(T,k) ≪ |U| = (T−2)k for every T ≥ 8, k ≥ 2 (margin
> growing linearly in k), the cascade exhausts the fiber budget: every
> admissible family has no covering placement valid on all of U, hence*
>
> **(iii) the cell (1K, T−3U) at p² is empty; and by the lossless
> reformulation equivalence, LRC holds at that rung.**
>
> *Measured range of the calibration: T = 8, 9, 10 at p ≤ 73 (the
> reviewer's three rungs; the committed GT cells p ≤ 61 at T = 7/8). The
> statement for p beyond the calibration is conditional on H1g/H2m holding
> there, which is exactly what the hypotheses purchase.*

**What is conditional and what is not.** The reduction chain
(i)→(ii)→(iii) is proved from the committed machinery (Lemma 1's
dilate-intersection reduction, the drift law, the lossless reformulation
equivalence — all machine-verified; see the trace). The arithmetical
content — that the cut beats the pool within C(T,k) fibers — is proved
 GIVEN the five inputs. None of the five inputs is silently established;
that is the point of this page.

## The hypothesis ledger

| # | hypothesis (precise form) | status | evidence / route |
|---|---|---|---|
| **NF** | Census normal-form completeness: every covering family at an uncovered fiber is census-admissible (the normal form with normalizer difference 1, canonical differences, two-tier sizes loses no coverings) | **proved in structure, machine-validated** | the footprint/arc derivations are T-free and verified 957/957 against direct Z_{p²} enumeration; the T=7 census-restricted search equals brute force at the boundary prime p=5 (the 8-covering control); the T=8 zoo census itself is the measured object |
| **H0** | Sequential composition: the per-fiber cuts compose, Cᵢ ≈ γ̄ⁱ·\|Sol\|^{i+1}/p^{(m−1)i} (no tangent-drift trajectory stays inside the solution family for many fibers) | **measured benign in range** | depths ≤ 4 at the committed cells (764-pair record, T-19); ≤ 11 across the reviewer's three rungs; exit-time ≤ 2 exhaustive over admissible families at T=7/T=8c; route = tangency codimension ≥ 1 |
| **H1g** | Distinct-family volume law: \|Sol(f, sz)\| ≤ c₁·ov^{m−1} for every ±-distinct family f and admissible size tuple, with c₁ = O(1) uniform in (T, k, f) | **measured, with margin growing in k; chain class proved** | T=8 this session: max ratio \|Sol\|/ov⁴ = 0.214 / 0.083 / 0.013 at p = 17 / 19 / 29 (chain-class benchmark: 2.17 / 3.61 / 2.03); chain class proved (B1′/B2, T-19); remaining distance = the strand census (see the companion note) |
| **H2m** | Mid-region dilate-sparsity: survivors ≤ γ₀·\|Sol₀\|·\|Sol_j\|/p^{m−1} for every fiber pair with λ̄ > Λ₀, γ₀ = O(1) uniform | **measured; proved for the chain class (loose constants)** | 764 pairs at 5 committed cells: mid-zone γ median 0 (instant kills), max 40.6; chain class proved via Theorem B3 with constants (m!)²·4m·s_max — loose by orders of magnitude; route = the Lemma-E/three-distance pair-line count |
| **H2r** | Resonance budget: the weak-cut fibers (λ̄ ≤ Λ₀ or small-rational λ) number ≤ 2Λ₀(T−2)/T + 1, k-free; the λ = −1 fiber is never a cut | **derived given Λ₀; the λ=−1 clause PROVED (Lemma R)** | Lemma R proved (the x ↦ −x involution; survival fraction exactly 1 at every committed family/cell); the count is the Gauss-circle-type small-rational dilation count; **Λ₀ ≈ 4 is data-chosen** |
| **EQ** | Lossless reformulation: the cell emptiness at the rung implies LRC there | **committed theorem** | the pair-collapse equality M(V) = max over pairs of τ (machine-verified, 2,870/2,870 independent re-verification); the fibration threshold-compatibility |

## The free-parameter disclosure (standing practice, per the reviewer)

Every constant entering Theorem SC, with its provenance:

- **c₁ (H1g):** bounded by data only — ≤ 0.214 at T=8 over three primes
  (p = 17, 19, 29), and the bound tightens with k; **not proved**. The
  chain-class value (the only proved case) is 2.0–3.6.
- **Λ₀ ≈ 4 (H2r, and the resonance tail of C(T)):** data-chosen — it is
  the λ̄ zone where the measured per-fiber cuts at the committed cells
  first exceed ~10× the random model. N_res = 1 + 2Λ₀(T−2)/T is therefore
  a **derived formula containing one fitted constant** — partially
  derived, not derived outright.
- **γ₀ (H2m):** bounded by data (≤ 41 mid-zone over the committed cells);
  proved only for the chain class, with loose constants.
- **R (the T-18 junction-volume scale, fitted 4.2):** **retired.** For the
  chain class it is replaced by the derived overlap budget (T−6); for the
  thin class it is superseded by c₁ under H1g. No fitted R remains in
  Theorem SC.
- **The depth formula's own constants** are otherwise derived: the leading
  term from the overlap budget T−6 and the modulus T (both structural),
  the tail from Lemma R (proved) + the Λ₀-count above.

## The trace of the derivation

1. **T-17:** the general-T fiber machinery (three-tier arc law, drift law
   — proved, 3,528/3,528); the T=8 zoo census: 840 ±-distinct families at
   p=17, 840 = 7·6·5·4 (then believed exact).
2. **T-18:** Lemma 1 (the cascade step IS a dilate-intersection
   λ·Sol₀ + τ); Lemma 2 (pool ≤ (D−1)_{T−4}); Lemma 3 (spacing:
   N(λ,t) ≤ L/b_min + 1); Prop 4 (k-cancellation: cut = k·T·(T/R)^{T−5});
   Theorem 5 (conditional comparison; the interpolating depth law). R and
   N_res fitted.
3. **T-19:** the pool law corrected — the falling-factorial identity is a
   p=17 accident (saturation 1.00 / 0.571 / 0.357 / 0.133 / 0.015 at
   p = 17 / 19 / 23 / 29 / 31, all GT-verified; T=9 measured 90.3% /
   49.0% by stride sampling); Lemma R proved (the reflection fiber);
   B1′/B2 proved (chain-box parameterization, overlap-sum identity);
   Theorem B3 proved (chain-class dilate reduction — H2 for the interval
   family); N_res derived modulo Λ₀; the H1 exponent corrected to m−1 =
   T−4 for the chain class.
4. **T-20 (this session):** Lemma V proved (renormalization invariance —
   the volume is a function of the dilation class of the difference set;
   verified at p=17: 7 orbits, and p=19: 14 orbits in 9 volume classes);
   the full distinct-family volume census at T=8 (840 families at p=17,
   GT'd; 1,680 at p=19, GT'd; targeted at p=29, GT'd); the volume law
   measured with the margin profile above; the top-class mechanism
   identified (the binary-cascade orbits); H1g reduced to the **strand
   census** — a named, bounded combinatorial obligation (companion note).

## The sampled numbers, with intervals (standing practice)

- T=9 pool saturation, p=29: **0.903, 95% Wilson CI [0.889, 0.915]**
  (n = 2,005 stride-sampled tuples, full decisions); projected pool
  **1.39×10⁵, CI [1.37, 1.41]×10⁵** — *projected from the falling
  factorial times the sampled fraction, at one prime*.
- T=9 pool saturation, p=31: **0.490, CI [0.468, 0.511]** (n = 2,002);
  projected pool **1.18×10⁵, CI [1.12, 1.23]×10⁵** — same caveat.
- The rung-monotone reading (T=9 pools ≈ the distinct-tuple pool) rests
  on **two primes, ε ∈ {2, 4}**; the ε-modulation is large (0.90 vs 0.49)
  and unmodeled. T=10 saturation (predicted ≥ 0.9 at the 40% budget) is
  **not measured**.
- All T=8 saturation figures (1.00 / 0.571 / 0.357 / 0.133 / 0.015) are
  exhaustive censuses, not samples — no intervals needed.
- Cascade-depth figures: committed cells are exact runs (median 2, max 4
  at the zoo random-40); the 9 / 11 / 10 maxima are the reviewer's
  measured plateau across k ≤ 9 / 8 / 3.

## The honest closing

The program is either one proof away (H1g + the sharp H2m, both reduced
to concrete objects: the strand census and the pair-line count) or one
systematic obstruction away (a family class whose junction structure
escapes the ov^{m−1} volume law — for which the T=8 evidence shows the
opposite: the measured margin **grows** with k, and no distinct family
reaches even 1/10 of the chain-class volume at equal ov beyond k = 2).
The distinction will be settled by whether the strand census closes at
T=8 with uniform constants, not by further measurement; the attempt
record and its verdict are in the companion note. The pattern the
reviewer flagged — each closure revealing the next gap at the same scale —
is alive in this structure: the strand census at T=8 decomposes by
resonance class, and uniformity across classes is itself a lemma. That is
the cost of the conditional form stated above; it is also exactly what
the hypotheses of Theorem SC now name explicitly, which the scattered
prior conditionals did not.
