# Audit Adjudication Addendum — The m=10 Resolution (v3 Revision Record)

**Date:** 2026-10-08 (second certification session).
**Applies to:** `flagship_audit_adjudication.md` (points O-4/O-5 and the
"search evidence thin" observation) and `companion_audit_adjudication.md`
(the zonotope-strand ladder points).
**Trigger:** the reviewer's directive — "m=10 remains uncertified, and the
audit showed random local search finds deeper holes than expected — worth
one more certification attempt."

**Outcome up front: ρ(1,2,3,4,10) = 7/22 EXACTLY, CERTIFIED.** The last open
row of the n=5 deepest-hole ladder is closed. New artifacts:
`scripts/certify_m10_v2.py` → `scripts/out_certify_m10_v2.json` (+ the
checkpoint state file), shipped in the v3 paper bundles.

---

## 1. The two audit points, and their disposition

### The "search evidence thin" point (opus, flaghship.txt:38)

- **opus's claim:** 40 random local searches sufficed to beat the published
  harmonic lower bound (D = 4183/12288 > 607/1792), so the 288-start search
  behind m=10's "consistent" label looked thin.
- **Adjudication at v2:** ACCEPT (the point stood; v2 hedged the label).
- **Disposition now: RESOLVED BY EXECUTION.** The second attempt ran the
  search at ~50× opus's demonstrated budget in starts and ~700× in raw
  evaluations: 4,194,304 uniform random points, the full {k/16}^4 (65,536)
  and {k/32}^4 (1,048,576) dyadic grids, a 500,000-point torsion
  neighborhood, 2,048 batched random uphill climbs (1,802,240 evaluations),
  and 64 Nelder–Mead polishes — all screened on a fixed s-grid over
  S = [0,r] ∪ [1−r,1), which is *sound* (the grid minimum of G(c,·) is an
  upper bound for D(c), so no point above the threshold can be missed), with
  exact rational re-evaluation of the 114 deepest candidates plus the 50
  persisted first-attempt survivor centers and the 24 deepest dyadic-16
  points (188 exact checks in total).
- **Result: no excess.** The deepest verified point in the entire search is
  exactly 7/22, at the torsion center (1/2, 0, 1/2, 1/2) itself; the deepest
  first-attempt survivor center is 1135/3584; nothing exceeds 7/22 anywhere.
  The audit's worry was legitimate and is now answered at scale.

### O-5: m=10 "quantized minimax loss, not an excess" (opus, flaghship.txt:39)

- **opus's claim:** boxes with U > r are by definition the unresolved ones;
  residual volume < 6×10⁻⁷ does not exclude a deeper hole; the
  characterization was unproved.
- **Adjudication at v2:** ACCEPT (row kept open, wording hedged).
- **Disposition now: MOOTED BY CERTIFICATION.** The new exact **Lipschitz box
  certificate** proves what the first attempt's minimax bound could not:
  D is 1-Lipschitz in the sup norm (for fixed σ, each term ‖c_j − v_jσ‖
  moves by at most |c_j − c′_j| while the σ-term is fixed; minimize over σ),
  so for a dyadic box C with center c₀ and half-width h,
  **max_C D ≤ D(c₀) + h** with D(c₀) computed exactly in rational arithmetic.
  At the first attempt's residual width (1/512 per side, h = 1/1024):
  1135/3584 + 1/1024 = 2277/7168 ≈ 0.31766 < 7/22 ≈ 0.31818 — the observed
  gap closes exactly. The re-run branch and bound (cheap grid →
  grid-Lipschitz filter with sampled exact confirmation → exact minimax U →
  split) **emptied the tree**: 190,295 box evaluations (9,608 cheap +
  42,502 grid-Lipschitz + 43,038 exact-minimax prunes; 95,147 splits;
  max depth 34), zero survivors, 1,286 s on one core. Together with the
  torsion center attaining 7/22 (hand check: at σ = 7/22 the five terms are
  7/22, 3/22, 1/22, 5/22, 7/22, with the binding pair at opposite slopes
  +1 and −10), this certifies equality.

## 2. Soundness profile of the certificate

- **Negative control:** the harmonic instance v = (1,2,3,4,5) at r = 1/3
  again REFUSES to certify (84,986 box evaluations, 38,117 survivors
  persisting around its true deep holes, the deepest known 4183/12288 —
  opus's own witness, now guarding the control).
- **Determinism check:** the n=4 certificates re-certify under the new
  cascade with identical box counts (875 and 443 pops) — the pop sequence
  is fully determined by the heap priority, so the certificate is
  checkpoint-slice-independent.
- **Exact re-confirmation:** 846/846 sampled float-stage prunes of the
  target run (1-in-64 sampling plus the first 32 of each stage) were
  re-confirmed in exact rational arithmetic; likewise 295/295 for the
  refused control and 36/36, 34/34 for the n=4 re-certifications. Zero
  failures.
- **Adversarial unit tests:** the 1/512-boxes around the three known deeper
  holes — 607/1792 and 4183/12288 (harmonic, r = 1/3) and 6/19
  (v = (1,2,3,4,15), r = 5/16) — all REFUSE the Lipschitz test exactly
  (center depths 347/1024, 1045/3072, 6133/19456; each plus 1/1024 exceeds
  its r). The certificate cannot kill a true deeper hole.
- **Witness assertions:** all five canonical witnesses (7/22; 607/1792;
  4183/12288; 97/288; 6/19) re-asserted exactly at run start; D
  cross-validated against the independent lattice-DFS implementation on 80
  random rational points across both instances.

## 3. What changes in the papers (v3)

- **Flagship (v3):** abstract and §5 updated — the n=5 modular-law row now
  reads: m ∈ {5,15} fail exactly, m ∈ {10,20,25,30} certified, all 20
  non-multiples fail exactly; §5 records the second attempt in full; §5's
  B&B subsection gains "The second attempt's cascade" paragraph; §6 drops
  the "deeper point not excluded" caveat; Appendix B's certificate table
  gains the v3 row (95,148 leaves, 1,286 s, 0 survivors) with the cascade
  footnote and the fifth reproduction command; Appendix C gains the
  sup-norm Lipschitz lemma, the box-certificate corollary, the
  2277/7168 worked example, the adversarial refusals, and the 7/22 hand
  check. A "What changed in this revision (v3)" paragraph is added to §1.
- **Companion (v3):** no theorem, lemma, or hypothesis changes; §6 and
  §sec:scale records updated (m=10 closed at 7/22; the harmonic lower-bound
  record harmonized with the audited 4183/12288 witness); version marker
  recorded.
- **Zonotope report (v2 file):** the ladder, the §3 table row, the verdict,
  and the validation ledger all updated; no "consistent" label remains
  anywhere in the record.

## 4. Bottom line

The audit improved the record twice: first by finding the deeper harmonic
witness (4183/12288, adopted in v2), and now by pressing the two weak points
of the m=10 row until the only fully satisfactory answer — an exact
certificate — was produced. The n=5 modular-law story is complete:
**the odd multiples 5 and 15 fail exactly; the even and higher multiples
10, 20, 25, 30 are certified; the law is a small-n phenomenon.**
