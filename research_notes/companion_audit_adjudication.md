# Companion Monograph — Audit Adjudication (v2 Revision Record)

**Paper:** *The Lonely Runner Problem on Pair-Sum Lattices — a lossless finitization
and a conditional scaling-closure theorem* (v1: 65 pp, Oct 8, including the T-26
final revision).

**Audits received** (`audits/companion.txt`): **opus** (mathematical-gap audit with
a proposed repair), **gemini** (readability audit).

**Adjudication standard.** Each audited claim was checked against the paper's actual
statements (secsc.tex Theorem SC + ledger; secscale.tex Prop. depth, spacing lemma,
volume laws; sec2.tex equality theorem; sec7.tex soundness record) and the persisted
program records. The central mathematical objection (the missing multi-fiber bridge)
was accepted and repaired exactly along the auditor's proposed line.

---

## 1. The central finding — ACCEPTED and REPAIRED

**opus:** "Theorem SC is not currently derivable from the stated ledger. … H2m
supplies a *pairwise* intersection bound. Pairwise sparsity does not control
higher-order correlations: every pair can look random while a large common
intersection survives. … You cannot prove the paper's higher-order cascade bound
from H2m as currently stated."

**Verification against the paper:** confirmed. The ledger's H0 row is
"C_i ≈ γ̄^i |Sol|^{i+1}/p^{(m−1)i}; no tangent-drift trajectory stays inside the
solution family for many fibers" — an unquantified recurrence for the
*independent-extension model* with no constants, no quantifiers, and no permitted
error. Proposition `prop:depth` is explicitly stated "in the independent-extension
composition." The actual cascade is history-dependent (repeated dilations, resonance
fibers, changing solution sets), and no hypothesis as stated bounds *conditional*
intersections |S_{i−1} ∩ A_i| after prior filters. H2m is exactly the i = 1
(unconditional) case.

**The repair (opus's, adopted):** a new quantified hypothesis replaces H0 —

> **(Hq)** For every admissible cascade history S_{i−1} and every next fiber i,
> |S_{i−1} ∩ A_i| ≤ Γ·|S_{i−1}|·|Sol_i| / p^{m−1}, where A_i = λ_i^{-1}(Sol_i − τ_i),
> with Γ uniform over families, size tuples, dilations, translations, and prior
> surviving fibers.

Then the derivation becomes a genuine induction: with d = m−1, B = c₁·ov^d and
r = ΓB/p^d, |S_i| ≤ r|S_{i−1}| gives |S_q| ≤ B·r^q, so the cascade is empty whenever
q > log B / ((m−1)log(p/ov) − log(c₁Γ)), provided the denominator is uniformly
positive. This is stated in v2 as a theorem (the Hq-induction), and the depth
denominator is written exactly — (m−1)log(p/ov) − log(c₁Γ) — with the paper's
T/(T−6) plateau identified as the k → ∞ form of p/ov with the O(1/k) correction
absorbed into the O(1) tail.

**v2 two-layer restructure (opus's "most important repair"):**
1. **Theorem SC (layer 1, rigorous):** conditional on the quantified hypotheses
   {NF, Hq, H1g(η), H2r, EQ} — with H1g sharpened to carry the margin (see §2) and
   H2m retained as Hq's one-step base case (measured; proved for the chain class).
   The derivation from the ledger is now complete (the induction above).
2. **Calibration (layer 2, measurements):** H0's former content (depths ≤ 4 at
   committed cells, ≤ 11 across the calibration rungs, exit times ≤ 2, the 764-pair
   record) is reported as *measurements supporting Hq*, with no language suggesting
   they nearly constitute proof.

**Consequences propagated:** the abstract's "under two named open hypotheses" becomes
"under three named open hypotheses — the multi-fiber survivor inequality Hq, the
distinct-family volume law H1g (with its uniform margin), and the resonance budget
H2r"; the ledger table gains the Hq row (status: OPEN); the epistemic status is
updated accordingly.

## 2. The remaining mathematical points (opus)

| # | Point | Verdict | v2 action |
|---|---|---|---|
| 1 | H0 "not a mathematical hypothesis" | ACCEPT | Replaced by Hq (quantified); H0's measurements moved to the calibration layer |
| 2 | H1g "not equivalent to the threshold claim" — need c₁ < (T/(T−6))^{m−1} with quantified margin η uniform in (T,k) | ACCEPT | H1g restated: \|Sol\| ≤ c₁·ov^{m−1} with c₁ ≤ (T/(T−6))^{m−1}·e^{−η}, some η > 0 uniform in (T, k, family, size tuple); Proposition `threshold` re-described as locating the vanishing value of the depth denominator, not as an equivalence |
| 3 | H2m quantifiers incomplete (conditional intersections) | ACCEPT | Superseded by Hq; H2m's row reframed as the one-step base case + chain-class theorem |
| 4 | H2r "not proved from the definition; Λ₀ data-chosen; small-rational λ undefined" | ACCEPT | Row relabeled OPEN (measured support; formula shape derived from the proved reflection lemma + a Gauss-circle-type heuristic); the resonance set defined explicitly (fibers with dilation ratio λ of small rational height \|λ\| ≤ Λ₀ in the height function fixed in §scale); "derived given Λ₀" wording removed |
| 5 | Prop 9.5 → exhaustion underproved (history-dependent cascade) | ACCEPT | The exhaustion is now Theorem-level via the Hq-induction; Prop `depth` retained, explicitly relabeled as the independent-extension *model* + calibration, its (a)–(d) conclusions restated as properties of the model/the formula |
| 6 | "C(T,k) ≪ \|U\| is not enough by itself" | ACCEPT | Theorem (ii) now cites the Hq-contraction explicitly (every surviving history loses ≥ 1 − r of its mass per fiber, r = Γc₁(ov/p)^{m−1} ≤ e^{−η}); the "margin growing linearly in k" phrase tagged as the fiber-budget surplus, conditional on the ledger |
| 7 | "Margin growing linearly in k is not established" | ACCEPT | Same tagging: all margin language in theorem statements is explicitly conditional; the abstract already says "under … hypotheses"; strengthened where loose |
| 8 | p = Tk + ε underspecified (ε never given a range; k non-unique) | ACCEPT | Theorem SC now: "write p = Tk + ε with k ≥ 2 and 0 ≤ ε < T (k = ⌊p/T⌋ canonical)"; the two-value size window and ov-formula footnoted with the same convention |
| 9 | Lemma 9.3 → cut heuristic jump (b_min = 1 only gives N ≤ L+1, not N ≈ (L+1)/λ̄) | ACCEPT | The jump is now marked explicitly as the heuristic step and the missing arithmetic count is stated as part of H2m's remaining obligation (the pair-line count); the lemma's own statement unchanged (it was correctly stated as a spacing bound only) |
| 10 | Lemma 3.5(a) (attainment) perturbation "too compressed" — different slopes, cusp-crossing non-binders | ACCEPT | Proof expanded: the non-binder estimate is the *global* Lipschitz property f_p(t±s) ≥ f_p(t) − v_p\|s\| (valid across cusps — this is why cusp crossings are immaterial); the binder linearity window is guaranteed per-binder by the min with the factor ½(½−m)/v_p, which is what handles differing slopes |
| 11 | M = ½ case should be its own lemma; all-equal patch late | ACCEPT | Split out as its own lemma ("the half case") with the gcd argument as its proof and the degenerate all-ones case handled inside it |
| 12 | Distinct vs non-distinct speeds | ACCEPT (clarified, no math change) | §2 already proves the theorem for possibly-non-distinct speeds and says so; a remark now states which version LRC needs (distinct, WLOG) and that nothing downstream uses distinctness |
| 13 | "New as stated" over-claims computational verification's reach | ACCEPT | Replaced by "not located in two literature rounds (56 queries, arXiv 2024–2026 enumerated); the correct prior, given how actively the spectrum line computes exact optima [Cordella], is that pieces may exist in another language" |

## 3. Internal inconsistencies (opus) and readability (gemini)

| # | Point | Verdict | v2 action |
|---|---|---|---|
| 14 | Ledger statuses look like evidence ("measured, margin growing in k") | ACCEPT | Statuses rewritten in caps-locked form: OPEN (measured support: …); the table caption now says "conditional inputs, not evidence toward validity" |
| 15 | "Closed by theorem" for the transfer-matrix route too broad (§9.7 correctly scopes to bijective-drift automata) | ACCEPT | All broad instances (abstract, intro, §SC) narrowed to "automata built on the cascade's bijective drift"; the epistemic status's scope statement retained verbatim |
| 16 | "p² side closed" vs "frontier cell closed" vs "rigid-zone lifts open" — scope confusion | ACCEPT | A scope paragraph in §SC fixes the three levels: *cells* (T = 7, 8 critical cells closed within tested primes p ≤ 61), *rungs* (conditional on the ledger), *rigid-zone finiteness* (open, empirical) |
| 17 | Table 8 T=9 stride samples read as structural | ACCEPT | Every projected-pool citation now carries "sampled/projected (two primes, ε ∈ {2,4}, unmodeled modulation)" |
| 18 | "All computations are exact integer computations" | ACCEPT | Notation §: "exact integer/rational arithmetic wherever a theorem or table entry asserts an exact value; the Wilson intervals, the fitted Λ₀, the projected pools, and displayed asymptotics are statistical or asymptotic, labeled as such" |
| 19 | §9.8 chain volume law \|Sol\| ≈ ov^{m−1} vs Theorem 9.9's s_max-factors | ACCEPT | Passage now states: the ov^{m−1} shape with O(1) constant is H1g's empirical content; the proved chain-class bound carries (m!)²·4m·s_max factors, loose by orders of magnitude; nothing in between is claimed |
| 20 | "T=10 saturation predicted ≥ 0.9" not a prediction with error control | ACCEPT | "extrapolated from the rung-monotone pattern (two primes, unmodeled ε-modulation); no error control; not measured" |
| 21 | Abstract megasentence (108 words, gemini) | ACCEPT | Split into declarative steps (gemini's proposed text adopted nearly verbatim) |
| 22 | "Honest/honestly" motif (>15 occurrences, gemini) | ACCEPT | Reduced to the load-bearing four (abstract's "stated honestly as a conditional theorem", the epistemic prior, the honest boundary, one §cal remark); others replaced by structural signposting |
| 23 | §8.2 placement-count dump (gemini) | ACCEPT | The eleven-family census paragraph keeps the family list; the ε-class placement counts moved to a compact two-row display table |
| 24 | h = ⌊(N−1)/T⌋ missing from §1.4 notation (gemini) | ACCEPT | Added to the notation paragraph |
| 25 | Bibliography (same phantom entries as the flagship) | ACCEPT | Lam–Leung split into the two real papers (Monthly 103 (1996) 562–564; J. Algebra 224 (2000) 91–109 — Crossref-verified); the phantom "de Bruijn–Schoenberg, On the two vanishing identities of Euler" replaced by the Rédei–de Bruijn–Schoenberg attribution carried by Lam–Leung's account; Mirsky–Newman kept with the unpublished caveat; Cusick 1973 added and cited at the two attribution sites |

**Rejected / no change:** opus's framing "the equality theorem is the strongest part,
but the computational verification cannot support the novelty claim" — the novelty
*hedge* was already present (the "spectrum line" caveat); only the "new as stated"
phrase needed replacing (done, #13). No audited numerical claim of the companion was
found false. The auditor's own arithmetic (256 = 4⁴ etc.) was verified exactly.

## 4. What does NOT change in v2

- Theorem SC's truth content: v1's statement was conditional on an under-quantified
  ledger; v2's is conditional on the quantified ledger. No lemma, proposition, or
  machine record of v1 is withdrawn; the Hq repair *names* the missing input rather
  than claiming it.
- All persisted records, the 225-item provenance ledger, the T-16..T-26 session
  record, and the paper's empirical findings (rigid-zone classifications, decay
  curve, census laws) — unchanged.
- The T-26 third-bridge disposition and the two-paper cross-referencing — unchanged
  except the flagship's new title in the cross-reference.

## 5. v2 deliverable

- `lonely_runner_pair_sum_lattices_conditional_theorem_paper_v2.pdf` — the revised
  monograph.
- `paper_sources_v2/` — complete LaTeX sources of v2 (v1 files untouched).
- Verification artifact: `scripts/audit_verify_v2.json` (shared with the flagship
  adjudication; includes the threshold-constant exactification used in §H1g).
