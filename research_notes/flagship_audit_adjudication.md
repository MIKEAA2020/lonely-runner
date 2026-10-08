# Flagship Paper — Audit Adjudication (v2 Revision Record)

**Paper:** *Every Natural Strengthening Is False: Type Mismatch in Lossless
Reformulations of the Lonely Runner Conjecture* (v1: 20 pp, Oct 8).

**Audits received** (`audits/flaghship.txt`): **opus** (arithmetic verification +
logic critique), **grok** (verification + defense), **gemini** (line-level audit with
drop-in LaTeX), **gemini2** (prose audit).

**Adjudication standard.** Every checkable claim was re-verified against the persisted
computational record (the cross-validated exact library `lrc_zono_lib.py`, the run
JSONs, and a new Sturm-sequence run: `scripts/audit_verify_v2.py` →
`scripts/out_audit_verify_v2.json`). Opposing points between auditors were resolved
on evidence, not rhetoric. Verdicts: **ACCEPT** (change made), **PARTIAL** (qualified
change), **REJECT** (no change, reason given).

---

## 1. Verification results (all auditors' arithmetic)

| Claim | Verifier | Result |
|---|---|---|
| All five witness depths (97/288 at σ=31/288, 29/80, 4/11, 607/1792), λ/δ values, Table 4 rows, 3,712,576 corpus count, 571/1792 = 7/22 + 9/19712 | opus, gemini | **CONFIRMED** (independently by both; gemini's hand calculation of 97/288 and h\*(1)=d!·Σvⱼ reproduced) |
| **opus's new witness: D = 4183/12288 ≈ 0.34041 at c = (465/4096, 605/1024, 1019/4096, 315/4096)** | this run (library D_exact, the implementation cross-validated 40 pts/instance) | **CONFIRMED EXACTLY — and it beats 607/1792**. Adopted in v2 as the harmonic lower-bound record with credit to the audit. |
| Spearman p-values "≈0.31 and 0.46" | this run (scipy on the persisted arrays) | **CONFIRMED**: ρ=+0.3162, p=0.3167 (N=12); ρ=−0.3117, p=0.4523 (N=8). |
| h\* real-rootedness at 45 dps is numerical, not exact | opus | **CONFIRMED as a wording gap; remedy executed**: exact Sturm-sequence verification run on all 16 rows (square-free part, sign-change count = degree) — **all 16 real-rooted, now exact**. v2 claims upgraded from "45 dps" to "exact Sturm verification". |
| Newton floor ≈ 1.5 at n=10; rung-2 n=3 nS = 60.5 | opus, gemini | **CONFIRMED**: floor at argmin k = 3/2 at n=10; rung-2 nS(3) = 60.5. |
| Table 2 last row: D = 5/14 = thr₆ exactly, no coordinates given | opus, grok | **CONFIRMED**: record has witness (1/16, 11/16, 1/8, 7/16, 7/8) with D = 5/14 = thr exactly. Boundary case, not a strict refutation. |
| m=15 is a multiple of 5 (sentence contradiction) | opus | **CONFIRMED**: the JSON record shows m=15 is a *multiple* whose law-equality clause fails (6/19 > 5/16 = δ); the sentence conflated the 20 non-multiple witnesses with the 2 multiple refutations. |
| "6/19 < 1/3 and 71/208 < 5/14" — the three further families don't refute the ρ-form | opus | **CONFIRMED**: they refute the *modular deepest-hole law* (ρ > δ), not the ρ-form (ρ > thr). Strict ρ-form refutations number exactly three: 97/288, 29/80, 4/11. |
| New finding of this adjudication | this run | The record also carries **boundary rows not in v1**: n=5 m=6,7,8 have ρ ≥ 1/3 = thr **exactly** (witness value 1/3, δ = 3/10). Added to the v2 table as zero-slack boundary instances. |
| Threshold constants 256/243/244 = (T/(T−6))^{T−4} at T=8,9,10 | opus's reverse-engineering | **CONFIRMED** (exact: 256, 243, 15625/64 = 244.14). The missing definition: **m = T−3** (arc count, defined in the companion's Theorem SC; absent from the flagship text). |
| U(C) displayed formula drops ‖σ‖ | opus, gemini | **CONFIRMED as textual**: the implementation (`certify_rung2_n5.py::_R_values`) *does* include ‖σ‖ as the first term in the max; only the paper's displayed equation omitted it. Formula fixed to match code and eq. (1). |
| m=10 control row: leaves "—", time in wrong column | opus | **CONFIRMED**: the persisted run shows 17,520 certified leaves for the refused control. Table corrected (leaves 17,520 (421); surv 45,355). |

## 2. Adjudicated oppositions

### O-1: The headline claim — opus ("cut it down") vs grok ("keep, qualify")
- **opus:** "Every natural strengthening is false" only holds for Instance 4 and part of
  Instance 3; Instances 1–2 show absence/closure, not falsity; the four instances are
  not independent (three sit on the pair-sum cascade); "realistically 1.5 settings."
- **grok:** the headline is an observational summary, not a meta-theorem; exhibited
  four times with checkable witnesses; "valuable and honestly scoped"; proposes
  "every toolkit-opening strengthening examined".
- **VERDICT: PARTIAL — grok's substance, opus's rigor.** opus is right that the
  unrestricted "every" over-claims (Instances 1–2 are audit/closure findings, not
  falsifications) and that "independent" misdescribes the shared cascade. grok is
  right that the response is qualification, not amputation — the paper never claims a
  meta-theorem and the archive is the contribution.
- **v2 actions:** title gains "Examined"; the abstract/intro/§3/§6 replace "four
  independent settings" with "four settings across two reformulation worlds" (three
  manifestations on the pair-sum cascade, one on the zonotope); Observation 3's
  "provably" split by instance (theorem / audit / measured); Table 1's Status column
  distinguishes *false (witness)* / *closed by theorem* / *zero in proved inventory
  (audit)*.

### O-2: The mechanism — opus ("tautology or false") vs grok/gemini2 ("compelling")
- **opus:** "as hard as" undefined; "no instrument invariant under the rewriting can
  distinguish" is true by definition; equivalent rewrites *can* be easier (the paper's
  own finitization is one); "strictly stronger wherever false" is circular.
- **VERDICT: ACCEPT opus on the letter; the substance survives restatement.**
- **v2 actions:** Observation 1 restated: losslessness transfers *truth* exactly
  (proof of the rewrite = proof of LRC and conversely); rewrites may and do change
  tractability — the paper's own finitization bought computability — what losslessness
  forbids is the free lunch of a strictly-stronger-yet-true statement; the tautological
  "no invariant instrument can distinguish" sentence removed. The circular sentence
  in §1.1 rewritten ("strictly stronger (they imply LRC; the converse fails); at the
  extremal instances the demanded margin does not exist").

### O-3: Spearman "no correlation" — opus vs grok
- **opus:** ρ=+0.32 (N=12), ρ=−0.31 (N=8) give p≈0.31/0.46 — cannot support "no
  correlation."
- **grok:** arrays are in the provenance run, checkable.
- **VERDICT: ACCEPT opus** (p-values verified exactly). grok's point (checkability)
  is true but orthogonal: the data are checkable *and* underpowered.
- **v2 actions:** "no measurable correlation" → "no statistically distinguishable
  correlation at the achievable sample sizes (p = 0.32, 0.45); the sign flip between
  adjacent rungs is consistent with noise; these arrays cannot establish correlation
  or its absence."

### O-4: Self-containedness — opus ("none of it can be checked here") vs grok/gemini
- **opus:** the positive core rests on an unpublished companion monograph.
- **grok/gemini:** acceptable for a companion paper; add a one-page self-contained
  verification core (the δ = 1/2 − λ identity, the runner formulation with the
  Lipschitz certificate, τ(u,w) written out, the one-line hand checks).
- **VERDICT: ACCEPT the grok/gemini remedy.**
- **v2 actions:** new **Appendix C (self-contained verification core)** with: the
  four-forms dictionary incl. the elementary identity proof; the explicit τ(u,w)
  formula (gemini's drop-in, adapted); the linear-form/Lipschitz facts; hand
  evaluations of 97/288 (σ = 31/288), 607/1792, and the audit's 4183/12288; and the
  statement of Theorem 4.1 in fully explicit form.

### O-5: m=10 "quantized minimax loss, not an excess" — opus vs the paper's hedge
- **opus:** boxes with U > r are by definition unresolved; residual volume < 6×10⁻⁷
  doesn't exclude a deeper hole; "characterized as pure minimax loss" (§6.1) is
  unproved.
- **VERDICT: ACCEPT opus.** The *diagnostic* (all top U-values exactly
  571/1792 = 7/22 + 9/19712 — a fixed quantization offset) is sound; the
  *characterization* is not.
- **v2 actions:** §5.2/§6.1 reworded to "consistent with quantized minimax loss … a
  deeper point within the unresolved volume is not excluded; the row stays open."

### O-6: Θ(1/n) strictness decay — opus vs gemini's verified table
- **opus:** S ≥ 1 always, so S ~ 30/n must break by n ≈ 26–31; the Newton floor
  (≈1+4/d) accounts for most of the decay — any real-rooted family does this; track
  S − 1 (or S/floor); "nS ≈ 37 at n=3,4" ignores rung-2's 60.5.
- **gemini:** Table 4 arithmetic exact; wants the family label made explicit.
- **VERDICT: ACCEPT opus's mathematics; gemini's table stands.** Computed: S/floor
  declines 3.06→1.99 (harmonic), 5.04→2.00 (rung-2); the 1/n form must break near
  n ≈ 26 (nS ≈ 30 vs Newton floor nS ≥ n + 4).
- **v2 actions:** abstract/§3.5/§A.2 reworded to the precise statement: S ≈ 30/n **in
  the measured range n ≤ 10**; log-concavity forces S ≥ 1 and Newton forces
  S ≥ floor(k) ≈ 1 + 4/d at the middle coefficients, so the 1/n behavior cannot
  persist (it would cross the floor near n ≈ 26); much of the measured decline is the
  generic Newton-floor decline, and the family-specific excess S/floor declines
  3.06→1.99 / 5.04→2.00. "≈37 at n=3,4" labeled harmonic-only; rung-2 n=3 = 60.5
  stated. Table 4 gains Newton-floor and S/floor columns.

### O-7: Instance 3's numbers — opus
- The ~250× requirement vs the first measured gain 95×; exponent m undefined.
- **VERDICT: ACCEPT.** v2 defines m = T−3 at first use and states the k=2 shortfall
  explicitly: the measured thinning factors (95×/1,207×/3,803× at k=2,3,4) *start
  below* the threshold constant (244–256) — the margin must both grow in k and clear
  the constant, which is exactly H1g's open content; the saturation record shows the
  envelope tight where the margin is needed.

## 3. Remaining points, adjudicated

| # | Point (auditor) | Verdict | v2 action |
|---|---|---|---|
| 1 | "Provably" in abstract vs §3.2's hypothesis admission (opus) | ACCEPT | Abstract: "whose proved inventory — audited statement-by-statement — contains only…"; §3.1 splits by instance |
| 2 | 97/288 superseded by 607/1792 in lead positions (opus, grok) | ACCEPT | Abstract and Table 2 now lead with 607/1792, note 4183/12288 (audit); 97/288 kept as the hand-checkable witness |
| 3 | Table 2: add "strictly > thr?" column; 607/1792 row; (1,2,3,4,5,7) δ=1/3 is computed-δ, flag it; emphasize non-extremal failure (grok, gemini) | ACCEPT | All done; boundary rows (m=6,7,8 at n=5; m=9, m=12 at n=6) grouped with a yes/no strict column |
| 4 | "For n ≤ 4 the two statements coincide" mixes center-is-deepest with ρ ≤ thr (opus) | ACCEPT | Rewritten: "the ρ-form *holds* on every instance tested at n ≤ 4 (certified n=3 m ≤ 24, n=4 m ≤ 16); the center-is-deepest property holds exactly at the multiples" |
| 5 | "Rung-2 instances" for m=20,25,30 — only m=10 is rung-2 (opus) | ACCEPT | "the consistent n=5 multiples" |
| 6 | "Three scripts" vs four listed; "two ways" vs three (opus, gemini) | ACCEPT | "Four scripts"; §2.3 "three independent ways" |
| 7 | "Normalized volume Σvⱼ" wrong — it is the relative volume (opus) | ACCEPT | "relative volume Σvⱼ (normalized volume d!·Σvⱼ)"; h\*(1) = d!·Σvⱼ stated with it |
| 8 | gcd-empty-set convention never used (opus) | ACCEPT | Clause removed |
| 9 | "known for at most seven runners" next to "8–15 computer-assisted"; n speeds = n+1 runners never said; n ≥ 2 vs n ≥ 3 (opus) | ACCEPT | Rewritten: "proved for n ≤ 7; for 8 ≤ n ≤ 15 only computational verifications of special families"; observer note added; n ≥ 2/3 harmonized with a remark |
| 10 | Table 4 h\*(t) vs h\*(z) (opus) | ACCEPT | z throughout |
| 11 | Table 5 control row (opus) | ACCEPT | leaves 17,520 (421); surv 45,355 |
| 12 | Bibliography: Cusick missing; Shephard/Frobenius uncited; Lam–Leung metadata; de Bruijn–Schoenberg phantom; Mirsky–Newman unpublished caveat (opus, gemini) | ACCEPT | Cusick 1973 added and cited; Shephard cited (zonotope volume), Frobenius cited (Eulerian); **Crossref verification**: Lam–Leung is *Monthly* 103 (1996) 562–564 **and** *J. Algebra* 224 (2000) 91–109 (the v1 entry merged the two into a phantom 107 (2000)); "de Bruijn–Schoenberg, On the two vanishing identities of Euler" not found in Crossref → replaced by the Rédei–de Bruijn–Schoenberg attribution as carried by Lam–Leung's J. Algebra account; Mirsky–Newman kept with the unpublished caveat |
| 13 | Notation collisions: V speed vector vs quotient space; w pair-member vs bad-set index; T vs n vs 𝕋 (gemini) | ACCEPT | Quotient space renamed 𝒳_v with ‖·‖_{𝒳_v}; bad sets re-indexed B_p, p ∉ {u,w}; T kept (companion's rung notation) with "T = n+1 runners" parenthetical at first use |
| 14 | τ(u,w) under-specified (gemini) | ACCEPT | Explicit formula in §4.1 and Appendix C, incl. the u/w reflection identity |
| 15 | Abstract megasentence (gemini2); "honest" motif (gemini2); §3.5 inline dump → cross-ref (gemini2); §3.2 transition (gemini2) | ACCEPT | Abstract restructured (four short setting-sentences); "deliberately honest" and self-referential audit prose neutralized (2 remaining uses are load-bearing); §3.5 defers to Table 2; transition adopted |
| 16 | Conjecture A.1 formatting (gemini strategic rec) | ACCEPT | Explicit Conjecture (arithmetic h\*-real-rootedness of projected-basis zonotopes), stated decoupled |
| 17 | nS plot (grok) | ACCEPT | Small two-panel figure added to Appendix A (nS and S/floor vs n, both families) |
| 18 | One-line reproduction commands (grok) | ACCEPT | Appendix B: `python3 hstar_n9_n10.py` etc. for all four scripts |
| 19 | Cite the 225-item ledger + negative control in the abstract (grok) | ACCEPT | One clause added |

**Rejected / no change:** none — every audited point was either accepted or accepted
in part. opus's "arithmetic holds up" and gemini's full-table verification stand as
the record's confirmation; no numerical claim of v1 was found false, and one audited
number (the harmonic lower bound) was *improved* by the audit itself.

## 4. v2 deliverable

- `lonely_runner_type_mismatch_paper_v2.pdf` — revised paper (title: *Every Natural
  Strengthening **Examined** Is False: …*).
- `paper2_sources_v2/` — complete LaTeX sources of v2 (v1 files untouched).
- New artifacts: `scripts/audit_verify_v2.py` → `scripts/out_audit_verify_v2.json`
  (the adjudication's verification run: opus-witness two-implementation check,
  Spearman+p, 16-row exact Sturm verification, Newton-floor analysis, h\*(1)
  identity, threshold constants, table-record extraction).
