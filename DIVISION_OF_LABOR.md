# Division of Labor & Implementation Plan — Main Text vs. Supplementary Information

Guiding principle (user's): **main text = a complete, compelling, self-contained
story for the general reader; SI = resource for the specialist (verification,
depth, reproducibility).** Essential-to-the-conclusion → main; supportive/
detailed/tangential → SI. Every SI item must be cited in the main text. No
"data dump"; no key findings hidden in SI.

## Current page budgets (measured from the compiled PDFs)

**Flagship** (31 pp): Intro 279 ln; Preliminaries 151; type-mismatch 337;
positive core 137; witness archive 555 ln (≈8 pp); Conclusion 107; App A
(h*) 203 ln; App B (provenance) 230 ln; App C (verification core) 180 ln.

**Companion** (69 pp): Intro 387; main theorem (conditional) 381; equality 291;
covering 322; rigid zone 798 ln (≈9 pp); group-ring 549; pq 288; higher rungs
253; **scaling-closure 1290 ln (≈15 pp)**; measured coverage 222; open problems
232; reproducibility 312; App A (h*).

## Division of labor — the decision rule applied

### FLAGSHIP (negative-results paper)
- **Main text (the story):** Intro; Preliminaries; §3 type-mismatch mechanism +
  the four instances (the core claim); §4 positive core (bounded statement);
  §5 **condensed** witness table (only the headline exact witnesses that drive
  the narrative: 16/47, 17/47, 29/80, 4/11, 1793/5000, 7/22); Conclusion.
- **SI (specialist):** full witness archive (the 555-line §5 detail, the
  per-family certificates, the negative-control refusals); the h*-observations
  (App A); the computational provenance (App B); the self-contained
  verification core (App C).

### COMPANION (conditional-program paper)
- **Main text (the story):** Intro; §1 main theorem (conditional) + ledger;
  §2 pair-sum equality (the unconditional core, with proof); §3 covering
  framework; §4 rigid-zone classification (the K-theorems — the unconditional
  classifications, with proofs); §5 measured coverage + the decay law (the
  empirical narrative, incl. the new descent corollary's consequence); §6 open
  problems.
- **SI (specialist):** the **scaling-closure derivation** (the 1290-line
  `secscale` — this is the "specialist" heart: the cascade model, pool/cut/
  volume laws, the strand census, the transfer-matrix analysis, the calibration
  record — none of it is needed to *state* the theorem, only to *derive* it);
  the group-ring translation proofs (the computational details); the pq tensor
  proofs; the higher-rung census tables; the reproducibility/soundness record
  (the 312-line ledger); the h* appendix.

## Why this split (essential vs. supportive)

- **Companion scaling-closure (`secscale`, 15 pp) → SI.** It is the *derivation*
  of the theorem's hypotheses, not the theorem. A general reader needs the
  theorem's statement, the ledger's shape, and the named open inputs; the
  cascade/pool/cut/volume calibration is specialist verification. This is the
  single biggest page saving and the clearest "essential vs. supportive" call.
- **Companion rigid-zone (§4) stays in main.** The K-classifications are
  *unconditional proved results* and the paper's positive content; they are the
  story, not support.
- **Flagship witness archive (§5) → SI.** The headline witnesses are the
  narrative; the per-family certificate detail is verification. (Already
  partly appendices; the remaining 555-line §5 detail moves to SI.)
- **Both: reproducibility/soundness ledger → SI.** It is verification
  infrastructure, not results.
- **Both: h* appendix → SI.** Tangential to the core conclusion.

## Implementation plan (one paper per turn)

### Turn 1 — Companion
1. **Extract SI:** create `sources/companion_v5/si/` with `si_scale.tex`
   (the whole `secscale`, minus the theorem-statement pointer), `si_gr.tex`
   (group-ring proof details), `si_pq.tex`, `si_higher.tex`, `si_sound.tex`,
   `si_hstar.tex`. Each `\input`ed from a new `main_si.tex`.
2. **Main text:** replace `\input{secscale}` with a short "the closure theorem
   and its hypotheses" section (the theorem statement, the ledger table, and a
   one-paragraph pointer to SI §S1 for the derivation). Keep the descent lemma
   + small-p corollary in main (they are short and are proved results).
3. **Every SI item cited:** add "(see SI §S1)" etc. at each main-text claim that
   rests on a moved proof.
4. **Rebuild:** `main.tex` (main, target ≤ ~45 pp) and `main_si.tex` (SI);
   verify refs resolve across the two documents (use `xr` or a shared `.aux`,
   or restate the few cross-document labels).

### Turn 2 — Flagship
1. **Extract SI:** `si_witness.tex` (the 555-line §5 detail), `si_prov.tex`
   (App B), `si_self.tex` (App C), `si_hstar.tex` (App A).
2. **Main text:** §5 becomes the condensed witness table + the headline
  certificates; appendices become a "Supplementary Information" document.
3. **Cite every SI item** in the main text.
4. **Rebuild** both; verify.

### Cross-cutting (both turns)
- Cross-document references: either (a) compile main+SI as one PDF with
  `\appendix` (simplest, keeps refs working — "SI" as appendices), or (b) use
  the `xr` package for cross-document `\ref`. Recommend (a) for robustness:
  main = §1–§N, then `\appendix` + SI sections, all in one build, with the SI
  clearly demarcated and cited from the main text.
- Keep every proof in full (user's standing instruction): moving to SI ≠
  condensing; the SI must contain the *complete* proofs, not summaries.

## Decisions to confirm with the user
1. Deliver main + SI as **one PDF** (main, then `\appendix` SI) or **two PDFs**
  (main + separate SI)? (Recommend one PDF with a clear SI boundary — simplest
  and keeps cross-refs working.)
2. Companion: is the **scaling-closure derivation** (15 pp) acceptable to move
  wholesale to SI, given it is the derivation of the paper's central
  conditional theorem? (Recommend yes — it is the clearest "supportive vs.
  essential" case.)
3. Flagship: move the full witness archive to SI, keeping only the headline
  witnesses in main? (Recommend yes.)