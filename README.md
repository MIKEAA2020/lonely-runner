# Lonely Runner Conjecture — Pair-Sum Lattices, Lossless Reformulations, and a Type-Mismatch Obstruction

Research program on the Lonely Runner Conjecture (LRC) via lossless reformulations:
pair-sum lattices, constrained transfer cascades, character-rigidity systems, and the
Lonely-Runner zonotope / Ehrhart-theoretic formulation.

**Honest scope statement (stated first, as in the papers).** No progress on the Lonely
Runner Conjecture itself is claimed: no new cases of LRC are proved anywhere in this
repository. What the program establishes is (i) a systematic negative result — every
natural *strengthening* of the reformulation framework proposed as a route to LRC is
false, each refuted by an exact rational witness — (ii) a bounded positive core stated
with exact theorems, and (iii) a machine-reasserted provenance record for every cited
number (225 items, zero failures).

## The papers (v4 harmonic-certification revision — current)

| Paper | File | Contents |
|---|---|---|
| **Flagship v4** — *Every Natural Strengthening Examined Is False: Type Mismatch in Lossless Reformulations of the Lonely Runner Conjecture* | `papers/lonely_runner_type_mismatch_paper_v4.pdf` (30 pages) | **The harmonic instance is closed exactly: rho(1,2,3,4,5) = 16/47.** The third certification session found a *balance family* of witnesses — c(k) = (2k+160, 3k-13, 4k+91, 5k)/235, k = 20..35 — on which the three pair-dips of the runner triple {3,4,5} tie at 16/47, forced by the weighted identity 32 r_(3,5) + 27 r_(4,5) + 35 r_(3,4) = 32 (the c-terms cancel). The branch and bound at r = 16/47 emptied only after two new certificates were added: a vertex-checked *sliding-sigma gap certificate* (for a fixed lift vector m, every c in the box admits a valid sigma iff an interval of affine forms intersects S; the gap is concave in c, so 2^d vertex checks certify the whole box — sigma free to slide along the band) and a *triple-identity certificate* (three pair-dips with weights w_(a,b) = v_c(v_a+v_b) whose weighted values sum to a constant <= (sum w) r over the box), together with lattice-aligned splitting at the dip lattice q = 470. Tree: 782,255 box evaluations, 0 survivors, 5,935 s; 4,841 sampled prunes re-confirmed exactly. The two n=6 boundary rows are resolved as strict refutations: rho >= 17/47 > 5/14 at m=12 (triple {2,5,12}: identity 84 r_(2,5) + 70 r_(2,12) + 34 r_(5,12) = 68; witness (0, 27/47, 0, 25/94, 22/47)) and rho >= 455/1269 > 5/14 at m=9 (witness (13/423, 22/47, 58/423, 695/846, 401/846)). The n=5 boundary rows m=6,7,8 remain open. Sources: `sources/flagship_v4/` (includes certify_harmonic_v4.py and out_certify_harmonic_v4.json). |
| **Companion v4** — *Pair-Sum Lattices and a Conditional Scaling-Closure Theorem* | `papers/lonely_runner_pair_sum_lattices_conditional_theorem_paper_v4.pdf` (67 pages) | No theorem, lemma, or hypothesis changed; the zonotope-strand record updated (the harmonic instance closed exactly at 16/47; the two n=6 boundary rows resolved as strict refutations). Sources: `sources/companion_v4/`. |

## The papers (v3 m=10-certification revision — superseded, kept for the record)

| Paper | File | Contents |
|---|---|---|
| **Flagship v3** — *Every Natural Strengthening Is False: Type Mismatch in Lossless Reformulations of the Lonely Runner Conjecture* | `papers/lonely_runner_type_mismatch_paper_v3.pdf` (28 pages) | **The m=10 row is certified: rho(1,2,3,4,10) = 7/22 exactly.** The second, audit-driven certification attempt paired a sound random-local-search screen (~50x the audit's demonstrated search budget in starts: 4,194,304 random points, full dyadic-16/32 grids, 500k torsion points, 2,048 climbs, 64 Nelder–Mead polishes, 188 exact re-evaluations — nothing above 7/22; the deepest verified point is the torsion center itself) with a new exact Lipschitz box certificate (D is 1-Lipschitz in the sup norm; max_C D <= D(center) + half-width; at the first attempt's residual width, 1135/3584 + 1/1024 = 2277/7168 < 7/22). The re-run branch and bound emptied the tree (190,295 box evaluations, 0 survivors), with the negative control refusing again, the n=4 certificates reproducing identically, 846/846 sampled prunes exact-confirmed, and three adversarial deeper-hole boxes refusing the Lipschitz test. Sup-norm Lipschitz lemma + box-certificate corollary + the 7/22 hand check added to Appendix C. Adjudication addendum: `papers/audit_adjudication_addendum_m10.md`. |
| **Companion v3** — *Pair-Sum Lattices and a Conditional Scaling-Closure Theorem* | `papers/lonely_runner_pair_sum_lattices_conditional_theorem_paper_v3.pdf` (67 pages) | No theorem, lemma, or hypothesis changed; the zonotope-strand record updated (m=10 closed at 7/22; the harmonic lower-bound record harmonized with the audited 4183/12288 witness); version marker recorded. |

## The papers (v2 audit revision — superseded, kept for the record)

| Paper | File | Contents |
|---|---|---|
| **Flagship v2** — *Every Natural Strengthening Examined Is False: Type Mismatch in Lossless Reformulations of the Lonely Runner Conjecture* | `papers/lonely_runner_type_mismatch_paper_v2.pdf` (26 pages) | The audit revision: headline qualified ("examined"); two-worlds framing replaces "four independent"; boundary rows separated from strict refutations with a ">thr?" column; the audit's improved harmonic witness 4183/12288 adopted; exact Sturm verification of h\* real-rootedness; Spearman claims carry p-values; U(C) formula matches the implementation; self-contained verification core (Appendix C); notation collisions fixed; bibliography repaired (Cusick added; Lam–Leung split into its two real papers; phantom de Bruijn–Schoenberg entry replaced). Adjudication: `papers/flagship_audit_adjudication.md`. |
| **Companion v2** — *Pair-Sum Lattices and a Conditional Scaling-Closure Theorem* | `papers/lonely_runner_pair_sum_lattices_conditional_theorem_paper_v2.pdf` (67 pages) | The Hq repair: Theorem SC restated under three quantified open hypotheses (the multi-fiber survivor inequality Hq — new, replacing the unquantified H0 and subsuming H2m — H1g with its uniform margin, H2r as explicit open); the derivation is now the explicit Hq-induction (Proposition); the independent-extension model relabeled as calibration; scope levels (cells / rungs / rigid-zone) fixed; attainment proof expanded; the half case split into its own lemma; statuses labeled OPEN; bibliography repaired. Adjudication: `papers/companion_audit_adjudication.md`. |

## The two papers (v1 — superseded, kept for the record)

| Paper | File | Contents |
|---|---|---|
| **Flagship** — *Every Natural Strengthening Is False: Type Mismatch in Lossless Reformulations of the Lonely Runner Conjecture* | `papers/lonely_runner_type_mismatch_paper.pdf` (20 pages) | The type-mismatch mechanism; four falsified strengthening families with exact witnesses (covering/audit, transfer cascade closed by Lemma M, character rigidity falsified/saturated, the zonotope rho-form false from n=5 onward); the bounded positive core stated exactly; the witness archive with the certification ladder; Appendix A: h\* real-rootedness through n=10 with strictness decay and non-correlation; Appendix B: the 225-item provenance ledger and reproduction map. |
| **Companion monograph** — *Pair-Sum Lattices and a Conditional Scaling-Closure Theorem* | `papers/lonely_runner_pair_sum_lattices_conditional_theorem_paper.pdf` (65 pages) | The full development behind the positive core: the sum-only equality theorem; prime, prime-square and pq structure; the rigid-zone classification (families K4–K37); the complete derivation-and-calibration record; the transfer-cascade monomial dichotomy (Lemma M); the six-entry hypothesis ledger; and the conditional scaling-closure theorem SC with its free-parameter disclosure and projected labels. The final revision records the disposition of the third bridge: the zonotope–Ehrhart instrument engaged and closed by the same type mismatch (h\* real-rooted through n=10; deepest-hole reading false from n=5 with exact witnesses). |

The flagship paper is **not** an abridgement of the monograph. It was written to a
separate, compact specification (mechanism + four instances + bounded positive core +
witness archive + conclusion, plus two appendices) and states the positive-core results
exactly, cross-referencing the monograph for full proofs (Theorem 1.1, Lemmas 2.2–2.7,
Theorem SC of the companion). The monograph in turn cross-references the flagship as
the companion paper where the negative half is worked out. The two are designed to be
submitted together as main paper + companion.

## Repository layout

```
papers/                  Final PDFs (flagship + companion monograph)
sources/flagship/        v1 LaTeX sources + scripts + persisted outputs
sources/flagship_v2/     v2 LaTeX sources (audit revision) + outputs
sources/flagship_v3/     v3 LaTeX sources (m=10 certified) + outputs (incl. certify_m10_v2.py + out_certify_m10_v2.json)
sources/companion/       v1 LaTeX sources + figures + scripts + outputs
sources/companion_v2/    v2 LaTeX sources (Hq repair) + outputs
sources/companion_v3/    v3 LaTeX sources (zonotope-strand record update) + outputs
research_notes/          The 20 dated research notes of the program (T-series record)
scripts/                 All research/verification scripts and their out_*.json
                         outputs (the canonical provenance record)
UPSTREAM.md              Note on the upstream computational-audits repository
LICENSE                  MIT
```

## Key machine-checkable results

- **Refuted (exact rational witnesses):** the zonotope "deepest hole = center" rho-form
  fails at n=5 harmonic m=5 (rho ≥ 97/288, improved to 607/1792) and at three further
  instances through n=6, including a non-extremal instance (n=6, v=(1,2,3,4,5,7),
  witness (1/16, 9/16, 0, 3/8, 0), D = 4/11 > 5/14).
- **Certified (exact branch-and-bound minimax certificates):** n=3 cores; n=4 m ≤ 16
  (m=16 at 9/34, 222 leaves); n=5 m = 20, 25, 30 (13/42, 4/13, 19/62) — with a negative
  control that correctly refuses the harmonic n=5 instance.
- **Open (honestly recorded):** n=5 m=10 — one serious certification attempt
  (61,038 exact leaves, residual volume < 6e-7, minimax-loss signature
  1135/3584 < 7/22 at centers, top bound U = 571/1792), recorded as open.
- **Transfer cascade:** every constrained transfer matrix of the lossless cascade is a
  partial permutation with spectral radius exactly 0 or 1 (Lemma M) — the margin is
  provably external to the lossless dynamics.
- **h\* polynomials:** real-rooted, log-concave, Newton through n=10 for both families
  (45-dps exact-integer arithmetic); strictness decay ~33/n; Spearman rank correlation
  between strictness and gap ≈ +0.32/−0.31 (noise).

## Reproducing

- Papers: `cd sources/flagship && tectonic main.tex` (same for `sources/companion`;
  covers are pre-rendered HTML → PDF and merged; merge scripts included).
- Scripts: plain Python 3 (stdlib + `numpy`/`sympy`/`mpmath` where noted). Each script
  writes an `out_*.json` next to itself; the committed `out_*.json` files are the
  persisted runs cited by the papers. The flagship's Appendix B maps every cited
  number to script + output file.
- Provenance re-assertion: `python3 sources/flagship/paper_provenance.py` re-checks all
  225 inherited cited numbers against the persisted JSONs (expected: 0 failures).
- Audit adjudication verification (v2): `python3 scripts/audit_verify_v2.py` re-checks
  the auditors' claims (the 4183/12288 witness, Spearman p-values, 16-row exact Sturm
  real-rootedness, Newton floors, threshold constants).

## Publishing this repository

The full internal session history (33 commits, complete worklog, every intermediate
state) is preserved in the working repository and exported as a git bundle
(`lrc_full_provenance.bundle`) alongside this release. To publish this curated release:

```bash
git remote add origin https://github.com/<user>/<repo>.git
git push -u origin main
```

## License

MIT — see `LICENSE`.
