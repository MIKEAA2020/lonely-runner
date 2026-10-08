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

## The two papers

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
sources/flagship/        LaTeX sources, cover, and the certification/provenance
                         scripts with their persisted outputs (builds with tectonic)
sources/companion/       LaTeX sources, figures, cover, and the derivation scripts
                         with their persisted outputs (builds with tectonic)
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
