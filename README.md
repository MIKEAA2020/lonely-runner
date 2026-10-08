# Lonely Runner

Computational audits and research on the **Lonely Runner Conjecture**.

## The papers

The research program culminated in two papers, revised through adjudicated
independent audits (v2), completed by a second certification attempt (v3),
and then closed at the harmonic instance itself (v4, current):

1. **Every Natural Strengthening Examined Is False: Type Mismatch in
   Lossless Reformulations of the Lonely Runner Conjecture** (v4, 30 pp) —
   the flagship.
2. **Pair-Sum Lattices and a Conditional Scaling-Closure Theorem**
   (v4, 67 pp) — the companion monograph (the Hq repair: Theorem SC under
   three quantified open hypotheses; v4 changes no theorem).

The v4 revision certifies the archive's oldest open row:
**rho(1,2,3,4,5) = 16/47 exactly** — attained by a balance family of
witnesses whose three tied pair-dips obey the weighted identity
32 r_(3,5) + 27 r_(4,5) + 35 r_(3,4) = 32, and certified by branch and
bound with two new certificates (a vertex-checked sliding-sigma gap
certificate and a triple-identity certificate) plus lattice-aligned
splitting: 782,255 box evaluations, 0 survivors. The two n=6 boundary
rows are resolved as strict refutations (rho >= 17/47 at m=12,
rho >= 455/1269 at m=9). The v3 revision had certified the last open row
of the n=5 deepest-hole ladder: rho(1,2,3,4,10) = 7/22 exactly.

Full sources, scripts, and machine-checked provenance live in the
[`papers` branch](https://github.com/MIKEAA2020/lonely-runner/tree/papers);
releases:
[v4.0-harmonic-certified](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v4.0-harmonic-certified)
(current — the harmonic instance closed exactly, the n=6 boundary rows
resolved),
[v3.0-m10-certified](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v3.0-m10-certified)
(the m=10 row certified — includes the audit adjudication addendum),
[v2.0-audit-revision](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v2.0-audit-revision)
(the audit revision — includes the audit adjudications), and
[v1.0-papers](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v1.0-papers)
(original).

Honest scope: no progress on the Lonely Runner Conjecture itself is claimed;
the papers establish a systematic negative result (each natural strengthening
of the reformulation framework examined is false, with exact rational
witnesses), a bounded positive core, and a 225-item zero-failure provenance
ledger.

## This branch (`main`)

The upstream computational-audit layer: the `chalf` exhaustive checks
(N=2..7), the `diag2` diagnostic battery, and supporting scripts
(`audits/computational-check/`), plus the auditor reports
(`audits/flaghship.txt`, `audits/companion.txt`) that drove the v2 revision.

License: MIT.
