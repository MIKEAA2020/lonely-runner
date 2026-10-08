# Lonely Runner

Computational audits and research on the **Lonely Runner Conjecture**.

## The papers

The research program culminated in two papers, revised through adjudicated
independent audits (v2, current):

1. **Every Natural Strengthening Examined Is False: Type Mismatch in
   Lossless Reformulations of the Lonely Runner Conjecture** (v2, 26 pp) —
   the flagship.
2. **Pair-Sum Lattices and a Conditional Scaling-Closure Theorem**
   (v2, 67 pp) — the companion monograph (the Hq repair: Theorem SC under
   three quantified open hypotheses).

Full sources, scripts, and machine-checked provenance live in the
[`papers` branch](https://github.com/MIKEAA2020/lonely-runner/tree/papers);
releases:
[v2.0-audit-revision](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v2.0-audit-revision)
(current — includes the audit adjudications) and
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
