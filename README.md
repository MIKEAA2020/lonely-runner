# Lonely Runner

Computational audits and research on the **Lonely Runner Conjecture**.

## The papers

The research program culminated in two papers, revised through adjudicated
independent audits (v2) and then completed by a second certification
attempt (v3, current):

1. **Every Natural Strengthening Examined Is False: Type Mismatch in
   Lossless Reformulations of the Lonely Runner Conjecture** (v3, 28 pp) —
   the flagship.
2. **Pair-Sum Lattices and a Conditional Scaling-Closure Theorem**
   (v3, 67 pp) — the companion monograph (the Hq repair: Theorem SC under
   three quantified open hypotheses; v3 changes no theorem).

The v3 revision certifies the last open row of the n=5 deepest-hole ladder:
**rho(1,2,3,4,10) = 7/22 exactly** (an audit-driven second attempt: a sound
random-local-search screen at 4.7M evaluations found no excess, and a new
exact Lipschitz box certificate emptied the branch-and-bound tree).

Full sources, scripts, and machine-checked provenance live in the
[`papers` branch](https://github.com/MIKEAA2020/lonely-runner/tree/papers);
releases:
[v3.0-m10-certified](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v3.0-m10-certified)
(current — includes the audit adjudication addendum),
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
