# Lonely Runner

Computational audits and research on the **Lonely Runner Conjecture**.

## The papers

The research program culminated in two papers, revised through adjudicated
independent audits (v2), completed by two certification closures (v3, v4),
and unified in v5 (current) with the weight identity promoted to a proved
lemma:

1. **Every Natural Strengthening Examined Is False: Type Mismatch in
   Lossless Reformulations of the Lonely Runner Conjecture** (v5, 32 pp) —
   the flagship.
2. **Pair-Sum Lattices and a Conditional Scaling-Closure Theorem**
   (v4, 67 pp) — the companion monograph (the Hq repair: Theorem SC under
   three quantified open hypotheses; unchanged in v5).

The v5 revision unifies the record and proves the arithmetic behind it:
the pair-character identity c*chi_ab + a*chi_bc = b*chi_ac (a one-line
lemma) forces the denominators of the exact values through two binding
archetypes — triple binding (rho = N/(2(ab+ac+bc)); the harmonic
**rho(1,2,3,4,5) = 16/47** and the rung-2 17/47) and pair bisection
(rho = N/(2(u+v)); the m=7 witness 4/11) — and all four certified ladder
values (7/22, 13/42, 4/13, 19/62) are pair-bisecting at the center class
on the pair (1,m): the two archetypes govern every exact value in the
archive. The harmonic closure's triple-identity certificate is the lemma
in box form. m=9 is a strict refutation (1793/5000 > 5/14) with the exact
value unknown; the global equality at m=12 and m=6 is honestly open.
Every witness-level claim of the two closure runs was re-asserted against
the independent implementation before incorporation (11/11 checks).

Full sources, scripts, and machine-checked provenance live in the
[`papers` branch](https://github.com/MIKEAA2020/lonely-runner/tree/papers);
releases:
[v5.0-weight-lemma](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v5.0-weight-lemma)
(current — the unified revision: the weight lemma proved, both closures
incorporated, the archetype governance of every exact value),
[v4.0-harmonic-certified](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v4.0-harmonic-certified)
(the harmonic instance closed exactly, the n=6 boundary rows
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
