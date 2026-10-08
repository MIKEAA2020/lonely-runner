# Lonely Runner

Computational audits and research on the **Lonely Runner Conjecture**.

## The papers

The research program culminated in two papers (full sources, scripts, and
machine-checked provenance in the [`papers` branch](https://github.com/MIKEAA2020/lonely-runner/tree/papers)
and the [v1.0-papers release](https://github.com/MIKEAA2020/lonely-runner/releases/tag/v1.0-papers)):

1. **Every Natural Strengthening Is False: Type Mismatch in Lossless
   Reformulations of the Lonely Runner Conjecture** (20 pp) — the flagship.
2. **Pair-Sum Lattices and a Conditional Scaling-Closure Theorem** (65 pp) —
   the companion monograph.

Honest scope: no progress on the Lonely Runner Conjecture itself is claimed;
the papers establish a systematic negative result (each natural strengthening
of the reformulation framework is false, with exact rational witnesses), a
bounded positive core, and a 225-item zero-failure provenance ledger.

## This branch (`main`)

The upstream computational-audit layer: the `chalf` exhaustive checks
(N=2..7), the `diag2` diagnostic battery, and supporting scripts
(`audits/computational-check/`).

License: MIT.
