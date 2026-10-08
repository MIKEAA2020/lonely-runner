# Upstream repository

This program builds on the upstream computational-audits repository
(https://github.com/MIKEAA2020/..., "Add files via upload" history, MIT license),
which contains the original `audits/computational-check` material:
the `chalf` exhaustive checks (N=2..7), the `diag2` diagnostic battery
(bridge, casebook, climb, critical-hunt, gap), and the supporting
`lrc_chalf_patch.py` / `lrc_chalf_probe.py` / `lrc_exhaust_check.py` scripts.

That repository is owned and published separately by its author and is **not**
duplicated here. The program's papers cite it as the upstream audit layer; all
results in this repository were produced by the scripts under `scripts/` and
`sources/*/`, whose persisted `out_*.json` outputs constitute the provenance
record (see the flagship paper's Appendix B).

A full-fidelity copy of the upstream tree as it existed during the program
(171 tracked files, including every audit report) is preserved in the internal
working repository and in the exported `lrc_full_provenance.bundle`.
