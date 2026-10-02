Computational check of the Lonely Runner polytope integrality (n <= 5, v_n <= 50)
Extended to n = 6 (v_n <= 50) and n = 7 (v_n <= 38)
================================================================================

Contents
  lrc_ilp_check.c   exhaustive solver (exact integer arithmetic) with two modes:
                      run   n V prefix   full sweep, writes dump + tight list + summary
                      brute n            literal k-box enumeration, vectors on stdin
                    Build: gcc -O2 -o lrc_ilp_check lrc_ilp_check.c
  lrc_verify.py     validation driver for the n <= 5 check: re-verifies every
                    witness, runs an independent reference implementation,
                    brute-force and HiGHS MILP cross-checks, verifies the source
                    document's claims, and writes stats.json. Run: python3
                    lrc_verify.py (expects lrc_ilp_check in the same directory
                    and out/ data).
  lrc_extend67.py   extension driver for n = 6 / n = 7: full witness
                    re-verification (L1), independent reference (L2), literal
                    brute enumeration (L3), HiGHS MILP (L4), census statistics
                    (tight sets, second-tightest gap, gap distribution, modular
                    witness checks, targeted family lookups) -> stats67.json.
  lrc_theory_experiments.py
                    experiments behind "integrality analysis.txt" and
                    "ladder analysis.txt": exact vertex enumeration of P(v) via
                    oriented spanning trees of the constraint network (E1),
                    census witness structure incl. fiber widths (E2), modular
                    rung census on the n = 4..7 dumps (E3), Smith normal forms
                    of the pair matrix (E4) -> theory_experiments.json.
  stats.json        machine-readable aggregate results of the n <= 5 validation.
  stats67.json      machine-readable aggregate results of the n = 6 / n = 7 run.
  theory_experiments.json   machine-readable E1-E4 data.
  tight_n*.txt      complete tight-set census (gap = 1/(n+1)) with witnesses,
                    n = 1..5 (v_n <= 50), n = 6 (v_n <= 50), n = 7 (v_n <= 38).

Reproducing the full sweep (about 30 seconds, single core, n <= 5):
  gcc -O2 -o lrc_ilp_check lrc_ilp_check.c
  mkdir -p out
  for n in 1 2 3 4 5; do ./lrc_ilp_check run $n 50 out/n${n}_V50; done
  # dumps (74 MB total) land in out/; SHA-256 in "computational check results.txt"
  python3 lrc_verify.py     # re-runs every validation layer -> stats.json

Reproducing the extension (about 11 minutes, single core):
  ./lrc_ilp_check run 6 50 out/n6_V50     # 15,890,700 vectors, ~5.5 min
  ./lrc_ilp_check run 7 38 out/n7_V38     # 12,620,256 vectors, ~5.5 min
  python3 lrc_extend67.py                 # L1-L4 + census -> stats67.json
  python3 lrc_theory_experiments.py       # E1-E4 -> theory_experiments.json
  # dumps (508 + 404 MB) are not committed; SHA-256 in "computational check
  # results.txt" section 9. V = 38 for n = 7 is the largest bound completing
  # within the original environment's per-process wall-time limit.

Note: lrc_verify.py and lrc_extend67.py as committed read/write absolute paths
under /home/z/my-project/scripts/out; adjust OUT/BIN at the top when running
elsewhere.
