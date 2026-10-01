Computational check of the Lonely Runner polytope integrality (n <= 5, v_n <= 50)
================================================================================

Contents
  lrc_ilp_check.c   exhaustive solver (exact integer arithmetic) with two modes:
                      run   n V prefix   full sweep, writes dump + tight list + summary
                      brute n            literal k-box enumeration, vectors on stdin
                    Build: gcc -O2 -o lrc_ilp_check lrc_ilp_check.c
  lrc_verify.py     validation driver: re-verifies every witness, runs an
                    independent reference implementation, brute-force and HiGHS
                    MILP cross-checks, verifies the source document's claims,
                    and writes stats.json. Run: python3 lrc_verify.py
                    (expects lrc_ilp_check in the same directory and out/ data).
  stats.json        machine-readable aggregate results of the validation run.
  tight_n*.txt      complete tight-set census (gap = 1/(n+1)) with witnesses.

Reproducing the full sweep (about 30 seconds, single core):
  gcc -O2 -o lrc_ilp_check lrc_ilp_check.c
  mkdir -p out
  for n in 1 2 3 4 5; do ./lrc_ilp_check run $n 50 out/n${n}_V50; done
  # dumps (74 MB total) land in out/; SHA-256 in "computational check results.txt"
  python3 lrc_verify.py     # re-runs every validation layer -> stats.json

Note: lrc_verify.py as committed reads/writes absolute paths under
/home/z/my-project/scripts/out; adjust OUT/BIN at the top when running elsewhere.
