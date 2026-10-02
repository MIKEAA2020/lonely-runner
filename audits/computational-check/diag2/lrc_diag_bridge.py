"""Attribution experiment: what kills would-be deficit-law violations?

Bridge region for size n: m(n+1) < b <= n(2m+1)+1  (deficit law wants
b <= m(n+1); the covering/orbit lemma only gives b <= n(2m+1)+1; beyond
the covering wall the union-bound kills trivially, so the bridge is where
the factor-2 obstruction lives).

A bridge CONFIG is a residue set R = {m} U interiors U {b-m} at time a/b.
Kill levels:
  L2 (orbit / b-grid): some j in [1, b-1] has min_i dist(j r_i, 0) > m,
      i.e. the b-grid time j/b beats the config -- a necessary condition
      of the argmax (f(j a/b) <= gap). Includes the divisor lemma
      (j = b/d) and unit kills.
  L3 (off-grid crossings): for every unit s mod b, the canonical lift
      speeds c_i = s r_i mod b (in [1, b-1]) has gap > m/b -- the config
      is killed by pair-sum crossings with denominators not tied to b.
  REALIZED: some lift has gap == m/b exactly -> an actual deficit-law
      violation (bug or counterexample); re-verified with the exact
      reference implementation.

The L2 share as a function of n is the tractability gradient of the
residue reformulation: it measures how much of the bridge the b-grid
arithmetic (the "clothing") closes on its own.

Usage: python3 lrc_diag_bridge.py out.json
"""
import sys
import json
import random
from collections import Counter
from fractions import Fraction
from math import gcd, comb
import time

import numpy as np

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_diag_gap import GapEngine, reduce_candidate, gap_int

CAP = 80          # max denominator considered
CFG_CAP = 500     # max configs enumerated/sampled per (n, m, b)
CFG_CAP_BIG = 150  # for n >= 8
UNIT_CAP_BIG = 12  # sampled units for n >= 8


def orbit_fail(R, b, m):
    """All j in [1, b-1] with min_i dist(j r_i mod b, 0) > m (empty=pas s)."""
    Ra = np.array(R, dtype=np.int64)
    j = np.arange(1, b, dtype=np.int64)
    M = (j[:, None] * Ra[None, :]) % b
    dist = np.minimum(M, b - M)
    mins = dist.min(axis=1)
    return j[mins > m]


def bridge_scan(eng):
    results = {}
    for n in range(3, 11):
        m_list = list(range(1, 9)) if n == 4 else [1, 2, 3]
        agg = dict(configs=0, L2=0, L3=0, realized=0,
                   L2_divisor=0, L2_unit_only=0,
                   realized_detail=[], winner_bgrid=0, winner_offgrid=0,
                   winner_below_b=0, winner_above_b=0,
                   winner_deficits=Counter(), per_mb=[])
        cfg_cap = CFG_CAP if n < 8 else CFG_CAP_BIG
        for m in m_list:
            lo, hi = m * (n + 1) + 1, min(n * (2 * m + 1) + 1, CAP)
            for b in range(lo, hi + 1):
                pool = list(range(m + 1, b - m))
                need = n - 2
                if need < 0:
                    continue
                total_cfgs = comb(len(pool), need) if need <= len(pool) else 0
                if total_cfgs == 0:
                    continue
                rng = random.Random(12345 + 97 * n + 13 * m + b)
                if total_cfgs <= cfg_cap:
                    cfgs = [list(c) for c in __import__("itertools").combinations(pool, need)]
                else:
                    cfgs = [sorted(rng.sample(pool, need)) for _ in range(cfg_cap)]
                units = [s for s in range(1, b) if gcd(s, b) == 1]
                if n >= 8 and len(units) > UNIT_CAP_BIG:
                    units = rng.sample(units, UNIT_CAP_BIG)
                mb = dict(n=n, m=m, b=b, configs=len(cfgs),
                          sampled=total_cfgs > cfg_cap,
                          L2=0, L3=0, realized=0)
                target = Fraction(m, b)
                for interiors in cfgs:
                    R = [m] + interiors + [b - m]
                    agg["configs"] += 1
                    fails = orbit_fail(R, b, m)
                    if len(fails):
                        agg["L2"] += 1
                        mb["L2"] += 1
                        div = any(int(gcd(int(j), b)) > 1 for j in fails)
                        if div:
                            agg["L2_divisor"] += 1
                        else:
                            agg["L2_unit_only"] += 1
                        continue
                    # L3: canonical lifts over units
                    realized = False
                    win = None
                    for s in units:
                        c = sorted((s * r) % b for r in R)
                        g = eng.gap(c)
                        if g == target:
                            realized = True
                            gref, _ = gap_int(c)
                            agg["realized_detail"].append(
                                dict(n=n, m=m, b=b, R=R, s=s, speeds=c,
                                     gap=str(g), gap_ref=str(gref)))
                            break
                        if s == 1 and win is None:
                            g1, (mp, D, jp) = eng.gap(c, with_argmax=True)
                            a2, b2, m2 = reduce_candidate(mp, D, jp)
                            win = (D, b2, m2, m2 * (n + 1) - b2)
                    if realized:
                        agg["realized"] += 1
                        mb["realized"] += 1
                        continue
                    agg["L3"] += 1
                    mb["L3"] += 1
                    if win:
                        D, b2, m2, dwin = win
                        if D % b == 0:
                            agg["winner_bgrid"] += 1
                        else:
                            agg["winner_offgrid"] += 1
                        if D < b:
                            agg["winner_below_b"] += 1
                        else:
                            agg["winner_above_b"] += 1
                        agg["winner_deficits"][dwin] += 1
                agg["per_mb"].append(mb)
        wd = {str(k): v for k, v in sorted(agg["winner_deficits"].items())}
        results[n] = dict(
            configs=agg["configs"], L2=agg["L2"], L3=agg["L3"],
            realized=agg["realized"],
            L2_share=round(agg["L2"] / agg["configs"], 4) if agg["configs"] else None,
            L2_divisor=agg["L2_divisor"], L2_unit_only=agg["L2_unit_only"],
            winner_bgrid=agg["winner_bgrid"], winner_offgrid=agg["winner_offgrid"],
            winner_below_b=agg["winner_below_b"], winner_above_b=agg["winner_above_b"],
            winner_deficit_hist=wd,
            realized_detail=agg["realized_detail"],
            per_mb=[{k: v for k, v in mb.items()} for mb in agg["per_mb"]],
        )
        print(f"n={n}: configs={agg['configs']} L2={agg['L2']} "
              f"({results[n]['L2_share']}) L3={agg['L3']} realized={agg['realized']} "
              f"[divisor-kills {agg['L2_divisor']}, unit-only {agg['L2_unit_only']}] "
              f"winners: b-grid {agg['winner_bgrid']} / off-grid {agg['winner_offgrid']}")
    return results


if __name__ == "__main__":
    eng = GapEngine(CAP)
    t0 = time.time()
    res = bridge_scan(eng)
    with open(sys.argv[1], "w") as f:
        json.dump(res, f, indent=1)
    print(f"total {round(time.time()-t0,1)}s")
