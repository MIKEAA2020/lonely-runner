"""Top-pair rescue law: empirical coverage test across n=4..10.

Claim under test: for a bridge config (m(n+1) < b <= n(2m+1)+1) that
survives the orbit lemma, the canonical-family lift with speeds
c_i = s*r_i mod b is killed by the pair-sum crossing family of the TWO
LARGEST speeds of that lift (the residues nearest the down-pole b-m).

Coverage is measured per (config, unit-lift): does the top pair's best
crossing value exceed m/b?  When it fails, record which pairs do win.

Usage: python3 lrc_diag_rescue.py out.json
"""
import sys
import json
import random
import itertools
from collections import Counter
from fractions import Fraction
from math import gcd, comb
import time

import numpy as np

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_diag_gap import GapEngine
from lrc_diag_bridge import orbit_fail

CAP = 80


def pair_best_value(V, u, w, eng):
    """Best envelope value on the crossing family t = j/(u+w)."""
    D = u + w
    best = -1
    for j in range(D):
        mm = min(min((v * j) % D, D - (v * j) % D) for v in V)
        if mm > best:
            best = mm
    return Fraction(best, D)


def rescue_scan(eng):
    results = {}
    for n in range(4, 11):
        m_list = list(range(1, 9)) if n == 4 else [1, 2, 3]
        cfg_cap = 500 if n < 8 else 150
        unit_cap = None if n < 8 else 12
        tot = topp = 0
        fails = []
        winners_when_fail = Counter()
        t0 = time.time()
        for m in m_list:
            lo, hi = m * (n + 1) + 1, min(n * (2 * m + 1) + 1, CAP)
            for b in range(lo, hi + 1):
                pool = list(range(m + 1, b - m))
                need = n - 2
                if need > len(pool):
                    continue
                rng = random.Random(999 + 31 * n + m)
                ncfg = comb(len(pool), need)
                if ncfg <= cfg_cap:
                    cfgs = list(itertools.combinations(pool, need))
                else:
                    cfgs = [tuple(sorted(rng.sample(pool, need)))
                            for _ in range(cfg_cap)]
                units = [s for s in range(1, b) if gcd(s, b) == 1]
                if unit_cap and len(units) > unit_cap:
                    units = rng.sample(units, unit_cap)
                target = Fraction(m, b)
                for interiors in cfgs:
                    R = [m] + list(interiors) + [b - m]
                    if len(orbit_fail(R, b, m)):
                        continue  # orbit already kills; not a rescue case
                    for s in units:
                        c = sorted((s * r) % b for r in R)
                        val = pair_best_value(c, c[-2], c[-1], eng)
                        tot += 1
                        if val > target:
                            topp += 1
                        else:
                            wins = []
                            for u, w in itertools.combinations(c, 2):
                                if pair_best_value(c, u, w, eng) > target:
                                    wins.append((u, w))
                            winners_when_fail[len(wins)] += 1
                            if len(fails) < 25:
                                fails.append(dict(
                                    n=n, m=m, b=b, R=R, s=s, lift=c,
                                    top_val=str(val), target=str(target),
                                    n_winning_pairs=len(wins),
                                    winning_pairs=[list(p) for p in wins[:6]]))
        results[n] = dict(
            config_lifts=tot, top_pair_rescues=topp,
            coverage=round(topp / tot, 4) if tot else None,
            fail_winner_pair_counts=dict(winners_when_fail),
            fails=fails,
            seconds=round(time.time() - t0, 1),
        )
        print(f"n={n}: orbit-surviving (config, lift) pairs={tot}; "
              f"top-pair rescue coverage={results[n]['coverage']}")
    return results


if __name__ == "__main__":
    eng = GapEngine(CAP)
    res = rescue_scan(eng)
    with open(sys.argv[1], "w") as f:
        json.dump(res, f, indent=1)
