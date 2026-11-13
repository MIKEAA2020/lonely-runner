#!/usr/bin/env python3
"""
Higher-rung step 3: confirmation experiments.

(X1) T=7 census at p=113 (eps=1, k=16) -- the tail-decay law at scale.
(X2) DIRECT sampled ground truth for the T=8 (1K,5U) cell at p in
     {17, 19}: random unit 5-tuples (distinct bad sets), count
     consecutive good fibers -- no census restriction (the m=5 census is
     a zoo; this measures the drift-exit mechanism directly).
(X3) DIRECT sampled ground truth for T=7 (1K,4U) at p in {23, 29} --
     cross-check of the census-restricted GT.

Output: scripts/out_higher_rung_step3.json
"""
import json
import random
import sys
import time
from collections import Counter

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import (profile_sizes, ball_mask, census_m4,
                               family_stats)

OUT = {}


def unit_footprints(p, T, v):
    """footprint[v][j] = bitmask of B_v cap F_j in t-coords."""
    N = p * p
    c = (N - 1) // T
    tab = []
    for j in range(p):
        m = 0
        for t in range(p):
            x = j + p * t
            wx = (v * x) % N
            if min(wx, N - wx) <= c:
                m |= 1 << t
        tab.append(m)
    return tab


def direct_sampled_gt(p, T, m_units, nsamp, seed=20261005):
    """Random normalized families {B_{pa'}, B_1, B_{v2..v_m}}: count
    consecutive good fibers from the first uncovered fiber."""
    N = p * p
    c = (N - 1) // T
    k_eff = c // p
    bm = ball_mask(p, k_eff)
    full = (1 << p) - 1
    U_all = [j for j in range(p) if not (bm >> j) & 1]
    rng = random.Random(seed)
    units = [v for v in range(2, N) if v % p and v != N - 1]
    best = 0
    best_desc = None
    n_good_total = 0
    n_fams = 0
    t0 = time.time()
    foot_cache = {}

    def foot(v):
        if v not in foot_cache:
            foot_cache[v] = unit_footprints(p, T, v)
        return foot_cache[v]

    for _ in range(nsamp):
        ap = rng.randrange(1, p)
        ainv = pow(ap, -1, p)
        U = [ainv * j % p for j in U_all]
        # sample m distinct bad-set units
        vs = []
        guard = 0
        while len(vs) < m_units and guard < 200:
            v = rng.choice(units)
            ok = all((v - w) % N not in (0,) and (v + w) % N != 0
                     for w in [1] + vs)
            if (v - 1) % N == 0:
                ok = False
            if ok:
                vs.append(v)
            guard += 1
        if len(vs) < m_units:
            continue
        n_fams += 1
        tabs = [foot(1)] + [foot(v) for v in vs]
        ng = 0
        for j in U:
            mm = 0
            for tb in tabs:
                mm |= tb[j]
            if mm == full:
                ng += 1
            else:
                break
        n_good_total += ng
        if ng > best:
            best = ng
            best_desc = (ap, tuple(vs))
    return {'p': p, 'T': T, 'm': m_units, 'n_fams': n_fams,
            'best_run': best, 'best_desc': best_desc,
            'mean_run': n_good_total / max(1, n_fams),
            'n_U': len(U_all), 'elapsed': time.time() - t0}


def main():
    t0 = time.time()
    print('=== higher_rung_step3 ===')
    # (X2) T=8 (1K,5U) sampled
    for p in (17, 19):
        r = direct_sampled_gt(p, 8, 5, 20000)
        OUT['samp_T8_1K5U_p%d' % p] = r
        print('  [X2] p=%d: %d fams, best run %d/%d, mean %.2f (%.0fs)' %
              (p, r['n_fams'], r['best_run'], r['n_U'], r['mean_run'],
               r['elapsed']))
        sys.stdout.flush()
    # (X3) T=7 (1K,4U) sampled cross-check
    for p in (23, 29):
        r = direct_sampled_gt(p, 7, 4, 20000)
        OUT['samp_T7_1K4U_p%d' % p] = r
        print('  [X3] p=%d: %d fams, best run %d/%d, mean %.2f (%.0fs)' %
              (p, r['n_fams'], r['best_run'], r['n_U'], r['mean_run'],
               r['elapsed']))
        sys.stdout.flush()
    OUT['elapsed_s'] = time.time() - t0
    json.dump(OUT, open('/home/z/my-project/scripts/'
                        'out_higher_rung_step3.json', 'w'), indent=1)
    print('DONE %.0fs' % OUT['elapsed_s'])


if __name__ == '__main__':
    main()
