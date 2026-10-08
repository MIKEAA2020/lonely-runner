#!/usr/bin/env python3
"""
h2_pair_check.py -- zone-resolved check of H2 (dilate-sparsity) on the
COMMITTED cells/primes only (same scope as scaling_closure_verify.py;
per-pair recording added).

H2 (precise form under test): away from resonant dilations, the
two-fiber survivor count
    surv = |{a in Sol0 : lambda*a + tau in Sol_j}|
satisfies surv <= gamma0 * |Sol0|*|Sol_j|/p^(m-1) with gamma0 = O(1)
uniform over cells/families/fibers.  Resonance classifier: b_min(lambda,
L) = min{b >= 1 : ||lambda*b||_p <= L} with L = 2k+2 (the arc-size
scale = the per-coordinate box scale of the chain pullback).  Zones:
  small  : lambda_bar <= 4   (the classic small-dilate resonance)
  large  : lambda_bar >= p-4
  mid    : otherwise         (expected gamma ~ O(1) or instant kills)
Also records per-pair b_min at L = 2k+2 for the resonance-budget law:
#{fibers with b_min <= B} <= B*(2L+1), and the weak set lambda_bar <= 2.

Cells (committed): T=7 (1K,4U) p=17,29 top-10 census fams;
T=8-critical (2K,4U) p=17 top-10; T=8 zoo (1K,5U) p=17 top-10 +
random-40 distinct (seed 20261006 -- same as T-18).

Output: scripts/out_h2_pair_check.json + stdout log.
"""
import ast
import json
import random
import sys
import time
from collections import defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import (ball_mask, profile_sizes, size_tuples,
                               census_m4, k_of)
from higher_rung_step2 import sigma0, fiber_data
from scaling_closure_verify import (solutions_m5, fiber_data_m5,
                                    m4_sizedicts, sign_patterns_for)

OUT = {}


def b_min_of(lam, p, L):
    """min{b >= 1 : ||lam*b||_p <= L}; p if none (b <= p-1 scanned)."""
    x = 0
    for b in range(1, p):
        x += lam
        if x >= p:
            x -= p
        r = min(x, p - x)
        if r <= L:
            return b
    return p


def run_pairs(p, T, res_tuple, sizedict, U):
    """Anchor at first admissible fiber; per-second-fiber records."""
    fdat = fiber_data if len(res_tuple) == 3 else fiber_data_m5
    fibers = []
    for j in U:
        sz, cc = fdat(p, T, res_tuple, j)
        if sz in sizedict:
            fibers.append((j, sz, cc))
    if not fibers:
        return None
    j0, sz0, c0 = fibers[0]
    S0 = sizedict[sz0]
    inv_j0 = pow(j0, -1, p)
    recs = []
    for (j, sz, cc) in fibers[1:]:
        Solj = sizedict[sz]
        l = (j * inv_j0) % p
        tau = tuple((cc[i] - l * c0[i]) % p for i in range(len(c0)))
        surv = 0
        for a in S0:
            t = tuple((l * a[i] + tau[i]) % p for i in range(len(a)))
            if t in Solj:
                surv += 1
        lb = l if l <= p // 2 else p - l
        k = k_of(p, T)
        recs.append({'j': j, 'lambda': l, 'lambda_bar': lb,
                     'b_min': b_min_of(l, p, 2 * k + 2),
                     'n_sol0': len(S0), 'n_solj': len(Solj),
                     'survivors': surv})
    # cascade depth (same protocol as T-18)
    cands = {tuple((c0[i] - a[i]) * inv_j0 % p
                   for i in range(len(a))) for a in S0}
    alive = list(cands)
    depth = len(fibers)
    d = 1
    for (j, sz, cc) in fibers[1:]:
        Solj = sizedict[sz]
        alive = [mu for mu in alive
                 if tuple((cc[i] - mu[i] * j) % p
                          for i in range(len(mu))) in Solj]
        d += 1
        if not alive:
            depth = d
            break
    return {'depth': depth, 'n_fibers_admissible': len(fibers), 'pairs': recs}


def zone_of(lb, p):
    if lb <= 4:
        return 'small'
    if lb >= p - 4:
        return 'large'
    return 'mid'


def analyze(p, T, fam_list, m, U, sizedicts, label=''):
    print('  --- %s (T=%d, p=%d, %d fams, m=%d) ---' %
          (label, T, p, len(fam_list), m))
    sys.stdout.flush()
    t0 = time.time()
    pairs = []
    depths = []
    for fam in fam_list:
        sizedict = sizedicts[fam]
        for pat in sign_patterns_for(p, fam)[:1]:   # canonical sign
            r = run_pairs(p, T, pat, sizedict, U)
            if r is None:
                continue
            depths.append(r['depth'])
            n_coord = m - 1
            for pr in r['pairs']:
                rand = pr['n_sol0'] * pr['n_solj'] / float(p) ** n_coord
                g = (pr['survivors'] / rand) if rand > 0 else \
                    (float('inf') if pr['survivors'] else 0.0)
                z = zone_of(pr['lambda_bar'], p)
                pairs.append({'family': str(fam), 'zone': z,
                              'gamma': (g if g < 1e9 else None),
                              'gamma_raw': g,
                              'surv_frac': pr['survivors'] / float(pr['n_sol0'])
                              if pr['n_sol0'] else 0.0,
                              **{kk: pr[kk] for kk in
                                 ('j', 'lambda', 'lambda_bar', 'b_min',
                                  'n_sol0', 'n_solj', 'survivors')}})
    by_zone = defaultdict(list)
    for pr in pairs:
        by_zone[pr['zone']].append(pr['gamma_raw'] if pr['gamma'] is not None
                                   else 1e9)
    out = {'label': label, 'T': T, 'p': p, 'm': m, 'k': k_of(p, T),
           'n_families': len(fam_list),
           'depths': sorted(depths),
           'n_pairs': len(pairs),
           'zone_stats': {}, 'pairs': pairs}
    for z, gs in by_zone.items():
        gs = sorted(gs)
        n = len(gs)
        fin = [g for g in gs if g < 1e9]
        out['zone_stats'][z] = {
            'n': n,
            'gamma_median': gs[n // 2],
            'gamma_p90': gs[int(n * 0.9)],
            'gamma_max': gs[-1],
            'frac_zero': sum(1 for g in gs if g == 0) / float(n),
            'frac_inf': sum(1 for g in gs if g >= 1e9) / float(n),
            'gamma_max_finite': (fin[-1] if fin else None),
        }
    # resonance budget: fibers with small b_min among ALL pair lambdas
    lams = sorted({pr['lambda'] for pr in pairs})
    L = 2 * k_of(p, T) + 2
    for B in (1, 2, 4, 8):
        cnt = sum(1 for l in lams if b_min_of(l, p, L) <= B)
        out['res_set_b%d' % B] = {'count': cnt, 'bound': min(p - 1, B * (2 * L + 1))}
    out['elapsed_s'] = round(time.time() - t0, 1)
    print('    depths: %s | pairs: %d | zones: %s' %
          (out['depths'], len(pairs),
           {z: (round(v['gamma_median'], 2), round(v['gamma_p90'], 2),
                round(v['gamma_max'], 1), v['n'])
            for z, v in out['zone_stats'].items()}))
    sys.stdout.flush()
    return out


def main():
    print('=== H2 zone-resolved pair check (committed cells) ===')
    # ---------------- T=7 (1K,4U): p = 17, 29 ----------------
    for p in (17, 29):
        T = 7
        bm = ball_mask(p, k_of(p, T))
        U = [j for j in range(p) if not (bm >> j) & 1]
        fams_c, _ = census_m4(p, T, verbose=False)
        from collections import Counter
        agg = Counter()
        for (s1, s2, s3, s4, d2, d3, d4), st in fams_c.items():
            agg[(d2, d3, d4)] += st['n']
        top = [f for f, _ in agg.most_common(10)]
        sds = m4_sizedicts(p, T, top)
        OUT['T7_p%d' % p] = analyze(p, T, top, 4, U, sds,
                                    label='T=7 (1K,4U) top-10')
    # ---------------- T=8-critical (2K,4U): p = 17 ----------------
    p, T = 17, 8
    bm = ball_mask(p, k_of(p, T))
    U = [j for j in range(p) if not (bm >> j) & 1]
    fams_c, _ = census_m4(p, T, verbose=False)
    from collections import Counter
    agg = Counter()
    for (s1, s2, s3, s4, d2, d3, d4), st in fams_c.items():
        agg[(d2, d3, d4)] += st['n']
    top = [f for f, _ in agg.most_common(10)]
    sds = m4_sizedicts(p, T, top)
    OUT['T8c_p17'] = analyze(p, T, top, 4, U, sds, label='T=8-critical top-10')
    # ---------------- T=8 zoo (1K,5U): p = 17 ----------------
    print('  --- building m=5 solutions for selected zoo families ---')
    sys.stdout.flush()
    with open('/home/z/my-project/scripts/out_higher_rung_step1.json') as f:
        step1 = json.load(f)
    zoo = step1['T8_m5_census']['p17']
    top_zoo = [tuple(fl) for fl, n in zoo['top_families'][:10]]
    rng = random.Random(20261006)
    distinct = [tuple(ast.literal_eval(s)) for s in zoo['families'].keys()]
    distinct = [f for f in distinct if len(set(f)) == 4 and 1 not in f]
    sample = rng.sample(distinct, 40)
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(5, on, off, p)
    sols_by_fam = {}
    for fam in top_zoo + sample:
        keys = [szt + fam for szt in szts]
        sols = solutions_m5(p, T, keys)
        sizedict = defaultdict(set)
        for key, lst in sols.items():
            sizedict[key[:5]].update(lst)
        sols_by_fam[fam] = dict(sizedict)
    OUT['T8zoo_p17_top'] = analyze(p, T, top_zoo, 5, U, sols_by_fam,
                                   label='zoo top-10')
    OUT['T8zoo_p17_sample'] = analyze(p, T, sample, 5, U, sols_by_fam,
                                      label='zoo random-40 distinct')
    with open('/home/z/my-project/scripts/out_h2_pair_check.json', 'w') as f:
        json.dump(OUT, f, indent=1)
    print('DONE; wrote scripts/out_h2_pair_check.json')


if __name__ == '__main__':
    main()
