#!/usr/bin/env python3
"""
Scaling-closure lemma attempt -- targeted verification of the proof's
load-bearing structural claims, on COMMITTED cells and primes only.

The proof attempt (download/scaling_closure_lemma_attempt.md) derives the
cascade as a bounded-junction dilate-intersection and needs:

  [J]  the junction-volume law: |Sol(family, size)| <= (R(T) k)^(m-2)
       (m = number of footprints; m-2 = T-5 at the frontier cells);
  [C]  the two-fiber dilate-cut law: the anchor solutions surviving the
       second fiber number
           surv = |{a in Sol0 : lambda*a + tau in Sol_j}|,
       tau = c0(j) - lambda*c0(j0),  lambda = j/j0 (mod p),
       and the INDEPENDENT-EXTENSION model predicts surv ~ gamma *
       |Sol0|*|Sol_j| / p^(m-1) with gamma = O(1) (the correlation
       factor -- the proof's deepest input);
  [R]  the resonance structure: weak cuts concentrate at dilation ratios
       lambda whose cyclic representative is small / near p / a small
       rational (the b_min resonance of Lemma D);
  [D]  the cascade depth per family ~ the formula
       (T-4) log(Tk) / (log(Tk) + (T-5) log(T/R)) + tail.

Verification scope (committed):
  - T=7 (1K,4U)  at p = 17, 29   (committed census cells, committed code);
  - T=8 (2K,4U)  at p = 17       (committed critical census);
  - T=8 (1K,5U)  at p = 17       (the ZOO cell -- committed census counts;
    solutions for selected keys rebuilt with the census_m4_solutions
    pattern, cross-checked against the committed counts).

NO new primes, NO T=9/T=10/T=11, NO cell-emptiness claims beyond what
the committed GT already established.  This is lemma verification under
persistent rule (iii).

Output: scripts/out_scaling_closure_verify.json + stdout log.
"""
import ast
import json
import random
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import (value_arc, ball_mask, profile_sizes,
                               ap_base, rot, fit_placements, size_tuples,
                               census_m4, census_m4_solutions, k_of)
from higher_rung_step2 import sigma0, fiber_data

OUT = {}


def inv_table(p):
    return [0] + [pow(d, -1, p) for d in range(1, p)]


# ----------------------------------------------------------------------
# m=5 selected-key solutions (the census_m4_solutions pattern, one AP up)
# ----------------------------------------------------------------------

def solutions_m5(p, T, keys):
    """{key: [(a2,a3,a4,a5), ...]} for selected
    (s1,s2,s3,s4,s5,d2,d3,d4,d5) keys.  Pruning mirrors census_m5."""
    full = (1 << p) - 1
    inv = inv_table(p)
    want = set(keys)
    got = defaultdict(list)
    for key in sorted(want):
        (s1, s2, s3, s4, s5, d2, d3, d4, d5) = key
        ov = s1 + s2 + s3 + s4 + s5 - p
        m1 = (1 << s1) - 1
        b2 = ap_base(d2, s2, p)
        b3 = ap_base(d3, s3, p)
        b4 = ap_base(d4, s4, p)
        for a2 in range(p):
            m2 = rot(b2, a2, p, full)
            if (m1 & m2).bit_count() > ov:
                continue
            m12 = m1 | m2
            if (full & ~m12).bit_count() > s3 + s4 + s5:
                continue
            for a3 in range(p):
                m3 = rot(b3, a3, p, full)
                m123 = m12 | m3
                if m123.bit_count() < p - s4 - s5:
                    continue
                if (m12 & m3).bit_count() > ov:
                    continue
                for a4 in range(p):
                    m4 = rot(b4, a4, p, full)
                    m1234 = m123 | m4
                    if m1234.bit_count() < p - s5:
                        continue
                    if (m123 & m4).bit_count() > ov:
                        continue
                    R = full & ~m1234
                    nr = R.bit_count()
                    if nr == 0:
                        for a5 in range(p):
                            got[key].append((a2, a3, a4, a5))
                        continue
                    if nr > s5:
                        continue
                    Rl = [i for i in range(p) if (R >> i) & 1]
                    for a5 in fit_placements(Rl, p, s5, d5, inv[d5]):
                        got[key].append((a2, a3, a4, a5))
    return dict(got)


# ----------------------------------------------------------------------
# m=5 fiber data (the fiber_data pattern, one unit up)
# ----------------------------------------------------------------------

def fiber_data_m5(p, T, res_tuple, j):
    """For units (1, r2, r3, r4, r5) at fiber j: sizes (5-tuple) and
    canonical relative starts (4-tuple)."""
    k = k_of(p, T)
    eps = p % T
    s1 = sigma0(p, T, 1, j)
    out_s, out_c = [], []
    for r in res_tuple:
        ubar = pow(r, -1, p)
        d = min(ubar, p - ubar)
        neg = (ubar != d)
        s = sigma0(p, T, r, j)
        a0 = (r * j) % p
        _, L = value_arc(p, T, k, eps, a0)
        c = (s - s1) % p
        if neg:
            c = (c - (L - 1) * d) % p
        out_s.append(L)
        out_c.append(c)
    _, L1 = value_arc(p, T, k, eps, j % p)
    return (L1,) + tuple(out_s), tuple(out_c)


# ----------------------------------------------------------------------
# the cascade + the two-fiber dilate analysis
# ----------------------------------------------------------------------

def run_cascade(p, T, res_tuple, sizedict, U, pair_analysis=False):
    """Anchor at the first admissible fiber; cascade through U.

    sizedict: {size-tuple: set of a-tuples}.
    Returns depth, |Sol0|, and (if pair_analysis) per-second-fiber
    survivor counts + lambda metadata."""
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
    cands = {tuple((c0[i] - a[i]) * inv_j0 % p
                   for i in range(len(a))) for a in S0}
    alive = list(cands)
    depth = len(fibers)      # survived all admissible fibers
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
    res = {'depth': depth, 'n_anchor_cands': len(cands),
           'n_fibers_admissible': len(fibers)}
    if pair_analysis:
        pairs = []
        for (j, sz, cc) in fibers[1:]:
            Solj = sizedict[sz]
            l = (j * inv_j0) % p
            surv = 0
            for a in S0:
                # a survives fiber j iff lambda*a + tau in Solj with
                # tau = cc - lambda*c0
                t = tuple((l * a[i] + cc[i] - l * c0[i]) % p
                          for i in range(len(a)))
                if t in Solj:
                    surv += 1
            lb = l if l <= p // 2 else p - l
            pairs.append({'j': j, 'lambda': l, 'lambda_bar': lb,
                          'n_sol0': len(S0), 'n_solj': len(Solj),
                          'survivors': surv})
        res['pairs'] = pairs
    return res


def sign_patterns_for(p, fam):
    """All residue tuples whose canonical differences are fam
    (2^(m-1) patterns; d=1 entries give r in {1, p-1})."""
    inv = inv_table(p)
    pats = []
    for code in range(1 << len(fam)):
        pat = []
        for i, d in enumerate(fam):
            base = inv[d] if d > 1 else 1
            pat.append(base if (code >> i) & 1 else p - base)
        pats.append(tuple(pat))
    return pats


def analyze(p, T, fam_list, m, U, sizedicts, sign_mode='canonical',
            pair_analysis=False, label=''):
    """sizedicts: {family: {size-tuple: set of a-tuples}}."""
    print('  --- %s (T=%d, p=%d, %d families, m=%d, signs=%s) ---'
          % (label, T, p, len(fam_list), m, sign_mode))
    sys.stdout.flush()
    depths = []
    gammas = []
    resonance = defaultdict(list)
    sol_sizes = []
    t0 = time.time()
    per_fam = {}
    for fam in fam_list:
        sizedict = sizedicts[fam]
        if sign_mode == 'canonical':
            patterns = [sign_patterns_for(p, fam)[0]]
        else:
            patterns = sign_patterns_for(p, fam)
        fam_depths = []
        for pat in patterns:
            r = run_cascade(p, T, pat, sizedict, U,
                            pair_analysis=pair_analysis)
            if r is None:
                continue
            fam_depths.append(r['depth'])
            depths.append(r['depth'])
            if pair_analysis:
                for pr in r.get('pairs', []):
                    n_coord = m - 1
                    rand = pr['n_sol0'] * pr['n_solj'] / float(p) ** n_coord
                    g = pr['survivors'] / rand if rand > 0 else \
                        (float('inf') if pr['survivors'] else 0.0)
                    if g < 1e9:
                        gammas.append(g)
                    zone = ('small' if pr['lambda_bar'] <= 4 else
                            'large' if pr['lambda_bar'] >= p - 4 else 'mid')
                    frac = (pr['survivors'] / float(pr['n_sol0'])
                            if pr['n_sol0'] else 0.0)
                    resonance[zone].append(frac)
        for S in sizedict.values():
            sol_sizes.append(len(S))
        per_fam[str(fam)] = {'depths': fam_depths,
                             'n_sizes': len(sizedict)}
    depths.sort()
    gammas.sort()
    sol_sizes.sort()
    out = {
        'label': label, 'T': T, 'p': p, 'm': m, 'k': k_of(p, T),
        'n_families': len(fam_list), 'sign_mode': sign_mode,
        'depths': depths,
        'depth_median': depths[len(depths) // 2] if depths else None,
        'depth_max': depths[-1] if depths else None,
        'gamma_median': gammas[len(gammas) // 2] if gammas else None,
        'gamma_p90': gammas[int(len(gammas) * 0.9)] if gammas else None,
        'gamma_max': gammas[-1] if gammas else None,
        'gamma_min': gammas[0] if gammas else None,
        'n_gamma': len(gammas),
        'sol_max': sol_sizes[-1] if sol_sizes else None,
        'sol_median': sol_sizes[len(sol_sizes) // 2] if sol_sizes else None,
        'sol_n': len(sol_sizes),
        'resonance_mean_survival': {z: (sum(v) / len(v) if v else None)
                                    for z, v in resonance.items()},
        'resonance_n': {z: len(v) for z, v in resonance.items()},
        'per_fam': per_fam,
        'elapsed_s': round(time.time() - t0, 1),
    }
    print('    depths: median %s max %s (n=%d) | gamma: median %.2f '
          'p90 %.2f max %.2f (n=%d) | |Sol|: median %s max %s (n=%d) | %.1fs'
          % (out['depth_median'], out['depth_max'], len(depths),
             out['gamma_median'] if out['gamma_median'] is not None else -1,
             out['gamma_p90'] if out['gamma_p90'] is not None else -1,
             out['gamma_max'] if out['gamma_max'] is not None else -1,
             len(gammas), out['sol_median'], out['sol_max'],
             len(sol_sizes), out['elapsed_s']))
    if out['resonance_mean_survival']:
        print('    mean survival fraction by lambda-zone: %s'
              % {z: (round(v, 4) if v is not None else None)
                 for z, v in out['resonance_mean_survival'].items()})
    sys.stdout.flush()
    return out


# ----------------------------------------------------------------------

def m4_sizedicts(p, T, fams):
    """{family: {size-tuple: set of solutions}} from the committed
    census code (census_m4 + census_m4_solutions)."""
    fams_c, _ = census_m4(p, T, verbose=False)
    keys_by_fam = defaultdict(list)
    for key in fams_c:
        keys_by_fam[(key[4], key[5], key[6])].append(key)
    sols = census_m4_solutions(p, T, [k for f in fams for k in
                                       keys_by_fam.get(f, [])])
    out = defaultdict(lambda: defaultdict(set))
    for key, lst in sols.items():
        out[(key[4], key[5], key[6])][key[:4]].update(lst)
    return {f: dict(out[f]) for f in fams}


def main():
    print('=== scaling-closure verification (committed cells/primes) ===')
    # ---------------- T=7 (1K,4U): p = 17, 29 ----------------
    for p in (17, 29):
        T = 7
        k = k_of(p, T)
        bm = ball_mask(p, k)
        U = [j for j in range(p) if not (bm >> j) & 1]
        fams_c, _ = census_m4(p, T, verbose=False)
        agg = Counter()
        for (s1, s2, s3, s4, d2, d3, d4), st in fams_c.items():
            agg[(d2, d3, d4)] += st['n']
        top = [f for f, _ in agg.most_common(10)]
        sds = m4_sizedicts(p, T, top)
        OUT['T7_p%d' % p] = analyze(p, T, top, 4, U, sds,
                                    pair_analysis=True,
                                    label='T=7 (1K,4U) top-10 census fams')
        OUT['T7_p%d_signs' % p] = analyze(p, T, top[:3], 4, U, sds,
                                          sign_mode='all',
                                          label='T=7 top-3, all 8 signs')
    # ---------------- T=8-critical (2K,4U): p = 17 ----------------
    p, T = 17, 8
    k = k_of(p, T)
    bm = ball_mask(p, k)
    U = [j for j in range(p) if not (bm >> j) & 1]
    fams_c, _ = census_m4(p, T, verbose=False)
    agg = Counter()
    for (s1, s2, s3, s4, d2, d3, d4), st in fams_c.items():
        agg[(d2, d3, d4)] += st['n']
    top = [f for f, _ in agg.most_common(10)]
    sds = m4_sizedicts(p, T, top)
    OUT['T8c_p17'] = analyze(p, T, top, 4, U, sds, pair_analysis=True,
                             label='T=8-critical (2K,4U) top-10')
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
    fam_list = top_zoo + sample
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(5, on, off, p)
    t0 = time.time()
    sols_by_fam = {}
    for fam in fam_list:
        keys = [szt + fam for szt in szts]
        sols = solutions_m5(p, T, keys)
        sizedict = defaultdict(set)
        for key, lst in sols.items():
            sizedict[key[:5]].update(lst)
        sols_by_fam[fam] = dict(sizedict)
    print('    m=5 solutions built in %.1fs (%d families, %d size tuples)'
          % (time.time() - t0, len(sols_by_fam), len(szts)))
    # self-check vs committed counts for the top zoo family
    chk = sum(len(v) for v in sols_by_fam[top_zoo[0]].values())
    committed = next(n for fl, n in zoo['top_families']
                     if tuple(fl) == top_zoo[0])
    print('    self-check %s: rebuilt %d vs committed %d '
          '(difference = the nr==0 bystander placements: census_m5 '
          'EXCLUDES first-4-already-cover + free-a5, the covering '
          'condition INCLUDES them; delta = %d)'
          % (top_zoo[0], chk, committed, chk - committed))
    OUT['zoo_selfcheck'] = {'family': str(top_zoo[0]),
                            'rebuilt': chk, 'committed': committed,
                            'note': 'census_m5 excludes nr==0 bystander '
                                    'placements; solutions include them '
                                    '(covering condition, m=4 GT '
                                    'convention)'}
    OUT['T8zoo_p17_top'] = analyze(p, T, top_zoo, 5, U, sols_by_fam,
                                   pair_analysis=True,
                                   label='zoo top-10')
    OUT['T8zoo_p17_sample'] = analyze(p, T, sample, 5, U, sols_by_fam,
                                      pair_analysis=True,
                                      label='zoo random-40 distinct')
    OUT['T8zoo_p17_signs'] = analyze(p, T, top_zoo[:2], 5, U, sols_by_fam,
                                     sign_mode='all',
                                     label='zoo top-2, all 16 signs')
    with open('/home/z/my-project/scripts/out_scaling_closure_verify.json',
              'w') as f:
        json.dump(OUT, f, indent=1)
    print('DONE; wrote scripts/out_scaling_closure_verify.json')


if __name__ == '__main__':
    main()
