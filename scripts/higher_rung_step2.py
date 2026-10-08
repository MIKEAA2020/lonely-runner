#!/usr/bin/env python3
"""
Higher-rung step 2: the shear side at T=7/T=8.

(A) DRIFT LAW (general T, proved + verified).  For a unit v = r + p m
    (r = v mod p, m the lift): the footprint start on fiber j obeys
        sigma_v(j) = sigma_r^0(j) - ubar(r) * m * j  (mod p),
    with sigma_r^0(j) = ubar*(S2(a0) - floor(rj/p)), a0 = rj mod p.
    For a colliding pair w = u + p m (shared residue) the RELATIVE drift
    is exactly linear: sigma_w(j) - sigma_u(j) = -ubar*m*j.  [S2-GEN]

    CANONICALIZATION: an AP with difference -d equals the AP with
    difference +d started (s-1)*d earlier; the census uses canonical
    d in [1,(p-1)/2], so a fiber footprint with ubar_i = -d_i converts
    to the canonical start sigma_i - (s_i - 1) d_i.

(B) GROUND TRUTH (the decisive experiment): does a (1K,4U) family
    {B_{pa'}, B_1, B_{v2}, B_{v3}, B_{v4}} cover Z_{p^2} at T=7?
      - brute force (all unit triples, distinct bad sets = v_i =/= +-v_j
        mod p^2, INCLUDING same-residue colliding pairs) at p in {5,7,11};
      - census-restricted enumeration (residue triples admissible by the
        fiber census; lifts via the drift-vector candidates fixed by the
        first fiber) at p in {5,7,11,13,17}.
    [GT], [GT-BF]

(C) EXIT TIME: for the top census families, the maximum number of
    consecutive good fibers over all drift vectors mu (the analogue of
    T=6's "max 6 good rows").  [EXIT]

Output: scripts/out_higher_rung_step2.json + stdout log.
"""
import json
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import (value_arc, profile_sizes, ball_mask,
                               census_m4, census_m4_solutions)

OUT = {}


def k_eps(p, T):
    return (p - p % T) // T, p % T


def arc_start(V, p):
    """The cyclic-AP start (step 1) of the value arc V (a set)."""
    v = sorted(V)
    m = len(v)
    if m == 1:
        return v[0]
    gaps = [(v[(i + 1) % m] - v[i]) % p for i in range(m)]
    gmax = max(gaps)
    gi = gaps.index(gmax)
    return v[(gi + 1) % m]


def sigma0(p, T, r, j):
    """Residue-part start of the footprint of residue r at fiber j."""
    ubar = pow(r, -1, p)
    a0 = (r * j) % p
    a1 = (r * j) // p
    V, _ = value_arc(p, T, *k_eps(p, T), a0)
    return (ubar * (arc_start(V, p) - a1)) % p


def sigma_full(p, T, r, m, j):
    """Start of the footprint of unit v = r + p m at fiber j
    (drift law form)."""
    ubar = pow(r, -1, p)
    return (sigma0(p, T, r, j) - ubar * m * j) % p


# ----------------------------------------------------------------------
# (A) drift law
# ----------------------------------------------------------------------

def drift_law_check():
    res = {'checks': 0, 'bad': 0}
    for T, plist in ((7, (17, 29, 43)), (8, (17, 41))):
        for p in plist:
            for r in (2, 3, 5, 7):
                if r % p == 0:
                    continue
                ubar = pow(r, -1, p)
                for m in range(0, min(p, 6)):
                    v = r + p * m
                    for j in range(p):
                        # direct: a1 from the full unit
                        a0 = (v * j) % p
                        a1 = (v * j) // p
                        V, _ = value_arc(p, T, *k_eps(p, T), a0)
                        s_dir = (ubar * (arc_start(V, p) - a1)) % p
                        s_law = sigma_full(p, T, r, m, j)
                        res['checks'] += 1
                        if s_dir != s_law:
                            res['bad'] += 1
                            if res['bad'] < 5:
                                print('DRIFT FAIL T=%d p=%d r=%d m=%d j=%d'
                                      % (T, p, r, m, j))
    print('[S2-GEN] drift law: %d/%d pass' %
          (res['checks'] - res['bad'], res['checks']))
    return res


# ----------------------------------------------------------------------
# census solution sets + fiber data
# ----------------------------------------------------------------------

def census_solution_sets(p, T):
    """{(d2,d3,d4) -> {(s1,s2,s3,s4) -> set of (a2,a3,a4)}}."""
    fams, meta = census_m4(p, T, verbose=False)
    keys = list(fams.keys())
    sols = census_m4_solutions(p, T, keys)
    by_fam = defaultdict(dict)
    for (s1, s2, s3, s4, d2, d3, d4), lst in sols.items():
        d = by_fam.setdefault((d2, d3, d4), {})
        d.setdefault((s1, s2, s3, s4), set()).update(lst)
    return by_fam, meta


def fiber_data(p, T, res_tuple, j):
    """For units (1, r2, r3, r4) at fiber j: the canonical relative starts
    (c2, c3, c4) (census-normalized, sign-converted) and the sizes."""
    k, eps = k_eps(p, T)
    s1 = sigma0(p, T, 1, j)
    out_s = []
    out_c = []
    for r in res_tuple:
        ubar = pow(r, -1, p)
        d = min(ubar, p - ubar)          # canonical difference
        neg = (ubar != d)                # footprint runs with -d
        s = sigma0(p, T, r, j)
        a0 = (r * j) % p
        _, L = value_arc(p, T, k, eps, a0)
        c = (s - s1) % p
        if neg:
            c = (c - (L - 1) * d) % p
        out_s.append(L)
        out_c.append(c)
    a0 = (1 * j) % p
    _, L1 = value_arc(p, T, k, eps, a0)
    return (L1,) + tuple(out_s), tuple(out_c)


# ----------------------------------------------------------------------
# (B) ground truth
# ----------------------------------------------------------------------

def ground_truth_bruteforce(p, T=7, verbose=True):
    """Direct enumeration of ALL normalized families (distinct bad sets,
    including colliding same-residue pairs).  Small p only.
    Ball radius k_eff = floor(c/p) (robust at the degenerate p=T)."""
    N = p * p
    c = (N - 1) // T
    k_eff = c // p
    units = [v for v in range(1, N) if v % p and v <= N - v]
    foot = {}
    for v in units:
        tab = []
        for j in range(p):
            m = 0
            for t in range(p):
                x = j + p * t
                wx = (v * x) % N
                if min(wx, N - wx) <= c:
                    m |= 1 << t
            tab.append(m)
        foot[v] = tab
    bm = ball_mask(p, k_eff)
    full = (1 << p) - 1
    n_cov = 0
    n_chk = 0
    t0 = time.time()
    for ap in range(1, p):
        ainv = pow(ap, -1, p)
        U = [ainv * j % p for j in range(p) if not (bm >> j) & 1]
        nu = len(units)
        for i2 in range(nu):
            v2 = units[i2]
            f2 = foot[v2]
            for i3 in range(i2 + 1, nu):
                v3 = units[i3]
                f3 = foot[v3]
                for i4 in range(i3 + 1, nu):
                    v4 = units[i4]
                    # distinct bad sets: v_i =/= +- v_j mod N
                    if (v3 - v2) % N in (0, N - 2 * v2 + 0) or \
                       (v4 - v2) % N in (0, N - 2 * v2) or \
                       (v4 - v3) % N in (0, N - 2 * v3):
                        continue
                    # quicker explicit check:
                    if any((b - a) % N == 0 or (b + a) % N == 0
                           for (a, b) in ((v2, v3), (v2, v4), (v3, v4))):
                        continue
                    n_chk += 1
                    good = True
                    for j in U:
                        if foot[1][j] | f2[j] | f3[j] | foot[v4][j] != full:
                            good = False
                            break
                    if good:
                        n_cov += 1
                        print('  [GT-BF] COVERING p=%d a\'=%d v=(1,%d,%d,%d)'
                              % (p, ap, v2, v3, v4))
                        sys.stdout.flush()
    if verbose:
        print('  [GT-BF] p=%d: %d triples checked, %d coverings, %.1fs' %
              (p, n_chk, n_cov, time.time() - t0))
        sys.stdout.flush()
    return {'p': p, 'n_checked': n_chk, 'n_coverings': n_cov}


def ground_truth_census_restricted(p, T=7, verbose=True):
    """Census-restricted (1K,4U) enumeration via drift candidates."""
    by_fam, meta = census_solution_sets(p, T)
    k, eps = k_eps(p, T)
    bm = ball_mask(p, k)
    U_all = [j for j in range(p) if not (bm >> j) & 1]
    # admissible residue triples (with the family they certify)
    res_options = {}
    for (d2, d3, d4) in by_fam:
        i2 = pow(d2, -1, p)
        i3 = pow(d3, -1, p)
        i4 = pow(d4, -1, p)
        for s2 in (i2, p - i2):
            for s3 in (i3, p - i3):
                for s4 in (i4, p - i4):
                    res_options[(s2 % p, s3 % p, s4 % p)] = (d2, d3, d4)
    if verbose:
        print('  [GT] p=%d k=%d: %d census families, %d residue triples' %
              (p, k, len(by_fam), len(res_options)))
        sys.stdout.flush()
    n_cov = 0
    n_combo = 0
    best = 0
    best_desc = None
    cov_seen = set()
    t0 = time.time()
    for ap in range(1, p):
        ainv = pow(ap, -1, p)
        U = [ainv * j % p for j in U_all]
        if not U:
            continue
        for (r2, r3, r4), fam in res_options.items():
            sizedict = by_fam[fam]
            fibers = []
            ok = True
            for j in U:
                sz, cc = fiber_data(p, T, (r2, r3, r4), j)
                if sz not in sizedict:
                    ok = False
                    break
                fibers.append((j, sz, cc))
            if not ok:
                continue
            n_combo += 1
            j0, sz0, (b2, b3, b4) = fibers[0]
            S0 = sizedict[sz0]
            inv_j0 = pow(j0, -1, p)
            cands = set()
            for (a2, a3, a4) in S0:
                cands.add((((b2 - a2) * inv_j0) % p,
                           ((b3 - a3) * inv_j0) % p,
                           ((b4 - a4) * inv_j0) % p))
            for (mu2, mu3, mu4) in cands:
                ng = 1
                for (j, sz, (c2, c3, c4)) in fibers[1:]:
                    if ((c2 - mu2 * j) % p, (c3 - mu3 * j) % p,
                            (c4 - mu4 * j) % p) in sizedict[sz]:
                        ng += 1
                    else:
                        break
                if ng > best:
                    best = ng
                    best_desc = (ap, r2, r3, r4, mu2, mu3, mu4)
                if ng == len(fibers):
                    # reconstruct the units and enforce distinct bad sets
                    m2 = (mu2 * r2) % p
                    m3 = (mu3 * r3) % p
                    m4 = (mu4 * r4) % p
                    v2 = r2 + p * m2
                    v3 = r3 + p * m3
                    v4 = r4 + p * m4
                    N = p * p
                    vs = (v2, v3, v4)
                    ok_dist = True
                    for ia in range(3):
                        for ib in range(ia + 1, 3):
                            if (vs[ib] - vs[ia]) % N == 0 or \
                               (vs[ib] + vs[ia]) % N == 0:
                                ok_dist = False
                    for v in vs:
                        if v % N in (1, N - 1):
                            ok_dist = False
                    if ok_dist:
                        sig = (ap,) + tuple(sorted(
                            min(v, N - v) for v in vs))
                        if sig not in cov_seen:
                            cov_seen.add(sig)
                            n_cov += 1
                            print('  [GT] COVERING p=%d a\'=%d '
                                  'r=(%d,%d,%d) mu=(%d,%d,%d) v=(%d,%d,%d)'
                                  % (p, ap, r2, r3, r4, mu2, mu3, mu4,
                                     v2, v3, v4))
                            sys.stdout.flush()
    if verbose:
        print('  [GT] p=%d: %d combos, %d coverings, best exit %d/%d (%s), '
              '%.1fs' % (p, n_combo, n_cov, best, len(U_all), best_desc,
                         time.time() - t0))
        sys.stdout.flush()
    return {'p': p, 'n_combos': n_combo, 'n_coverings': n_cov,
            'best_exit': best, 'n_U': len(U_all), 'best_desc': best_desc}


# ----------------------------------------------------------------------
# (C) exit time
# ----------------------------------------------------------------------

def exit_time_top(p, T=7, topn=10, verbose=True):
    fams, meta = census_m4(p, T, verbose=False)
    fam_agg = Counter()
    for (s1, s2, s3, s4, d2, d3, d4), st in fams.items():
        fam_agg[(d2, d3, d4)] += st['n']
    top = [f for f, _ in fam_agg.most_common(topn)]
    keys = [key for key in fams if (key[4], key[5], key[6]) in set(top)]
    sols = census_m4_solutions(p, T, keys)
    by_fam = defaultdict(dict)
    for (s1, s2, s3, s4, d2, d3, d4), lst in sols.items():
        d = by_fam.setdefault((d2, d3, d4), {})
        d.setdefault((s1, s2, s3, s4), set()).update(lst)
    k, eps = k_eps(p, T)
    bm = ball_mask(p, k)
    U_all = [j for j in range(p) if not (bm >> j) & 1]
    res = {}
    for fam in top:
        sizedict = by_fam.get(fam)
        if not sizedict:
            continue
        r2 = pow(fam[0], -1, p)
        r3 = pow(fam[1], -1, p)
        r4 = pow(fam[2], -1, p)
        fibers = []
        for j in U_all:
            sz, cc = fiber_data(p, T, (r2, r3, r4), j)
            if sz in sizedict:
                fibers.append((j, sz, cc))
        if len(fibers) < 2:
            res[str(fam)] = {'n_fibers': len(fibers),
                             'best': len(fibers)}
            continue
        j0, sz0, (b2, b3, b4) = fibers[0]
        inv_j0 = pow(j0, -1, p)
        cands = set()
        for (a2, a3, a4) in sizedict[sz0]:
            cands.add((((b2 - a2) * inv_j0) % p,
                       ((b3 - a3) * inv_j0) % p,
                       ((b4 - a4) * inv_j0) % p))
        best = 0
        for (mu2, mu3, mu4) in cands:
            ng = 1
            for (j, sz, (c2, c3, c4)) in fibers[1:]:
                if ((c2 - mu2 * j) % p, (c3 - mu3 * j) % p,
                        (c4 - mu4 * j) % p) in sizedict[sz]:
                    ng += 1
                else:
                    break
            best = max(best, ng)
        res[str(fam)] = {'n_fibers': len(fibers), 'best': best,
                         'n_cand_mu': len(cands)}
        if verbose:
            print('  [EXIT] p=%d fam=%s: %d/%d fibers, %d mu-cands' %
                  (p, fam, best, len(fibers), len(cands)))
            sys.stdout.flush()
    return res


# ----------------------------------------------------------------------

def main():
    t0 = time.time()
    print('=== higher_rung_step2: shear side at T=7/T=8 ===')
    OUT['drift_law'] = drift_law_check()
    print('--- [GT-BF] brute-force control ---')
    for p in (5, 7, 11):
        OUT['gt_bf_p%d' % p] = ground_truth_bruteforce(p, 7)
    print('--- [GT] census-restricted ground truth ---')
    for p in (5, 11, 13, 17):
        OUT['gt_p%d' % p] = ground_truth_census_restricted(p, 7)
    print('--- [EXIT] exit-time at top families ---')
    for p in (17, 29):
        OUT['exit_p%d' % p] = exit_time_top(p, 7)
    OUT['elapsed_s'] = time.time() - t0
    with open('/home/z/my-project/scripts/out_higher_rung_step2.json',
              'w') as f:
        json.dump(OUT, f, indent=1)
    print('DONE in %.1fs' % OUT['elapsed_s'])


if __name__ == '__main__':
    main()


# ----------------------------------------------------------------------
# (B') T=8 (2K,4U) critical cell ground truth
# ----------------------------------------------------------------------

def ground_truth_T8_2K4U(p, T=8, verbose=True):
    """The (2K,4U) cell at N=p^2, T=8: two kernel sets (a', a'' distinct
    mod +-p) + four units.  U = fibers uncovered by both kernels.  On each
    U-fiber the four unit footprints must cover Z_p (the T=8-critical
    census problem).  Census-restricted + distinctness-filtered."""
    by_fam, meta = census_solution_sets(p, T)
    k, eps = k_eps(p, T)
    bm = ball_mask(p, k)
    U_all = [j for j in range(p) if not (bm >> j) & 1]
    res_options = {}
    for (d2, d3, d4) in by_fam:
        i2 = pow(d2, -1, p)
        i3 = pow(d3, -1, p)
        i4 = pow(d4, -1, p)
        for s2 in (i2, p - i2):
            for s3 in (i3, p - i3):
                for s4 in (i4, p - i4):
                    res_options[(s2 % p, s3 % p, s4 % p)] = (d2, d3, d4)
    if verbose:
        print('  [GT8] p=%d k=%d: %d census families, %d residue triples' %
              (p, k, len(by_fam), len(res_options)))
        sys.stdout.flush()
    N = p * p
    n_cov = 0
    n_combo = 0
    best = 0
    best_desc = None
    cov_seen = set()
    t0 = time.time()
    # precompute per-residue-triple fiber data ONCE (a'=independent)
    fdata = {}
    for (r2, r3, r4), fam in res_options.items():
        sizedict = by_fam[fam]
        rows = []
        for j in range(p):
            sz, cc = fiber_data(p, T, (r2, r3, r4), j)
            if sz in sizedict:
                rows.append((j, sz, cc))
        fdata[(r2, r3, r4)] = (sizedict, rows)
    for ap in range(1, p):
        ainv = pow(ap, -1, p)
        Ua = {ainv * j % p for j in U_all}
        for app in range(ap + 1, p):
            if (app - ap) % p == 0 or (app + ap) % p == 0:
                continue  # kernel bad sets must be distinct
            ainvp = pow(app, -1, p)
            Ub = {ainvp * j % p for j in U_all}
            U = sorted(Ua & Ub)
            if not U:
                continue
            for (r2, r3, r4), (sizedict, rows) in fdata.items():
                rowmap = {j: (sz, cc) for (j, sz, cc) in rows}
                fibers = [(j, rowmap[j][0], rowmap[j][1]) for j in U
                          if j in rowmap]
                if len(fibers) < len(U):
                    continue
                n_combo += 1
                j0, sz0, (b2, b3, b4) = fibers[0]
                inv_j0 = pow(j0, -1, p)
                cands = set()
                for (a2, a3, a4) in sizedict[sz0]:
                    cands.add((((b2 - a2) * inv_j0) % p,
                               ((b3 - a3) * inv_j0) % p,
                               ((b4 - a4) * inv_j0) % p))
                for (mu2, mu3, mu4) in cands:
                    ng = 1
                    for (j, sz, (c2, c3, c4)) in fibers[1:]:
                        if ((c2 - mu2 * j) % p, (c3 - mu3 * j) % p,
                                (c4 - mu4 * j) % p) in sizedict[sz]:
                            ng += 1
                        else:
                            break
                    if ng > best:
                        best = ng
                        best_desc = (ap, app, r2, r3, r4, mu2, mu3, mu4)
                    if ng == len(fibers):
                        m2 = (mu2 * r2) % p
                        m3 = (mu3 * r3) % p
                        m4 = (mu4 * r4) % p
                        vs = (r2 + p * m2, r3 + p * m3, r4 + p * m4)
                        ok_dist = True
                        for ia in range(3):
                            for ib in range(ia + 1, 3):
                                if (vs[ib] - vs[ia]) % N == 0 or \
                                   (vs[ib] + vs[ia]) % N == 0:
                                    ok_dist = False
                        for v in vs:
                            if v % N in (1, N - 1):
                                ok_dist = False
                        if ok_dist:
                            sig = (ap, app) + tuple(sorted(
                                min(v, N - v) for v in vs))
                            if sig not in cov_seen:
                                cov_seen.add(sig)
                                n_cov += 1
                                print('  [GT8] COVERING p=%d a\'=(%d,%d) '
                                      'r=(%d,%d,%d) mu=(%d,%d,%d)' %
                                      (p, ap, app, r2, r3, r4, mu2, mu3,
                                       mu4))
                                sys.stdout.flush()
    if verbose:
        print('  [GT8] p=%d: %d combos, %d coverings, best exit %d/%d, '
              '%.1fs' % (p, n_combo, n_cov, best,
                         len(U_all), time.time() - t0))
        sys.stdout.flush()
    return {'p': p, 'n_combos': n_combo, 'n_coverings': n_cov,
            'best_exit': best, 'best_desc': best_desc}
