#!/usr/bin/env python3
"""
strand_rigidity_debug.py -- T-22 companion: (a) localize the Lemma V'
projection-invariance failure (which orbit members differ, and how);
(b) the fallback coordinate probe: marginals of the arc-normalized starts
b_i = a_i * d_i^{-1} mod p (the per-arc phase characters) and of the
pairwise differences a_j - a_i mod p, on the top-volume families at each
committed prime; (c) identify the true top orbit per prime.
"""
import json
import sys
from collections import Counter

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples
from scaling_closure_verify import solutions_m5

T, M = 8, 5


def canon(x, p):
    x %= p
    return min(x, p - x)


def renorm_orbit(S, p):
    D = sorted(set([1] + list(S)))
    out = []
    for delta in D:
        inv = pow(delta, -1, p)
        out.append(tuple(sorted(canon(d * inv, p) for d in D if d != delta)))
    return out


def enum(p, fam, szts):
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    out = {}
    for key in sorted(sols):
        out[key[:5]] = sols[key]
    return out


def marg(vals):
    return len(set(vals))


def phase_probe(p, fam, szt, sols):
    """Marginals in the fallback coordinates for one (family, szt)."""
    inv = [pow(d, -1, p) for d in fam]
    b = [[(a[i] * inv[i]) % p for a in sols] for i in range(4)]
    return {
        'n': len(sols),
        'raw': [marg([a[i] for a in sols]) for i in range(4)],
        'bnorm': [marg(b[i]) for i in range(4)],
        'pairdiff': [marg([(a[j] - a[i]) % p for a in sols])
                     for i in range(4) for j in range(i + 1, 4)],
        'bnorm_diff': [marg([(b[j][t] - b[i][t]) % p
                             for t in range(len(sols))])
                       for i in range(4) for j in range(i + 1, 4)],
    }


def main():
    OUT = {}
    for p in (17, 19, 29):
        k, (on, off), cls, eps = profile_sizes(p, T)
        szts = size_tuples(M, on, off, p)
        D = (p - 1) // 2
        from itertools import combinations
        allsets = [tuple(sorted(c)) for c in combinations(range(2, D + 1), 4)]
        seen = set()
        orbits = []
        for s in allsets:
            if s in seen:
                continue
            orb = sorted(set(renorm_orbit(s, p)))
            seen |= set(orb)
            orbits.append(orb)
        # per-set totals (Lemma V: equal across members)
        tot = {}
        for orb in orbits:
            for s in orb:
                rec = enum(p, s, szts)
                tot[s] = sum(len(v) for v in rec.values())
        mx = max(tot.values())
        top_sets = [s for s in tot if tot[s] == mx]
        print('=== p=%d: %d orbits; per-set total max %d at %d set(s) %s ==='
              % (p, len(orbits), mx, len(top_sets), top_sets[:3]))
        OUT['p%d' % p] = {'top_sets': [list(s) for s in top_sets],
                          'top_total': mx}

        # (a) V' localization: the orbit containing the top set
        orb = [o for o in orbits if top_sets[0] in o][0]
        print('  orbit of %s: %s' % (top_sets[0], orb))
        profs = {}
        for s in orb:
            rec = enum(p, s, szts)
            profs[s] = Counter(tuple(sorted([marg([a[i] for a in lst])
                                             for i in range(4)]))
                               for lst in rec.values() if lst)
            profs[s]['TOTALS'] = sum(len(v) for v in rec.values())
        base = profs[orb[0]]
        for s in orb[1:]:
            if profs[s] != base:
                d1 = base - profs[s]
                d2 = profs[s] - base
                print('  V\' MISMATCH %s vs %s:' % (orb[0], s))
                print('    profiles only in rep: %s' % dict(d1))
                print('    profiles only in member: %s' % dict(d2))

        # (b) fallback coordinates on the top set (and runner-up orbit)
        s0 = top_sets[0]
        rec = enum(p, s0, szts)
        best = max(rec.items(), key=lambda kv: len(kv[1]))
        szt, sols = best
        pr = phase_probe(p, s0, szt, sols)
        print('  top set %s @ szt %s: n=%d' % (s0, szt, pr['n']))
        print('    raw  marginals: %s' % pr['raw'])
        print('    bnorm marginals: %s' % pr['bnorm'])
        print('    pairdiff marginals: %s' % pr['pairdiff'])
        print('    bnorm_diff marginals: %s' % pr['bnorm_diff'])
        OUT['p%d' % p]['top_detail'] = {'set': list(s0), 'szt': list(szt),
                                        **{k_: v for k_, v in
                                           pr.items()}}

        # runner-up: second-highest total orbit, for contrast
        uniq = sorted(set(tot.values()), reverse=True)
        if len(uniq) > 1:
            s1 = [s for s in tot if tot[s] == uniq[1]][0]
            rec1 = enum(p, s1, szts)
            best1 = max(rec1.items(), key=lambda kv: len(kv[1]))
            szt1, sols1 = best1
            pr1 = phase_probe(p, s1, szt1, sols1)
            print('  runner set %s @ szt %s: n=%d' % (s1, szt1, pr1['n']))
            print('    raw %s | bnorm %s | pairdiff %s | bdiff %s'
                  % (pr1['raw'], pr1['bnorm'], pr1['pairdiff'],
                     pr1['bnorm_diff']))
            OUT['p%d' % p]['runner_detail'] = {'set': list(s1),
                                               'szt': list(szt1),
                                               **{k_: v for k_, v in
                                                  pr1.items()}}
        # min positive orbit
        pos = [s for s in tot if tot[s] > 0]
        mn = min(tot[s] for s in pos)
        s2 = [s for s in pos if tot[s] == mn][0]
        rec2 = enum(p, s2, szts)
        nz = [(szt, lst) for szt, lst in rec2.items() if lst]
        if nz:
            szt2, sols2 = max(nz, key=lambda kv: len(kv[1]))
            pr2 = phase_probe(p, s2, szt2, sols2)
            print('  min set %s @ szt %s: n=%d' % (s2, szt2, pr2['n']))
            print('    raw %s | bnorm %s | pairdiff %s | bdiff %s'
                  % (pr2['raw'], pr2['bnorm'], pr2['pairdiff'],
                     pr2['bnorm_diff']))
            OUT['p%d' % p]['min_detail'] = {'set': list(s2),
                                            'szt': list(szt2),
                                            **{k_: v for k_, v in
                                               pr2.items()}}
        print()

    # p=29 targeted top set phase probe (committed top {2,4,7,8})
    p = 29
    k, (on, off), cls, eps = profile_sizes(p, T)
    szts = size_tuples(M, on, off, p)
    for s0 in [(2, 4, 7, 8), (2, 4, 8, 13)]:
        rec = enum(p, s0, szts)
        nz = [(szt, lst) for szt, lst in rec.items() if lst]
        if not nz:
            print('p=29 set %s: inadmissible in this scope' % (s0,))
            continue
        szt, sols = max(nz, key=lambda kv: len(kv[1]))
        pr = phase_probe(p, s0, szt, sols)
        print('p=29 set %s @ szt %s: n=%d' % (s0, szt, pr['n']))
        print('  raw %s | bnorm %s | pairdiff %s | bdiff %s'
              % (pr['raw'], pr['bnorm'], pr['pairdiff'], pr['bnorm_diff']))
        OUT.setdefault('p29_targeted', {})[str(list(s0))] = {
            'set': list(s0), 'szt': list(szt),
            **{k_: v for k_, v in pr.items()}}
        print()

    with open('/home/z/my-project/scripts/out_strand_rigidity_debug.json',
              'w') as f:
        json.dump(OUT, f, indent=1, default=int)
    print('wrote scripts/out_strand_rigidity_debug.json')


if __name__ == '__main__':
    main()
