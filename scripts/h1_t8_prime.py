#!/usr/bin/env python3
"""
h1_t8_prime.py -- the H1-general T=8 volume census at further zoo primes.

p=19 (k=2, eps=3): FULL distinct-family census (1,680 families) -- the
  orbit law and the eps-class dependence at fixed k.
p=29 (k=3, eps=5): TARGETED census -- the binary-cascade / resonant-anchor
  families (the p=17 top class + the anchor candidates) + a random sample
  for the median trend.  The k-scaling adjudication.

Same covering convention and cascade as h1_t8_census.py (Stage 1); GT via
the direct brute force on sample keys at each prime.

Resumable: state in scripts/h1_t8_prime_state.json.
Output: scripts/out_h1_t8_prime.json.
"""
import ast
import json
import random
import sys
import time
from collections import defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import (profile_sizes, size_tuples, ap_base, rot,
                               fit_placements)

STATE = '/home/z/my-project/scripts/h1_t8_prime_state.json'
OUT = '/home/z/my-project/scripts/out_h1_t8_prime.json'
T, M = 8, 5


def inv_table(p):
    return [0] + [pow(d, -1, p) for d in range(1, p)]


def count_sols(p, inv, key):
    (s1, s2, s3, s4, s5, d2, d3, d4, d5) = key
    ov = s1 + s2 + s3 + s4 + s5 - p
    full = (1 << p) - 1
    m1 = (1 << s1) - 1
    b2 = ap_base(d2, s2, p)
    b3 = ap_base(d3, s3, p)
    b4 = ap_base(d4, s4, p)
    n = 0
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
                    n += p
                    continue
                if nr > s5:
                    continue
                Rl = [i for i in range(p) if (R >> i) & 1]
                n += len(fit_placements(Rl, p, s5, d5, inv[d5]))
    return n


def brute_sols(p, key):
    (s1, s2, s3, s4, s5, d2, d3, d4, d5) = key
    full = (1 << p) - 1
    b1 = (1 << s1) - 1
    bs = [ap_base(d, s, p) for (d, s) in
          ((d2, s2), (d3, s3), (d4, s4), (d5, s5))]
    n = 0
    for a2 in range(p):
        m2 = rot(bs[0], a2, p, full)
        for a3 in range(p):
            m3 = rot(bs[1], a3, p, full)
            for a4 in range(p):
                m4 = rot(bs[2], a4, p, full)
                m1234 = b1 | m2 | m3 | m4
                for a5 in range(p):
                    if (m1234 | rot(bs[3], a5, p, full)) == full:
                        n += 1
    return n


def fams_for(p, mode):
    D = (p - 1) // 2
    if mode == 'full':
        vals = [d for d in range(2, D + 1)]
    elif mode == 'targeted':
        # anchors: doubling chain, small steps, t = floor((p+1)/3) ~ 3^-1,
        # q = (p-1)/2 antipode, and the p=17 top-class members
        t = (p + 1) // 3
        q = D
        vals = sorted(set([2, 3, 4, 5, 6, 7, 8, t, q, 2 * t % p]))
        vals = [d for d in vals if 2 <= d <= D]
    elif mode == 'cascade':
        # the binary-cascade measurement: all sets containing {2,4,8}
        vals = sorted(set([2, 4, 8] + [x for x in range(3, (p - 1) // 2 + 1)
                                       if x not in (2, 4, 8)]))
        vals = [d for d in vals if 2 <= d <= D]
    import itertools
    fams = []
    for comb in itertools.combinations(sorted(set(vals)), 4):
        if mode == 'cascade' and not set(comb) >= {2, 4, 8}:
            continue
        for perm in itertools.permutations(comb):
            fams.append(perm)
    return fams


def run_prime(p, mode, sample_extra=0, gt_n=3):
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(M, on, off, p)
    inv = inv_table(p)
    print('=== p=%d (k=%d, eps=%d, %s): sizes (%d,%d), %d size tuples, '
          'ov range [%d,%d] ===' % (p, k, eps, cls, on, off, len(szts),
                                    sum(szts[0]) - p, sum(szts[-1]) - p))
    fams = fams_for(p, mode)
    if sample_extra:
        rng = random.Random(20261007 + p)
        pool = []
        import itertools as it
        D = (p - 1) // 2
        allsets = list(it.combinations(range(2, D + 1), 4))
        chosen = rng.sample(allsets, min(sample_extra, len(allsets)))
        have = {tuple(sorted(f)) for f in fams}
        for c in chosen:
            if c not in have:
                pool.append(c)
        fams = fams + pool
    print('families to census: %d (%s)' % (len(fams), mode))
    try:
        with open(STATE) as f:
            st = json.load(f)
    except Exception:
        st = {}
    key_ = 'p%d' % p
    st.setdefault(key_, {'done': {}, 'gt': []})
    counts = st[key_]['done']
    todo = [f for f in fams if str(list(f)) not in counts]
    t0 = time.time()
    for idx, fam in enumerate(todo):
        per_size = {}
        for szt in szts:
            per_size[str(szt)] = count_sols(p, inv, tuple(szt) + fam)
        counts[str(list(fam))] = per_size
        if idx % 20 == 0 or idx == len(todo) - 1:
            with open(STATE, 'w') as f:
                json.dump(st, f)
            el = time.time() - t0
            if idx % 100 == 0 or idx == len(todo) - 1:
                print('  [%4d/%4d] %s total=%d (%.1fs, %.2fs/fam)'
                      % (idx + 1, len(todo), fam,
                         sum(per_size.values()), el, el / (idx + 1)))
                sys.stdout.flush()
    with open(STATE, 'w') as f:
        json.dump(st, f)
    # stats
    tot = {tuple(json.loads(kk)): sum(v.values())
           for kk, v in counts.items()}
    vols = sorted(tot.values())
    n_adm = sum(1 for v in vols if v > 0)
    print('  census: %d families, admissible %d, totals median %d p90 %d '
          'max %d' % (len(tot), n_adm, vols[len(vols) // 2],
                      vols[int(len(vols) * .9)], vols[-1]))
    top = sorted(tot.items(), key=lambda kv: -kv[1])[:12]
    for f, v in top:
        per = counts[str(list(f))]
        bs = max(per.items(), key=lambda kv: kv[1])
        s = ast.literal_eval(bs[0])
        print('    %s: total %7d, max/size %6d at %s (ov=%d)'
              % (f, v, bs[1], bs[0], sum(s) - p))
    # orbit law check (full mode only)
    if mode == 'full':
        by_set = defaultdict(list)
        for f, v in tot.items():
            by_set[tuple(sorted(f))].append(v)
        nclasses = len({vols_i for vols_i in
                        (s[0] for s in by_set.values())})
        allord = all(len(set(vl)) == 1 for vl in by_set.values())
        print('  orbit law: %d sets, %d volume classes, all 24 orders '
              'equal: %s' % (len(by_set), nclasses, allord))
    # GT
    if gt_n and len(top):
        rng = random.Random(p)
        gtrec = []
        for f, v in top[:gt_n - 1]:
            per = counts[str(list(f))]
            bs = max(per.items(), key=lambda kv: kv[1])
            key = tuple(ast.literal_eval(bs[0])) + f
            nc = per[bs[0]]
            t1 = time.time()
            nb = brute_sols(p, key)
            gtrec.append({'key': list(key), 'brute': nb, 'cascade': nc,
                          'match': nb == nc,
                          'brute_s': round(time.time() - t1, 1)})
            print('    GT %s: brute %d vs cascade %d -> %s'
                  % (key, nb, nc, 'MATCH' if nb == nc else 'MISMATCH'))
        st[key_]['gt'] = gtrec
        with open(STATE, 'w') as f:
            json.dump(st, f)
    return tot


def main():
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = {}
    if which in ('all', 'p19'):
        out['p19'] = {str(list(f)): v for f, v in
                      run_prime(19, 'full', gt_n=2).items()}
    if which in ('all', 'cascade29'):
        out['p29_cascade'] = {str(list(f)): v for f, v in
                              run_prime(29, 'cascade', sample_extra=120,
                                        gt_n=2).items()}
    if which in ('all', 'targeted29'):
        out['p29'] = {str(list(f)): v for f, v in
                      run_prime(29, 'targeted', sample_extra=150,
                                gt_n=2).items()}
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)
    print('wrote %s' % OUT)


if __name__ == '__main__':
    main()
