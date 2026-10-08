#!/usr/bin/env python3
"""
h1_t8_census.py -- the H1-general attack at T=8, Stage 1 (T-20 directive).

Rebuilds |Sol| for ALL 840 +/-distinct-difference families at the committed
zoo cell (T=8 (1K,5U), p=17), over all 32 admissible size tuples, in the
COVERING convention (the cascade's Sol: includes the first-4-already-cover
bystanders -- the m=4 GT convention of T-18/T-19).

Counting-only pass (no solution lists stored; lists are Stage 2 for the top
families).  The counting cascade mirrors solutions_m5 pruning exactly.

GT: independent direct brute force over all p^4 start tuples (fresh code
path, no cascade, no pruning) on a sample of keys spanning admissible and
empty ones, compared count-by-count.

Also records, per family: total volume, max per-size volume, and the
chain-law benchmark 24*#{o>=0, sum=ov, o_i<=s} for the top tuples (the
c_family calibration).

Output: scripts/out_h1_t8_census.json (+ stdout).  Resumable via state file
scripts/h1_t8_state.json (this environment reaps background jobs).
"""
import json
import sys
import time
from collections import defaultdict
from math import comb

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import (profile_sizes, size_tuples, ap_base, rot,
                               fit_placements)

P, T, M = 17, 8, 5
STATE = '/home/z/my-project/scripts/h1_t8_state.json'
OUT = '/home/z/my-project/scripts/out_h1_t8_census.json'


def inv_table(p):
    return [0] + [pow(d, -1, p) for d in range(1, p)]


def count_sols(p, key):
    """Count covering placements (a2..a5) for key=(s1..s5,d2..d5).
    Mirrors solutions_m5 exactly (same pruning), counts only."""
    (s1, s2, s3, s4, s5, d2, d3, d4, d5) = key
    ov = s1 + s2 + s3 + s4 + s5 - p
    full = (1 << p) - 1
    inv = INV
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
                    n += p            # free a5 (bystanders included)
                    continue
                if nr > s5:
                    continue
                Rl = [i for i in range(p) if (R >> i) & 1]
                n += len(fit_placements(Rl, p, s5, d5, inv[d5]))
    return n


def brute_sols(p, key):
    """Independent GT: direct enumeration of ALL p^4 start tuples, no
    cascade, no pruning, no fit_placements -- plain mask OR test."""
    (s1, s2, s3, s4, s5, d2, d3, d4, d5) = key
    full = (1 << p) - 1
    b1 = (1 << s1) - 1
    b2 = ap_base(d2, s2, p)
    b3 = ap_base(d3, s3, p)
    b4 = ap_base(d4, s4, p)
    b5 = ap_base(d5, s5, p)
    n = 0
    for a2 in range(p):
        m2 = rot(b2, a2, p, full)
        for a3 in range(p):
            m3 = rot(b3, a3, p, full)
            m23 = m2 | m3
            for a4 in range(p):
                m4 = rot(b4, a4, p, full)
                m1234 = b1 | m23 | m4
                for a5 in range(p):
                    m5 = rot(b5, a5, p, full)
                    if (m1234 | m5) == full:
                        n += 1
    return n


def chain_benchmark(s, ov, m=M):
    """24 * #{o >= 0 in Z^m : sum o = ov, o_i <= s_i} (cyclic orders with
    the normalizer first: (m-1)! = 24), caps by the multiset s."""
    import itertools
    cnt = 0
    for perm in itertools.permutations(range(m)):
        if perm[0] != 0:
            continue
        caps = [s[perm[t]] for t in range(m)]
        dp = [1] + [0] * ov
        for cap in caps:
            ndp = [0] * (ov + 1)
            run = 0
            for v in range(ov + 1):
                run += dp[v]
                if v - cap - 1 >= 0:
                    run -= dp[v - cap - 1]
                ndp[v] = run
            dp = ndp
        cnt += dp[ov]
    return cnt


def main():
    global INV
    p = P
    INV = inv_table(p)
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(M, on, off, p)
    print('cell T=%d (1K,%dU) p=%d: k=%d eps=%d class=%s sizes=(%d,%d), '
          '%d size tuples' % (T, M - 1, p, k, eps, cls, on, off, len(szts)))
    print('size-tuple sums: %s' % sorted({sum(s) for s in szts}))

    # the 840 distinct families from the committed census
    with open('/home/z/my-project/scripts/out_higher_rung_step1.json') as f:
        step1 = json.load(f)
    zoo = step1['T8_m5_census']['p17']
    import ast
    distinct = sorted(tuple(ast.literal_eval(s)) for s in zoo['families']
                      if len(set(ast.literal_eval(s))) == 4
                      and 1 not in ast.literal_eval(s))
    print('distinct families: %d (committed n_distinct_families=%d)'
          % (len(distinct), zoo['n_distinct_families']))

    # resumable state
    try:
        with open(STATE) as f:
            st = json.load(f)
    except Exception:
        st = {'done': {}, 'gt': []}
    t0 = time.time()
    counts = {k_: v for k_, v in st['done'].items()}
    keys_done = {tuple(json.loads(k_)) for k_ in counts}
    todo = [f for f in distinct if f not in keys_done]
    print('resuming: %d done, %d todo' % (len(counts), len(todo)))
    for idx, fam in enumerate(todo):
        per_size = {}
        for szt in szts:
            key = tuple(szt) + fam
            n = count_sols(p, key)
            per_size[str(szt)] = n
        counts[str(list(fam))] = per_size
        st['done'] = counts
        if idx % 40 == 0 or idx == len(todo) - 1:
            with open(STATE, 'w') as f:
                json.dump(st, f)
            el = time.time() - t0
            print('  [%4d/%4d] %s total=%d  (%.1fs elapsed, %.2fs/fam)'
                  % (idx + 1, len(todo), fam,
                     sum(per_size.values()), el, el / (idx + 1)))
            sys.stdout.flush()
    with open(STATE, 'w') as f:
        json.dump(st, f)

    # ---------------- aggregate stats ----------------
    tot_by_fam = {tuple(json.loads(k)): sum(v.values())
                  for k, v in counts.items()}
    max_by_fam = {tuple(json.loads(k)): max(v.values())
                  for k, v in counts.items()}
    vols = sorted(tot_by_fam.values())
    maxes = sorted(max_by_fam.values())
    n_adm = sum(1 for v in vols if v > 0)
    print()
    print('=== H1 volume census (covering convention), 840 distinct '
          'families ===')
    print('per-family totals: median %d, p90 %d, max %d; admissible %d/%d'
          % (vols[len(vols) // 2], vols[int(len(vols) * .9)], vols[-1],
             n_adm, len(vols)))
    print('per-(family,size) maxes: median %d, p90 %d, max %d'
          % (maxes[len(maxes) // 2], maxes[int(len(maxes) * .9)],
             maxes[-1]))
    # ov benchmark
    print('ov by size tuple: %s' %
          {str(szt): sum(szt) - p for szt in szts[:4]})
    top = sorted(tot_by_fam.items(), key=lambda kv: -kv[1])[:20]
    print('top-20 families by total volume:')
    for f, v in top:
        ms = max_by_fam[f]
        # best size tuple for this family
        per = counts[str(list(f))]
        bs = max(per.items(), key=lambda kv: kv[1])
        s = ast.literal_eval(bs[0])
        ov = sum(s) - p
        bm = chain_benchmark(s, ov)
        print('  %s: total %6d, max/size %5d at %s (ov=%d, chain-bench %d,'
               ' ratio %.2f)' % (f, v, ms, bs[0], ov, bm, ms / float(bm)))

    # ---------------- GT sample ----------------
    print()
    print('=== GT: independent brute force (direct p^4 enumeration) ===')
    rng_vols = vols[:]
    gt_keys = []
    # top-3 families' best tuples + a couple of mid/empty keys
    for f, v in top[:3]:
        per = counts[str(list(f))]
        bs = max(per.items(), key=lambda kv: kv[1])
        gt_keys.append(tuple(ast.literal_eval(bs[0])) + f)
    # a mid family and the minimum-volume one (all 840 are admissible at
    # p=17 -- the committed pool identity -- so no zero family exists)
    mid = vols[len(vols) // 2]
    midf = next(f for f, v in tot_by_fam.items() if v == mid)
    per = counts[str(list(midf))]
    bs = max(per.items(), key=lambda kv: kv[1])
    gt_keys.append(tuple(ast.literal_eval(bs[0])) + midf)
    minf = min(tot_by_fam.items(), key=lambda kv: kv[1])[0]
    gt_keys.append(tuple(szts[0]) + minf)
    # also the committed chain family top tuple as a control
    gt_keys.append((5, 5, 5, 5, 5) + (1, 1, 1, 1))
    gt_rec = []
    for key in gt_keys:
        t1 = time.time()
        nb = brute_sols(p, key)
        t2 = time.time()
        nc = count_sols(p, key)
        ok = (nb == nc)
        gt_rec.append({'key': list(key), 'brute': nb, 'cascade': nc,
                       'match': ok, 'brute_s': round(t2 - t1, 1)})
        print('  key %s: brute %d vs cascade %d -> %s (%.1fs)'
              % (key, nb, nc, 'MATCH' if ok else 'MISMATCH', t2 - t1))
        sys.stdout.flush()
    st['gt'] = gt_rec
    with open(STATE, 'w') as f:
        json.dump(st, f)

    out = {
        'cell': {'T': T, 'p': p, 'k': k, 'eps': eps, 'cls': cls,
                 'on': on, 'off': off, 'm': M},
        'n_distinct': len(distinct),
        'n_admissible': n_adm,
        'totals': {str(list(f)): v for f, v in tot_by_fam.items()},
        'per_size': counts,
        'stats': {
            'total_median': vols[len(vols) // 2],
            'total_p90': vols[int(len(vols) * .9)],
            'total_max': vols[-1],
            'max_per_size_median': maxes[len(maxes) // 2],
            'max_per_size_p90': maxes[int(len(maxes) * .9)],
            'max_per_size_max': maxes[-1],
        },
        'gt': gt_rec,
    }
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)
    print()
    print('wrote %s' % OUT)


if __name__ == '__main__':
    main()
