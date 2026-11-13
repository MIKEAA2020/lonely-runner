#!/usr/bin/env python3
"""
(2,2)-family structural check: every census parity-family solution has the
integer-parity parametrization (non-wrapping step-2 APs = one parity class
of I1 (possibly truncated) + at most ov spill points in A1), and the counts
per multiset equal the p-independent table.  Also: the (1,q)/(q,1) two-block
spans.  Machine check of the p-independence mechanism.
"""
import sys
from collections import Counter
from math import gcd


def pmrep(d, p):
    r = d % p
    return min(r, p - r)


def census_all(p, k, cls):
    sizes = (2 * k + 1, 2 * k + 2) if cls == 5 else (2 * k, 2 * k + 1, 2 * k + 2)
    ovs = (0, 1) if cls == 5 else (0, 1, 2)
    sols = []
    full = (1 << p) - 1
    inv = [0] * p
    for d in range(1, p):
        inv[d] = pow(d, -1, p)

    def fit(R, s3, d3, iv):
        if not R:
            return list(range(p))
        v = sorted((iv * r) % p for r in R)
        m = len(v)
        if m > s3:
            return []
        gaps = [(v[(i + 1) % m] - v[i]) % p for i in range(m)]
        gmax = max(gaps)
        if m == 1:
            arc_len, arc_start = 1, v[0]
        else:
            gi = gaps.index(gmax)
            arc_start = v[(gi + 1) % m]
            arc_len = p - gmax + 1
        if arc_len > s3:
            return []
        return [(((arc_start - u) % p) * d3) % p
                for u in range(s3 - arc_len + 1)]

    for s1 in sizes:
        for s2 in sizes:
            for s3 in sizes:
                ov = s1 + s2 + s3 - p
                if ov not in ovs:
                    continue
                m1 = 0
                for j in range(s1):
                    m1 |= 1 << j
                for d2 in range(1, p):
                    base2 = 0
                    for j in range(s2):
                        base2 |= 1 << ((d2 * j) % p)
                    for a2 in range(p):
                        m2 = ((base2 << a2) | (base2 >> (p - a2))) & full \
                            if a2 else base2
                        if (m1 & m2).bit_count() > ov:
                            continue
                        Rmask = full & ~(m1 | m2)
                        R = [i for i in range(p) if (Rmask >> i) & 1]
                        if len(R) > s3:
                            continue
                        for d3 in range(1, p):
                            for a3 in fit(R, s3, d3, inv[d3]):
                                sols.append(((0, 1, s1), (a2, d2, s2),
                                             (a3, d3, s3)))
    return sols


def main():
    print('=== (2,2) parity structure check ===')
    for p in (23, 29, 41, 31, 37, 43):
        cls = p % 6
        k = (p - cls) // 6
        sols = census_all(p, k, cls)
        par_ok = par_bad = 0
        ms_counts = Counter()
        tb_ok = tb_bad = 0
        for cfg in sols:
            (_, _, s1), (a2, d2, s2), (a3, d3, s3) = cfg
            f = (pmrep(d2, p), pmrep(d3, p))
            ms = tuple(sorted((s1, s2, s3)))
            ms_counts[(ms, f)] += 1
            if f == (2, 2):
                ok = True
                for (a, d, s) in ((a2, d2, s2), (a3, d3, s3)):
                    A = [(a + d * j) % p for j in range(s)]
                    inI1 = [x for x in A if s1 <= x <= p - 1]
                    out = [x for x in A if not (s1 <= x <= p - 1)]
                    # one integer parity, non-wrapping, spills <= ov in A1
                    pars = set(x % 2 for x in inI1)
                    if len(pars) > 1:
                        ok = False
                    if any(not (0 <= x <= s1 - 1) for x in out):
                        ok = False
                    if len(out) > (s1 + s2 + s3 - p):
                        ok = False
                    # non-wrapping: the sorted values are consecutive step 2
                    sv = sorted(inI1)
                    if any(b - a != 2 for a, b in zip(sv, sv[1:])):
                        ok = False
                par_ok += int(ok)
                par_bad += int(not ok)
            if f in ((1, (p - 1) // 2), ((p - 1) // 2, 1)):
                # two-block: blocks of ceil/floor consecutive ints with
                # right endpoints differing by exactly q (mod +-)
                q = (p - 1) // 2
                for (a, d, s) in ((a2, d2, s2), (a3, d3, s3)):
                    if pmrep(d, p) == q:
                        A = sorted((a + d * j) % p for j in range(s))
                        # block structure: consecutive runs
                        runs = []
                        cur = [A[0]]
                        for x in A[1:]:
                            if x == cur[-1] + 1:
                                cur.append(x)
                            else:
                                runs.append(cur)
                                cur = [x]
                        runs.append(cur)
                        if len(runs) == 2:
                            r1, r2 = runs
                            if (r2[-1] - r1[-1]) % p in (q, p - q) and \
                                    len(r1) + len(r2) == s and \
                                    abs(len(r1) - len(r2)) <= 1:
                                tb_ok += 1
                            else:
                                tb_bad += 1
                        else:
                            tb_bad += 1
        print('p=%2d (cls %d): (2,2) parity-structure %d ok / %d bad; '
              'two-block placements %d ok / %d bad'
              % (p, cls, par_ok, par_bad, tb_ok, tb_bad))
        # p-independence of the per-(multiset,family) counts
        tbl = {}
        for (ms, f), c in sorted(ms_counts.items()):
            tbl['%s|%s' % (ms, f)] = c
        print('   per-(multiset,family): %s' % tbl)
        sys.stdout.flush()


if __name__ == '__main__':
    main()
