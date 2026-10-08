#!/usr/bin/env python3
"""
pq side, part 2: the size-restricted row census + the pinning/drift analysis.

(1) ROWCEN-RESTRICTED: at p = 1 mod 6 the row profile is {2k, 2k+1}
    (cap 2k+1, slack 2): the row-relevant census is s1 = 2k+1 with
    (s2,s3) in {2k, 2k+1}^2, ov in {0,1,2} -- per-ordered-triple family
    counts.  At p = 5 mod 6: s1 = 2k+1, (s2,s3) = (2k+2, 2k+2) (partition).

(2) PINNING: for N = pq, over all (a, b, u1 < u2):
      - rows y outside the kernel-p AP Y_a; good row := X_b^c subset of
        f_1(y) u f_2(y);
      - every good row's normalized family (pmrep(b*d1), pmrep(b*d2))
        must be a census family at the row modulus, and its normalized
        footprint SETS must equal a census solution (the pinned set);
      - max #good rows over all families (the shear margin);
      - drift rates of the normalized starts.
"""
import json
import sys
from collections import Counter
from math import gcd

OUT = {}


def pmrep(d, p):
    r = d % p
    return min(r, p - r)


def norm(v, N):
    v %= N
    return min(v, N - v)


# ----------------------------------------------------------------------
# (1) size-restricted row census (reuses the pq_side census, restricted)
# ----------------------------------------------------------------------

def census_s1(p, k, cls, s1_fixed, allowed_sizes, ovs):
    sizes = allowed_sizes
    sols = []
    full = (1 << p) - 1
    inv = [0] * p
    for d in range(1, p):
        inv[d] = pow(d, -1, p)

    def fit(R, p, s3, d3, iv):
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

    for s2 in sizes:
        for s3 in sizes:
            ov = s1_fixed + s2 + s3 - p
            if ov not in ovs:
                continue
            m1 = 0
            for j in range(s1_fixed):
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
                        for a3 in fit(R, p, s3, d3, inv[d3]):
                            sols.append(((0, 1, s1_fixed), (a2, d2, s2),
                                         (a3, d3, s3)))
    return sols


def row_census_restricted():
    res = {}
    q_of = lambda p: (p - 1) // 2
    print('--- row census, size-restricted ---')
    for p in (11, 17, 23, 29, 41, 47, 53, 59):
        k = (p - 5) // 6
        s1 = 2 * k + 1
        sols = census_s1(p, k, 5, s1, (2 * k + 1, 2 * k + 2), {0})
        fam = Counter()
        for cfg in sols:
            (_, _, _), (_, d2, _), (_, d3, _) = cfg
            fam[(pmrep(d2, p), pmrep(d3, p))] += 1
        res['eps5_%d' % p] = {'n': len(sols),
                               'families': {str(t): c for t, c in
                                            sorted(fam.items())}}
        print('p=%2d (5 mod 6): %d coverings, families %s'
              % (p, len(sols), dict(sorted(fam.items()))))
    for p in (19, 31, 37, 43, 61):
        k = (p - 1) // 6
        s1 = 2 * k + 1
        sols = census_s1(p, k, 1, s1, (2 * k, 2 * k + 1), {0, 1, 2})
        fam = Counter()
        fam_ms = Counter()
        for cfg in sols:
            (_, _, _), (_, d2, s2), (_, d3, s3) = cfg
            fam[(pmrep(d2, p), pmrep(d3, p))] += 1
            fam_ms[(pmrep(d2, p), pmrep(d3, p), s2, s3)] += 1
        res['eps1_%d' % p] = {'n': len(sols),
                              'families': {str(t): c for t, c in
                                           sorted(fam.items())},
                              'per_sizes': {'%d_%d' % (s2, s3): v
                                            for (f2, f3, s2, s3), v in
                                            sorted(fam_ms.items())}}
        print('p=%2d (1 mod 6): %d coverings, families %s'
              % (p, len(sols), dict(sorted(fam.items()))))
        print('   per (s2,s3): %s' % res['eps1_%d' % p]['per_sizes'])
    sys.stdout.flush()
    return res


# ----------------------------------------------------------------------
# (2) pinning / drift analysis
# ----------------------------------------------------------------------

class PQ:
    def __init__(self, p, q):
        self.p, self.q, self.N = p, q, p * q
        self.c = (self.N - 1) // 6
        self.rho_p = self.c // q
        self.rho_q = self.c // p
        self.units_pm = [u for u in range(1, self.N)
                         if gcd(u, self.N) == 1 and u < self.N - u]

    def footprint_formula(self, u, y):
        p, q, c = self.p, self.q, self.c
        alpha = pow(u, -1, p)
        sigma = (u * y) % q
        jlo = -((c + sigma) // q)
        jhi = (c - sigma) // q
        pts = frozenset((alpha * (sigma + q * jp)) % p
                        for jp in range(jlo, jhi + 1))
        start = (alpha * (sigma + q * jlo)) % p
        return pts, start, jhi - jlo + 1


def pinning():
    res = {}
    for (p, q) in [(5, 7), (5, 11), (5, 19), (7, 11), (7, 13), (7, 17),
                   (11, 13), (5, 23)]:
        W = PQ(p, q)
        p_, q_, c = W.p, W.q, W.c
        rho = W.rho_p
        # the pinned census set at the row modulus (normal form)
        k = (p_ - p_ % 6) // 6
        if p_ % 6 == 5:
            s1 = 2 * k + 1
            pinned = census_s1(p_, k, 5, s1, (2 * k + 1, 2 * k + 2), {0})
            fam_ok = {(1, 1), (1, (p_ - 1) // 2), (2, 2),
                      ((p_ - 1) // 2, 1)}
        elif p_ >= 19:
            s1 = 2 * k + 1
            pinned = census_s1(p_, k, 1, s1, (2 * k, 2 * k + 1), {0, 1, 2})
            fam_ok = {(1, 1), (1, (p_ - 1) // 2), (2, 2),
                      ((p_ - 1) // 2, 1)}
        else:
            s1 = 2 * k + 1
            pinned = census_s1(p_, k, 1, s1, (2 * k, 2 * k + 1, 2 * k + 2),
                               {0, 1, 2})
            fam_ok = None   # boundary primes: record, don't assert
        pinned_sets = [(frozenset((a2 + d2 * j) % p_ for j in range(s2)),
                        frozenset((a3 + d3 * j) % p_ for j in range(s3)),
                        (pmrep(d2, p_), pmrep(d3, p_)))
                       for (_, _, _), (a2, d2, s2), (a3, d3, s3) in pinned]
        pinned_pairs = set((A2, A3) for A2, A3, _ in pinned_sets)

        X = frozenset((pow(1, 1, p_) * t) % p_ for t in [])  # placeholder
        # column AP X_b and row AP Y_a
        maxgood = (None, -1)
        fam_good = Counter()
        n_in_family = n_pinned = n_good_rows = 0
        drift_counter = Counter()
        collapse_cases = 0
        n_fam = 0
        for a in [a for a in range(1, q_) if gcd(a, q_) == 1 and a <= q_ - a]:
            Ymask = [norm(a * y, q_) <= W.rho_q for y in range(q_)]
            Yc = [y for y in range(q_) if not Ymask[y]]
            for b in [b for b in range(1, p_) if gcd(b, p_) == 1
                      and b <= p_ - b]:
                Xb = frozenset(x for x in range(p_) if norm(b * x, p_) <= rho)
                Xc = frozenset(x for x in range(p_) if x not in Xb)
                # normalization map: x -> b*x + rho  (sends X_b -> [0,2rho])
                def nrm(S):
                    return frozenset((b * x + rho) % p_ for x in S)
                for i1 in range(len(W.units_pm)):
                    u1 = W.units_pm[i1]
                    d1 = (pow(u1, -1, p_) * q_) % p_
                    bd1 = (b * d1) % p_
                    F1 = {y: W.footprint_formula(u1, y)[0] for y in Yc}
                    for i2 in range(i1 + 1, len(W.units_pm)):
                        u2 = W.units_pm[i2]
                        d2 = (pow(u2, -1, p_) * q_) % p_
                        bd2 = (b * d2) % p_
                        n_fam += 1
                        good = []
                        for y in Yc:
                            f2, _, _ = W.footprint_formula(u2, y)
                            if Xc <= (F1[y] | f2):
                                good.append(y)
                                n_good_rows += 1
                                fam = (pmrep(bd1, p_), pmrep(bd2, p_))
                                in_fam = (fam_ok is None or fam in fam_ok)
                                n_in_family += int(in_fam)
                                fam_good[fam] += 1
                                # pinned-set check (either assignment order)
                                A1n, A2n = nrm(F1[y]), nrm(f2)
                                pin = (A1n, A2n) in pinned_pairs or \
                                      (A2n, A1n) in pinned_pairs
                                n_pinned += int(pin)
                        if len(good) > maxgood[1]:
                            maxgood = ((a, b, u1, u2), len(good))
                        if len(good) >= 2:
                            # drift of the normalized starts
                            s_prev = None
                            for y in good:
                                _, st1, _ = W.footprint_formula(u1, y)
                                _, st2, _ = W.footprint_formula(u2, y)
                                cur = ((b * st1 + rho) % p_,
                                       (b * st2 + rho) % p_)
                                if s_prev is not None and y == s_prev[0] + 1:
                                    drift_counter[(cur[0] - s_prev[1][0])
                                                  % p_,
                                                  (cur[1] - s_prev[1][1])
                                                  % p_] += 1
                                s_prev = (y, cur)
                            # collapse detection: equal drift on both units
                            rates1 = set()
                            prev = None
                            for y in range(min(Yc), max(Yc) + 1):
                                _, st1, _ = W.footprint_formula(u1, y)
                                _, st2, _ = W.footprint_formula(u2, y)
                                cur = ((b * st1 + rho) % p_,
                                       (b * st2 + rho) % p_)
                                if prev is not None:
                                    rates1.add((cur[0] - prev[0]) % p_)
                                prev = cur
                            if len(rates1) == 1:
                                collapse_cases += 1
        res['%d_%d' % (p_, q_)] = {
            'families_scanned': n_fam,
            'good_rows_total': n_good_rows,
            'good_rows_in_census_family': n_in_family,
            'good_rows_pinned': n_pinned,
            'max_good_rows': maxgood[1], 'argmax': maxgood[0],
            'family_histogram': {str(t): c for t, c in
                                 sorted(fam_good.items())},
            'drift_pairs': {str(t): c for t, c in
                            sorted(drift_counter.items())},
            'collapse_rate_families': collapse_cases,
        }
        print('[PIN] N=%3d: %d families, %d good rows; in-family %d/%d, '
              'pinned %d/%d; max good rows %d (argmin a,b,u=%s); collapse '
              'families %d'
              % (W.N, n_fam, n_good_rows, n_in_family, n_good_rows,
                 n_pinned, n_good_rows, maxgood[1], maxgood[0],
                 collapse_cases))
        sys.stdout.flush()
    return res


def main():
    OUT['row_census_restricted'] = row_census_restricted()
    OUT['pinning'] = pinning()
    with open('out_pq_side2.json', 'w') as f:
        json.dump(OUT, f, indent=1, default=str)
    print('written: out_pq_side2.json')


if __name__ == '__main__':
    main()
