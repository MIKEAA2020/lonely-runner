#!/usr/bin/env python3
"""
pq side: the tensor translation (Theorem R) + cell verification.

THEOREM R (row translation; proved on paper, machine-checked here).
At N = pq (5 <= p < q primes), c = floor((N-1)/6), rho_p = floor(c/q),
rho_q = floor(c/p):
  - kernel-p member (speed p*a, a in Z_q^*) = the union of the
    2*rho_q+1 rows Y_a = {y in Z_q : ||a y||_q <= rho_q};
    kernel-q member (speed q*b, b in Z_p^*) = the union of the
    2*rho_p+1 columns X_b = {x in Z_p : ||b x||_p <= rho_p}.
  - a unit u meets the row y (the pi_q-fiber {k = y + qj}) in an AP in Z_p:
    f_u(y) = { (u^{-1} mod p)(sigma + q j') mod p : j' in J(y) },
    sigma = (u*y mod q) in [0,q),  J(y) = [ceil((-c-sigma)/q), floor((c-sigma)/q)],
    an interval of n_u(y) = floor((c-sigma)/q)+floor((c+sigma)/q)+1
    consecutive integers, n_u(y) in {2*rho_p, 2*rho_p+1, 2*rho_p+2}.
  - (1,1,2) covering <=> for every row y in Y_a^c: X_b^c subset of
    f_u1(y) u f_u2(y).  At p = 5 mod 6, rho_p = (p-5)/6 and
    |X_b^c| = p - 2*rho_p - 1 = 2*(2*rho_p+2): both footprints must have
    the maximal size and partition X_b^c exactly -- i.e. (f_1, f_2, X_b)
    is a critical three-AP PARTITION of Z_p of the B-side multiset
    {2k+2, 2k+2, 2k+1} with the third AP the fixed interval X_b.

Experiments:
  [CONST]   the closed-form constants (rho_p = floor(p/6) etc.).
  [STRUCT]  footprint formula: exact set equality, all units u, all rows y,
            several N; size profile; drift rates of the AP starts.
  [ROWCEN]  row-level census at the row modulus p: critical coverings with
            s1 = |X| = 2*rho_p+1 fixed (the X-restricted census) -- the
            four families with an interval member; counts constant in p.
  [CELL112] (1,1,2) exhaustive at N in {35,55,65,77,85,91,95,115,119,143}:
            0 coverings; margins; T-10 reproduction at 77/91.
  [CELL13U] (1K,3U), both (kappa_p,kappa_q) in {(1,0),(0,1)}: exhaustive
            unit triples: 0 coverings; margins.
  [SHEAR]   for the closest (1,1,2) families: per-row deficits, the good
            rows' membership in the census families, the drift rates.

Output: scripts/out_pq_side.json + stdout log.
"""
import json
import sys
from collections import Counter
from math import comb, gcd

import numpy as np

OUT = {}


def pmrep(d, p):
    r = d % p
    return min(r, p - r)


def norm(v, N):
    v %= N
    return min(v, N - v)


# ----------------------------------------------------------------------
# structures at N = pq
# ----------------------------------------------------------------------

class PQ:
    def __init__(self, p, q):
        assert p < q and p >= 5
        self.p, self.q, self.N = p, q, p * q
        self.c = (self.N - 1) // 6
        self.rho_p = self.c // q
        self.rho_q = self.c // p
        self.units = [u for u in range(1, self.N) if gcd(u, self.N) == 1]
        self.units_pm = [u for u in self.units if u < self.N - u]

    def B_unit(self, u):
        k = np.arange(self.N)
        uk = (u * k) % self.N
        return uk, (np.minimum(uk, self.N - uk) <= self.c)

    def rowset(self, a):
        """Y_a = rows covered by kernel-p member of speed p*a."""
        y = np.arange(self.q)
        ay = (a * y) % self.q
        return (np.minimum(ay, self.q - ay) <= self.rho_q)

    def colset(self, b):
        x = np.arange(self.p)
        bx = (b * x) % self.p
        return (np.minimum(bx, self.p - bx) <= self.rho_p)

    def B_kernelp(self, a):
        Y = self.rowset(a)
        k = np.arange(self.N)
        return Y[k % self.q]

    def B_kernelq(self, b):
        X = self.colset(b)
        k = np.arange(self.N)
        return X[k % self.p]

    def footprint(self, u, y):
        """direct: {x : (x,y) in B_u}; returns sorted list."""
        p, q, N, c = self.p, self.q, self.N, self.c
        out = []
        for j in range(p):
            k = (y + q * j) % N
            if norm(u * k, N) <= c:
                out.append(k % p)
        return sorted(set(out))

    def footprint_formula(self, u, y):
        """the AP formula of Theorem R; returns (points, j_interval)."""
        p, q, c = self.p, self.q, self.c
        alpha = pow(u, -1, p)
        sigma = (u * y) % q
        jlo = -((c + sigma) // q)          # ceil((-c-sigma)/q)
        jhi = (c - sigma) // q             # floor((c-sigma)/q)
        pts = [(alpha * (sigma + q * jp)) % p for jp in range(jlo, jhi + 1)]
        return sorted(set(pts)), (jlo, jhi)


# ----------------------------------------------------------------------
# [CONST] + [STRUCT]
# ----------------------------------------------------------------------

def const_and_struct():
    res = {}
    for (p, q) in [(5, 7), (5, 11), (7, 11), (7, 13), (11, 13), (5, 19),
                   (7, 17), (5, 23)]:
        W = PQ(p, q)
        ok_const = (W.rho_p == p // 6 and W.rho_q == q // 6)
        eqn = (p - 2 * W.rho_p - 1) == 2 * (2 * W.rho_p + 2) if p % 6 == 5 \
            else None
        # STRUCT: set equality + size profile + drift
        n_ok = n_eq = 0
        sizes = Counter()
        drifts = Counter()
        starts = {}
        for u in W.units_pm:
            alpha = pow(u, -1, p)
            for y in range(q):
                f1 = W.footprint(u, y)
                f2, J = W.footprint_formula(u, y)
                n_eq += int(f1 == f2)
                n_ok += 1
                sizes[len(f1)] += 1
                if f1:
                    starts[(u, y)] = f1[0] if len(f1) > 1 else f1[0]
        # drift of the *set* start is ill-defined under rotation; use the
        # formula start a(y) = alpha*(sigma + q*jlo)
        for u in W.units_pm[:12]:
            alpha = pow(u, -1, p)
            for y in range(q - 1):
                s0 = (alpha * ((u * y) % q + q * (-((W.c + (u * y) % q) // q)))) % p
                s1 = (alpha * ((u * (y + 1)) % q + q * (-((W.c + (u * (y + 1)) % q) // q)))) % p
                drifts[(s1 - s0) % p] += 1
        res['%d_%d' % (p, q)] = {
            'rho_p': W.rho_p, 'rho_q': W.rho_q,
            'rho_p_is_floor_p6': ok_const,
            'X_size': 2 * W.rho_p + 1, 'Y_size': 2 * W.rho_q + 1,
            'Xc_size': p - 2 * W.rho_p - 1, 'Yc_size': q - 2 * W.rho_q - 1,
            'partition_equality_p5': eqn,
            'struct_checks': n_eq, 'struct_total': n_ok,
            'struct_ok': n_eq == n_ok,
            'size_profile': {str(k): v for k, v in sizes.items()},
            'max_footprint': max(sizes),
            'drift_values': {str(k): v for k, v in drifts.items()},
        }
        print('[CONST] N=%3d: rho_p=%d(=%d?) rho_q=%d |X|=%d |X^c|=%d '
              'equality(p5)=%s ; [STRUCT] %d/%d set-equalities, sizes %s'
              % (W.N, W.rho_p, p // 6, W.rho_q, 2 * W.rho_p + 1,
                 p - 2 * W.rho_p - 1, eqn, n_eq, n_ok, dict(sizes)))
        sys.stdout.flush()
    return res


# ----------------------------------------------------------------------
# [ROWCEN] X-restricted census at the row modulus
# ----------------------------------------------------------------------

def row_census():
    """census of critical coverings at modulus p with s1 = |X| = 2*rho+1
    fixed, over all (d2,d3); returns family counts."""
    res = {}

    def fit(R, p, s3, d3, inv):
        if not R:
            return list(range(p))
        v = sorted((inv * r) % p for r in R)
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

    def census(p, k, cls, s1_fixed, ovs):
        sizes = (2 * k + 1, 2 * k + 2) if cls == 5 else \
            (2 * k, 2 * k + 1, 2 * k + 2)
        sols = []
        full = (1 << p) - 1
        inv = [0] * p
        for d in range(1, p):
            inv[d] = pow(d, -1, p)
        for s2 in sizes:
            for s3 in sizes:
                if (s1_fixed + s2 + s3 - p) not in ovs:
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
                        if (m1 & m2).bit_count() > (s1_fixed + s2 + s3 - p):
                            continue
                        Rmask = full & ~(m1 | m2)
                        R = [i for i in range(p) if (Rmask >> i) & 1]
                        if len(R) > s3:
                            continue
                        for d3 in range(1, p):
                            for a3 in fit(R, p, s3, d3, inv[d3]):
                                sols.append(((0, 1, s1_fixed),
                                             (a2, d2, s2), (a3, d3, s3)))
        return sols

    for p in (5, 11, 17, 23, 29, 41, 47, 53, 59):
        k = (p - 5) // 6
        rho = p // 6
        s1 = 2 * rho + 1
        sols = census(p, k, 5, s1, {0, 1})
        fam = Counter()
        for cfg in sols:
            (_, _, _), (_, d2, _), (_, d3, _) = cfg
            fam[(pmrep(d2, p), pmrep(d3, p))] += 1
        res['%d' % p] = {'s1': s1, 'n': len(sols),
                         'families': {str(t): c for t, c in
                                      sorted(fam.items())}}
        print('[ROWCEN] row p=%2d (5 mod 6): s1=%d, %d coverings, '
              'families %s' % (p, s1, len(sols), dict(sorted(fam.items()))))
        sys.stdout.flush()

    for p in (7, 13, 19, 31, 37, 43):
        k = (p - 1) // 6
        rho = p // 6
        s1 = 2 * rho + 1
        # census scope ov in {0,1,2}; row-condition scope ov in {0..4}
        sols = census(p, k, 1, s1, {0, 1, 2})
        fam = Counter()
        for cfg in sols:
            (_, _, _), (_, d2, _), (_, d3, _) = cfg
            fam[(pmrep(d2, p), pmrep(d3, p))] += 1
        sols4 = census(p, k, 1, s1, {3, 4})
        fam4 = Counter()
        for cfg in sols4:
            (_, _, _), (_, d2, _), (_, d3, _) = cfg
            fam4[(pmrep(d2, p), pmrep(d3, p))] += 1
        res['%d' % p] = {'s1': s1, 'n': len(sols),
                         'families_ov012': {str(t): c for t, c in
                                            sorted(fam.items())},
                         'n_ov34': len(sols4),
                         'families_ov34': {str(t): c for t, c in
                                           sorted(fam4.items())}}
        print('[ROWCEN] row p=%2d (1 mod 6): s1=%d, ov0-2: %d coverings '
              'families %s; ov3-4: %d coverings families %s'
              % (p, s1, len(sols), dict(sorted(fam.items())), len(sols4),
                 dict(sorted(fam4.items()))))
        sys.stdout.flush()
    return res


# ----------------------------------------------------------------------
# [CELL112] + [CELL13U]
# ----------------------------------------------------------------------

def cells():
    res = {}
    Ns = [(5, 7), (5, 11), (5, 13), (7, 11), (5, 17), (7, 13), (5, 19),
          (5, 23), (7, 17), (11, 13)]
    for (p, q) in Ns:
        W = PQ(p, q)
        N = W.N
        # unit bad sets as bitmasks
        U = {}
        for u in W.units_pm:
            _, mask = W.B_unit(u)
            m = 0
            for k in np.nonzero(mask)[0]:
                m |= 1 << int(k)
            U[u] = m
        us = sorted(U)
        full = (1 << N) - 1

        # ---- (1,1,2): kernel-p (a) x kernel-q (b) x unit pairs
        best112 = (None, N + 1)   # (family, uncovered)
        n_cov = 0
        n_checked = 0
        for a in [a for a in range(1, q) if gcd(a, q) == 1 and a <= q - a]:
            Bp = W.B_kernelp(a)
            mp = 0
            for k in np.nonzero(Bp)[0]:
                mp |= 1 << int(k)
            for b in [b for b in range(1, p) if gcd(b, p) == 1 and b <= p - b]:
                Bq = W.B_kernelq(b)
                mq = 0
                for k in np.nonzero(Bq)[0]:
                    mq |= 1 << int(k)
                R = full & ~(mp | mq)
                r = R
                for i1 in range(len(us)):
                    R1 = r & ~U[us[i1]]
                    if R1 == 0:
                        # single unit covers the rectangle
                        n_cov += 1
                        continue
                    for i2 in range(i1 + 1, len(us)):
                        n_checked += 1
                        miss = R1 & ~U[us[i2]]
                        if miss == 0:
                            n_cov += 1
                        else:
                            cnt = miss.bit_count()
                            if cnt < best112[1]:
                                best112 = ((a, b, us[i1], us[i2]), cnt)
        # ---- (1K,3U): both orientations
        best13 = (None, N + 1)
        n_cov13 = 0
        n_checked13 = 0
        n_distinct_killed = 0
        for orient in ('p', 'q'):
            kern = range(1, q if orient == 'p' else p)
            kk = q if orient == 'p' else p
            for a in [a for a in kern if gcd(a, kk) == 1 and a <= kk - a]:
                if orient == 'p':
                    Mk = W.B_kernelp(a)
                else:
                    Mk = W.B_kernelq(a)
                    kk = p
                m = 0
                for k in np.nonzero(Mk)[0]:
                    m |= 1 << int(k)
                R = full & ~m
                nu = len(us)
                for i1 in range(nu):
                    R1 = R & ~U[us[i1]]
                    if R1 == 0:
                        n_cov13 += 1
                        continue
                    for i2 in range(i1 + 1, nu):
                        R2 = R1 & ~U[us[i2]]
                        if R2 == 0:
                            n_cov13 += 1
                            continue
                        for i3 in range(i2 + 1, nu):
                            n_checked13 += 1
                            miss = R2 & ~U[us[i3]]
                            if miss == 0:
                                n_cov13 += 1
                            else:
                                cnt = miss.bit_count()
                                if cnt < best13[1]:
                                    best13 = ((orient, a, us[i1], us[i2],
                                               us[i3]), cnt)
        # step-class statistics of the closest families
        res['%d' % N] = {
            'p': p, 'q': q, 'c': W.c, 'rho_p': W.rho_p, 'rho_q': W.rho_q,
            'cell112': {'coverings': n_cov, 'checked_pairs': n_checked,
                        'min_uncovered': best112[1],
                        'argmin': best112[0]},
            'cell13U': {'coverings': n_cov13, 'checked_triples':
                        n_checked13, 'min_uncovered': best13[1],
                        'argmin': best13[0]},
        }
        print('[CELL] N=%3d: (1,1,2) %d coverings / %d pairs checked, '
              'min uncover %d; (1K,3U) %d coverings / %d triples, '
              'min uncover %d'
              % (N, n_cov, n_checked, best112[1], n_cov13, n_checked13,
                 best13[1]))
        sys.stdout.flush()
    return res


# ----------------------------------------------------------------------
# [SHEAR] per-row analysis of the closest (1,1,2) families
# ----------------------------------------------------------------------

def shear():
    res = {}
    for (p, q) in [(5, 7), (5, 11), (7, 11), (7, 13)]:
        W = PQ(p, q)
        Xc_size = p - 2 * W.rho_p - 1
        # find the closest families by margin over a reduced grid
        best = (None, p * q + 1)
        for a in [a for a in range(1, q) if gcd(a, q) == 1 and a <= q - a]:
            Y = W.rowset(a)
            Yc = [y for y in range(q) if not Y[y]]
            for b in [b for b in range(1, p) if gcd(b, p) == 1 and b <= p - b]:
                X = W.colset(b)
                Xc = [x for x in range(p) if not X[x]]
                for i1 in range(len(W.units_pm)):
                    u1 = W.units_pm[i1]
                    for i2 in range(i1 + 1, len(W.units_pm)):
                        u2 = W.units_pm[i2]
                        tot = 0
                        for y in Yc:
                            f1 = set(W.footprint(u1, y))
                            f2 = set(W.footprint(u2, y))
                            tot += len([x for x in Xc if x not in f1
                                        and x not in f2])
                        if tot < best[1]:
                            best = ((a, b, u1, u2), tot)
        (a, b, u1, u2), tot = best
        Y = W.rowset(a)
        X = W.colset(b)
        Xc = [x for x in range(p) if not X[x]]
        rows = []
        d1 = (pow(u1, -1, p) * q) % p
        d2 = (pow(u2, -1, p) * q) % p
        for y in range(q):
            if Y[y]:
                continue
            f1 = W.footprint(u1, y)
            f2 = W.footprint(u2, y)
            miss = len([x for x in Xc if x not in f1 and x not in f2])
            rows.append({'y': y, 'n1': len(f1), 'n2': len(f2),
                         'miss': miss,
                         'steps_pm': (pmrep(d1, p), pmrep(d2, p))})
        res['%d' % W.N] = {
            'argmin': {'a': a, 'b': b, 'u1': u1, 'u2': u2},
            'total_miss': tot, 'Xc_size': Xc_size,
            'steps_pmrep': [pmrep(d1, p), pmrep(d2, p)],
            'rows': rows,
        }
        ngood = sum(1 for r in rows if r['miss'] == 0)
        print('[SHEAR] N=%3d argmin a=%d b=%d u=(%d,%d) steps=%s: '
              'total miss %d, good rows %d/%d'
              % (W.N, a, b, u1, u2, (pmrep(d1, p), pmrep(d2, p)), tot,
                 ngood, len(rows)))
        sys.stdout.flush()
    return res


def main():
    print('=== pq side: tensor translation verification ===')
    OUT['const_struct'] = const_and_struct()
    OUT['row_census'] = row_census()
    OUT['cells'] = cells()
    OUT['shear'] = shear()
    with open('out_pq_side.json', 'w') as f:
        json.dump(OUT, f, indent=1, default=str)
    print('written: out_pq_side.json')


if __name__ == '__main__':
    main()
