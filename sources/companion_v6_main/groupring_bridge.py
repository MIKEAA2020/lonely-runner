#!/usr/bin/env python3
"""
Group-ring bridge verification (Priority 1-2 of the reviewer's programme).

THEOREM GT (proved on paper, machine-checked here).  For ANY three APs
A_i = {a_i + j d_i mod p, 0 <= j < s_i} in Z_p with d_i != 0:

    Q(z) := sum_i z^{a_i} (1 - z^{d_i s_i}) prod_{j != i} (1 - z^{d_j})
    P(z) := prod_i (1 - z^{d_i})
    E(z) := sum_x ( counts(x) - 1 ) z^x        (signed excess)

satisfy   Q = E * P   in  R = Z[z]/(z^p - 1)   (the 24-term identity),
UNCONDITIONALLY.  Covering of Z_p  <=>  E >= 0 pointwise; with
sum s_i - p in {0,1} covering <=> E in {0, z^e}.  The solution of Q = X P
with the right mass X(1) = sum s_i - p is unique (kernel of mult-by-P is
Z.Phi_p, and Phi_p has mass p).

Experiments (independent code path from the prior harness ap_lemmas_verify):
  [CTRL-A] the 52 +-distinct exceptions (Lemma A failures at p=7: 50, p=13: 2,
           reconstructed from out_ap_lemmas.json) satisfy the identity.
  [CENSUS] exhaustive census at Lemma B primes (p=5 mod 6, 11..59) and
           Lemma A primes (p=1 mod 6, 19..61) over ALL (d2,d3) in (Z_p^*)^2:
           which +-rep families (pmrep(d2),pmrep(d3)) cover; identity check
           on every covering found; per-family per-multiset counts.
  [NEG]    +-distinct exhaustive re-verification (0 coverings) same primes.
  [SAMPLE] random configs: sharp form of the identity (signed E always),
           covering <=> E >= 0, kernel contains Z.Phi_p.
  [DEFECT] distance-to-solvability statistics on the +-distinct side.

Output: scripts/out_groupring_bridge.json + stdout log.
"""
import json
import random
import sys
from collections import Counter, defaultdict

import numpy as np

OUT = {}


# ----------------------------------------------------------------------
# core helpers
# ----------------------------------------------------------------------

def pmrep(d, p):
    r = d % p
    return min(r, p - r)


def ap_points(a, d, s, p):
    return (a + d * np.arange(s)) % p


def counts_of(cfg, p):
    c = np.zeros(p, dtype=np.int64)
    for (a, d, s) in cfg:
        c[ap_points(a, d, s, p)] += 1
    return c


def Q_array(cfg, p):
    (a1, d1, s1), (a2, d2, s2), (a3, d3, s3) = cfg
    Q = np.zeros(p, dtype=np.int64)
    for i in range(3):
        ai, di, si = cfg[i]
        dj = cfg[(i + 1) % 3][1]
        dk = cfg[(i + 2) % 3][1]
        cip = (di * si) % p
        for b1 in (0, 1):
            for b2 in (0, 1):
                for b3 in (0, 1):
                    e = (ai + b1 * cip + b2 * dj + b3 * dk) % p
                    Q[e] += (-1) ** (b1 + b2 + b3)
    return Q


def P_array(cfg, p):
    D = np.zeros(p, dtype=np.int64)
    D[0] = 1
    for (a, d, s) in cfg:
        D = D - np.roll(D, d)
    return D


def circ_mul(E, P, p):
    supp = np.nonzero(E)[0]
    acc = np.zeros(p, dtype=np.int64)
    for y in supp:
        acc += E[y] * np.roll(P, int(y))
    return acc


def check_identity(cfg, p):
    c = counts_of(cfg, p)
    E = c - 1
    Q = Q_array(cfg, p)
    P = P_array(cfg, p)
    ok = np.array_equal(Q, circ_mul(E, P, p))
    return ok, E, Q, P


def covering(cfg, p):
    return bool(counts_of(cfg, p).min() >= 1)


# ----------------------------------------------------------------------
# bitmask + residual-fit enumeration (exhaustive over all (d2,d3))
# ----------------------------------------------------------------------

def critical_ordered_sizes(k, cls):
    p = 6 * k + cls
    if cls == 5:
        allowed = (2 * k + 1, 2 * k + 2)
        ovs = (0, 1)
    else:
        allowed = (2 * k, 2 * k + 1, 2 * k + 2)
        ovs = (0, 1, 2)
    out = []
    for s1 in allowed:
        for s2 in allowed:
            for s3 in allowed:
                if (s1 + s2 + s3 - p) in ovs:
                    out.append((s1, s2, s3))
    return sorted(set(out))


def fit_placements(R, p, s3, d3, inv):
    """R = sorted list of residues to be covered by A3 = a3 + [0,s3) d3.
    Returns the list of valid a3.  A window [t, t+s3) in inv-coordinates
    must contain inv*R; the minimal cyclic arc containing inv*R has length
    l = p - maxgap; valid window starts t in [arcend - s3, arcstart]
    (cyclically), i.e. s3 - l + 1 placements if l <= s3."""
    if not R:
        return list(range(p))  # any placement covers empty R (still coverings
        # require nothing else) -- but s3>=1 always has points; R empty means
        # A1 u A2 = Z_p already, impossible when s3 >= 1 unless overlap huge
    v = sorted((inv * r) % p for r in R)
    m = len(v)
    if m > s3:
        return []
    # cyclic gaps
    gaps = [(v[(i + 1) % m] - v[i]) % p for i in range(m)]
    gmax = max(gaps)
    if m == 1:
        arc_len = 1  # single point
        arc_start = v[0]
    else:
        gi = gaps.index(gmax)
        arc_start = v[(gi + 1) % m]  # right after the largest gap
        # residues in the minimal covering arc (complement of the open gap)
        arc_len = p - gmax + 1
    if arc_len > s3:
        return []
    out = []
    nplace = s3 - arc_len + 1
    for u in range(nplace):
        t = (arc_start - u) % p
        out.append((t * d3) % p)
    return out


def enumerate_coverings_fast(p, k, cls, mode='census'):
    """All normalized coverings (a1=0, d1=1), sizes critical, over:
    mode='census'  : all d2, d3 in Z_p^*
    mode='distinct': d2, d3 with pairwise +-distinctness from 1 and each other
    Returns list of configs ((0,1,s1),(a2,d2,s2),(a3,d3,s3))."""
    sizes = critical_ordered_sizes(k, cls)
    sols = []
    full = (1 << p) - 1
    inv = [0] * p
    for d in range(1, p):
        inv[d] = pow(d, -1, p)
    for (s1, s2, s3) in sizes:
        ov = s1 + s2 + s3 - p
        m1 = 0
        for j in range(s1):
            m1 |= 1 << j
        for d2 in range(1, p):
            if mode == 'distinct' and pmrep(d2, p) == 1:
                continue
            # A2 masks for all a2
            base2 = 0
            for j in range(s2):
                base2 |= 1 << ((d2 * j) % p)
            for a2 in range(p):
                sh = a2
                m2 = ((base2 << sh) | (base2 >> (p - sh))) & full \
                    if sh else base2
                h = (m1 & m2).bit_count()
                if h > ov:
                    continue
                m12 = m1 | m2
                Rmask = full & ~m12
                R = [i for i in range(p) if (Rmask >> i) & 1]
                if len(R) > s3:
                    continue
                for d3 in range(1, p):
                    if mode == 'distinct' and (
                            pmrep(d3, p) == 1
                            or pmrep(d3, p) == pmrep(d2, p)):
                        continue
                    for a3 in fit_placements(R, p, s3, d3, inv[d3]):
                        sols.append(((0, 1, s1), (a2, d2, s2),
                                     (a3, d3, s3)))
    return sols


# ----------------------------------------------------------------------
# experiments
# ----------------------------------------------------------------------

def ctrl_A_exceptions():
    res = {}
    prior = json.load(open('out_ap_lemmas.json'))
    for pk in ('7', '13'):
        sols = prior['lemmaA'][pk]['solutions']
        p = int(pk)
        ok_n = 0
        ov_counter = Counter()
        prof_counter = Counter()
        for s in sols:
            s1, s2, s3 = s['assign']
            cfg = ((0, 1, s1), (s['a2'], s['d2'], s2), (s['a3'], s['d3'], s3))
            ok, E, Q, P = check_identity(cfg, p)
            assert E.min() >= 0 and covering(cfg, p)
            ov_counter[int(E.sum())] += 1
            prof_counter[(int(E.sum()), int((E != 0).sum()),
                          int(E.max()))] += 1
            ok_n += int(ok)
            if not ok:
                print('IDENTITY FAIL', pk, s)
        res[pk] = {'n': len(sols), 'identity_ok': ok_n,
                   'ov_distribution': dict(ov_counter),
                   'E_profiles_(sum,support,max)': {
                       str(t): c for t, c in prof_counter.items()}}
        print('[CTRL-A] p=%s: %d configs, identity holds on %d; '
              'ov dist %s; (sum,support,max) %s' %
              (pk, len(sols), ok_n, dict(ov_counter), dict(prof_counter)))
    return res


def census_and_negatives():
    res = {}
    primes_B = [11, 17, 23, 29, 41, 47, 53, 59]
    primes_A = [19, 31, 37, 43, 61]
    for p in primes_B + primes_A:
        k = (p - 5) // 6 if p % 6 == 5 else (p - 1) // 6
        cls = 5 if p % 6 == 5 else 1
        sols = enumerate_coverings_fast(p, k, cls, mode='census')
        sols_d = enumerate_coverings_fast(p, k, cls, mode='distinct')
        fam = Counter()
        fam_ms = Counter()
        idok = ncheck = 0
        for cfg in sols:
            (_, _, s1), (_, d2, s2), (_, d3, s3) = cfg
            fam[(pmrep(d2, p), pmrep(d3, p))] += 1
            fam_ms[(pmrep(d2, p), pmrep(d3, p),
                    tuple(sorted((s1, s2, s3))))] += 1
        for cfg in sols[:3000]:
            ok, E, Q, P = check_identity(cfg, p)
            assert E.min() >= 0 and covering(cfg, p)
            idok += int(ok)
            ncheck += 1
            if not ok:
                print('IDENTITY FAIL (census)', p, cfg)
        # per-multiset family counts for the p-independence table
        ms_table = defaultdict(dict)
        for (f2, f3, ms), c in fam_ms.items():
            ms_table[str(ms)][str((f2, f3))] = c
        res[p] = {'k': k, 'cls': cls,
                  'census_total': len(sols),
                  'families_(pmrep_d2,pmrep_d3)': {str(t): c
                                                   for t, c in
                                                   sorted(fam.items())},
                  'per_multiset': dict(ms_table),
                  'identity_checked': ncheck, 'identity_ok': idok,
                  'distinct_coverings': len(sols_d)}
        print('[CENSUS] p=%2d k=%2d cls=%d: total=%d families=%s '
              'distinct=%d identity %d/%d' %
              (p, k, cls, len(sols), dict(sorted(fam.items())), len(sols_d),
               idok, ncheck))
        sys.stdout.flush()
    return res


def sample_checks():
    rng = random.Random(20261004)
    res = {}
    phi = {}
    for p in (11, 13, 17, 19, 23, 29, 31, 37):
        k = (p - 5) // 6 if p % 6 == 5 else (p - 1) // 6
        cls = 5 if p % 6 == 5 else 1
        sizes = critical_ordered_sizes(k, cls)
        ones = np.ones(p, dtype=np.int64)
        n_id = n_agree = n_ker = 0
        N = 3000
        for _ in range(N):
            s1, s2, s3 = rng.choice(sizes)
            d2 = rng.randrange(1, p)
            d3 = rng.randrange(1, p)
            cfg = ((0, 1, s1), (rng.randrange(p), d2, s2),
                   (rng.randrange(p), d3, s3))
            ok, E, Q, P = check_identity(cfg, p)
            n_id += int(ok)
            if not ok:
                print('IDENTITY FAIL (sample)', p, cfg)
            if (E.min() >= 0) == covering(cfg, p):
                n_agree += 1
            else:
                print('SIGN DISAGREEMENT', p, cfg)
            c = rng.choice([-1, 1])
            if np.array_equal(circ_mul(E + c * ones, P, p), Q):
                n_ker += 1
        res[p] = {'n': N, 'identity_ok': n_id, 'sign_agree': n_agree,
                  'kernel_ok': n_ker}
        print('[SAMPLE] p=%2d: identity %d/%d, sign<->covering %d/%d, '
              'kernel(Z.Phi_p) %d/%d' % (p, n_id, N, n_agree, N, n_ker, N))
        sys.stdout.flush()
    return res


def defect_stats():
    rng = random.Random(777)
    res = {}
    for p in (11, 17, 23, 29, 41):
        k = (p - 5) // 6 if p % 6 == 5 else (p - 1) // 6
        cls = 5 if p % 6 == 5 else 1
        sizes = critical_ordered_sizes(k, cls)
        # distinct pairs list (for sampling)
        pairs = []
        for d2 in range(2, p):
            if pmrep(d2, p) == 1:
                continue
            for d3 in range(2, p):
                if pmrep(d3, p) == 1 or pmrep(d3, p) == pmrep(d2, p):
                    continue
                pairs.append((d2, d3))
        best, best_cfg = None, None
        vals = []
        N = 3000
        for _ in range(N):
            s1, s2, s3 = rng.choice(sizes)
            d2, d3 = rng.choice(pairs)
            cfg = ((0, 1, s1), (rng.randrange(p), d2, s2),
                   (rng.randrange(p), d3, s3))
            if covering(cfg, p):
                continue
            Q = Q_array(cfg, p)
            if s1 + s2 + s3 == p:
                dv = int(np.count_nonzero(Q))
            else:
                P = P_array(cfg, p)
                C = np.stack([np.roll(P, e) for e in range(p)])
                dv = int(min(np.count_nonzero(Q - C[e]) for e in range(p)))
            vals.append(dv)
            if best is None or dv < best:
                best, best_cfg = dv, (s1, s2, s3, d2, d3)
        vals.sort()
        res[p] = {'n': len(vals), 'min': best,
                  'p05': vals[max(0, len(vals) // 20)],
                  'median': vals[len(vals) // 2],
                  'argmin': best_cfg}
        print('[DEFECT] p=%2d: n=%d min=%d p05=%d median=%d argmin=%s' %
              (p, len(vals), best, vals[max(0, len(vals) // 20)],
               vals[len(vals) // 2], best_cfg))
        sys.stdout.flush()
    return res


def main():
    print('=== Group-ring bridge verification ===')
    OUT['ctrl_A'] = ctrl_A_exceptions()
    OUT['census'] = census_and_negatives()
    OUT['sample'] = sample_checks()
    OUT['defect'] = defect_stats()
    with open('out_groupring_bridge.json', 'w') as f:
        json.dump(OUT, f, indent=1, default=str)
    print('written: out_groupring_bridge.json')


if __name__ == '__main__':
    main()
