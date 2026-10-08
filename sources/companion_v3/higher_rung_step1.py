#!/usr/bin/env python3
"""
Higher-rung step 1: the general-T fiber machinery + the T=7 / T=8-critical
fiber censuses (4-difference level).  Targets (a)/(b) substrate.

SETTINGS.  Rung T (= n+1): composite modulus N = p^2, c = floor((N-1)/T),
bad sets B_w = {x in Z_N : ||w x||_N <= c}.  Kernel member pa: covers
complete fibers; the fiber ball B_k = {r : min(r,p-r) <= k}, k = (p-eps)/T,
|B_k| = 2k+1, |U| = p - 2k - 1 = (T-2)k + eps - 1.
Unit u on fiber F_j = {j + pt}: B_u cap F_j in the t-coordinate is the AP
with difference ubar = u^{-1} mod p, size L(a0) where a0 = uj mod p:

    L(a0) = 2k+1  if a0 in B_k ("on"),
    L(a0) = 2k+2  if a0 off-ball and 2(k*eps+d0) >= p-1   ("off+", eps > T/2)
    L(a0) = 2k    if a0 off-ball and 2(k*eps+d0) <= p-2   ("off-", eps < T/2)
    d0 = floor((eps^2-1)/T).

(This reproduces the committed T=6 tables of the colliding-units note,
boundary for boundary, both eps classes -- checked below.)

FIBER CENSUS (the object of targets (a)/(b)).  On an uncovered fiber the
T-2-j unit footprints must cover Z_p.  At T=7 the frontier cell is
(1K,4U): four footprints, sizes in {on, off}, structural overlap ~k
(linear slack -- the "growing pinned sets" regime).  At T=8 the frontier
cells are (2K,4U) -- four footprints, overlap O(1): the CRITICAL 4-difference
census, the tight analogue of the T=6 three-difference census -- and
(1K,5U) (five footprints, linear slack; deferred to step 2).

Normal form: a1 = 0, d1 = 1 (the normalized unit's footprint is an
interval), remaining differences canonical in [1,(p-1)/2] (an AP-set has a
unique +-canonical difference).  Classification key: the ordered tuple
(d2,...,dm) itself (pmrep = d for canonical d).

Checks:
  [M-CTRL] T=6 arc/size table vs the colliding note, p in {11,13,17,19,23,
           29,37}, both eps classes; plus direct Z_{p^2} footprint
           verification (t-set formula vs brute force) at T=6 and T=7.
  [M-GEN]  on = B_k exactly, at T=7/T=8 for every eps class (all a0).
  [C7]     T=7 4-difference census at 16 primes (all six eps classes):
           per-(d2,d3,d4) placement counts; family lists; p-independence;
           the +-distinct question (SL clause (i) analogue).
  [C8]     T=8-critical 4-difference census ((2K,4U) fiber problem) at
           8 primes (eps in {1,3,5,7}): same record.

Output: scripts/out_higher_rung_step1.json + stdout log.
"""
import json
import sys
import time
from collections import Counter, defaultdict

OUT = {}


# ----------------------------------------------------------------------
# general-T fiber machinery
# ----------------------------------------------------------------------

def eps_of(p, T):
    return p % T


def k_of(p, T):
    return (p - p % T) // T


def ball_mask(p, k):
    """{r : min(r, p-r) <= k} as a python int bitmask."""
    m = 0
    for r in list(range(0, k + 1)) + list(range(p - k, p)):
        m |= 1 << r
    return m


def value_arc(p, T, k, eps, a0):
    """V(a0) = {s in Z_p : min(a0+p s, p^2-a0-p s) <= c}, c=floor((p^2-1)/T).
    Returns (set_of_s, L)."""
    c = (p * p - 1) // T
    good = []
    for s in range(p):
        v = a0 + p * s
        if min(v, p * p - v) <= c:
            good.append(s)
    return set(good), len(good)


def arc_three_tier(p, T, k, eps, a0):
    """Three-tier arc law (derived; machine-verified below).
    M = k*eps + d0, d0 = floor((eps^2-1)/T).
    L = 2k+1  on the ball B_k AND on the shallow off bands;
    L = 2k+2  on the deep band (p-1-M, M]      when 2M >= p-1 (off+ class);
    L = 2k    on the deep band (M, p-1-M]      when 2M <= p-2 (off- class).
    At eps = +-1 mod T the deep band IS the whole off-ball (the T=6
    phenomenon); at intermediate eps a shallow band of width |M - k|
    appears with the on-size."""
    d0 = (eps * eps - 1) // T
    M = k * eps + d0
    if a0 <= min(M, p - 1 - M) or a0 > max(M, p - 1 - M):
        return 2 * k + 1
    if 2 * M >= p - 1:
        return 2 * k + 2
    return 2 * k


def machinery_control():
    """[M-CTRL] + [M-GEN]."""
    res = {}
    # --- T=6 control against the committed note ---
    ok_all = True
    for p in (11, 13, 17, 19, 23, 29, 37):
        T = 6
        eps = p % 6
        k = (p - eps) // 6
        # note's table: on = [1,k] U [p-k,p-1] plus 0; off in between
        bad = []
        for a0 in range(p):
            _, L = value_arc(p, T, k, eps, a0)
            Lf = arc_three_tier(p, T, k, eps, a0)
            inball = min(a0, p - a0) <= k
            # expected from the note
            if eps == 5:
                exp = 2 * k + 1 if inball else 2 * k + 2
            else:
                exp = 2 * k + 1 if inball else 2 * k
            if L != exp or Lf != L:
                bad.append((a0, L, Lf, exp, inball))
        res['T6_p%d' % p] = {'eps': eps, 'k': k, 'bad': bad}
        ok_all &= not bad
    print('[M-CTRL] T=6 arc table: %s' % ('ALL PASS' if ok_all else 'FAIL'))
    # --- on = B_k at T=7 / T=8, every eps class ---
    ok_gen = True
    for T, plist in ((7, (11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)),
                     (8, (17, 19, 29, 31, 41, 43, 47))):
        for p in plist:
            eps = p % T
            k = (p - eps) // T
            bad = []
            for a0 in range(p):
                _, L = value_arc(p, T, k, eps, a0)
                exp = arc_three_tier(p, T, k, eps, a0)
                if L != exp:
                    bad.append((a0, L, exp))
            res['T%d_p%d' % (T, p)] = {'eps': eps, 'k': k, 'bad': bad}
            ok_gen &= not bad
    print('[M-GEN] three-tier arc law at T=7/T=8: %s' %
          ('ALL PASS' if ok_gen else 'FAIL'))
    res['m_ctrl_pass'] = bool(ok_all)
    res['m_gen_pass'] = bool(ok_gen)
    return res


def footprint_direct(p, T, u, j):
    """B_u cap F_j by brute force over Z_{p^2}: returns the t-set."""
    c = (p * p - 1) // T
    N = p * p
    ts = set()
    for t in range(p):
        x = j + p * t
        wx = (u * x) % N
        d = min(wx, N - wx)
        if d <= c:
            ts.add(t)
    return ts


def footprint_formula(p, T, u, j):
    """t-set = ubar * (V(a0) - a1) mod p, a0 = uj mod p, a1 = floor(uj/p)."""
    ubar = pow(u, -1, p)
    a0 = (u * j) % p
    a1 = (u * j) // p
    V, _ = value_arc(p, T, (p - p % T) // T, p % T, a0)
    return {(ubar * (s - a1)) % p for s in V}


def footprint_control():
    """Direct vs formula footprints at T=6 (control) and T=7/T=8."""
    res = {}
    nbad = 0
    ntot = 0
    for T, plist in ((6, (11, 13, 17, 19)), (7, (29, 43)), (8, (17, 41))):
        for p in plist:
            for u in (1, 2, 3, 5, 7, 11, p - 2, p - 3):
                if u % p == 0:
                    continue
                for j in range(0, p, max(1, p // 12)):
                    d = footprint_direct(p, T, u, j)
                    f = footprint_formula(p, T, u, j)
                    ntot += 1
                    if d != f:
                        nbad += 1
                        print('FOOTPRINT MISMATCH T=%d p=%d u=%d j=%d' %
                              (T, p, u, j))
    res = {'checks': ntot, 'bad': nbad}
    print('[M-FOOT] footprint formula vs direct: %d/%d pass' %
          (ntot - nbad, ntot))
    return res


# ----------------------------------------------------------------------
# m-AP covering census (normal form), bitmask + arc-fit cascade
# ----------------------------------------------------------------------

def profile_sizes(p, T):
    """(on, off) sizes for the rung; also the class label."""
    eps = p % T
    k = (p - eps) // T
    d0 = (eps * eps - 1) // T
    M = k * eps + d0
    if 2 * M >= p - 1:
        return k, (2 * k + 1, 2 * k + 2), 'off+', eps
    return k, (2 * k, 2 * k + 1), 'off-', eps


def ap_base(d, s, p):
    """mask of {d*j mod p : j < s}."""
    m = 0
    for j in range(s):
        m |= 1 << ((d * j) % p)
    return m


def rot(base, a, p, full):
    if a == 0:
        return base
    return ((base << a) | (base >> (p - a))) & full


def fit_placements(R, p, s, d, inv):
    """Valid starts a for the AP a + [0,s) d containing all of R
    (R = sorted residues).  Returns list of a."""
    if not R:
        return list(range(p))
    v = sorted((inv * r) % p for r in R)
    m = len(v)
    if m > s:
        return []
    if m == 1:
        arc_len, arc_start = 1, v[0]
    else:
        gaps = [(v[(i + 1) % m] - v[i]) % p for i in range(m)]
        gmax = max(gaps)
        gi = gaps.index(gmax)
        arc_start = v[(gi + 1) % m]
        arc_len = p - gmax + 1
    if arc_len > s:
        return []
    out = []
    for u in range(s - arc_len + 1):
        t = (arc_start - u) % p
        out.append((t * d) % p)
    return out


def size_tuples(m, on, off, p):
    """ordered m-tuples of sizes from {on, off} with sum >= p
    (necessary for covering)."""
    out = []
    for code in range(1 << m):
        sz = tuple(off if (code >> i) & 1 else on for i in range(m))
        if sum(sz) >= p:
            out.append(sz)
    return out


def census_m4(p, T, verbose=True, projections=True):
    """4-AP covering census in normal form (0,1,s1),(a2,d2,s2),(a3,d3,s3),
    (a4,d4,s4).  d2,d3,d4 canonical in [1,(p-1)/2].  COMPACT: per
    (s1,s2,s3,s4,d2,d3,d4) records count + projection sets (a2, a3, a4,
    a3-a2, a4-a3, a4-a2) -- exactly the pinning measurements."""
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    full = (1 << p) - 1
    inv = [0] * p
    for d in range(1, p):
        inv[d] = pow(d, -1, p)
    dcan = list(range(1, (p + 1) // 2))
    fams = defaultdict(lambda: {'n': 0})
    t0 = time.time()
    n_inner = 0
    n_cfg = 0

    def rec(key, a2, a3, a4):
        f = fams[key]
        f['n'] += 1
        if projections:
            for nm, v in (('a2', a2), ('a3', a3), ('a4', a4),
                          ('d32', (a3 - a2) % p), ('d43', (a4 - a3) % p),
                          ('d42', (a4 - a2) % p)):
                if nm not in f:
                    f[nm] = set()
                f[nm].add(v)

    for (s1, s2, s3, s4) in size_tuples(4, on, off, p):
        ov = s1 + s2 + s3 + s4 - p
        m1 = (1 << s1) - 1
        base2 = {d: ap_base(d, s2, p) for d in dcan}
        masks3 = {d: [rot(ap_base(d, s3, p), a, p, full)
                      for a in range(p)] for d in dcan}
        for d2 in dcan:
            b2 = base2[d2]
            for a2 in range(p):
                m2 = rot(b2, a2, p, full)
                if (m1 & m2).bit_count() > ov:
                    continue
                m12 = m1 | m2
                if (full & ~m12).bit_count() > s3 + s4:
                    continue
                for d3 in dcan:
                    row = masks3[d3]
                    for a3 in range(p):
                        n_inner += 1
                        m123 = m12 | row[a3]
                        if m123.bit_count() < p - s4:
                            continue
                        if (m12 & row[a3]).bit_count() > ov:
                            continue
                        R = full & ~m123
                        nr = R.bit_count()
                        if nr == 0:
                            for d4 in dcan:
                                key = (s1, s2, s3, s4, d2, d3, d4)
                                for a4 in range(p):
                                    rec(key, a2, a3, a4)
                                    n_cfg += 1
                            continue
                        if nr > s4:
                            continue
                        Rl = [i for i in range(p) if (R >> i) & 1]
                        for d4 in dcan:
                            key = (s1, s2, s3, s4, d2, d3, d4)
                            for a4 in fit_placements(Rl, p, s4, d4,
                                                     inv[d4]):
                                rec(key, a2, a3, a4)
                                n_cfg += 1
    if verbose:
        print('  census p=%d T=%d eps=%d cls=%s: %d fam-size keys, %d '
              'configs, %d inner, %.1fs' %
              (p, T, eps, cls, len(fams), n_cfg, n_inner, time.time() - t0))
        sys.stdout.flush()
    out = {k: {'n': v['n'],
               **{nm: len(s) for nm, s in v.items() if nm != 'n'}}
           for k, v in fams.items()}
    return out, dict(k=k, on=on, off=off, cls=cls, eps=eps, n_cfg=n_cfg)


def census_m4_solutions(p, T, keys, verbose=False):
    """Recompute the full solution list for SELECTED (s1,s2,s3,s4,d2,d3,d4)
    keys (for the shear analysis).  Returns {key: [(a2,a3,a4), ...]}."""
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    full = (1 << p) - 1
    inv = [0] * p
    for d in range(1, p):
        inv[d] = pow(d, -1, p)
    dcan = list(range(1, (p + 1) // 2))
    want = set(keys)
    got = defaultdict(list)
    for key in want:
        (s1, s2, s3, s4, d2, d3, d4) = key
        ov = s1 + s2 + s3 + s4 - p
        m1 = (1 << s1) - 1
        m2 = rot(ap_base(d2, s2, p), 0, p, full)  # placeholder; loop below
        for a2 in range(p):
            m2 = rot(ap_base(d2, s2, p), a2, p, full)
            if (m1 & m2).bit_count() > ov:
                continue
            m12 = m1 | m2
            for a3 in range(p):
                m3 = rot(ap_base(d3, s3, p), a3, p, full)
                if (m12 | m3).bit_count() < p - s4:
                    continue
                if (m12 & m3).bit_count() > ov:
                    continue
                R = full & ~(m12 | m3)
                nr = R.bit_count()
                if nr == 0:
                    for a4 in range(p):
                        got[key].append((a2, a3, a4))
                    continue
                if nr > s4:
                    continue
                Rl = [i for i in range(p) if (R >> i) & 1]
                for a4 in fit_placements(Rl, p, s4, d4, inv[d4]):
                    got[key].append((a2, a3, a4))
    return dict(got)


# ----------------------------------------------------------------------
# experiments
# ----------------------------------------------------------------------

def family_stats(fams):
    """Aggregate compact census -> per-(d2,d3,d4) family stats."""
    pat = Counter()
    pat_distinct = Counter()
    maxproj = defaultdict(lambda: defaultdict(int))
    for (s1, s2, s3, s4, d2, d3, d4), st in fams.items():
        key = (d2, d3, d4)
        pat[key] += st['n']
        for nm in ('a2', 'a3', 'a4', 'd32', 'd43', 'd42'):
            if nm in st:
                maxproj[key][nm] = max(maxproj[key][nm], st[nm])
        if len({d2, d3, d4}) == 3 and 1 not in (d2, d3, d4):
            pat_distinct[key] += st['n']
    vals = Counter()
    for (d2, d3, d4) in pat:
        for d in (d2, d3, d4):
            vals[d] += 1
    return dict(pat=pat, pat_distinct=pat_distinct, vals=vals,
                maxproj=maxproj)


def run_census_batch(primes, T, m=4):
    res = {}
    for p in primes:
        fams, meta = census_m4(p, T)
        st = family_stats(fams)
        pat = st['pat']
        patd = st['pat_distinct']
        maxproj = st['maxproj']
        res['p%d' % p] = {
            'eps': meta['eps'], 'k': meta['k'], 'cls': meta['cls'],
            'n_families': len(pat),
            'n_configs': sum(pat.values()),
            'n_distinct_families': len(patd),
            'n_distinct_configs': sum(patd.values()),
            'd_value_inventory': {str(v): c for v, c in
                                  sorted(st['vals'].items())},
            'top_families': [[list(k), v] for k, v in pat.most_common(15)],
            'max_proj_per_family': {str(k): dict(v) for k, v in
                                    maxproj.items()},
            'families': {str(k): v for k, v in sorted(pat.items())},
            'distinct_families': {str(k): v for k, v in
                                  sorted(patd.items())},
            'famsize_detail': {str(k): {'n': v['n'],
                                        **{nm: s for nm, s in v.items()
                                           if nm != 'n'}}
                               for k, v in fams.items()},
        }
        print('  p=%3d eps=%d k=%d: %4d fams, %6d cfgs; distinct: %4d fams '
              '%6d cfgs; max proj sizes: %s' %
              (p, meta['eps'], meta['k'], len(pat), sum(pat.values()),
               len(patd), sum(patd.values()),
               {nm: max(mp[nm] for mp in maxproj.values() if nm in mp)
                for nm in ('a2', 'd32', 'd43', 'd42')}))
        sys.stdout.flush()
        # incremental flush
        OUT['partial'] = res
        with open('/home/z/my-project/scripts/out_higher_rung_step1.json',
                  'w') as f:
            json.dump(OUT, f, indent=1)
    return res


def run_T7_census():
    """[C7] the T=7 4-difference census (the (1K,4U) fiber problem)."""
    print('--- [C7] T=7 4-difference census (linear-slack regime) ---')
    primes = [29, 43, 71,           # eps=1
              23, 37, 79,           # eps=2
              17, 31, 59,           # eps=3
              11, 53,               # eps=4
              19, 47, 61,           # eps=5
              13, 41, 83]           # eps=6
    return run_census_batch(primes, 7, 'T7')


def run_T8critical_census():
    """[C8] the T=8-critical 4-difference census ((2K,4U) fiber problem)."""
    print('--- [C8] T=8-critical 4-difference census (O(1) overlap) ---')
    primes = [17, 41, 73,           # eps=1
              19, 43,               # eps=3
              29,                   # eps=5
              31, 47]               # eps=7
    return run_census_batch(primes, 8, 'T8c')


def census_m5(p, T, verbose=True):
    """5-AP covering census (the T=8 (1K,5U) fiber problem), COUNT-ONLY
    (memory-safe): per (s1,d2,d3,d4,d5) the config count."""
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    full = (1 << p) - 1
    inv = [0] * p
    for d in range(1, p):
        inv[d] = pow(d, -1, p)
    dcan = list(range(1, (p + 1) // 2))
    counts = Counter()
    t0 = time.time()
    n_inner = 0
    for (s1, s2, s3, s4, s5) in size_tuples(5, on, off, p):
        ov = s1 + s2 + s3 + s4 + s5 - p
        m1 = (1 << s1) - 1
        base2 = {d: ap_base(d, s2, p) for d in dcan}
        masks3 = {d: [rot(ap_base(d, s3, p), a, p, full)
                      for a in range(p)] for d in dcan}
        masks4 = {d: [rot(ap_base(d, s4, p), a, p, full)
                      for a in range(p)] for d in dcan}
        for d2 in dcan:
            b2 = base2[d2]
            for a2 in range(p):
                m2 = rot(b2, a2, p, full)
                if (m1 & m2).bit_count() > ov:
                    continue
                m12 = m1 | m2
                if (full & ~m12).bit_count() > s3 + s4 + s5:
                    continue
                for d3 in dcan:
                    row3 = masks3[d3]
                    for a3 in range(p):
                        m123 = m12 | row3[a3]
                        if m123.bit_count() < p - s4 - s5:
                            continue
                        if (m12 & row3[a3]).bit_count() > ov:
                            continue
                        for d4 in dcan:
                            row4 = masks4[d4]
                            for a4 in range(p):
                                n_inner += 1
                                m1234 = m123 | row4[a4]
                                if m1234.bit_count() < p - s5:
                                    continue
                                if (m123 & row4[a4]).bit_count() > ov:
                                    continue
                                R = full & ~m1234
                                nr = R.bit_count()
                                if nr == 0 or nr > s5:
                                    continue
                                Rl = [i for i in range(p) if (R >> i) & 1]
                                for d5 in dcan:
                                    n = len(fit_placements(Rl, p, s5, d5,
                                                           inv[d5]))
                                    if n:
                                        counts[(s1, d2, d3, d4, d5)] += n
    if verbose:
        print('  census m5 p=%d T=%d eps=%d cls=%s: %d families, %d configs,'
              ' %d inner, %.1fs' %
              (p, T, eps, cls, len(counts), sum(counts.values()), n_inner,
               time.time() - t0))
        sys.stdout.flush()
    return dict(counts), dict(k=k, on=on, off=off, cls=cls, eps=eps)


def run_T8_m5_census():
    """[C8b] the T=8 5-difference census ((1K,5U) fiber problem)."""
    print('--- [C8b] T=8 5-difference census (linear slack) ---')
    res = {}
    for p in (17, 19):
        counts, meta = census_m5(p, 8)
        pat = Counter()
        patd = Counter()
        for (s1, d2, d3, d4, d5), n in counts.items():
            key = (d2, d3, d4, d5)
            pat[key] += n
            if len({d2, d3, d4, d5}) == 4 and 1 not in (d2, d3, d4, d5):
                patd[key] += n
        vals = Counter()
        for (d2, d3, d4, d5) in pat:
            for d in (d2, d3, d4, d5):
                vals[d] += 1
        res['p%d' % p] = {
            'eps': meta['eps'], 'k': meta['k'], 'cls': meta['cls'],
            'n_families': len(pat), 'n_configs': sum(pat.values()),
            'n_distinct_families': len(patd),
            'n_distinct_configs': sum(patd.values()),
            'd_value_inventory': {str(v): c for v, c in
                                  sorted(vals.items())},
            'top_families': [[list(k), v] for k, v in pat.most_common(15)],
            'families': {str(k): v for k, v in sorted(pat.items())},
        }
        print('  p=%d: %d fams %d cfgs; distinct %d fams %d cfgs' %
              (p, len(pat), sum(pat.values()), len(patd),
               sum(patd.values())))
        sys.stdout.flush()
        OUT['partial'] = res
    return res


def main():
    t0 = time.time()
    print('=== higher_rung_step1: general-T machinery + T=7/T=8 censuses ===')
    OUT['machinery_T6_control'] = machinery_control()
    OUT['footprint_control'] = footprint_control()
    OUT['T7_census'] = run_T7_census()
    OUT['T8critical_census'] = run_T8critical_census()
    OUT['T8_m5_census'] = run_T8_m5_census()
    OUT['elapsed_s'] = time.time() - t0
    with open('/home/z/my-project/scripts/out_higher_rung_step1.json', 'w') as f:
        json.dump(OUT, f, indent=1)
    print('DONE in %.1fs; output out_higher_rung_step1.json' % OUT['elapsed_s'])


if __name__ == '__main__':
    main()
