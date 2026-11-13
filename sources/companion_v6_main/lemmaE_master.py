#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lemmaE_master.py -- machine verification of the Lemma E architecture
(the three-distance foot count; discharges Input E).

Setting: p = 6k+1 prime, m in {2k+1, 2k+2} (the two consumer set sizes;
the consumed form uses m = 2k+1 at arc lengths s1 in {2k+1, 2k+2}), d in
[2, 3k] (d = q = 3k included to see the endpoint exceptions).

X = S_d^(m) = {d*j mod p : j < m} sorted; cyclic gaps; G2 = max adjacent
pair sum = max over x in X of next(x) - prev(x).

STRUCTURE (verified at every (p, d, m), p <= STRUCT_P_MAX):
  A  Runs and columns.  q_j = floor(d*j/p); run w = {d*j - w*p} over the
     j with q_j = w; its column is c_w = (-w*p) mod d; the run starts at
     the BOTTOM c_w of its column (value in [0,d)), steps by d, and for
     w < W is FULL: height F(c) = m0 + [c < r] (m0 = floor(p/d),
     r = p mod d), top = p - d + c_{w+1}.  Run W is a prefix of height
     h_W with top H = d*(m-1) - W*p.  W = floor(d*(m-1)/p).
  B  Bands.  Band i (points x with x div d = i) has residue set R_i =
     C := {c_0..c_W} for i < min(h_W, m0), C_full := C minus c_W for
     h_W <= i <= m0-1, and top band T = (C_full cap [0,r)) plus c_W if
     h_W > m0.  The gap sequence of X is the concatenation of the cyclic
     gap blocks of R_0..R_{m0-1} (on Z_d, ending in the wrap gap
     d - max R_i) followed by the top block of T on [0,r) (ending in
     r - max T).
  C  Master identity.  G2(X) = max over the candidates
       Gamma2(C), Gamma2(C_full) [if C_full occurs as a band set],
       Gamma2(T on Z_r) [if |T| >= 2], span0, spanTop,
     where span0 = (r - max T) + firstgap(C)  [span around 0],
     spanTop = (d - max R_{m0-1}) + (firstgap(T) or r)  [span around
     m0*d], plus (W = 1 case) the singleton-band spans
     (d - max R_{i-1}) + d and 2d.
  D  Aux lemmas.  W-inequalities (m=2k+1: W >= d-2k+1 for d <= 3k-2,
     W = d-2k at d = 3k-1, W = d-2k-1 at d = 3k; m=2k+2: W >= d-2k+1
     for d <= 3k-2, W = d-2k at d in {3k-1, 3k}); Gamma2 <= n - |R| + 2
     for |R| >= 2 (pigeonhole + averaging); max C <= max T + sigma
     (sigma = d - r); r <= 2k-1 for d <= 3k-1; span0 <= Gamma2(C);
     spanTop <= d - W + 2 (W >= 2); singleton bands iff W = 1.
  E  Theorem.  p >= 31, d in [3, 3k-1] => G2 <= 2k+1 (both m).
     p in {13, 19} dirty (exhibits).  Endpoints: d = 2 gives G2 = 2k+3
     (m = 2k+1); d = 3k gives G2 = 2k+2 (both m); d = 3k-1 gives 6.
  F  Consumer equivalence.  min_a feet(A(a,d,s)) >= 2  <=>  G2 <= s1
     (sample, brute force); full viable-set cross-check against
     out_foot_census.json (p <= 419).
"""
import json
import sys
import time

import numpy as np

STRUCT_P_MAX = 500      # full structural verification (A-D) + identity (C)
DIRECT_P_MAX = 2000     # direct G2 theorem check (E)
FAST_P_MAX = 4000       # candidate-path theorem check (E), numpy
CENSUS_JSON = '/home/z/my-project/scripts/out_foot_census.json'
OUT_JSON = '/home/z/my-project/scripts/out_lemmaE.json'

FAILS = []
STATS = {}


def fail(tag, info):
    FAILS.append((tag, info))
    if len(FAILS) < 40:
        print('FAIL[%s]: %r' % (tag, info))
    elif len(FAILS) == 40:
        print('... suppressing further failure prints')


def primes_6k1(n):
    sieve = np.ones(n + 1, bool)
    sieve[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return [int(x) for x in np.nonzero(sieve)[0] if x % 6 == 1]


def gaps_of(S, p):
    g = [S[i + 1] - S[i] for i in range(len(S) - 1)]
    g.append(p - S[-1] + S[0])
    return g


def G2_of(S, p):
    g = gaps_of(S, p)
    n = len(g)
    return max(g[i] + g[(i + 1) % n] for i in range(n))


def gamma2_cyc(R, n):
    """Cyclic 2-step gap max of sorted residues R (R[0]=0) on Z_n.
    Singleton convention: 2n (a lone column point has both gaps = n
    when the band below is also a singleton; the exact per-band span is
    handled separately in the W=1 branch)."""
    if len(R) == 1:
        return 2 * n
    g = [R[i + 1] - R[i] for i in range(len(R) - 1)] + [n - R[-1]]
    return max(g[i] + g[(i + 1) % len(g)] for i in range(len(g)))


def gamma2_bound(R, n):
    """Aux lemma: Gamma2 <= n - |R| + 2 for |R| >= 2."""
    if len(R) < 2:
        return True
    return gamma2_cyc(R, n) <= n - len(R) + 2


def check_structure(p, k, m, d):
    """A-D verification for one (p, m, d). Appends to FAILS on error.
    Returns a result dict (for aggregation)."""
    S = sorted((d * j) % p for j in range(m))
    G2 = G2_of(S, p)
    W = (d * (m - 1)) // p
    r = p % d
    sig = d - r
    m0 = p // d
    res = dict(p=p, m=m, d=d, W=W, r=r, sig=sig, m0=m0, G2=G2)

    # ---------- W = 0 : the AP case ----------
    if W == 0:
        if S != [d * i for i in range(m)]:
            fail('W0-AP', (p, m, d))
        formula = max(2 * d, p - d * (m - 2))
        if G2 != formula:
            fail('W0-G2', (p, m, d, G2, formula))
        res['case'] = 'W0'
        STATS['W0'] = STATS.get('W0', 0) + 1
        return res

    # ---------- A: runs and columns ----------
    orbit = [(-w * p) % d for w in range(W + 1)]
    if orbit[0] != 0 or len(set(orbit)) != W + 1:
        fail('orbit-distinct', (p, m, d))
    C = sorted(orbit)
    C_full = sorted(orbit[:-1])
    cW = orbit[-1]
    H = d * (m - 1) - W * p
    h_W = None
    for w in range(W + 1):
        j0 = -((-w * p) // d)                 # ceil(w*p/d)
        j1 = min(-((-((w + 1) * p) // d)) - 1, m - 1)
        pts = [d * j - w * p for j in range(j0, j1 + 1)]
        if not pts or pts[0] != orbit[w]:
            fail('run-bottom', (p, m, d, w))
        if any(pts[i + 1] - pts[i] != d for i in range(len(pts) - 1)):
            fail('run-step', (p, m, d, w))
        if w < W:
            if pts[-1] != p - d + orbit[w + 1]:
                fail('run-top', (p, m, d, w, pts[-1], p - d + orbit[w + 1]))
            if len(pts) != m0 + (1 if orbit[w] < r else 0):
                fail('run-full-height', (p, m, d, w))
        else:
            if pts[-1] != H:
                fail('prefix-top', (p, m, d, pts[-1], H))
            if len(pts) != m - j0:
                fail('prefix-height', (p, m, d))
            h_W = len(pts)
    F_cW = m0 + (1 if cW < r else 0)
    if not (1 <= h_W <= F_cW):
        fail('h_W-range', (p, m, d, h_W, F_cW))

    # ---------- B: bands and the block identity ----------
    T = sorted(set(x for x in C_full if x < r) | ({cW} if h_W > m0 else set()))
    pred = [C if i < h_W else C_full for i in range(m0)] + [T]
    for i in range(m0 + 1):
        band = sorted(set(x % d for x in S if x // d == i))
        if band != pred[i]:
            fail('band-set', (p, m, d, i, band, pred[i]))
    blocks = []
    for i in range(m0):
        R = pred[i]
        blocks.append([R[j + 1] - R[j] for j in range(len(R) - 1)] + [d - R[-1]])
    topb = [T[j + 1] - T[j] for j in range(len(T) - 1)] + [r - T[-1]]
    pgaps = [x for b in blocks for x in b] + topb
    if pgaps != gaps_of(S, p):
        fail('block-identity', (p, m, d))

    # ---------- C: master identity ----------
    span0 = (r - T[-1]) + C[1]
    if W >= 2:
        R_last = C_full if h_W <= m0 - 1 else C
        gC = gamma2_cyc(C, d)
        gCf = gamma2_cyc(C_full, d) if h_W < m0 else None
        gT = gamma2_cyc(T, r) if len(T) >= 2 else None
        spanTop = (d - R_last[-1]) + (T[1] if len(T) >= 2 else r)
        cands = [gC, span0, spanTop]
        if gCf is not None:
            cands.append(gCf)
        if gT is not None:
            cands.append(gT)
        if max(cands) != G2:
            fail('master-identity', (p, m, d, cands, G2))
        res.update(case='W2+', gC=gC, gCf=gCf, gT=gT,
                   span0=span0, spanTop=spanTop)
        # ---------- D ----------
        if gC > d - W + 1:
            fail('D-gC-bound', (p, m, d, gC, d - W + 1))
        if gCf is not None and gCf > d - W + 2:
            fail('D-gCf-bound', (p, m, d, gCf, d - W + 2))
        if gT is not None and gT > r:
            fail('D-gT-bound', (p, m, d, gT, r))
        if max(C) > T[-1] + sig:
            fail('D-maxC-maxT-sigma', (p, m, d, max(C), T[-1], sig))
        if span0 > gC:
            fail('D-span0-gC', (p, m, d, span0, gC))
        if spanTop > d - W + 2:
            fail('D-spanTop-bound', (p, m, d, spanTop, d - W + 2))
        if not gamma2_bound(C, d) or not gamma2_bound(C_full, d):
            fail('D-gamma2-general', (p, m, d))
        if len(T) >= 2 and not gamma2_bound(T, r):
            fail('D-gamma2-T', (p, m, d))
    else:
        # W == 1: explicit singleton-band spans
        cands = [d, span0]                       # Gamma2(C) = d
        if h_W <= m0 - 1:
            cands.append((d - C[-1]) + d)        # first singleton band
            if m0 - 1 >= h_W + 1:
                cands.append(2 * d)              # consecutive singletons
        R_last = C_full if h_W <= m0 - 1 else C  # C_full = {0}
        spanTop = (d - R_last[-1]) + (T[1] if len(T) >= 2 else r)
        cands.append(spanTop)
        if len(T) >= 2:
            cands.append(gamma2_cyc(T, r))
        if max(cands) != G2:
            fail('master-identity-W1', (p, m, d, cands, G2))
        res.update(case='W1', span0=span0, spanTop=spanTop)
        if spanTop > 2 * d:
            fail('D-spanTop-W1', (p, m, d, spanTop))
    if C_full == [0] and W != 1:
        fail('D-singleton-iff-W1', (p, m, d))
    if d <= 3 * k - 1 and r > 2 * k - 1:
        fail('D-r-bound', (p, m, d, r))

    # W-inequality
    if m == 2 * k + 1:
        if d <= 3 * k - 2 and W < d - 2 * k + 1:
            fail('D-Wineq-m5', (p, m, d, W, d - 2 * k + 1))
        if d == 3 * k - 1 and W != d - 2 * k:
            fail('D-Wineq-m5-top', (p, m, d, W))
        if d == 3 * k and W != d - 2 * k - 1:
            fail('D-Wineq-m5-3k', (p, m, d, W))
    else:
        if 3 <= d <= 3 * k - 2 and W < d - 2 * k + 1:
            fail('D-Wineq-m4', (p, m, d, W, d - 2 * k + 1))
        if d == 3 * k - 1 and W != d - 2 * k:
            fail('D-Wineq-m4-top', (p, m, d, W))
        if d == 3 * k and W != d - 2 * k:
            fail('D-Wineq-m4-3k', (p, m, d, W))

    STATS[res['case']] = STATS.get(res['case'], 0) + 1
    # soft bonus: induced rotation on Z_r (C cap [0,r) as a segment of
    # step sigma mod r) -- counted, not asserted
    Cr = sorted(x for x in C if x < r)
    j = len(Cr)
    if j >= 1:
        step = sig % r
        seg = sorted(set((step * w) % r for w in range(j)))
        if len(seg) == j and seg == Cr:
            STATS['induced-rot-ok'] = STATS.get('induced-rot-ok', 0) + 1
        else:
            STATS['induced-rot-bad'] = STATS.get('induced-rot-bad', 0) + 1
    return res


def g2_direct_np(p, m, d):
    js = (np.arange(m, dtype=np.int64) * d) % p
    S = np.sort(js)
    g = np.diff(S)
    g = np.append(g, p - int(S[-1]) + int(S[0]))
    return int((g + np.roll(g, -1)).max())


def candidates_fast(p, m, d):
    """Candidate-path G2 for W >= 2 (identity verified in stage 1);
    direct for W <= 1."""
    W = (d * (m - 1)) // p
    if W <= 1:
        return g2_direct_np(p, m, d)
    r = p % d
    sig = d - r
    m0 = p // d
    orbit = (-p * np.arange(W + 1, dtype=np.int64)) % d
    cW = int(orbit[-1])
    C = np.sort(orbit)
    C_full = np.sort(orbit[:-1])
    h_W = m - (-((-W * p) // d))          # m - ceil(W p / d)
    T = np.sort(C_full[C_full < r])
    if h_W > m0:
        T = np.sort(np.append(T, cW))
        T = np.unique(T)
    span0 = (r - int(T[-1])) + int(C[1])
    R_last = C_full if h_W <= m0 - 1 else C
    if len(T) >= 2:
        spanTop = (d - int(R_last[-1])) + int(T[1])
    else:
        spanTop = (d - int(R_last[-1])) + r
    cands = [span0, spanTop]
    for R, n in ((C, d), (C_full, d), (T, r)):
        if n == d and len(R) == len(C_full) and h_W >= m0:
            continue          # C_full not an occurring band set
        if len(R) < 2:
            continue
        g = np.diff(R)
        g = np.append(g, n - int(R[-1]))
        cands.append(int((g + np.roll(g, -1)).max()))
    return max(cands)


def main():
    t0 = time.time()
    pr_struct = [p for p in primes_6k1(STRUCT_P_MAX) if p >= 7]
    print('=== Stage 1: structural verification, %d primes p <= %d ==='
          % (len(pr_struct), STRUCT_P_MAX))
    n_struct = 0
    for p in pr_struct:
        k = (p - 1) // 6
        for m in (2 * k + 1, 2 * k + 2):
            for d in range(2, 3 * k + 1):
                check_structure(p, k, m, d)
                n_struct += 1
        sys.stdout.flush()
    print('  %d (p,m,d) triples checked; fails so far: %d'
          % (n_struct, len(FAILS)))

    # ---------------- Stage 2: direct theorem check ----------------
    print('=== Stage 2: direct G2 theorem check, p <= %d ===' % DIRECT_P_MAX)
    n_thm = 0
    worst = None
    maxg2_by_p = {}
    for p in [x for x in primes_6k1(DIRECT_P_MAX) if x >= 31]:
        k = (p - 1) // 6
        best = 0
        best_at = None
        for m in (2 * k + 1, 2 * k + 2):
            for d in range(3, 3 * k):
                G2 = g2_direct_np(p, m, d)
                n_thm += 1
                if G2 > best:
                    best, best_at = G2, (m, d)
                slack = 2 * k + 1 - G2
                if worst is None or slack < worst[0]:
                    worst = (slack, p, m, d, G2)
                if G2 > 2 * k + 1:
                    fail('E-theorem', (p, m, d, G2, 2 * k + 1))
        maxg2_by_p[p] = (best, best_at)
    print('  %d checks; min slack = %r; max G2 per p: first=%s last=%s'
          % (n_thm, worst,
             {p: maxg2_by_p[p] for p in sorted(maxg2_by_p)[:4]},
             {p: maxg2_by_p[p] for p in sorted(maxg2_by_p)[-4:]}))

    # dirty exhibits and endpoint exceptions (all 6k+1 primes <= DIRECT)
    print('=== Stage 2b: dirty primes and endpoint exceptions ===')
    dirty_ex = {}
    endpt = {}
    for p in [x for x in primes_6k1(DIRECT_P_MAX) if 7 <= x]:
        k = (p - 1) // 6
        if p in (13, 19):
            dirty_ex[p] = [d for d in range(3, 3 * k)
                           if g2_direct_np(p, 2 * k + 1, d) > 2 * k + 1]
            print('  p=%d (k=%d): dirty d (m=2k+1): %s'
                  % (p, k, dirty_ex[p]))
        if 2 * k + 1 >= 3:
            v2 = g2_direct_np(p, 2 * k + 1, 2)
            if v2 != 2 * k + 3:
                fail('E-endpoint-d2', (p, v2, 2 * k + 3))
        if k >= 2:
            for m, exp in ((2 * k + 1, 2 * k + 2), (2 * k + 2, 2 * k + 2)):
                v = g2_direct_np(p, m, 3 * k)
                if v != exp:
                    fail('E-endpoint-3k', (p, m, v, exp))
        if 3 * k - 1 >= 3:
            v = g2_direct_np(p, 2 * k + 1, 3 * k - 1)
            if v != 6:
                fail('E-endpoint-3km1', (p, v, 6))
            endpt[p] = v
    print('  endpoint exceptions verified (d=2: 2k+3; d=3k: 2k+2; '
          'd=3k-1: 6)')

    # ---------------- Stage 3: candidate-path check ----------------
    print('=== Stage 3: candidate-path theorem check, '
          '%d < p <= %d ===' % (DIRECT_P_MAX, FAST_P_MAX))
    n_fast = 0
    n_dis = 0
    for p in [x for x in primes_6k1(FAST_P_MAX) if x > DIRECT_P_MAX]:
        k = (p - 1) // 6
        for m in (2 * k + 1, 2 * k + 2):
            for d in range(3, 3 * k):
                G2 = candidates_fast(p, m, d)
                n_fast += 1
                if G2 > 2 * k + 1:
                    fail('E-fast', (p, m, d, G2, 2 * k + 1))
        # spot equality with direct
        for d in range(3, 3 * k, max(1, k // 3)):
            if candidates_fast(p, 2 * k + 1, d) != g2_direct_np(p, 2 * k + 1, d):
                fail('fast-vs-direct', (p, d))
            n_dis += 1
    print('  %d candidate checks + %d direct cross-checks' % (n_fast, n_dis))

    # ---------------- Stage 4: consumer equivalence ----------------
    print('=== Stage 4: consumer equivalence + census cross-check ===')
    # (a) brute-force feet vs G2, sample p <= 150
    n_eq = 0
    for p in [x for x in primes_6k1(150) if x >= 7]:
        k = (p - 1) // 6
        s = 2 * k + 1
        a_grid = np.arange(p)[:, None]
        for d in range(2, 3 * k + 1):
            js = (d * np.arange(s)) % p
            arr = (a_grid + js[None, :]) % p
            for s1 in (2 * k + 1, 2 * k + 2):
                tmax = int((arr >= s1).sum(axis=1).max())
                feet_min = s - tmax
                G2 = g2_direct_np(p, s, d)
                if (feet_min >= 2) != (G2 <= s1):
                    fail('F-equivalence', (p, d, s1, feet_min, G2))
                n_eq += 1
    print('  equivalence min_a feet >= 2 <=> G2 <= s1: %d checks' % n_eq)
    # (b) viable-set cross-check vs the census JSON
    with open(CENSUS_JSON) as f:
        census = json.load(f)
    n_v = 0
    v_ok = 0
    for ps, row in census['per_prime'].items():
        p = int(ps)
        if row.get('cls') != 1 or 'consumed_onefoot' not in row:
            continue
        k = (p - 1) // 6
        for (s1, key) in ((2 * k + 2, '%d_%d' % (2 * k + 2, 2 * k + 1)),
                          (2 * k + 1, '%d_%d' % (2 * k + 1, 2 * k + 1))):
            if key not in row['consumed_onefoot']:
                continue
            v_cen = set(row['consumed_onefoot'][key]['viable'])
            v_g2 = set(d for d in range(2, 3 * k + 1)
                       if g2_direct_np(p, 2 * k + 1, d) > s1)
            n_v += 1
            if v_cen == v_g2:
                v_ok += 1
            else:
                fail('F-viable', (p, key, sorted(v_cen), sorted(v_g2)))
    print('  viable-set cross-check: %d/%d match' % (v_ok, n_v))

    # ---------------- report ----------------
    out = {
        'fails': [list(map(str, f)) for f in FAILS],
        'n_fails': len(FAILS),
        'stats': STATS,
        'n_struct': n_struct,
        'n_theorem_direct': n_thm,
        'n_theorem_fast': n_fast,
        'worst_slack': worst,
        'dirty_exhibits': {str(p): v for p, v in dirty_ex.items()},
        'max_g2_by_p': {str(p): v for p, v in sorted(maxg2_by_p.items())
                        [-12:]},
        'seconds': round(time.time() - t0, 1),
    }
    with open(OUT_JSON, 'w') as f:
        json.dump(out, f, indent=1)
    print('\n=== SUMMARY ===')
    print('total fails: %d' % len(FAILS))
    print('stats: %r' % STATS)
    print('worst slack (2k+1 - G2): %r' % (worst,))
    print('time: %.1fs' % (time.time() - t0))
    print('written: %s' % OUT_JSON)
    return 0 if not FAILS else 1


if __name__ == '__main__':
    sys.exit(main())
