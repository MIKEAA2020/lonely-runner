#!/usr/bin/env python3
"""
EQ_independent_verify.py — INDEPENDENT verification of the equality theorem
    M(V) = max over pairs {u,w} of tau(u,w)        (pair-sum lattices)
by a from-scratch implementation sharing NO code with the battery pipeline
(the reviewer's pre-submission recommendation 5: "a from-scratch
implementation, different code path, on a random subset of sets").

Route A (unrestricted, time-domain, folklore candidates):
    M_A(V) = max of f(t) = min_i ||v_i t|| over the exact candidate set
        CAND = { j/(2 v_i) }                (breakpoints of the individual
                                             distance curves: v_i t in Z/2)
             U { c/(v_i+v_j), c=0..S-1 }     (pairwise crossings, SUMS)
             U { c/|v_i-v_j|, c=0..D-1 }     (pairwise crossings, DIFFERENCES)
    Justification (classical, folklore): between consecutive breakpoints
    each ||v_i t|| is affine, so f is concave piecewise-linear there; the
    max on such an interval is attained at an endpoint (a breakpoint) or at
    a kink of the min (a crossing of two affine pieces), and crossings of
    s_i v_i t - m_i = s_j v_j t - m_j have denominator v_i +- v_j.  Hence
    max over CAND = max over t in [0,1).  This route does NOT presuppose
    the pair-sum structure or the sum-only sharpening: difference
    denominators are included.

Route B (pair-sum lattices only — the theorem's right-hand side):
    M_B(V) = max over pairs {p,q}, N = v_p + v_q, k in Z_N of
             min_i ||v_i k / N||.

Checks per set (all exact integer arithmetic; no floating point):
    C1  M_A == M_B            — the equality theorem AND the sum-only
                                 sharpening (sum lattices alone suffice:
                                 the unrestricted candidate set, which
                                 strictly contains the sum grids and the
                                 breakpoints, never beats them).
    C2  M_A >= 1/(n+1)        — LRC ground-truth consistency.
    C3  difference-only max < or = M_B, reported (never >, implied by C1);
                                 the strict count is the strong form of
                                 the sharpening observation.

Corpora: random primitive sets at n=2..6 over mixed ranges, a random
subset of the paper's measurement domain (primitive 5-sets, v<=56),
a beyond-corpus probe (v<=200), and targeted anchors (regular sets,
known tight sets, the n=4 open set, the n=5 max-speed-pair
counterexample).  Deterministic (fixed seed).  Runtime ~1 min.
Run: python3 -u EQ_independent_verify.py
"""
from math import gcd
from itertools import combinations
import random
import time

SEED = 20261003


# ------------------------------------------------------------------ helpers

def dist_mod(v, a, b):
    """|| v*a / b || as an integer residue distance (0 <= d <= b/2)."""
    r = (v * a) % b
    return r if r <= b - r else b - r


def f_min(vs, a, b):
    """min_i ||v_i a / b||, returned as the exact fraction (c, b)."""
    c = b
    for v in vs:
        d = dist_mod(v, a, b)
        if d < c:
            c = d
    return (c, b)


def cmp_frac(p, q):
    """compare the rationals (c1,b1) vs (c2,b2): returns -1, 0, or +1."""
    l = p[0] * q[1]
    r = q[0] * p[1]
    return (l > r) - (l < r)


# ------------------------------------------------------- route A (folklore)

def route_A(vs):
    """Exact M(V) over the unrestricted folklore candidate set."""
    best = (0, 1)
    # breakpoints of the individual distance curves
    for v in vs:
        b = 2 * v
        for j in range(1, b):
            val = f_min(vs, j, b)
            if cmp_frac(val, best) > 0:
                best = val
    # pairwise crossings: sums and differences
    for (u, w) in combinations(vs, 2):
        S = u + w
        for c in range(1, S):
            val = f_min(vs, c, S)
            if cmp_frac(val, best) > 0:
                best = val
        D = abs(w - u)
        for c in range(1, D):
            val = f_min(vs, c, D)
            if cmp_frac(val, best) > 0:
                best = val
    return best


def route_A_diffonly(vs):
    """max of f over difference-denominator candidates only."""
    best = (0, 1)
    for (u, w) in combinations(vs, 2):
        D = abs(w - u)
        for c in range(1, D):
            val = f_min(vs, c, D)
            if cmp_frac(val, best) > 0:
                best = val
    return best


# ------------------------------------------------------ route B (pair-sums)

def route_B(vs):
    """max over pairs of max over the pair-sum lattice; returns
    (best_fraction, witness=(pair, k))."""
    best = (0, 1)
    wit = None
    for (p, q) in combinations(range(len(vs)), 2):
        N = vs[p] + vs[q]
        bc = 0
        bk = 0
        for k in range(1, N):
            m = N
            for v in vs:
                d = dist_mod(v, k, N)
                if d < m:
                    m = d
            if m > bc:
                bc = m
                bk = k
        if cmp_frac((bc, N), best) > 0:
            best = (bc, N)
            wit = ((vs[p], vs[q]), bk)
    return best, wit


# ---------------------------------------------------------------- checking

def check_set(vs):
    n = len(vs)
    ma = route_A(vs)
    mb, wit = route_B(vs)
    md = route_A_diffonly(vs)
    return {
        'n': n,
        'ma': ma,
        'mb': mb,
        'md': md,
        'wit': wit,
        'c1': cmp_frac(ma, mb) == 0,
        'c2': (n + 1) * ma[0] >= ma[1],
        'diff_strict': cmp_frac(md, mb) < 0,
    }


def rand_primitive_set(rng, n, vmax):
    while True:
        vs = tuple(sorted(rng.sample(range(1, vmax + 1), n)))
        g = 0
        for x in vs:
            g = gcd(g, x)
        if g == 1:
            return vs


GROUPS = [
    ('n=2, v<=60', 2, 60, 300),
    ('n=3, v<=50', 3, 50, 400),
    ('n=4, v<=48', 4, 48, 500),
    ('n=5, v<=56  (paper measurement domain)', 5, 56, 1200),
    ('n=6, v<=30', 6, 30, 400),
    ('n=5, v<=200 (beyond-corpus probe)', 5, 200, 60),
]

# (V, expected M as (num, den)) — committed values: regular sets tight at
# 1/(n+1); (1,2,3,5) tight at 1/4 (worklog T-3); non-regular tight sets
# (1,3,4,7) at 1/5 and (1,3,4,5,9) at 1/6 (worklog T-3 diagnostics).
ANCHORS = [
    ((1, 2), (1, 3)),
    ((1, 2, 3), (1, 4)),
    ((1, 2, 3, 4), (1, 5)),
    ((1, 2, 3, 4, 5), (1, 6)),
    ((1, 2, 3, 4, 5, 6), (1, 7)),
    ((1, 2, 3, 5), (1, 4)),
    ((1, 3, 4, 7), (1, 5)),
    ((1, 3, 4, 5, 9), (1, 6)),
]

# report-only anchors (no committed exact value asserted): the n=4 open
# set of the battery, and the n=5 max-speed-pair counterexample.
REPORT_ONLY = [(3, 5, 8, 13), (2, 6, 8, 10, 11)]


def fmt(frac):
    return "%d/%d" % (frac[0], frac[1])


def main():
    t0 = time.time()
    rng = random.Random(SEED)
    print("=" * 72)
    print("INDEPENDENT VERIFICATION OF THE EQUALITY THEOREM")
    print("M(V) = max_pair tau(u,w)   (from-scratch code path, exact")
    print("integer arithmetic, seed %d)" % SEED)
    print("=" * 72)

    fail = []
    total = 0
    eq_count = 0
    lrc_count = 0
    diff_eq = 0
    diff_examples = []

    for (label, n, vmax, count) in GROUPS:
        g_eq = g_lrc = g_diff_eq = 0
        for _ in range(count):
            vs = rand_primitive_set(rng, n, vmax)
            r = check_set(vs)
            total += 1
            if r['c1']:
                eq_count += 1
                g_eq += 1
            else:
                fail.append(('C1', vs, r))
            if r['c2']:
                lrc_count += 1
                g_lrc += 1
            else:
                fail.append(('C2', vs, r))
            if not r['diff_strict']:
                diff_eq += 1
                g_diff_eq += 1
                if len(diff_examples) < 3:
                    diff_examples.append((vs, fmt(r['md'])))
        print("  %-42s %5d sets | C1 equal %5d | C2 LRC %5d | diff==M %4d"
              % (label, count, g_eq, g_lrc, g_diff_eq))

    print("\n[anchors] (exact M asserted against committed values)")
    for (vs, exp) in ANCHORS:
        r = check_set(vs)
        ok = r['c1'] and r['c2'] and cmp_frac(r['ma'], exp) == 0
        total += 1
        if r['c1']:
            eq_count += 1
        if r['c2']:
            lrc_count += 1
        if not ok:
            fail.append(('ANCHOR', vs, r, exp))
        print("  V=%-22s M = %-8s expected %-8s witness pair/k %s  %s"
              % (str(vs), fmt(r['ma']), fmt(exp),
                 (str(r['wit'][0]) + " k=" + str(r['wit'][1]))
                 if r['wit'] else "-",
                 "PASS" if ok else "FAIL"))

    print("\n[report-only anchors]")
    for vs in REPORT_ONLY:
        r = check_set(vs)
        total += 1
        if r['c1']:
            eq_count += 1
        else:
            fail.append(('C1-report', vs, r))
        if r['c2']:
            lrc_count += 1
        else:
            fail.append(('C2-report', vs, r))
        print("  V=%-22s M = %-8s (LRC bound 1/%d)  C1 %s"
              % (str(vs), fmt(r['ma']), len(vs) + 1,
                 "PASS" if r['c1'] else "FAIL"))

    print("\n[summary]")
    print("  sets checked:            %d" % total)
    print("  C1 equality M_A == M_B:  %d/%d  (%.3f%%)"
          % (eq_count, total, 100.0 * eq_count / total))
    print("  C2 LRC bound holds:      %d/%d  (%.3f%%)"
          % (lrc_count, total, 100.0 * lrc_count / total))
    print("  difference-times also attaining M (allowed; never needed): %d"
          % diff_eq)
    if diff_examples:
        print("    examples: %s"
              % "; ".join("V=%s at %s" % (str(v), m)
                          for (v, m) in diff_examples))
    if fail:
        print("\n  FAILURES (%d):" % len(fail))
        for f in fail[:10]:
            print("    %s" % (f,))
        raise SystemExit("INDEPENDENT VERIFICATION FAILED")
    print("\n  INDEPENDENT VERIFICATION PASSED: the unrestricted folklore")
    print("  candidate scan (breakpoints + sums + differences, time domain)")
    print("  equals the pair-sum-lattice value on every set checked; sum")
    print("  lattices alone always suffice.  Elapsed %.1fs"
          % (time.time() - t0))


if __name__ == '__main__':
    main()
