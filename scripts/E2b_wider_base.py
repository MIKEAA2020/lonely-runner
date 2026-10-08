#!/usr/bin/env python3
"""
E2b_wider_base.py — (i) wider empirical base for (TAU-n) at n=5 (B=24) and
n=6 (B=18); (ii) identify and dump the sets certified ONLY by pairs not
containing the max speed (the counterexamples to the 'max-speed pair' rule);
(iii) full pair tables for those sets, for hand-verification.

Same conventions as E2_generalize_n5.py. tau_grid optimized: half-range
(k <-> N-k symmetry) + early pruning.
"""
from math import gcd
from itertools import combinations
from collections import Counter

def tau_grid(N, eff):
    """max over k in Z_N of min_w ||w*k||_N (integer units); k=0 excluded (m=0)."""
    best = 0
    half = N // 2
    for k in range(1, half + 1):
        m = N
        for w in eff:
            y = (w * k) % N
            d = y if y <= N - y else N - y
            if d < m:
                m = d
                if m <= best:
                    break
        if m > best:
            best = m
    return best

def jgon_label(N, eff, nmax):
    for j in range(2, nmax + 1):
        if N % j == 0 and all(w % j != 0 for w in eff):
            return j
    return None

def pair_table(V, T):
    n = len(V)
    rows = []
    for (p, q) in combinations(range(n), 2):
        N = V[p] + V[q]
        others = [V[i] for i in range(n) if i not in (p, q)]
        eff = (V[p],) + tuple(others)
        best, argk = tau_grid_full(N, eff)
        cert = T * best >= N
        mirror = any((others[a] + others[b]) % N == 0
                     for a in range(len(others)) for b in range(a + 1, len(others)))
        j = jgon_label(N, eff, T)
        rows.append((p, q, N, best, T * best - N, cert, mirror, j, argk, eff))
    return rows

def tau_grid_full(N, eff):
    """same as tau_grid but also returns argmax k."""
    best, argk = 0, 0
    half = N // 2
    for k in range(1, half + 1):
        m = N
        for w in eff:
            y = (w * k) % N
            d = y if y <= N - y else N - y
            if d < m:
                m = d
        if m > best:
            best, argk = m, k
    return best, argk

def battery_light(n, B, dump_only_nonmax=True):
    T = n + 1
    total = ok_any = ok_max = ok_top = ok_minmax = 0
    only_nonmax = []
    failures = []
    mech = Counter()
    for V in combinations(range(1, B + 1), n):
        g = 0
        for x in V:
            g = gcd(g, x)
        if g != 1:
            continue
        total += 1
        certs = []
        best_ratio, best_eff, best_N = 0, None, None
        for (p, q) in combinations(range(n), 2):
            N = V[p] + V[q]
            others = [V[i] for i in range(n) if i not in (p, q)]
            eff = (V[p],) + tuple(others)
            best = tau_grid(N, eff)
            if T * best >= N:
                certs.append((p, q))
            r = (T * best) / N
            if r > best_ratio:
                best_ratio, best_eff, best_N = r, eff, N
        if certs:
            ok_any += 1
            if any(n - 1 in pq for pq in certs):
                ok_max += 1
            else:
                only_nonmax.append(V)
            if (n - 2, n - 1) in certs:
                ok_top += 1
            if (0, n - 1) in certs:
                ok_minmax += 1
            j = jgon_label(best_N, best_eff, T)
            mech['jgon-j%d' % j if j else 'plain'] += 1
        else:
            failures.append(V)
    print(f"\n===== n = {n}, speeds <= {B}, target 1/{T} =====")
    print(f"sets: {total}   [any pair]: {ok_any}/{total}   "
          f"[max-speed pair]: {ok_max}/{total}   [top pair]: {ok_top}/{total}   "
          f"[(min,max)]: {ok_minmax}/{total}")
    print(f"mechanism census: {dict(mech)}")
    print(f"FAILURES: {len(failures)}" + (f"  e.g. {failures[:5]}" if failures else ""))
    print(f"sets certified ONLY by non-max pairs: {len(only_nonmax)}")
    for V in only_nonmax[:6]:
        print(f"  V = {V}")
        for (p, q, N, best, slack, cert, mirror, j, argk, eff) in pair_table(V, T):
            if cert:
                print(f"    CERT pair ({V[p]},{V[q]}) N={N} tau={best} "
                      f"slack={slack} mirror={mirror} jgon={j} argk={argk} eff={eff}")
    return only_nonmax

if __name__ == "__main__":
    battery_light(5, 16)
    battery_light(5, 24)
    battery_light(6, 18)
