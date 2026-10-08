#!/usr/bin/env python3
"""
E4_cross_check.py — cross-validation of the general-n Reduction Theorem.

For random n-speed sets (n = 5, 6), gcd 1:
  - M(V) computed EXACTLY by unrestricted breakpoint enumeration
    (candidates t = m/d, d in {2 v_i} U {v_i +- v_j}; the PL min-function
    attains its max at such breakpoints);
  - max over pairs of tau_{n-1}(p,q)/(v_p+v_q)  (grid computation);
  - assert Theorem 1:  M(V) >= max_pair tau/N   (exact fraction compare);
  - sanity: M(V) >= 1/(n+1) (LRC, known for these n).
Also reports the tightest Theorem-1 relations and whether M equals the
best-pair tau (grid attains the unrestricted optimum at the set level).
"""
from math import gcd
from itertools import combinations
import random

def exact_M(V):
    dens = set()
    n = len(V)
    for i in range(n):
        dens.add(2 * V[i])
        for j in range(i + 1, n):
            dens.add(V[i] + V[j])
            if V[i] != V[j]:
                dens.add(abs(V[i] - V[j]))
    dens.discard(0)
    best_num, best_den = 0, 1
    for d in dens:
        for m in range(d):
            val = min(min((w * m) % d, d - (w * m) % d) for w in V)
            if val * best_den > best_num * d:
                best_num, best_den = val, d
    return best_num, best_den

def tau_grid(N, eff):
    best = 0
    half = N // 2
    for k in range(1, half + 1):
        m = N
        for w in eff:
            y = (w * k) % N
            d = y if y <= N - y else N - y
            if d < m:
                m = d
        if m > best:
            best = m
    return best

def max_pair_tau(V):
    n = len(V)
    T = n + 1
    best = (0, 1)   # fraction tau/N
    for (p, q) in combinations(range(n), 2):
        N = V[p] + V[q]
        eff = (V[p],) + tuple(V[i] for i in range(n) if i not in (p, q))
        t = tau_grid(N, eff)
        if t * best[1] > best[0] * N:
            best = (t, N)
    return best   # best[0]/best[1]

def run(n, B, nsamp, seed):
    rng = random.Random(seed)
    tested = 0
    viol_T1 = 0
    viol_LRC = 0
    min_slack = None
    eq_count = 0
    while tested < nsamp:
        V = tuple(sorted(rng.sample(range(1, B + 1), n)))
        g = 0
        for x in V:
            g = gcd(g, x)
        if g != 1:
            continue
        tested += 1
        M_num, M_den = exact_M(V)
        t_num, t_den = max_pair_tau(V)
        # Theorem 1: M >= tau  <=>  M_num * t_den >= t_num * M_den
        if M_num * t_den < t_num * M_den:
            viol_T1 += 1
            print(f"  T1 VIOLATION: V={V} M={M_num}/{M_den} tau={t_num}/{t_den}")
        # slack = M - tau
        s = M_num / M_den - t_num / t_den
        if min_slack is None or s < min_slack[0]:
            min_slack = (s, V)
        if s == 0:
            eq_count += 1
        # LRC sanity: M >= 1/(n+1)  <=>  M_num * (n+1) >= M_den
        if M_num * (n + 1) < M_den:
            viol_LRC += 1
            print(f"  LRC VIOLATION: V={V} M={M_num}/{M_den}")
    print(f"n={n}, B<={B}: {nsamp} random gcd-1 sets | T1 violations: {viol_T1} | "
          f"LRC violations: {viol_LRC} | M == best-pair tau: {eq_count}/{nsamp}")
    if min_slack:
        print(f"   min slack (M - max_pair_tau): {min_slack[0]:.6f} at V={min_slack[1]}")

if __name__ == "__main__":
    run(5, 24, 150, 20261003)
    run(6, 18, 60, 20261003)
