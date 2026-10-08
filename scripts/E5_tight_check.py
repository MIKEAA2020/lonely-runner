#!/usr/bin/env python3
"""E5_tight_check.py — exact M(V) via breakpoint arithmetic for specific sets:
the n=5 tight candidates and the exceptional max-speed-pair counterexample."""
from itertools import combinations

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
    best_num, best_den, arg = 0, 1, None
    for d in sorted(dens):
        for m in range(d):
            val = min(min((w * m) % d, d - (w * m) % d) for w in V)
            if val * best_den > best_num * d:
                best_num, best_den, arg = val, d, (m, d)
    return best_num, best_den, arg

def tau_grid(N, eff):
    best, argk = 0, 0
    for k in range(1, N):
        m = N
        for w in eff:
            y = (w * k) % N
            d = y if y <= N - y else N - y
            if d < m:
                m = d
        if m > best:
            best, argk = m, k
    return best, argk

def report(V):
    n = len(V)
    M_num, M_den, arg = exact_M(V)
    best = (0, 1, None)
    for (p, q) in combinations(range(n), 2):
        N = V[p] + V[q]
        eff = (V[p],) + tuple(V[i] for i in range(n) if i not in (p, q))
        t, k = tau_grid(N, eff)
        if t * best[1] > best[0] * N:
            best = (t, N, (V[p], V[q], k))
    print(f"V={V}:  M = {M_num}/{M_den} (at t={arg[0]}/{arg[1]})   "
          f"max_pair tau = {best[0]}/{best[1]} via pair {best[2]}   "
          f"equality: {M_num * best[1] == best[0] * M_den}")

if __name__ == "__main__":
    report((2, 6, 8, 10, 11))     # exceptional max-speed-pair counterexample, n=5
    report((1, 3, 4, 5, 9))       # tight candidate from E2
    report((1, 2, 3, 4, 5))       # regular
    report((1, 3, 4, 5, 7))       # near-tight (ratio 1.0909)
    report((1, 4, 5, 6, 7, 11))   # exceptional set at n=6
