#!/usr/bin/env python3
"""N5_prime_scan_crosscheck.py — independent (set-based, no bitmask)
verification of the unit 4-tuple covering-capability scan at selected
primes: positive controls 17/19/37 (must come out capable) and extension
primes 53..71 (claimed safe by the mask-based scan)."""
from math import gcd
import time

def dist(x, N):
    r = x % N
    return r if r <= N - r else N - r

def bad_set(w, N, T=6):
    return frozenset(k for k in range(N) if T * dist(w * k, N) < N)

def capable_setbased(N, T=6):
    units = [w for w in range(1, N) if gcd(w, N) == 1]
    B = {w: bad_set(w, N, T) for w in units}
    B1 = B[1]
    full = frozenset(range(N))
    covers = []
    n = len(units)
    for i in range(n):
        b = units[i]
        sb = B1 | B[b]
        for jj in range(i, n):
            c = units[jj]
            sc = sb | B[c]
            for kk in range(jj, n):
                d = units[kk]
                if sc | B[d] == full:
                    covers.append((1, b, c, d))
    return covers

if __name__ == "__main__":
    for N in (17, 19, 37, 53, 59, 61, 67, 71, 73, 79):
        t0 = time.time()
        cov = capable_setbased(N)
        print("N=%3d: capable=%-5s  coverings=%5d  (%.1fs)"
              % (N, bool(cov), len(cov), time.time() - t0))
