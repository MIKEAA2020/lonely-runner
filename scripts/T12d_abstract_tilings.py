#!/usr/bin/env python3
"""T12d — Abstract (position-free) within-class covering solvability.

In class r of the (1K,3U) cell at p^2, the units must cover the fiber with
three arcs: after normalizing by u'_1, these are an interval I (diff 1) and
two APs with differences lambda, mu (the unit ratios).  Sizes:
  p = 6k+1: m=1 -> {2k+1, 2k, 2k}; m=2 -> {2k+1, 2k+1, 2k}; m=3 -> all 2k+1
  p = 6k+5: m=0 -> all 2k+2;        m=1 -> {2k+1, 2k+2, 2k+2}
  p = 3   : all 1.
Overlap budget = sum(sizes) - p.

Abstract question: do ANY (lambda, mu) and translations tile/cover Z_p?
Strict variant: lambda, mu not in {+-1}, lambda != +-mu (mirrors +/-
distinct units).  Relaxed variant: anything goes.

If the strict abstract problem is UNSOLVABLE for a given (p, m), then every
class with that m kills its configuration outright (level-b kill, no
coupling needed).  If solvable, the kill must come from residue selection
or coupling.

Also: quick check of the return-count floor at p = 5 mod 6.
"""
import json
from itertools import permutations

OUT = {}


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = b"\x00" * len(sieve[i * i:: i])
    return [i for i in range(2, n + 1) if sieve[i]]


def ap_set(lam, start, n, p):
    """AP with difference lam: {lam*s : s in [start, start+n)} as frozenset."""
    return frozenset((lam * s) % p for s in range(start, start + n))


def abstract_search(p, sizes, strict):
    """Existence of I(n1) U (lam J(n2) + t2) U (mu K(n3) + t3) = Z_p,
    overlap <= budget.  WLOG I = [0, n1), t2, t3 free."""
    n1, n2, n3 = sizes
    p_ = p
    budget = n1 + n2 + n3 - p
    if budget < 0:
        return None
    full = (1 << p) - 1
    Im = 0
    for q in range(n1):
        Im |= 1 << q
    units = list(range(1, p))
    count = 0
    example = None
    for lam in units:
        if strict and min(lam, p - lam) == 1:
            continue
        # J runs over all intervals of length n2 -> lam*J over all
        # lam-APs; equivalent to translating the base shape lam*[0,n2).
        B2 = frozenset((lam * s) % p for s in range(n2))
        for mu in units:
            if strict and (min(mu, p - mu) == 1 or
                           (lam * pow(mu, -1, p)) % p in (1, p - 1)):
                continue
            B3 = frozenset((mu * s) % p for s in range(n3))
            if len(B3) < n3:      # degenerate repeat (impossible: mu unit)
                continue
            B3list = sorted(B3)
            # precompute translate-masks are built on demand
            for t2 in range(p):
                B2m = 0
                for q in B2:
                    B2m |= 1 << ((q + t2) % p)
                D = Im | B2m
                comp = full & ~D
                if comp == 0:
                    continue       # no room for the third arc -> not a cover
                if comp.bit_count() > n3:
                    continue
                cand = full
                x = comp
                while x and cand:
                    lb = (x & -x).bit_length() - 1
                    m3 = 0
                    for q in B3list:
                        m3 |= 1 << ((lb - q) % p)
                    cand &= m3
                    x &= x - 1
                if cand:
                    # verify real covering with overlap budget
                    t = cand
                    while t:
                        lb = (t & -t).bit_length() - 1
                        K = frozenset((q + lb) % p for q in B3)
                        union = set(range(n1)) | \
                            frozenset((q + t2) % p for q in B2) | K
                        if len(union) == p:
                            count += 1
                            if example is None:
                                example = (lam, mu, t2, lb)
                        t &= t - 1
    return count, example


def size_multisets(p):
    k = p // 6
    if p == 3:
        return [(0, (1, 1, 1))]
    if p % 6 == 1:
        out = []
        base = {1: (2 * k + 1, 2 * k, 2 * k),
                2: (2 * k + 1, 2 * k + 1, 2 * k),
                3: (2 * k + 1, 2 * k + 1, 2 * k + 1)}
        for m, sizes in base.items():
            for perm in set(permutations(sizes)):
                out.append((m, perm))
        return out
    else:
        out = []
        base = {0: (2 * k + 2, 2 * k + 2, 2 * k + 2),
                1: (2 * k + 1, 2 * k + 2, 2 * k + 2)}
        for m, sizes in base.items():
            for perm in set(permutations(sizes)):
                out.append((m, perm))
        return out


PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 47]
rows = []
for p in PRIMES:
    k = p // 6
    for m, sizes in size_multisets(p):
        strict = abstract_search(p, sizes, strict=True)
        relaxed = abstract_search(p, sizes, strict=False)
        rows.append(dict(p=p, k=k, c6=p % 6, m=m, sizes=sizes,
                         strict=strict, relaxed=relaxed))
        if strict is None:
            s_str = "budget < 0 (impossible)"
        else:
            cnt, ex = strict
            s_str = (f"{cnt} solutions, e.g. lam={ex[0]}, mu={ex[1]}, "
                     f"t2={ex[2]}, t3={ex[3]}" if ex else "NONE")
        r_str = ("n/a" if relaxed is None else
                 (str(relaxed[0]) if relaxed else "0"))
        print(f"p={p:>3} ({p % 6} mod 6) m={m} sizes={sizes}: "
              f"strict: {s_str} | relaxed: {r_str}")

# ---- return-count floor at p = 5 mod 6 (quick) ----
print("\n=== floor check at p = 5 mod 6 (min_nu rho) ===")
floors5 = []
for p in primes_upto(300):
    if p % 6 != 5 or p < 29:
        continue
    k = p // 6
    best = None
    for nu in range(1, p):
        c = sum(1 for i in range(1, k + 1)
                if (i * nu) % p <= k or (i * nu) % p >= p - k)
        if best is None or c < best:
            best = c
    floors5.append((p, k, best, k // 5))
    print(f"p={p:>4} k={k:>3} min_rho={best:>3} floor(k/5)={k // 5} "
          f"{'OK' if best == k // 5 else 'MISMATCH'}")

OUT["abstract"] = rows
OUT["floors5"] = floors5
with open("/home/z/my-project/scripts/out_T12d.json", "w") as f:
    json.dump(OUT, f, indent=1, default=str)
print("saved scripts/out_T12d.json")
