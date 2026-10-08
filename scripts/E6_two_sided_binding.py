#!/usr/bin/env python3
"""
E6_two_sided_binding.py — n-genericity check of TWO-SIDED BINDING at the exact
global argmax (the ingredient that makes M(V) = max_pair tau an equality).

For a global argmax t* = m/d with value val/d (val = min residual):
  binders   = runners with residual val;
  two-sided = some binder has residue +val AND some binder has residue d-val
              (for val = d/2 the pole case: any two pole binders qualify).
The binder pair-sum identity (v_u+v_w)m = 0 mod d is then automatic.

Control: n = 4 (proved last turn — expect 100%).
Test:    n = 5, 6 random gcd-1 sets.
"""
from math import gcd
from itertools import combinations
import random

def exact_M_args(V):
    dens = set()
    n = len(V)
    for i in range(n):
        dens.add(2 * V[i])
        for j in range(i + 1, n):
            dens.add(V[i] + V[j])
            if V[i] != V[j]:
                dens.add(abs(V[i] - V[j]))
    dens.discard(0)
    best_num, best_den, args = 0, 1, []
    for d in sorted(dens):
        for m in range(d):
            val = min(min((w * m) % d, d - (w * m) % d) for w in V)
            if val * best_den > best_num * d:
                best_num, best_den, args = val, d, [(d, m, val)]
            elif val * best_den == best_num * d and val > 0:
                args.append((d, m, val))
    return best_num, best_den, args

def two_sided_at(V, d, m, val):
    rs = [(w * m) % d for w in V]
    cs = [min(r, d - r) for r in rs]
    binders = [i for i, c in enumerate(cs) if c == val]
    if val * 2 == d:   # pole case: any two binders at d/2
        return len(binders) >= 2
    pos = any(rs[i] == val for i in binders)
    neg = any(rs[i] == d - val for i in binders)
    return pos and neg

def run(n, B, nsamp, seed):
    rng = random.Random(seed)
    tested = 0
    ts_any = 0        # some argmax breakpoint is two-sided
    ts_all = 0        # every argmax breakpoint is two-sided
    singles = 0       # argmaxes with a single binder (would break the law)
    while tested < nsamp:
        V = tuple(sorted(rng.sample(range(1, B + 1), n)))
        g = 0
        for x in V:
            g = gcd(g, x)
        if g != 1:
            continue
        tested += 1
        M_num, M_den, args = exact_M_args(V)
        flags = [two_sided_at(V, d, m, val) for (d, m, val) in args]
        if any(flags):
            ts_any += 1
        if all(flags) and flags:
            ts_all += 1
        for (d, m, val) in args:
            rs = [(w * m) % d for w in V]
            cs = [min(r, d - r) for r in rs]
            if sum(1 for c in cs if c == val) == 1:
                singles += 1
    print(f"n={n}, B<={B}, {nsamp} sets: argmax two-sided (some): {ts_any}/{tested}, "
          f"(all): {ts_all}/{tested}; single-binder argmaxes: {singles}")

if __name__ == "__main__":
    run(4, 16, 200, 7)
    run(5, 24, 200, 7)
    run(6, 18, 80, 7)
