#!/usr/bin/env python3
"""T13b: dump the exact inventory of contained (h=0) and one-foot (h=1)
critical-size APs at representative large primes, for the proof write-up."""

import sys
sys.path.insert(0, "/home/z/my-project/scripts")
from ap_lemmas_verify import positions_by_hits, ap_set, k_of


def runs(A):
    """compress a sorted set of integers into maximal consecutive runs"""
    out = []
    for x in A:
        if out and out[-1][1] == x - 1:
            out[-1][1] = x
        else:
            out.append([x, x])
    return [f"[{a},{b}]" if a < b else f"[{a}]" for a, b in out]


def inventory(p, s1, s, d, hmax=1):
    hits = positions_by_hits(p, s1, s, d)
    res = []
    for a, h in enumerate(hits):
        if h <= hmax:
            A = sorted(ap_set(p, a, d, s))
            feet = [x for x in A if x < s1]
            res.append((h, a, feet, runs([x for x in A if x >= s1])))
    return res


for p, cls in ((239, 5), (241, 1)):
    k = k_of(p)
    dm = (p - 1) // 2
    if cls == 5:
        sizes = [2 * k + 1, 2 * k + 2]
    else:
        sizes = [2 * k, 2 * k + 1, 2 * k + 2]
    print(f"===== p={p} (6k+{cls}, k={k}): inventory, I1=[s1,{p-1}] =====")
    for s1 in sizes:
        for s in sizes:
            if s > p - s1:
                continue
            for d in (2, dm):
                inv = inventory(p, s1, s, d)
                if not inv:
                    continue
                print(f"--- s1={s1} s={s} d={d}: {len(inv)} APs")
                for h, a, feet, rr in inv:
                    print(f"    h={h} a={a} feet={feet} I1runs={rr}")
    print()
