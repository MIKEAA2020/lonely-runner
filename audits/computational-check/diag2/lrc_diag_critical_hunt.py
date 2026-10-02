"""Targeted hunt for NON-REGULAR critical sets at n=8,9,10 (gap == 1/(n+1)).

Mechanism (from the zoo/drop program): critical sets beyond the regular
family are floor-landing extensions of rung-2 / tight cores.  We grow them
recursively: start from critical and near-critical (rung-2) sets at small
n (from the exhaustive scans and the Theorem 6 family), add fillers
x <= 60, and keep any extension landing exactly on the floor 1/(n+1).

Usage: python3 lrc_diag_critical_hunt.py out.json
"""
import sys
import json
from fractions import Fraction
from math import gcd
import time

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_diag_gap import GapEngine, gap_int

M = 60


def normalize(V):
    g = gcd(*V)
    return tuple(sorted(x // g for x in V))


def rung2_family(k):
    """Theorem 6 family {1,...,k-1, 2k}: gap = 2/(2k+1) (rung 2 for size k)."""
    return list(range(1, k)) + [2 * k]


def hunt(eng):
    results = {}
    # seeds at each size: normalized critical + rung-2 sets, grown stepwise
    # size 3..7 seeds from the exhaustive scans (critical + d=1 examples)
    seeds = {n: set() for n in range(3, 11)}
    for n in range(3, 8):
        try:
            s = json.load(open(f"out/diag2_scan_n{n}_M18.json"))
        except FileNotFoundError:
            continue
        for V in s["critical_normalized_examples"]:
            seeds[n].add(tuple(V))
        for V, g, d in s["near_critical_examples"]:
            if Fraction(g) == Fraction(2, 2 * n + 1):  # rung-2 sets
                seeds[n].add(tuple(V))
    for k in range(3, 10):
        seeds[k].add(normalize(rung2_family(k)))
    seeds[6].add((1, 2, 3, 4, 5, 7))       # rung-2 core (census vector)
    seeds[7].add((1, 2, 3, 4, 5, 7, 12))   # non-regular critical (scan)

    for n in range(8, 11):
        bound = Fraction(1, n + 1)
        found = set()
        # extensions of size-(n-1) seeds by one filler
        pool = sorted(seeds[n - 1])
        checked = 0
        for S in pool:
            Sset = set(S)
            for x in range(2, M + 1):
                if x in Sset:
                    continue
                V = sorted(Sset | {x})
                checked += 1
                if eng.gap(V) == bound:
                    found.add(normalize(V))
        # also two-fillers from size-(n-2) seeds (bounded: only rung-2 seeds)
        pool2 = sorted(s for s in seeds[n - 2] if len(s) == n - 2)
        for S in pool2[:60]:
            Sset = set(S)
            for x in range(2, M + 1):
                if x in Sset:
                    continue
                V1 = sorted(Sset | {x})
                g1 = eng.gap(V1)
                if g1 < bound:
                    continue  # already below floor would be a violation; skip
                for y in range(x + 1, M + 1):
                    if y in Sset:
                        continue
                    V = sorted(Sset | {x, y})
                    checked += 1
                    if eng.gap(V) == bound:
                        found.add(normalize(V))
        # grow seeds for the next size: critical + rung-2 extensions
        for S in pool:
            Sset = set(S)
            for x in range(2, M + 1):
                if x in Sset:
                    continue
                V = sorted(Sset | {x})
                g = eng.gap(V)
                if g == bound:
                    seeds[n].add(normalize(V))
                elif g == Fraction(2, 2 * n + 1):
                    seeds[n].add(normalize(V))
        results[n] = dict(
            seeds_from=len(pool) + len(pool2[:60]),
            extensions_checked=checked,
            nonregular_critical=[list(t) for t in sorted(found)],
            nonregular_critical_count=len(found),
        )
    return results


if __name__ == "__main__":
    eng = GapEngine(M)
    t0 = time.time()
    res = hunt(eng)
    out = {str(n): r for n, r in res.items()}
    with open(sys.argv[1], "w") as f:
        json.dump(out, f, indent=1)
    for n, r in res.items():
        print(f"n={n}: non-regular critical sets (M<={M}): {r['nonregular_critical_count']}")
        for t in r["nonregular_critical"][:12]:
            print("   ", t)
    print(f"total {round(time.time()-t0,1)}s")
