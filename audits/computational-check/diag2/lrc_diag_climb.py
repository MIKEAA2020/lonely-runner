"""Test 2 (adversarial): gap-minimizing hill-climb at n=9,10, M<=60.

Random-restart hill climbing with plateau walks over n-subsets of [1..M],
minimizing the exact gap. Any set with gap < 1/(n+1) is a deficit-law
violation (bug or counterexample) and gets re-verified with the exact
reference implementation. Also censuses all critical sets (gap exactly
1/(n+1), gcd-normalized) encountered, flagging non-regular ones.

Usage: python3 lrc_diag_climb.py n M restarts steps out.json
"""
import sys
import json
import random
from collections import Counter
from fractions import Fraction
from math import gcd
import time

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_diag_gap import GapEngine, gap_int


def normalize(V):
    g = gcd(*V)
    return tuple(x // g for x in V)


def climb(n, M, restarts, steps, seed=0):
    rng = random.Random(seed)
    eng = GapEngine(M)
    bound = Fraction(1, n + 1)
    best_gap = None
    best_sets = set()
    critical = Counter()      # normalized sets with gap == bound
    critical_nonreg = set()   # those not containing {1..n-1}
    violations = []
    regular = set(range(1, n))
    evals = 0
    t0 = time.time()

    def gap_of(V):
        nonlocal evals
        evals += 1
        return eng.gap(V)

    def record(V, g):
        nonlocal best_gap
        if best_gap is None or g < best_gap:
            best_gap = g
        if g == bound:
            nv = normalize(V)
            critical[nv] += 1
            if not regular.issubset(set(nv)):
                critical_nonreg.add(nv)
        if g < bound:
            gref, _ = gap_int(V)
            violations.append((sorted(V), str(g), str(gref)))

    # seeded starts: regular and near-regular
    seeds = []
    seeds.append(list(range(1, n + 1)))              # regular critical
    seeds.append(list(range(1, n)) + [n + 2])        # filler variant
    seeds.append(list(range(1, n)) + [2 * n])        # Theorem 6 family
    seeds.append([1, 3, 4, 7] + [x for x in range(5, 5 + n - 4)] if n > 4 else [1, 3, 4, 7])
    for start in seeds:
        V = sorted(set(start))
        if len(V) == n and max(V) <= M:
            record(V, gap_of(V))

    for r in range(restarts):
        if r < len(seeds):
            V = sorted(set(seeds[r]))
            if len(V) != n or max(V) > M:
                V = sorted(rng.sample(range(1, M + 1), n))
        else:
            V = sorted(rng.sample(range(1, M + 1), n))
        cur = gap_of(V)
        record(V, cur)
        for s in range(steps):
            # propose a batch of single-speed replacements, take the best
            cands = []
            for _ in range(8):
                i = rng.randrange(n)
                x = rng.randint(1, M)
                if x in V:
                    continue
                W = list(V)
                W[i] = x
                W = sorted(W)
                cands.append((gap_of(W), W))
            if not cands:
                continue
            cands.sort(key=lambda cw: cw[0])
            gW, W = cands[0]
            if gW < cur or (gW == cur and rng.random() < 0.35):
                V, cur = W, gW
                record(V, cur)
            if cur < bound:
                break
        # keep elite sets
        if cur == best_gap or cur == bound:
            best_sets.add(normalize(V))
    stats = dict(
        n=n, M=M, restarts=restarts, steps=steps, evals=evals,
        bound=str(bound), best_gap=str(best_gap),
        violations=[(v, g, gr) for v, g, gr in violations],
        violation_count=len(violations),
        critical_count=len(critical),
        critical_examples=[list(k) for k in list(critical)[:60]],
        critical_nonregular=[list(k) for k in list(critical_nonreg)[:60]],
        critical_nonregular_count=len(critical_nonreg),
        elite_examples=[list(k) for k in list(best_sets)[:30]],
        seconds=round(time.time() - t0, 1),
    )
    return stats


if __name__ == "__main__":
    n, M = int(sys.argv[1]), int(sys.argv[2])
    restarts, steps = int(sys.argv[3]), int(sys.argv[4])
    outp = sys.argv[5]
    stats = climb(n, M, restarts, steps)
    with open(outp, "w") as f:
        json.dump(stats, f, indent=1)
    print(f"n={n} M={M}: evals={stats['evals']} in {stats['seconds']}s; "
          f"best_gap={stats['best_gap']} (bound {stats['bound']}); "
          f"violations={stats['violation_count']}; "
          f"critical(normalized)={stats['critical_count']} "
          f"(non-regular: {stats['critical_nonregular_count']})")
