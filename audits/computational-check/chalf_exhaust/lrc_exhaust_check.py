"""EXHAUSTIVENESS CHECK: is every set a seed or an extension of a productive core?

Advisor directive: "verify exhaustiveness separately.  Take small n, enumerate
all sets, and check that every set is either a seed or an extension of a
productive core.  If exhaustiveness fails, the whole program fails, and it is
better to know that now."

Universe: all subsets of [1..B] (B = 18), sizes 2..7 (every set, primitive or
not; scaling covariance noted).  For a k-set S (a potential core):

  P1(S) [tight]      X(S, 1/(k+2)) != {}  -- exact, unbounded filler
                     (General Filler Theorem, lrc_xs_theory).
  P2(S) [ladder]     some x <= XCAP (36) with gap(S+{x}) a LADDER value at
                     size k+1: tight 1/(k+1+1) or rung r/(r(k+1)+1), r >= 2.
                     (No-touch prefilter; cap noted honestly.)

  V (an N-set) is EXTENSION-like if some (N-1)-subset is productive (P1 or
  P2); otherwise SEED-like.  STRICT seed: additionally every removal strictly
  raises the gap (V is a drop from every core -- the reduction-relevant case).

Reported per N: counts, the gap distribution of seeds, MIN GAP OVER SEEDS vs
the ladder (floor 1/(N+1), rung-2 2/(2N+1)), seed lists, the fate of the
known inter-rung inhabitants, and (bonus) the GRID law and pole law at scale.

Stages: --size N runs core-productivity for size N-1 and classification for N.
Idempotent via JSON caches.  Run sizes ascending: 2 3 4 5 6 7.
"""
import argparse
import json
import os
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc_gap_lib import gap_int, _norm
from lrc_xs_theory import X_theory

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
B = 18
XCAP = 2 * B
GAPCACHE = os.path.join(OUT, "exhaust_gapcache.json")
PRODCACHE = os.path.join(OUT, "exhaust_prodcache.json")


def ladder_at(g, N):
    """tight 1/(N+1) or rung r/(rN+1), r>=2, at size N -> 'tight'/'rung'/None."""
    p, q = g.numerator, g.denominator
    if p == 1 and q == N + 1:
        return "tight"
    if p >= 2 and (q - 1) % p == 0 and (q - 1) // p == N:
        return "rung"
    return None


# ---------------------------------------------------------------- gap cache
_gapc = {}


def load_gapcache():
    if os.path.exists(GAPCACHE):
        with open(GAPCACHE) as f:
            raw = json.load(f)
        for k, v in raw.items():
            _gapc[tuple(int(x) for x in k.split(","))] = (
                Fraction(v[0], v[1]), v[2], v[3])
    return _gapc


def save_gapcache():
    raw = {",".join(map(str, k)): [v[0].numerator, v[0].denominator, v[1], v[2]]
           for k, v in _gapc.items()}
    tmp = GAPCACHE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(raw, f)
    os.replace(tmp, GAPCACHE)


def ginfo(S):
    """(gap, unified_grid, primitive) with cache."""
    key = tuple(S)
    if key not in _gapc:
        g, ts = gap_int(list(key))
        dens = {t.denominator for t in ts}
        uni = len(dens) == 1
        prim = gcd(*key) == 1 if len(key) > 1 else key[0] == 1
        _gapc[key] = (g, uni, prim)
    return _gapc[key]


# ---------------------------------------------------------------- prod cache
_prodc = {}


def load_prodcache():
    if os.path.exists(PRODCACHE):
        with open(PRODCACHE) as f:
            raw = json.load(f)
        for k, v in raw.items():
            _prodc[tuple(int(x) for x in k.split(","))] = v
    return _prodc


def save_prodcache():
    tmp = PRODCACHE + ".tmp"
    with open(tmp, "w") as f:
        json.dump({",".join(map(str, k)): v for k, v in _prodc.items()}, f)
    os.replace(tmp, PRODCACHE)


def argmaxes(S):
    g, ts = gap_int(list(S))
    return g, ts


def prod(S):
    """[P1, P2] for the k-set S."""
    key = tuple(sorted(S))
    if key in _prodc:
        return _prodc[key]
    k = len(key)
    g_S, ts = argmaxes(key)
    # P1: exact tight-extension existence (unbounded x)
    try:
        xs, R, tag = X_theory(list(key), Fraction(1, k + 2))
        p1 = len(xs) > 0
    except Exception:
        p1 = None
    # P2: ladder landing with x <= XCAP (no-touch prefilter)
    p2 = False
    own = ladder_at(g_S, k + 1)  # touch-extensions land on g_S itself
    if own is not None:
        p2 = True  # any touching x lands on a ladder value
    else:
        for x in range(1, XCAP + 1):
            if x in key:
                continue
            if any(_norm(x * t) >= g_S for t in ts):
                continue  # touch -> lands on g_S, not ladder
            g2, _ = gap_int(sorted(key + (x,)))
            if ladder_at(g2, k + 1) is not None:
                p2 = True
                break
    _prodc[key] = [p1, p2]
    return _prodc[key]


# ---------------------------------------------------------------- main stage
def run_size(N):
    R = []
    R.append(f"EXHAUSTIVENESS  universe=[1..{B}]  N={N}  (floor 1/{N+1}, "
             f"rung-2 2/{2*N+1})")
    # 1) productivity for all (N-1)-subsets
    k = N - 1
    subs = list(combinations(range(1, B + 1), k))
    todo = [S for S in subs if tuple(sorted(S)) not in _prodc]
    R.append(f"cores of size {k}: {len(subs)} total, {len(todo)} to compute")
    for i, S in enumerate(todo):
        prod(S)
        if (i + 1) % 2000 == 0:
            save_prodcache()
            print(f"  prod {i+1}/{len(todo)}", flush=True)
    save_prodcache()
    p1n = sum(1 for S in subs if _prodc[tuple(sorted(S))][0])
    p2n = sum(1 for S in subs if _prodc[tuple(sorted(S))][1])
    R.append(f"  P1(tight): {p1n}/{len(subs)}   P2(ladder, x<={XCAP}): {p2n}/{len(subs)}")

    # 2) classify all N-sets
    ext1 = ext2 = 0
    seeds1, seeds2 = [], []
    strict_seeds2 = []
    allsets = list(combinations(range(1, B + 1), N))
    gaps_of_seeds2 = []
    interrung_hits = []
    for V in allsets:
        gV, uni, prim = ginfo(V)
        best = None  # min gap over (N-1)-subsets
        all_strict = True
        for x in V:
            S = tuple(sorted(v for v in V if v != x))
            gS = ginfo(S)[0]
            if gS <= gV:
                all_strict = False
            p = _prodc[S]
            if p[0] or p[1]:
                best = "ext"
        has1 = any(_prodc[tuple(sorted(v for v in V if v != x))][0] for x in V)
        has2 = any(_prodc[tuple(sorted(v for v in V if v != x))][1] for x in V)
        if has1:
            ext1 += 1
        else:
            seeds1.append(V)
        if has2:
            ext2 += 1
        else:
            seeds2.append(V)
            gaps_of_seeds2.append((gV, V))
            if all_strict:
                strict_seeds2.append(V)
        # track the known inter-rung inhabitants
        if set(V) in ({1, 5, 6, 11, 16, 17}, {1, 2, 3, 4, 5, 7, 18},
                      {1, 3, 4, 5, 7, 13, 18}):
            interrung_hits.append(V)
    R.append(f"N-sets: {len(allsets)}  extension-like: P1 {ext1}, P2 {ext2} "
             f"| seed-like: P1 {len(seeds1)}, P2 {len(seeds2)} "
             f"(strict P2-seeds: {len(strict_seeds2)})")
    # 3) seed gap analysis
    if seeds2:
        gmin = min(g for g, V in gaps_of_seeds2)
        R.append(f"MIN GAP over P2-seed-like N-sets: {gmin} "
                 f"(floor 1/{N+1} = {Fraction(1, N+1)}, "
                 f"rung-2 2/{2*N+1} = {Fraction(2, 2*N+1)})")
        near = [(g, V) for g, V in gaps_of_seeds2
                if g <= Fraction(2, 2 * N + 1)]
        R.append(f"  seeds with gap <= rung-2: {len(near)}")
        for g, V in sorted(near)[:12]:
            R.append(f"    gap={g} V={V} "
                     f"{'STRICT' if V in [tuple(s) for s in strict_seeds2] else ''}")
        hist = {}
        for g, V in gaps_of_seeds2:
            hist[str(g)] = hist.get(str(g), 0) + 1
        R.append(f"  seed gap histogram: "
                 f"{sorted(hist.items(), key=lambda z: Fraction(z[0]))[:15]}")
        # primitive seeds list (small)
        pseeds = [V for V in seeds2 if ginfo(V)[2]]
        R.append(f"  primitive P2-seeds: {len(pseeds)}: {pseeds[:20]}")
    if seeds1 and not seeds2:
        R.append("(all P1-seeds are P2-extensions)")
    for V in interrung_hits:
        gV = ginfo(V)[0]
        has2 = any(_prodc[tuple(sorted(v for v in V if v != x))][1] for x in V)
        has1 = any(_prodc[tuple(sorted(v for v in V if v != x))][0] for x in V)
        R.append(f"  inter-rung inhabitant {V}: gap={gV} P1-ext={has1} P2-ext={has2}")
    # 4) grid law at scale (this size)
    uni_n = sum(1 for V in allsets if ginfo(V)[1])
    upri = sum(1 for V in allsets if ginfo(V)[1] and ginfo(V)[2])
    npri = sum(1 for V in allsets if ginfo(V)[2])
    R.append(f"[GRID at scale] N={N}: unified {uni_n}/{len(allsets)}; "
             f"primitive-only: {upri}/{npri}")
    report = "\n".join(R)
    print(report)
    with open(os.path.join(OUT, f"exhaust_N{N}.txt"), "w") as f:
        f.write(report + "\n")
    # persist details
    with open(os.path.join(OUT, f"exhaust_N{N}_detail.json"), "w") as f:
        json.dump(dict(
            N=N, B=B, XCAP=XCAP, total=len(allsets),
            ext1=ext1, ext2=ext2,
            seeds1=[list(map(int, V)) for V in seeds1],
            seeds2=[list(map(int, V)) for V in seeds2],
            strict_seeds2=[list(map(int, V)) for V in strict_seeds2],
            seed_gaps={",".join(map(str, V)): str(ginfo(V)[0])
                       for V in seeds2},
        ), f)
    save_gapcache()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", type=int, required=True)
    ap.add_argument("--cache-only", action="store_true")
    args = ap.parse_args()
    load_gapcache()
    load_prodcache()
    if args.cache_only:
        print(f"gapcache={len(_gapc)} prodcache={len(_prodc)}")
        return
    run_size(args.size)
    print(f"caches: gap={len(_gapc)} prod={len(_prodc)}")


if __name__ == "__main__":
    main()
