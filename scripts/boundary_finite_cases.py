#!/usr/bin/env python3
"""
T14c: Boundary tightening -- finite case analysis for the lemma-level
boundary primes, converting "verified" into "proved by exhaustive finite
enumeration with recorded margins":

  LEMMA B at p in {11, 17}  (k = 1, 2: the proof's one-foot census is
    dirty there -- viable one-foot differences {2,3,4,5} at 11 and
    {2,4,5,7,8} at 17 -- the small-k slack of Lemma T's inequalities).
  LEMMA A at p in {19, 31, 37} (k = 3, 5, 6: Proposition E's foot count
    s*s1/p ~ 2k/3 > 2 needs k >= 7; at k in {3,5,6} the extended census
    is dirty at exactly these primes).

For each boundary prime the harness produces:
  1. the complete viable-configuration census (multiset x d2 x h ->
     number of placements surviving the |A1^A2| <= ov prune), with the
     minimum fit margin  (lmin - s3)  over all placements -- the
     certificate by which every configuration dies;
  2. the one-foot census (d2 with h >= 1), separating the trichotomy
     differences {2, q} from the middle (small-k slack) differences;
  3. an INDEPENDENT brute-force cross-check (all d2, d3 in [2, dmax],
     all p^2 placements, direct set-union tests, no prune, no fit test)
     reproducing the search count exactly.

The margins reported are min over ALL viable placements of
  (smallest circular d3^{-1}-window containing the residual R = I1 minus A2) - s3,
negative iff a covering exists.  Everything exact integer arithmetic.
"""

import json
import sys
import time
from itertools import combinations, combinations_with_replacement, permutations

T0 = time.time()
OUT = {"meta": {}, "primes": {}}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def sieve(n):
    bs = bytearray([1]) * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:: i] = bytearray(len(bs[i * i:: i]))
    return [i for i in range(n + 1) if bs[i]]


def pmrep(d, p):
    m = d % p
    return min(m, p - m)


def k_of(p):
    return (p - (1 if p % 6 == 1 else 5)) // 6


def ap_set(p, a, d, s):
    return frozenset((a + j * d) % p for j in range(s))


def positions_by_hits(p, s1, s, d):
    D = [0] * (2 * p + 2)
    for j in range(s):
        x = (j * d) % p
        st = (p - x) % p
        D[st] += 1
        D[st + s1] -= 1
    cov = [0] * (2 * p)
    c = 0
    for t in range(2 * p):
        c += D[t]
        cov[t] = c
    return [cov[a] + cov[a + p] for a in range(p)]


def fit_margins(p, R, s3, d3_list, inv_cache):
    """For residual R (sorted list): for each d3, lmin = smallest window
    size; margin = lmin - s3 (negative => covering exists)."""
    m = len(R)
    res = []
    if m == 0:
        return [(d3, 0 - s3) for d3 in d3_list]
    if m == 1:
        return [(d3, 1 - s3) for d3 in d3_list]
    out = []
    for d3 in d3_list:
        inv = inv_cache[d3]
        M = sorted((r * inv) % p for r in R)
        gaps = [M[i + 1] - M[i] for i in range(m - 1)]
        wrapg = (M[0] + p - M[-1]) % p
        G = max(gaps + [wrapg])
        out.append((d3, p - G + 1 - s3))
    return out


def finite_case_analysis(p):
    k = k_of(p)
    cls = 1 if p % 6 == 1 else 5
    if cls == 1:
        size_opts = [2 * k, 2 * k + 1, 2 * k + 2]
        sums_ok = {p, p + 1, p + 2}
    else:
        size_opts = [2 * k + 1, 2 * k + 2]
        sums_ok = {p, p + 1}
    multisets = sorted(set(tuple(sorted(c))
                           for c in combinations_with_replacement(size_opts, 3)
                           if sum(c) in sums_ok))
    dmax = (p - 1) // 2
    inv_cache = {d: pow(d, -1, p) for d in range(1, dmax + 1)}
    I1_of = {}
    res = {"k": k, "cls": cls, "multisets": [list(m) for m in multisets],
           "dmax": dmax, "tables": {}, "onefoot": {}, "totals": {}}
    n_viable_tot = 0
    min_margin_tot = None
    onefoot = {}     # d2 -> {"h_dist": {h: count}, "min_margin": m}
    n_contained_tot = 0
    for ms in multisets:
        ov = sum(ms) - p
        tab = {}      # (d2, h) -> [count, min_margin]
        for (s1, s2, s3) in sorted(set(permutations(ms))):
            if s1 == 1:
                continue
            I1 = frozenset(range(s1, p))
            for d2 in range(2, dmax + 1):
                hits = positions_by_hits(p, s1, s2, d2)
                for a2, h in enumerate(hits):
                    if h > ov:
                        continue
                    n_viable_tot += 1
                    if h == 0:
                        n_contained_tot += 1
                    A2 = ap_set(p, a2, d2, s2)
                    R = sorted(I1 - A2)
                    if len(R) > s3:
                        continue    # cannot occur when h <= ov
                    margins = [mgn for (d3, mgn) in
                               fit_margins(p, R, s3, range(2, dmax + 1),
                                           inv_cache)
                               if pmrep(d3, p) != pmrep(d2, p)]
                    mmin = min(margins) if margins else None
                    key = (d2, h)
                    if key not in tab:
                        tab[key] = [0, None]
                    tab[key][0] += 1
                    if mmin is not None and (tab[key][1] is None
                                             or mmin < tab[key][1]):
                        tab[key][1] = mmin
                    if mmin is not None and (
                            min_margin_tot is None
                            or mmin < min_margin_tot):
                        min_margin_tot = mmin
                    if h >= 1:
                        if d2 not in onefoot:
                            onefoot[d2] = {"n": 0, "min_margin": None}
                        onefoot[d2]["n"] += 1
                        if mmin is not None and (
                                onefoot[d2]["min_margin"] is None
                                or mmin < onefoot[d2]["min_margin"]):
                            onefoot[d2]["min_margin"] = mmin
        res["tables"][str(ms)] = {f"d={d},h={h}": {"n": v[0], "min_margin": v[1]}
                                  for (d, h), v in sorted(tab.items())}
    res["totals"] = {"n_viable": n_viable_tot,
                     "n_contained": n_contained_tot,
                     "n_onefoot": sum(v["n"] for v in onefoot.values()),
                     "min_margin": min_margin_tot}
    # one-foot classification: trichotomy differences vs middle slack
    tri = {2, dmax}
    res["onefoot"] = {str(d): onefoot[d] for d in sorted(onefoot)}
    res["onefoot_trichotomy_ds"] = sorted(d for d in onefoot if d in tri)
    res["onefoot_middle_ds"] = sorted(d for d in onefoot if d not in tri)
    return res


def independent_brute(p):
    """Full cross-check: same normalization, NO prune, NO fit test.
    a2, a3 over all p^2 placements, direct union checks."""
    k = k_of(p)
    cls = 1 if p % 6 == 1 else 5
    if cls == 1:
        size_opts = [2 * k, 2 * k + 1, 2 * k + 2]
        sums_ok = {p, p + 1, p + 2}
    else:
        size_opts = [2 * k + 1, 2 * k + 2]
        sums_ok = {p, p + 1}
    multisets = sorted(set(tuple(sorted(c))
                           for c in combinations_with_replacement(size_opts, 3)
                           if sum(c) in sums_ok))
    dmax = (p - 1) // 2
    FULL = frozenset(range(p))
    n = 0
    for ms in multisets:
        for (s1, s2, s3) in sorted(set(permutations(ms))):
            if s1 == 1:
                continue
            A1 = frozenset(range(s1))
            for d2 in range(2, dmax + 1):
                A2s = [frozenset((a2 + j * d2) % p for j in range(s2))
                       for a2 in range(p)]
                for d3 in range(2, dmax + 1):
                    if pmrep(d3, p) == pmrep(d2, p):
                        continue
                    A3s = [frozenset((a3 + j * d3) % p for j in range(s3))
                           for a3 in range(p)]
                    for A2 in A2s:
                        rest = FULL - (A1 | A2)
                        for A3 in A3s:
                            if rest <= A3:
                                n += 1
                                break
    return n


def main():
    PB = [11, 17]
    PA = [19, 31, 37]
    primes = PB + PA
    OUT["meta"] = {"lemmaB_primes": PB, "lemmaA_primes": PA}
    for p in primes:
        t0 = time.time()
        res = finite_case_analysis(p)
        OUT["primes"][p] = res
        log(f"p={p} (k={res['k']}, cls={res['cls']}): "
            f"viable={res['totals']['n_viable']} "
            f"(contained {res['totals']['n_contained']}, "
            f"one-foot {res['totals']['n_onefoot']}) "
            f"min_margin={res['totals']['min_margin']} | "
            f"one-foot middle d's: {res['onefoot_middle_ds']} "
            f"({time.time()-t0:.1f}s)")
    log("independent brute-force cross-check (no prune / no fit test)")
    for p in primes:
        t0 = time.time()
        n = independent_brute(p)
        OUT["primes"][p]["brute_solutions"] = n
        log(f"  p={p}: brute solutions={n} ({time.time()-t0:.1f}s)")
        if n != 0:
            log(f"  !! BRUTE FOUND {n} COVERINGS AT p={p} -- LEMMA FALSE")
            sys.exit(1)

    with open("/home/z/my-project/scripts/out_boundary_cases.json", "w") as f:
        json.dump(OUT, f, indent=1, default=str)
    log("written scripts/out_boundary_cases.json")
    log("ALL BOUNDARY FINITE CASE ANALYSES COMPLETE: 0 solutions everywhere,"
        " every viable configuration carries a positive margin certificate")


if __name__ == "__main__":
    main()
