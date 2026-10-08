#!/usr/bin/env python3
"""
T13: Exhaustive verification of the AP lemmas (Lemma A / Lemma B), the
Trichotomy, and the reduction inputs.

Statements under test (collaboration record, audits/glm thought2.txt):

LEMMA A (p = 6k+1 prime): No three APs A1,A2,A3 in Z_p with pairwise
  +-distinct differences, sizes si in {2k, 2k+1, 2k+2}, sum(si) - p in
  {0,1,2}, and A1 U A2 A3 = Z_p.  [Human: true for k >= 3 (p >= 19);
  p = 7, 13 exhibit only non-+-distinct solutions.]

LEMMA B (p = 6k+5 prime): No three APs, pairwise +-distinct differences,
  sizes si in {2k+1, 2k+2}, sum(si) - p in {0,1} (multisets
  {2k+1,2k+2,2k+2} and {2k+2,2k+2,2k+2}), covering Z_p.

TRICHOTOMY: for k >= 2 (p = 6k+5) resp. k >= 3 (p = 6k+1): any AP of the
  critical sizes contained in the complement arc I1 = [s1, p-1] (s1 in the
  critical range) has difference d = 2 or d = (p-1)/2; the latter is
  exactly the two-block structure (blocks of ceil(s/2), floor(s/2)
  consecutive integers offset by (p-1)/2).

REDUCTION INPUTS: per-fiber unit bad-set sizes L(c) in {2k+1, 2k+2}
  (p = 6k+5) resp. {2k, 2k+1} (p = 6k+1), on/off-ball split at c in B_k;
  each fiber bad set is an AP with difference u^{-1} mod p.

CELL ANCHOR: (1K,3U) cells at p^2 empty - independent brute force over
  unit +-classes mod p^2 (small p).  RIGID ZONE: min number of
  +-distinct dilates of B_k covering Z_p (n=5 side), cross-check of the
  committed {7,13,17,19,37}.

Search algorithm (exact): normalize one AP (scale by d1^{-1}, translate)
to d1 = 1, A1 = [0, s1-1].  Exact prune: |A1 ^ A2| <= ov (inclusion-
exclusion: ov = p12+p13+p23 - triple >= p12).  Then R = I1 \\ A2 must fit
inside a d3-AP of size s3 (circular-arc fit).  Every recorded solution
is re-verified by direct set construction.

Anchors (cross-validation discipline): p=17 rigid family (1,2,2);
p=23 wrapping AP {5,13,21,6,14,22,7,15} in arc [5,22]; p=13 trichotomy
exception d=5; p=17 census {2,8}; parity families at 11/17; capacity
table; committed cell emptiness and rigid zone.
"""

import json
import sys
import time
from itertools import combinations, combinations_with_replacement, permutations
from math import gcd

import numpy as np

T0 = time.time()
OUT = {"meta": {}, "lemmaA": {}, "lemmaB": {}, "relaxed": {},
       "trichotomy": {}, "capacity": {}, "cells": {}, "rigid_zone": {},
       "anchors": {}, "near_miss": {}}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- utilities

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


def ap_set(p, a, d, s):
    return frozenset((a + j * d) % p for j in range(s))


def k_of(p):
    """k for p = 6k+1 or 6k+5."""
    return (p - (1 if p % 6 == 1 else 5)) // 6


def ball(p):
    """B_k = {r in [0,p): min(r, p-r) <= k}, the mod-p ball."""
    k = k_of(p)
    return frozenset(r for r in range(p) if min(r, p - r) <= k)


def positions_by_hits(p, s1, s, d):
    """hits[a] = #{j in [0,s): (a + j*d) mod p in [0, s1)} for a in [0,p).

    Difference-array over the circular forbidden arcs
    [(-x) mod p, (-x + s1 - 1) mod p], x = j*d mod p.  hits[a] = 0 means
    the AP starting at a is contained in I1 = [s1, p-1].
    """
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


def fit_all_d3(p, R, s3, d3_list, inv_cache):
    """For a fixed residual set R (sorted list), test over d3 in d3_list
    whether R fits inside some d3-AP of size s3 (circular-arc fit).

    Returns list of (d3, fits, lmin, a3_example) - a3 chosen as the tight
    window start (lattice coords) mapped back by *d3.
    """
    m = len(R)
    res = []
    if not d3_list:
        return res
    if m == 0:
        for d3 in d3_list:
            res.append((d3, True, 0, 0))
        return res
    if m == 1:
        for d3 in d3_list:
            res.append((d3, True, 1, R[0]))
        return res
    Rnp = np.array(R, dtype=np.int64)
    d3np = np.array(d3_list, dtype=np.int64)
    inv = np.array([inv_cache[d] for d in d3_list], dtype=np.int64)
    M = (Rnp[:, None] * inv[None, :]) % p
    Ms = np.sort(M, axis=0)
    gaps = np.diff(Ms, axis=0)
    wrapg = (Ms[0] + p - Ms[-1]) % p
    G = np.maximum(gaps.max(axis=0), wrapg)
    lmin = p - G + 1
    fits = lmin <= s3
    # window start in lattice coords = point following the max gap
    stack = np.vstack([gaps, wrapg[None, :]])   # (m, D): last row = wrap
    jstar = np.argmax(stack, axis=0)
    rows = np.minimum(jstar + 1, m - 1)
    cols = np.arange(stack.shape[1])
    nxt = Ms[rows, cols]              # point after gap i (per column)
    a_prime = np.where(jstar == m - 1, Ms[0], nxt)
    a3 = (d3np * a_prime) % p
    for i, d3 in enumerate(d3_list):
        res.append((d3, bool(fits[i]), int(lmin[i]), int(a3[i])))
    return res


# ------------------------------------------------------------ lemma search

def lemma_search(p, mode="distinct"):
    """Exhaustive search for 3-AP coverings of Z_p at the critical sizes.

    mode='distinct': pairwise +-distinct differences (the Lemma A/B claim).
    mode='relaxed' : differences may repeat / equal +-1 (rigid families).

    Returns dict with solutions, near-miss margin, and the viable census.
    """
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
    drange = list(range(2, dmax + 1)) if mode == "distinct" else list(range(1, dmax + 1))
    inv_cache = {d: pow(d, -1, p) for d in range(1, dmax + 1)}

    solutions = []
    best_margin = None   # min over configs of (lmin - s3); < 0 => solution
    census = []          # viable (s1,s2,s3,d2,a2,h) - the one-foot landscape
    census_dh = {}       # (d2, h) -> count, UNCAPPED (full landscape)
    hits_cache = {}

    def get_hits(s1, s2, d2):
        key = (s1, s2, d2)
        if key not in hits_cache:
            hits_cache[key] = positions_by_hits(p, s1, s2, d2)
        return hits_cache[key]

    n_viable = 0
    for ms in multisets:
        ov = sum(ms) - p
        for (s1, s2, s3) in sorted(set(permutations(ms))):
            if s1 == 1:
                # size-1 AP: difference undefined; covered by the other
                # normalizations (exhaustive - see header note on p=5).
                continue
            I1 = frozenset(range(s1, p))
            for d2 in drange:
                hits = get_hits(s1, s2, d2)
                for a2, h in enumerate(hits):
                    if h > ov:
                        continue
                    n_viable += 1
                    census_dh[(d2, h)] = census_dh.get((d2, h), 0) + 1
                    if len(census) < 4000:
                        census.append((s1, s2, s3, d2, a2, h))
                    A2 = ap_set(p, a2, d2, s2)
                    R = sorted(I1 - A2)
                    if len(R) > s3:
                        continue    # guard; cannot occur when h <= ov
                    d3list = [d for d in drange
                              if not (mode == "distinct"
                                      and pmrep(d, p) == pmrep(d2, p))]
                    for (d3, fits, lmin, a3) in fit_all_d3(p, R, s3, d3list, inv_cache):
                        marg = lmin - s3
                        if best_margin is None or marg < best_margin[0]:
                            best_margin = (marg, ms, (s1, s2, s3),
                                           d2, a2, d3, lmin)
                        if not fits:
                            continue
                        A1 = frozenset(range(s1))
                        A3 = ap_set(p, a3, d3, s3)
                        u = A1 | A2 | A3
                        if u != frozenset(range(p)) or len(A3) != s3:
                            # re-verification failure would indicate a bug
                            raise AssertionError(
                                f"fit-check bug p={p} {(s1,s2,s3)} "
                                f"d2={d2} a2={a2} d3={d3} a3={a3}")
                        # cross-ratio invariant (scaling-invariant key):
                        # with d1 = 1 the ratios are d2, d3, d2*d3^{-1}.
                        xr = sorted([pmrep(d2, p), pmrep(d3, p),
                                     pmrep(d2 * inv_cache[d3] % p, p)])
                        solutions.append({
                            "sizes": list(ms), "assign": [s1, s2, s3],
                            "d2": d2, "a2": a2, "d3": d3, "a3": a3,
                            "h2": h, "ov": ov,
                            "xr": xr,
                            "A1": sorted(range(s1)),
                            "A2": sorted(A2), "A3": sorted(A3),
                        })
    # size-multiset breakdown over ALL solutions (uncapped)
    bym = {}
    for s in solutions:
        key = tuple(s["sizes"])
        bym[key] = bym.get(key, 0) + 1
    red_sizes = ({2 * k, 2 * k + 1} if cls == 1 else {2 * k + 1, 2 * k + 2})
    n_red = sum(1 for s in solutions
                if all(sz in red_sizes for sz in s["sizes"]))
    return {"p": p, "k": k, "cls": cls, "mode": mode,
            "multisets": [list(m) for m in multisets],
            "n_viable": n_viable,
            "n_solutions": len(solutions),
            "n_solutions_reduction_sizes": n_red,
            "sol_by_sizes": {str(k_): v for k_, v in bym.items()},
            "solutions": solutions[:400],
            "best_margin": None if best_margin is None else {
                "margin": best_margin[0], "sizes": list(best_margin[1]),
                "assign": list(best_margin[2]), "d2": best_margin[3],
                "a2": best_margin[4], "d3": best_margin[5],
                "lmin": best_margin[6]},
            "census_d_set": sorted({d for (d, h) in census_dh}),
            "census_dh": {f"d={d},h={h}": c for (d, h), c in
                          sorted(census_dh.items())},
            "census_sample": census[:400]}


# ------------------------------------------------------- trichotomy census

def trichotomy_census(p):
    """For all critical (s1, s): which differences d admit an AP of size s
    contained in I1 = [s1, p-1]; verify d in {2, (p-1)/2} for k >= 2 (cls 5)
    / k >= 3 (cls 1); record exceptions and structure data."""
    k = k_of(p)
    cls = 1 if p % 6 == 1 else 5
    size_opts = [2 * k, 2 * k + 1, 2 * k + 2] if cls == 1 else [2 * k + 1, 2 * k + 2]
    dmax = (p - 1) // 2
    data = {}
    for s1 in size_opts:
        for s in size_opts:
            if s > p - s1:
                continue
            for d in range(2, dmax + 1):
                hits = positions_by_hits(p, s1, s, d)
                cont = [a for a in range(p) if hits[a] == 0]
                if cont:
                    data[f"s1={s1},s={s},d={d}"] = len(cont)
    # claim check
    claim_k = 2 if cls == 5 else 3
    exceptions = {}
    if k >= claim_k:
        for key in data:
            d = int(key.split("d=")[1])
            if d not in (2, dmax):
                exceptions[key] = data[key]
    # two-block structure verification for d = dmax:
    # from the two-lattice decomposition (p = 2*dmax + 1), the AP is
    # [a - ceil(s/2) + 1, a]  U  [a + dmax - floor(s/2) + 1, a + dmax]
    # (mod p): even-index block ends at a, odd-index block at a+dmax.
    twoblock_ok = True
    twoblock_detail = []
    for s1 in size_opts:
        for s in size_opts:
            if s < 2 or s > p - s1:
                continue
            hits = positions_by_hits(p, s1, s, dmax)
            for a in range(p):
                if hits[a] != 0:
                    continue
                A = ap_set(p, a, dmax, s)
                c1 = (s + 1) // 2
                c2 = s // 2
                arc1 = frozenset((a - m) % p for m in range(c1))
                arc2 = frozenset(((a + dmax) - m) % p for m in range(c2))
                if A != (arc1 | arc2):
                    twoblock_ok = False
                    twoblock_detail.append({"s1": s1, "s": s, "a": a,
                                            "A": sorted(A)})
    return {"k": k, "cls": cls, "data": data,
            "n_entries": len(data),
            "d_set": sorted({int(k_.split("d=")[1]) for k_ in data}),
            "claim_holds": (k < claim_k) or (not exceptions),
            "exceptions": exceptions,
            "twoblock_ok": twoblock_ok,
            "twoblock_bad": twoblock_detail[:10]}


# ------------------------------------------------------- capacity formula

def capacity_check(p):
    """Verify the per-fiber capacity formula and the AP structure of the
    fiber bad sets, by direct enumeration at modulus p^2.

    For u in [1,p) (u mod p determines everything; lifts share c and the
    difference u^{-1} mod p) and r' in [0,p): the bad set
      {j in [0,p): 6*||u*(r' + j*p) mod p^2|| < p^2}
    must be an AP with difference u^{-1} mod p and size
      L = 2k+1 if c = u r' mod p in B_k else (2k+2 if p%6==5 else 2k).
    """
    k = k_of(p)
    N = p * p
    Bk = ball(p)
    bad_size = 0
    bad_struct = 0
    checked = 0
    for u in range(1, p):
        winv = pow(u, -1, p)
        for rp in range(p):
            c = (u * rp) % p
            onball = min(c, p - c) <= k
            expected = (2 * k + 1) if onball else (
                2 * k + 2 if p % 6 == 5 else 2 * k)
            S = []
            for j in range(p):
                v = (u * (rp + j * p)) % N
                if 6 * min(v, N - v) < N:
                    S.append(j)
            checked += 1
            if len(S) != expected:
                bad_size += 1
                if bad_size <= 5:
                    log(f"  capacity SIZE mismatch p={p} u={u} r'={rp}: "
                        f"got {len(S)} expected {expected}")
            # AP structure with difference u^{-1}: scale by u^{-1}: the set
            # S*u must be a circular arc of consecutive residues (the
            # t-space bad set is an arc through 0, so it may wrap)
            sc = sorted((j * u) % p for j in S)
            if sc:
                gaps = [sc[i + 1] - sc[i] for i in range(len(sc) - 1)]
                wrapg = (sc[0] + p - sc[-1]) % p
                allg = gaps + [wrapg]
                big = [g for g in allg if g != 1]
                okstruct = (len(S) <= 1 or len(S) == p or
                            (len(big) == 1 and big[0] == p - len(S) + 1))
                if not okstruct:
                    bad_struct += 1
                    if bad_struct <= 5:
                        log(f"  capacity STRUCT mismatch p={p} u={u} "
                            f"r'={rp}: scaled={sc}")
    return {"p": p, "k": k, "checked": checked,
            "bad_size": bad_size, "bad_struct": bad_struct}


# ------------------------------------------------------------- cell anchor

def cell_check(p):
    """Independent brute-force: no (kernel=p, 3 units) covering of Z_{p^2}
    by the four bad sets.  Units enumerated up to the negation symmetry
    (B_{-u} = B_u since the norm is even); same-+-class unit pairs reduce
    to 3-set coverings, checked in the same pass via pairs."""
    N = p * p
    Bk = ball(p)
    FULL = (1 << N) - 1
    # kernel bad set: all fibers over B_k
    bp = 0
    for r in Bk:
        for j in range(p):
            bp |= 1 << (r + j * p)
    reps = [u for u in range(1, (N // 2) + 1) if gcd(u, N) == 1]
    bu = {}
    for u in reps:
        m = 0
        for kk_ in range(N):
            v = (u * kk_) % N
            if 6 * min(v, N - v) < N:
                m |= 1 << kk_
        bu[u] = m
    found = []
    # triples of distinct +-classes (covers all distinct-speed triples:
    # same-class pairs give identical bad sets -> 3-set covering, covered
    # by the pair loop below together with the triple loop's 2-subsets)
    for u1, u2, u3 in combinations(reps, 3):
        if (bu[u1] | bu[u2] | bu[u3] | bp) == FULL:
            found.append([u1, u2, u3])
            if len(found) > 20:
                break
    # pairs + kernel (i.e. a triple where two speeds share a +-class)
    pairs_found = []
    for u1, u2 in combinations(reps, 2):
        if (bu[u1] | bu[u2] | bp) == FULL:
            pairs_found.append([u1, u2])
            if len(pairs_found) > 20:
                break
    return {"p": p, "N": N, "n_reps": len(reps),
            "triple_coverings": found, "pair_coverings": pairs_found}


# ----------------------------------------------------------- rigid zone

def rigid_zone_check(p):
    """n=5 side: min number of +-distinct dilates of B_k covering Z_p
    (3 or 4 dilates).  Committed record: 3 at {7,13}; 4 at {17,19,37};
    >= 5 elsewhere."""
    Bk = ball(p)
    classes = list(range(1, (p + 1) // 2))
    balls = {u: frozenset((u * r) % p for r in Bk) for u in classes}
    full = frozenset(range(p))
    res = {"p": p, "ball_size": len(Bk), "min_j": None}
    for j in (3, 4):
        for comb in combinations(classes, j):
            cov = set()
            for u in comb:
                cov |= balls[u]
            if len(cov) == p:
                res["min_j"] = j
                res["witness"] = list(comb)
                return res
    return res


# ------------------------------------------- independent brute-force check

def independent_brute(p):
    """Full independent cross-check of the lemma search logic for small p:
    same normalization (scale d1 -> 1, translate a1 -> 0) but a2, a3
    enumerated brute-force over all p^2 placements with direct set-union
    checks (no prune, no fit check).  Confirms the prune+fit machinery."""
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
                        for A3 in A3s:
                            if (A1 | A2 | A3) == FULL:
                                n += 1     # one per (d2,a2,d3) config,
                                break      # matching the search's counting
    return {"p": p, "n_solutions": n}


# ------------------------------------------------------------------- main

def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 250
    primes = sieve(PMAX)
    pA = [p for p in primes if p % 6 == 1]
    pB = [p for p in primes if p % 6 == 5]
    OUT["meta"] = {"pmax": PMAX, "n_A": len(pA), "n_B": len(pB),
                   "pA": pA, "pB": pB,
                   "numpy": True}

    # ---- Lemma A / Lemma B (distinct mode) + relaxed mode ----
    for name, plist in (("lemmaA", pA), ("lemmaB", pB)):
        for p in plist:
            t0 = time.time()
            r = lemma_search(p, "distinct")
            OUT[name][p] = {kk_: v for kk_, v in r.items()
                            if kk_ != "census_sample"}
            OUT[name][p]["census_size"] = len(r["census_sample"])
            rr = lemma_search(p, "relaxed")
            OUT["relaxed"][p] = {"n_solutions": rr["n_solutions"],
                                 "n_viable": rr["n_viable"],
                                 "solutions": rr["solutions"][:200]}
            log(f"{name} p={p} (k={r['k']}): distinct "
                f"solutions={r['n_solutions']} viable={r['n_viable']} "
                f"margin={r['best_margin']['margin'] if r['best_margin'] else None} | "
                f"relaxed solutions={rr['n_solutions']} "
                f"({time.time()-t0:.1f}s)")

    # ---- Trichotomy ----
    for p in pA + pB:
        OUT["trichotomy"][p] = trichotomy_census(p)
        t = OUT["trichotomy"][p]
        log(f"trich p={p}: d_set={t['d_set']} claim={t['claim_holds']} "
            f"twoblock={t['twoblock_ok']} exceptions={len(t['exceptions'])}")

    # ---- Capacity ----
    for p in [q for q in pA + pB if q <= 37]:
        OUT["capacity"][p] = capacity_check(p)
        c = OUT["capacity"][p]
        log(f"capacity p={p}: checked={c['checked']} bad_size={c['bad_size']} "
            f"bad_struct={c['bad_struct']}")

    # ---- Cells ----
    for p in [q for q in pA + pB if q <= 31]:
        t0 = time.time()
        OUT["cells"][p] = cell_check(p)
        c = OUT["cells"][p]
        log(f"cell p={p}: reps={c['n_reps']} triples={len(c['triple_coverings'])} "
            f"pairs={len(c['pair_coverings'])} ({time.time()-t0:.1f}s)")

    # ---- Independent brute-force cross-validation of the search logic ----
    OUT["brute"] = {}
    for p in [q for q in pA + pB if 7 <= q <= 19]:
        t0 = time.time()
        b = independent_brute(p)
        src = OUT["lemmaA"].get(p) or OUT["lemmaB"].get(p)
        b["matches_search"] = (b["n_solutions"] == src["n_solutions"])
        OUT["brute"][p] = b
        log(f"brute p={p}: n={b['n_solutions']} "
            f"match={b['matches_search']} ({time.time()-t0:.1f}s)")

    # ---- Rigid zone (n=5 side) ----
    rz = {}
    for p in [q for q in pA + pB if q <= 61 and q >= 7]:
        rz[p] = rigid_zone_check(p)
    OUT["rigid_zone"] = {p: v for p, v in rz.items()}
    j3 = [p for p, v in rz.items() if v["min_j"] == 3]
    j4 = [p for p, v in rz.items() if v["min_j"] == 4]
    log(f"rigid zone: j=3 at {j3}, j=4 at {j4}")

    # ---- Anchors ----
    A = {}
    # (a) p=17 relaxed family (1,2,2), sizes (6,6,6), ov=1
    r17 = OUT["relaxed"].get(17, {})
    fam17 = [s for s in r17.get("solutions", [])
             if s["sizes"] == [6, 6, 6] and s["ov"] == 1
             and sorted(s["xr"]) == [1, 2, 2]]
    A["p17_rigid_family_found"] = len(fam17) > 0
    if fam17:
        s = fam17[0]
        A["p17_rigid_family_example"] = [s["A1"], s["A2"], s["A3"]]
    # (b) distinct searches: zero solutions everywhere
    A["lemmaA_zero_all"] = all(v["n_solutions"] == 0 for v in OUT["lemmaA"].values())
    A["lemmaB_zero_all"] = all(v["n_solutions"] == 0 for v in OUT["lemmaB"].values())
    # (c) p=11 parity family in relaxed
    r11 = OUT["relaxed"].get(11, {})
    A["p11_parity_family"] = any(sorted(s["xr"]) == [1, 2, 2]
                                 for s in r11.get("solutions", []))
    # (d) p=23 wrap counterexample arc [5,22], d=8, s=8
    hits23 = positions_by_hits(23, 5, 8, 8)
    A["p23_wrap_ap_contained"] = hits23[5] == 0
    A["p23_wrap_ap_set"] = sorted(ap_set(23, 5, 8, 8))
    # (e) p=13 trichotomy exception d=5 present
    A["p13_exception_d5"] = any(k_.startswith("s1=4,s=4,d=5")
                                for k_ in OUT["trichotomy"][13]["data"])
    # (f) p=17 census d in {2,8}
    t17 = OUT["trichotomy"][17]
    A["p17_census_d28"] = set(t17["d_set"]) <= {2, 8}
    # (g) capacity table anchors
    A["cap17"] = (OUT["capacity"][17]["bad_size"] == 0 and
                  OUT["capacity"][17]["bad_struct"] == 0)
    A["cap19"] = (OUT["capacity"][19]["bad_size"] == 0 and
                  OUT["capacity"][19]["bad_struct"] == 0)
    # (h) cells empty
    A["cells_empty"] = all(
        not v["triple_coverings"] and not v["pair_coverings"]
        for v in OUT["cells"].values())
    # (h2) brute force agrees with the pruned search
    A["brute_matches"] = all(v["matches_search"] for v in OUT["brute"].values())
    # (h3) trichotomy aggregate: d in {2,(p-1)/2} for k>=3 (cls 1), k>=2 (cls 5)
    A["trich_claim_all"] = all(v["claim_holds"] for v in OUT["trichotomy"].values())
    # (i) rigid zone
    A["rigid_j3"] = sorted(j3)
    A["rigid_j4"] = sorted(j4)
    OUT["anchors"] = A
    log(f"anchors: { {k_: v for k_, v in A.items() if k_ != 'p23_wrap_ap_set'} }")

    with open("/home/z/my-project/scripts/out_ap_lemmas.json", "w") as f:
        json.dump(OUT, f, indent=1, default=str)
    log("JSON written to scripts/out_ap_lemmas.json")


if __name__ == "__main__":
    main()
