#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kneser bridge — Step 1: the kernel-world map at the composite micro-cases.

Framework (byte-compatible with the committed battery N5..N9):
  n = 5 speeds, T = n+1 = 6; pair (p,q) with N = v_p + v_q; effective bad
  sets B_w = {k in Z_N : T*dist(w*k, N) < N} for the pair residue (B_p=B_q)
  and the three other speeds; the pair certifies tau(N) >= N/6 iff the FOUR
  bad sets do NOT cover Z_N.  "Kernel" speed: gcd(w, N) > 1.

This script delivers the strategy memo's requested first concrete step
(N = 49 stabilizer structure, generalised to p^2 and products), namely:
  1. ANCHORS: hand-verified covering facts + exact re-verification of the
     committed rigid-zone counts {7,13,17,19,37} = {20,68,16,64,32} and the
     committed composite unit-safety of 49/77/91 (pair-normalized tuples).
  2. Full covering maps (3-set and 4-set, over distinct bad sets) at
     N in {49, 77, 91, 121, 169}, classified by speed-type cell.
  3. Lemma bank (each proved on paper, verified exhaustively here):
     L-FIB  per-speed pullback: B_u^{(N)} = pi_{N/d}^{-1}(B_{u/d}^{(N/d)}),
            d = gcd(u,N)  (heterogeneous fibration; kernel bad sets are
            pullbacks of UNIT bad sets at the reduced modulus).
     L-HOM  homogeneous transfer (lift law, k-set version): if all speeds
            share divisor d, covering at N <=> covering at N/d.
     L-CRT  transverse covering at N = pq: a kernel-only family of p-speeds
            and q-speeds covers Z_{pq} iff the p-speeds' Z_q-images cover
            Z_q OR the q-speeds' Z_p-images cover Z_p.
     L-FC   fiber counting at N = p^2 (p >= 7): with >= 2 kernel speeds,
            covering <=> the kernel speeds' transferred sets cover Z_p;
            with exactly 1 kernel + 2 units (3-set): never covers.
  4. Stabilizer / Kneser tables: H(B_u + B_v) for every distinct-set pair,
     Kneser bound tightness, triple-sum stabilizers at N = 49.
  5. Multiplicative intersection spectrum |A_N ^ lam*A_N| over all units
     lam — the datum for the unit-world bridge (Lev/Konyagin/Ruzsa).
  6. L-FIB checks at N = 1001 = 7*11*13 (three-prime fibration tree).

Hand-verified anchors committed as assertions:
  B_7 U B_14 U B_21 = Z_49   (covering kernel triple exists at N=49)
  B_7 U B_14 U B_35 != Z_49  (classes {1,6},{2,5},{2,5} miss one class)
"""

from math import gcd
from itertools import combinations_with_replacement
from collections import defaultdict, Counter
import json
import time

T = 6  # n = 5 -> threshold 1/(n+1)

# ------------------------------------------------------------- conventions


def dist(x, N):
    r = x % N
    return r if r <= N - r else N - r


def bad_set(w, N):
    """Committed convention (byte-compatible with N5..N9): STRICT."""
    return frozenset(k for k in range(N) if T * dist(w * k, N) < N)


def arc(N):
    return bad_set(1, N)


def covers(sets, N):
    u = set()
    for s in sets:
        u |= s
    return len(u) == N


def stab(N, S):
    """Translation stabilizer {g : S+g = S} (a subgroup of Z_N)."""
    SS = set(S)
    return frozenset(g for g in range(N)
                     if {(x + g) % N for x in SS} == SS)


def sumset(X, Y, N):
    return frozenset((x + y) % N for x in X for y in Y)


def kneser(X, Y, N):
    """Kneser report for the pair (X, Y): |X+Y| >= |X+H| + |Y+H| - |H|."""
    S = sumset(X, Y, N)
    H = stab(N, S)

    def sat(Z):
        out = set()
        for x in Z:
            for h in H:
                out.add((x + h) % N)
        return frozenset(out)

    bound = len(sat(X)) + len(sat(Y)) - len(H)
    return dict(h=len(H), s=len(S), bound=bound,
                tight=(len(S) == bound), slack=len(S) - bound)


def factors(N):
    f = defaultdict(int)
    m, p = N, 2
    while p * p <= m:
        while m % p == 0:
            f[p] += 1
            m //= p
        p += 1
    if m > 1:
        f[m] += 1
    return dict(f)


def type_of(w, N):
    g = gcd(w, N)
    return "U" if g == 1 else "K%d" % g

# ----------------------------------------------------------------- anchors


def hand_anchors():
    B = {w: bad_set(w, 49) for w in (1, 7, 14, 21, 35)}
    assert covers([B[7], B[14], B[21]], 49), "(7,14,21) must cover Z_49"
    assert not covers([B[7], B[14], B[35]], 49), "(7,14,35) must NOT cover"
    assert len(B[7]) == 21 and len(B[1]) == 17
    assert B[7] == frozenset(t for t in range(49) if t % 7 in (0, 1, 6))
    assert B[14] == frozenset(t for t in range(49) if t % 7 in (0, 3, 4))
    assert B[21] == frozenset(t for t in range(49) if t % 7 in (0, 2, 5))


def normalized_unit_count(N):
    """# of pair-normalized unit 4-tuples (1,b,c,d), b<=c<=d, covering Z_N."""
    units = [w for w in range(1, N) if gcd(w, N) == 1]
    B = {w: bad_set(w, N) for w in units}
    B1 = B[1]
    c = 0
    n = len(units)
    for i in range(n):
        sb = B1 | B[units[i]]
        for j in range(i, n):
            sc = sb | B[units[j]]
            for l in range(j, n):
                if len(sc | B[units[l]]) == N:
                    c += 1
    return c


def committed_anchors():
    """Re-verify the committed record exactly (cross-validation discipline).

    Rigid zone at n=5 (T=6): unit-capable moduli {7,13,17,19,37} with
    pair-normalized covering counts 20/68/16/64/32; all other primes <= 47
    safe; composite micro-cases 49/77/91 unit-safe (verified, unproved).
    """
    for N, cnt in ((7, 20), (13, 68), (17, 16), (19, 64), (37, 32)):
        c = normalized_unit_count(N)
        assert c == cnt, "rigid-zone anchor N=%d: %d != %d" % (N, c, cnt)
    for N in (11, 23, 29, 31, 41, 43, 47):
        c = normalized_unit_count(N)
        assert c == 0, "safe-prime anchor N=%d: %d coverings" % (N, c)
    for N in (49, 77, 91):
        c = normalized_unit_count(N)
        assert c == 0, "composite unit-safety anchor N=%d: %d coverings" % (N, c)

# ------------------------------------------------------------ distinct sets


def build_distinct(N):
    """Distinct bad sets at N with representative speed, type, transferred
    (reduced-modulus) bad set for kernel speeds, and fiber-count cap."""
    speeds = range(1, N)
    Bs = {w: bad_set(w, N) for w in speeds}
    seen = {}
    for w in speeds:
        seen.setdefault(Bs[w], w)
    items = []
    for s, w in sorted(seen.items(), key=lambda kv: kv[1]):
        g = gcd(w, N)
        info = dict(speed=w, g=g, size=len(s))
        if g > 1:
            M = N // g
            info["M"] = M
            info["tset"] = bad_set(w // g, M)  # L-FIB transferred set
        items.append((w, s, type_of(w, N), info))
    # per-fiber cap c_p for unit sets (max #elements in any coset of the
    # relevant fibration) -- used by the L-FC condition.
    fac = factors(N)
    caps = {}
    for p in fac:
        M = N // p if len(fac) == 1 and fac[p] == 2 else None
        if M is None:
            continue
        c = 0
        for (w, s, t, info) in items:
            if t == "U":
                cnt = Counter(x % M for x in s)
                c = max(c, max(cnt.values()))
        caps[p] = c
    return items, caps

# -------------------------------------------------------------- predictions


def predict(N, types_, infos, k, caps):
    """(prediction, basis) for a k-multiset of distinct-set types at N.
    None prediction = open cell (measured, not proved here)."""
    nK = sum(1 for t in types_ if t != "U")
    if nK == 0:
        return None, "unit cell (open)"
    fac = factors(N)
    if len(fac) == 1 and list(fac.values()) == [2]:
        p = list(fac)[0]
        c = caps.get(p, 99)
        nU = k - nK
        if nK >= 2 and nU * c < p:
            ts = [info["tset"] for t, info in zip(types_, infos) if t != "U"]
            return covers(ts, p), "L-FC/L-HOM (>=2 kernel, p^2)"
        if nK == 1 and k == 3 and 2 * c < p:
            return False, "L-FC (1 kernel + 2 units, p^2)"
        return None, "open (cap condition fails)"
    if len(fac) == 2 and all(v == 1 for v in fac.values()):
        p, q = sorted(fac)
        if nK == k:
            side_q = covers([info["tset"] for t, info in zip(types_, infos)
                             if gcd(info["speed"], N) == p], q)
            side_p = covers([info["tset"] for t, info in zip(types_, infos)
                             if gcd(info["speed"], N) == q], p)
            return (side_q or side_p), "L-CRT (kernel-only, pq)"
        return None, "open (kernel+units at pq)"
    return None, "open (unhandled modulus shape)"


def cell_map(N, k, items, caps):
    """Covering map over k-multisets of DISTINCT bad sets, by type cell,
    with lemma predictions and mismatch accounting."""
    masks = [sum(1 << t for t in s) for (w, s, t, info) in items]
    full = (1 << N) - 1
    types = [t for (w, s, t, info) in items]
    infos = [info for (w, s, t, info) in items]
    stats = defaultdict(lambda: dict(total=0, covering=0, pred_ok=0,
                                     pred_mismatch=0, open=True, basis="",
                                     kern_expl=0, unit_assist=0))
    examples = defaultdict(list)
    incidents = []
    for combo in combinations_with_replacement(range(len(items)), k):
        m = 0
        for i in combo:
            m |= masks[i]
        cov = (m == full)
        ct = tuple(sorted(types[i] for i in combo))
        ci = [infos[i] for i in combo]
        cty = [types[i] for i in combo]
        pred, basis = predict(N, cty, ci, k, caps)
        st = stats[ct]
        st["total"] += 1
        st["open"] = pred is None
        st["basis"] = basis
        if cov:
            st["covering"] += 1
            ex = sorted(items[i][0] for i in combo)
            if len(examples[ct]) < 4:
                examples[ct].append(ex)
        if pred is None:
            if cov:
                # is the kernel subfamily alone responsible for the covering?
                kern = [masks[i] for i in combo if types[i] != "U"]
                if kern:
                    km = 0
                    for mm in kern:
                        km |= mm
                    if km == full:
                        st["kern_expl"] += 1
                    else:
                        st["unit_assist"] += 1
                        incidents.append(("UNIT-ASSISTED COVERING",
                                          sorted(items[i][0] for i in combo)))
                else:
                    st["unit_assist"] += 1
                    incidents.append(("PURE-UNIT COVERING",
                                      sorted(items[i][0] for i in combo)))
        elif pred == cov:
            st["pred_ok"] += 1
        else:
            st["pred_mismatch"] += 1
            incidents.append(("PREDICTION MISMATCH",
                              sorted(items[i][0] for i in combo),
                              "pred=%s" % pred, "actual=%s" % cov))
    return dict(stats={str(k): v for k, v in stats.items()},
                examples={str(k): v for k, v in examples.items()},
                incidents=incidents[:40], n_combos=sum(
                    v["total"] for v in stats.values()))

# ------------------------------------------------------ lemma verification


def verify_lfib(N):
    """L-FIB: B_u^{(N)} == pi_{N/d}^{-1}(B_{u/d}^{(N/d)}), d=gcd(u,N)>1."""
    checked = bad = 0
    for u in range(1, N):
        d = gcd(u, N)
        if d == 1:
            continue
        M = N // d
        BM = bad_set(u // d, M)
        lifted = frozenset(t for t in range(N) if (t % M) in BM)
        checked += 1
        if lifted != bad_set(u, N):
            bad += 1
            if bad <= 5:
                print("   L-FIB FAIL N=%d u=%d" % (N, u))
    return checked, bad


def verify_lhom(N, k):
    """L-HOM: (d*a_1..d*a_k) covers N iff (a_1..a_k) covers N/d."""
    bad = checked = 0
    divs = sorted(set(gcd(u, N) for u in range(1, N) if gcd(u, N) > 1))
    for d in divs:
        M = N // d
        if M < 2:
            continue
        for am in combinations_with_replacement(range(1, M), k):
            checked += 1
            covN = covers([bad_set(d * a, N) for a in am], N)
            covM = covers([bad_set(a, M) for a in am], M)
            if covN != covM:
                bad += 1
                if bad <= 5:
                    print("   L-HOM FAIL N=%d d=%d %s" % (N, d, am))
    return checked, bad


def verify_lcrt(N, k):
    """L-CRT at N=pq: kernel-only k-multisets cover iff p-side covers Z_q
    or q-side covers Z_p."""
    fac = sorted(factors(N))
    assert len(fac) == 2 and factors(N)[fac[0]] == 1 and factors(N)[fac[1]] == 1
    p, q = fac
    ps = set(u for u in range(1, N) if gcd(u, N) == p)
    qs = set(u for u in range(1, N) if gcd(u, N) == q)
    speeds = sorted(ps | qs)
    bad = checked = 0
    for sm in combinations_with_replacement(speeds, k):
        checked += 1
        cov = covers([bad_set(u, N) for u in sm], N)
        sideq = covers([bad_set(u // p, q) for u in sm if u in ps], q)
        sidep = covers([bad_set(u // q, p) for u in sm if u in qs], p)
        if cov != (sideq or sidep):
            bad += 1
            if bad <= 5:
                print("   L-CRT FAIL N=%d %s" % (N, sm))
    return checked, bad

# ----------------------------------------------------- Kneser / stabilizers


def kneser_pair_table(N, items):
    agg = defaultdict(Counter)
    for i in range(len(items)):
        for j in range(i, len(items)):
            (w1, s1, t1, _), (w2, s2, t2, _) = items[i], items[j]
            r = kneser(s1, s2, N)
            agg[tuple(sorted((t1, t2)))][(r["h"], r["s"], r["bound"],
                                         r["tight"])] += 1
    return {str(k): {str(kk): vv for kk, vv in sorted(v.items())}
            for k, v in agg.items()}


def kneser_triple_table(N, items):
    """Triple-sum stabilizers over all distinct-set 3-multisets (N=49)."""
    agg = defaultdict(Counter)
    n = len(items)
    for combo in combinations_with_replacement(range(n), 3):
        (w1, s1, t1, _), (w2, s2, t2, _), (w3, s3, t3, _) = (
            items[combo[0]], items[combo[1]], items[combo[2]])
        S12 = sumset(s1, s2, N)
        S = sumset(S12, s3, N)
        H = stab(N, S)
        cell = tuple(sorted((t1, t2, t3)))
        agg[cell][(len(H), len(S))] += 1
    return {str(k): {str(kk): vv for kk, vv in sorted(v.items())}
            for k, v in agg.items()}

# ------------------------------------------------------------------ spectra


def spectrum(N):
    A = arc(N)
    units = [w for w in range(1, N) if gcd(w, N) == 1]
    spec = {}
    for lam in units:
        spec[lam] = len(A & frozenset((lam * a) % N for a in A))
    return spec

# --------------------------------------------------------------------- main


def main():
    t0 = time.time()
    print("=== Kneser bridge, Step 1: kernel-world map ===")
    print("conventions: T=%d (n=5), strict bad-set inequality, "
          "4-set = tau-covering (pair + 3 others)\n" % T)

    hand_anchors()
    print("hand anchors: PASS  ((7,14,21) covers Z_49; (7,14,35) does not)")

    committed_anchors()
    print("committed anchors: PASS  (rigid-zone counts 20/68/16/64/32; "
          "primes 11,23,29,31,41,43,47 safe; 49/77/91 unit-safe)")

    report = {"T": T, "anchors": "pass",
              "rigid_zone_counts": {str(N): normalized_unit_count(N)
                                    for N in (7, 13, 17, 19, 37)},
              "moduli": {}}

    for N in (49, 77, 91, 121, 169):
        fac = factors(N)
        fstr = "*".join("%d^%d" % (p, e) if e > 1 else str(p)
                        for p, e in sorted(fac.items()))
        items, caps = build_distinct(N)
        nU = sum(1 for (w, s, t, i) in items if t == "U")
        print("\n=== N = %d = %s ===" % (N, fstr))
        print("distinct bad sets: %d (%d unit, %d kernel); |A| = %d; "
              "fiber caps %s" % (len(items), nU, len(items) - nU,
                                 len(arc(N)), caps))
        rep = dict(factorization=fstr, distinct=len(items),
                   unit_sets=nU, kernel_sets=len(items) - nU,
                   arc_size=len(arc(N)), fiber_caps={str(k): v
                                                     for k, v in caps.items()})
        for k in (3, 4):
            cm = cell_map(N, k, items, caps)
            mism = sum(v["pred_mismatch"] for v in cm["stats"].values())
            ke = sum(v["kern_expl"] for v in cm["stats"].values())
            ua = sum(v["unit_assist"] for v in cm["stats"].values())
            tot_cov = sum(v["covering"] for v in cm["stats"].values())
            print("  %d-set map: %d combos over distinct sets; coverings=%d "
                  "(kernel-explained=%d, unit-assisted=%d); prediction "
                  "mismatches: %d" % (k, cm["n_combos"], tot_cov, ke, ua, mism))
            for cell in sorted(cm["stats"]):
                v = cm["stats"][cell]
                tag = ("OPEN" if v["open"] else
                       ("PROVED" if v["pred_mismatch"] == 0 else "MISMATCH"))
                extra = ""
                if v["open"] and (v["kern_expl"] or v["unit_assist"]):
                    extra = " [kern_expl=%d unit_assist=%d]"
                    extra = extra % (v["kern_expl"], v["unit_assist"])
                print("    cell %-24s total=%-8d covering=%-6d [%s]%s"
                      % (cell, v["total"], v["covering"], tag, extra))
                if v["covering"] and cm["examples"].get(cell):
                    print("        e.g. %s" % (cm["examples"][cell][:2],))
            if ua:
                print("    !!! UNIT-ASSISTED COVERINGS EXIST: see incidents")
            if cm["incidents"]:
                for inc in cm["incidents"][:8]:
                    print("    !!! %s" % (inc,))
            rep["map_%dset" % k] = {"stats": cm["stats"],
                                    "examples": cm["examples"],
                                    "n_combos": cm["n_combos"],
                                    "prediction_mismatches": mism}
        lf = verify_lfib(N)
        lh3 = verify_lhom(N, 3)
        lh4 = verify_lhom(N, 4)
        print("  L-FIB: %d checks, %d failures | L-HOM(3): %d/%d | "
              "L-HOM(4): %d/%d" % (lf[0], lf[1], lh3[0] - lh3[1], lh3[0],
                                    lh4[0] - lh4[1], lh4[0]))
        rep["lfib"] = {"checked": lf[0], "bad": lf[1]}
        rep["lhom3"] = {"ok": lh3[0] - lh3[1], "total": lh3[0]}
        rep["lhom4"] = {"ok": lh4[0] - lh4[1], "total": lh4[0]}
        if len(fac) == 2:
            lc3 = verify_lcrt(N, 3)
            lc4 = verify_lcrt(N, 4)
            print("  L-CRT: 3-set %d/%d | 4-set %d/%d"
                  % (lc3[0] - lc3[1], lc3[0], lc4[0] - lc4[1], lc4[0]))
            rep["lcrt3"] = {"ok": lc3[0] - lc3[1], "total": lc3[0]}
            rep["lcrt4"] = {"ok": lc4[0] - lc4[1], "total": lc4[0]}
        spec = spectrum(N)
        smin = min(spec.values())
        argmin = sorted(l for l in spec if spec[l] == smin)
        print("  spectrum |A ^ lam*A|: min=%d at lam=%s (units=%d)"
              % (smin, argmin[:6], len(spec)))
        rep["spectrum_min"] = smin
        rep["spectrum_argmin"] = argmin[:12]
        rep["spectrum"] = {str(l): v for l, v in sorted(spec.items())}
        rep["kneser_pairs"] = kneser_pair_table(N, items)
        if N == 49:
            rep["kneser_triples"] = kneser_triple_table(N, items)
        report["moduli"][str(N)] = rep

    # three-prime tree: pullback identities at 1001 = 7*11*13
    lf1001 = verify_lfib(1001)
    print("\n=== N = 1001 = 7*11*13 (fibration tree) ===")
    print("  L-FIB: %d kernel speeds checked, %d failures" % lf1001)
    report["moduli"]["1001"] = {"lfib": {"checked": lf1001[0],
                                         "bad": lf1001[1]}}

    with open("/home/z/my-project/scripts/out_kneser_step1.json", "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)
    print("\nJSON report -> scripts/out_kneser_step1.json")
    print("total runtime %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
