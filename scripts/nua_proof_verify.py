#!/usr/bin/env python3
"""
nua_proof_verify.py -- targeted verification harness for the no-unit-
assistance theorem package (turn T-11, 2026-10-04).

NOT a census: no new moduli, no V-extension, no capability scan beyond the
committed record.  It verifies exactly the statements the proofs rely on:

  A. Boundary moduli N in {4, 9, 25} (small primes p in {2,3,5}):
     full 3-set and 4-set covering maps over distinct bad sets.  This is
     where the blanket "no unit assistance" statement is FALSE (N=9) and
     where the theorem's hypothesis p >= 7 is decided.
  B. Committed anchors: the T-10 covering maps at {49,77,91,121,169}
     reproduced exactly (per-cell covering counts, unit-assisted = 0).
  C. Lemma assertions:
     L1  fiber-cap closed form  cap(s) = 2*floor(s/6) + 1 + [s = 5 mod 6]
         holds EXACTLY (equality, not just <=) for every prime-indexed
         fibration of every composite modulus tested; plus the
         fiber-profile identity  B_u ^ F_j = u^{-1}(A ^ F_{u*j mod M}).
     L2  Dirichlet arc-intersection floor: |A ^ lam*A| >= 2 for every unit
         lam =/= +-1 at N in {49,77,91,121,169} (committed spectra
         re-verified exactly: 3/5/7/9/11).
     T1  no 3-multiset of unit bad sets covers Z_N (all eight moduli).
     T4  pure-unit reduction: {k in B_u : p | k} == p * B~_{u mod (N/p)}
         for every prime p | N and every unit u (all six composites);
         base cases Z_3, Z_5: every unit bad set is {0}.

Conventions byte-compatible with scripts/kneser_bridge_step1.py (T-10):
T = 6 strict; distinct bad sets deduped by frozenset, smallest-speed
representative; covering maps over combinations_with_replacement
(k-multisets of distinct sets -- the pair convention: the pair's two speeds
contribute the same bad set); type "U" (gcd(w,N)=1) / "K{g}".
Speeds range over 1..N-1 (the stuck speed w = 0 is excluded, as in T-10;
it is a kernel speed whose bad set is all of Z_N and is never needed).
"""

import json
import sys
from collections import defaultdict
from itertools import combinations_with_replacement
from math import gcd

T = 6
OUT = {}


def dist(x, N):
    r = x % N
    return min(r, N - r)


def bad_set(w, N):
    return frozenset(k for k in range(N) if T * dist(w * k, N) < N)


def arc(N):
    return bad_set(1, N)


def covers(sets, N):
    u = frozenset()
    for s in sets:
        u |= s
    return u == frozenset(range(N))


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


def cap_closed_form(s):
    """Uniform fiber cap for fibers of prime size s (Lemma 1)."""
    return 2 * (s // 6) + 1 + (1 if s % 6 == 5 else 0)


def build_distinct(N):
    speeds = range(1, N)
    Bs = {w: bad_set(w, N) for w in speeds}
    seen = {}
    for w in speeds:
        seen.setdefault(Bs[w], w)
    items = []
    for s, w in sorted(seen.items(), key=lambda kv: kv[1]):
        items.append((w, s, type_of(w, N)))
    return items


def inv(a, N):
    return pow(a, -1, N)


def fiber(N, M, j):
    return frozenset(k for k in range(N) if k % M == j)


def scale(X, lam, N):
    return frozenset((lam * x) % N for x in X)


# ------------------------------------------------------------ covering maps


def covering_map(N, items, k, full_stats=False):
    """k-multisets of distinct bad sets, by type cell.  Returns dict:
    cell -> {total, covering, kern_cover, unit_assist} (+ examples)."""
    masks = [sum(1 << t for t in s) for (w, s, t) in items]
    types = [t for (w, s, t) in items]
    full = (1 << N) - 1
    stats = defaultdict(lambda: dict(total=0, covering=0, kern_cover=0,
                                     unit_assist=0, examples=[]))
    n = len(items)
    n_combo = 0
    for combo in combinations_with_replacement(range(n), k):
        n_combo += 1
        if k == 3:
            m = masks[combo[0]] | masks[combo[1]] | masks[combo[2]]
        else:
            m = (masks[combo[0]] | masks[combo[1]] | masks[combo[2]]
                 | masks[combo[3]])
        if m == full:
            ct = tuple(sorted(types[i] for i in combo))
            st = stats[ct]
            st["total"] += 1
            st["covering"] += 1
            km = 0
            for i in combo:
                if types[i] != "U":
                    km |= masks[i]
            kc = (km == full)
            if kc:
                st["kern_cover"] += 1
            else:
                st["unit_assist"] += 1
            if len(st["examples"]) < 3:
                st["examples"].append([items[i][0] for i in combo])
        elif full_stats:
            ct = tuple(sorted(types[i] for i in combo))
            stats[ct]["total"] += 1
    return dict(stats), n_combo


# ---------------------------------------------------------------- moduli


def run_modulus(N, anchors=None, full_stats=False, report_cells=False):
    print("=" * 72)
    print("N = %d   factorization %s" % (N, factors(N)))
    items = build_distinct(N)
    nU = sum(1 for (w, s, t) in items if t == "U")
    nK = len(items) - nU
    print("distinct bad sets: %d  (%d U, %d K)" % (len(items), nU, nK))
    res = dict(distinct=len(items), nU=nU, nK=nK)
    OUT.setdefault(str(N), {}).update(res)

    # ---- covering maps (both arities) -----------------------------------
    for k in (3, 4):
        stats, n_combo = covering_map(N, items, k, full_stats=full_stats)
        cov = sum(v["covering"] for v in stats.values())
        ua = sum(v["unit_assist"] for v in stats.values())
        print("  %d-multisets: %d combos, coverings %d, unit-assisted %d"
              % (k, n_combo, cov, ua))
        for ct in sorted(stats):
            v = stats[ct]
            if v["covering"] or full_stats:
                print("    cell %-28s total=%-7d cov=%-5d kern=%-5d "
                      "uassist=%-4d ex=%s"
                      % (str(ct), v["total"], v["covering"], v["kern_cover"],
                         v["unit_assist"], v["examples"][:2]))
        if anchors and k in anchors:
            acov, acells = anchors[k]
            assert cov == acov, \
                "N=%d arity %d: covering %d != committed %d" % (N, k, cov, acov)
            for ct, n in acells.items():
                got = stats.get(ct, {}).get("covering", 0)
                assert got == n, "N=%d cell %s: %d != committed %d" % (
                    N, ct, got, n)
            print("    ANCHOR PASS (committed T-10 map reproduced exactly)")
        assert ua == 0 or N in (4, 9, 25), \
            "N=%d arity %d: unit-assisted %d (outside boundary moduli)" % (
                N, k, ua)
        res["map_%dset" % k] = {str(ct): v for ct, v in stats.items()}
        res["map_%dset_cov" % k] = cov
        res["map_%dset_ua" % k] = ua

    # ---- L1: fiber cap closed form + profile identity --------------------
    caps = {}
    for M in factors(N):
        s = N // M            # fiber size of the fibration Z_N -> Z_M
        if s < 2:
            continue
        mx = 0
        for (w, bs, t) in items:
            if t != "U":
                continue
            for j in range(M):
                mx = max(mx, len(bs & fiber(N, M, j)))
        caps[M] = mx
        want = cap_closed_form(s)
        assert mx == want, "N=%d fibration mod %d: cap %d != closed form %d" % (
            N, M, mx, want)
    if caps:
        print("  L1 fiber caps (max |B_u ^ fiber|) by fibration: %s  "
              "== closed form 2*floor(s/6)+1+[s=5 mod 6]  PASS"
              % {("mod %d (fiber size %d)" % (M, N // M)): v
                 for M, v in caps.items()})
    res["caps"] = {str(M): v for M, v in caps.items()}

    # profile identity (all units, all fibers, every fibration)
    A = arc(N)
    n_checks = 0
    for M in factors(N):
        if N // M < 2:
            continue
        for w in range(1, N):
            if gcd(w, N) != 1:
                continue
            Bu = bad_set(w, N)
            ui = inv(w, N)
            for j in range(M):
                lhs = Bu & fiber(N, M, j)
                rhs = scale(A & fiber(N, M, (w * j) % M), ui, N)
                assert lhs == rhs, "profile identity fails N=%d u=%d j=%d" % (
                    N, w, j)
                n_checks += 1
    if n_checks:
        print("  L1 profile identity B_u^F_j = u^{-1}(A^F_{uj}) : %d checks "
              "PASS" % n_checks)
    res["profile_checks"] = n_checks

    # ---- L2: spectrum |A ^ lam*A| ---------------------------------------
    if N >= 7:
        units = [w for w in range(1, N) if gcd(w, N) == 1]
        mn = None
        arg = None
        for lam in units:
            if lam % N in (1, N - 1):
                continue
            lamA = scale(A, lam, N)
            v = len(A & lamA)
            if mn is None or v < mn:
                mn, arg = v, lam
        print("  L2 spectrum: min |A ^ lam*A| over lam != +-1 = %d "
              "(argmin lam=%d)" % (mn, arg))
        res["spectrum_min"] = mn
        if N >= 31 and N % 6 != 0:
            assert mn >= 2, "N=%d: Dirichlet floor violated (%d)" % (N, mn)

    # ---- T4: pure-unit reduction through the p-multiples -----------------
    n_red = 0
    for p in factors(N):
        Mp = N // p
        if Mp < 2:
            continue
        for w in range(1, N):
            if gcd(w, N) != 1:
                continue
            lhs = frozenset(k for k in bad_set(w, N) if k % p == 0)
            rhs = frozenset(p * j for j in bad_set(w % Mp, Mp))
            assert lhs == rhs, "T4 reduction fails N=%d p=%d u=%d" % (
                N, p, w)
            n_red += 1
    if n_red:
        print("  T4 reduction {k in B_u : p|k} = p*B~_(u mod N/p) : %d checks "
              "PASS" % n_red)
    res["t4_checks"] = n_red
    return res


def main():
    print("hand anchor: (7,14,21) covers Z_49 ...", end=" ")
    B = {w: bad_set(w, 49) for w in (7, 14, 21, 35)}
    assert covers([B[7], B[14], B[21]], 49)
    assert not covers([B[7], B[14], B[35]], 49)
    print("PASS")
    OUT["hand_anchors"] = "pass"

    print("hand anchor: N=9 unit-assisted covering exists ...", end=" ")
    B9 = {w: bad_set(w, 9) for w in (1, 2, 3, 4)}
    assert B9[1] == frozenset({0, 1, 8})
    assert B9[2] == frozenset({0, 4, 5})
    assert B9[3] == frozenset({0, 3, 6})   # the only kernel-speed set
    assert B9[4] == frozenset({0, 2, 7})
    assert covers([B9[1], B9[2], B9[3], B9[4]], 9)
    assert not covers([B9[1], B9[2], B9[4]], 9)   # no 3-unit covering at 9
    print("PASS  (B_1^B_2^B_3^B_4 = Z_9, kernel subfamily {B_3} does not)")

    print("base cases for T4 at Z_3, Z_5: every unit bad set = {0} ...",
          end=" ")
    for p in (3, 5):
        for u in range(1, p):
            assert bad_set(u, p) == frozenset({0})
    print("PASS")

    print("\nDirichlet bound arithmetic (Lemma 2 prerequisite): "
          "floor(N/(c+1)) <= c ...")
    ok_from, fails = [], []
    for N in range(2, 400):
        c = (N - 1) // 6
        if N // (c + 1) <= c:
            ok_from.append(N)
        else:
            fails.append(N)
    for N in range(31, 400):
        if N % 6 != 0:
            c = (N - 1) // 6
            assert N // (c + 1) <= c, "Dirichlet case table broken at %d" % N
    print("  holds for all N in [31,400] with N = 0 mod 6 excluded; "
          "largest failures below 31: %s" % sorted(set(fails))[-8:])
    OUT["dirichlet_table"] = "pass"

    # ---------------- A. boundary moduli ---------------------------------
    print("\n" + "#" * 72)
    print("# A. Boundary moduli (small primes p = 2, 3, 5)")
    print("#" * 72)
    r4 = run_modulus(4, full_stats=True)
    assert r4["map_3set_cov"] == 0 and r4["map_4set_cov"] == 0
    r9 = run_modulus(9, full_stats=True)
    # the counterexample: exactly one covering, 4 sets, units essential
    m4 = r9["map_4set"]
    ua_cells = {ct: v for ct, v in m4.items() if v["unit_assist"]}
    assert len(ua_cells) == 1, ua_cells
    ct, v = list(ua_cells.items())[0]
    assert v["covering"] == 1 and v["unit_assist"] == 1
    assert r9["map_3set_cov"] == 0
    print("  >>> N=9: blanket no-unit-assistance is FALSE "
          "(cell %s: %s)" % (ct, v["examples"]))
    OUT["n9_counterexample"] = dict(cell=ct, example=v["examples"][0])

    r25 = run_modulus(25, full_stats=True)
    m3, m4 = r25["map_3set"], r25["map_4set"]
    # proved predictions at N=25:
    for ct, v in m3.items():
        assert v["covering"] == 0, ("N=25 3-set", ct, v)
    for ct, v in m4.items():
        if set(ct) == {"U"}:
            assert v["covering"] == 0, ("N=25 pure-unit 4-set", ct, v)
    k53u = m4.get("('K5', 'U', 'U', 'U')") or m4.get("('U', 'U', 'U', 'K5')")
    print("  >>> N=25 (1 kernel + 3 units) cell: %s"
          % (("covering=%d uassist=%d" % (k53u["covering"], k53u["unit_assist"]))
             if k53u else "absent"))

    # ---------------- B. committed anchors -------------------------------
    print("\n" + "#" * 72)
    print("# B. Committed T-10 anchors + lemma assertions "
          "(N = 49, 77, 91, 121, 169)")
    print("#" * 72)
    anchors = {
        49: {3: (1, {("K7", "K7", "K7"): 1}),
             4: (24, {("K7", "K7", "K7", "K7"): 3,
                      ("K7", "K7", "K7", "U"): 21})},
        77: {3: (1, {("K11", "K11", "K11"): 1}),
             4: (38, {("K11", "K11", "K11", "K11"): 3,
                      ("K11", "K11", "K11", "K7"): 5,
                      ("K11", "K11", "K11", "U"): 30})},
        91: {3: (3, {("K13", "K13", "K13"): 1, ("K7", "K7", "K7"): 2}),
             4: (138, {("K13", "K13", "K13", "K13"): 3,
                       ("K13", "K13", "K13", "K7"): 6,
                       ("K13", "K13", "K13", "U"): 36,
                       ("K13", "K7", "K7", "K7"): 6,
                       ("K7", "K7", "K7", "K7"): 15,
                       ("K7", "K7", "K7", "U"): 72})},
        121: {3: (0, {}), 4: (0, {})},
        169: {3: (2, {("K13", "K13", "K13"): 2}),
              4: (171, {("K13", "K13", "K13", "K13"): 15,
                        ("K13", "K13", "K13", "U"): 156})},
    }
    spectra = {49: 3, 77: 5, 91: 7, 121: 9, 169: 11}
    for N in (49, 77, 91, 121, 169):
        r = run_modulus(N, anchors=anchors[N])
        assert r["map_3set_ua"] == 0 and r["map_4set_ua"] == 0
        if N in spectra:
            assert r["spectrum_min"] == spectra[N], \
                "N=%d spectrum %s != committed %s" % (
                    N, r["spectrum_min"], spectra[N])
            print("    spectrum min matches committed record (%d) PASS"
                  % spectra[N])
        # T1: no 3-multiset of units covers
        for ct, v in r["map_3set"].items():
            if set(ct) == {"U"}:
                assert v["covering"] == 0, ("3-unit covering at", N, ct)

    print("\n" + "#" * 72)
    print("# ALL ASSERTIONS PASS")
    print("#" * 72)
    json.dump(OUT, open("/home/z/my-project/scripts/out_nua_verify.json",
                        "w"), indent=1, sort_keys=True)
    print("results -> scripts/out_nua_verify.json")


if __name__ == "__main__":
    main()
