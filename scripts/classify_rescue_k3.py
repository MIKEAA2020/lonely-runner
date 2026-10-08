#!/usr/bin/env python3
"""
classify_rescue_k3.py — follow-up to verify_top_pair_rescue_k3.py.

Classifies the rescue structure at two-sided failing witness points (n=4)
to pin the precise provable lemma and the obstruction taxonomy.

Tests:
  T1  Max-speed-pair rescue: does SOME pair containing the maximum-speed
      runner rescue?  (Such pairs are provably STUCK-FREE: N = v_max + v_x
      exceeds every other speed, so no remaining speed is 0 mod N.)
  T2  Top-speed pair (two largest speeds) rescue rate.
  T3  Mixed-pair rescue rate (control, expected 100%).
  T4  Pentagon lemma check: for pairs with 5 | N and no effective speed
      == 0 mod 5, the pentagon time k = N/5 must rescue (theoretical proof:
      residual = (N/5)*min(t,5-t) >= N/5 >= R). Verify on data.
  T5  Obstruction taxonomy for FAILING pairs: stuck (N | v_y), half-stuck
      (v_y == N/2), third-stuck (v_y in {N/3, 2N/3}), 5-blocked
      (5 | N and some effective speed == 0 mod 5), covering margin
      (6(R-1)+G >= N), small-N.
  T6  Mechanism: for winning rescues of max-speed pairs, record winning k
      and its pentagon alignment (5k mod N).
"""

from math import gcd
from itertools import combinations
from collections import Counter, defaultdict

VMAX = 12
BMAX = 50

def minres(v, x, B):
    r = B
    for vi in v:
        y = (vi * x) % B
        if y > B - y:
            y = B - y
        if y < r:
            r = y
    return r

def argmax_grid(v, B):
    best, args = -1, []
    for x in range(B):
        m = minres(v, x, B)
        if m > best:
            best, args = m, [x]
        elif m == best:
            args.append(x)
    return best, args

def rescue_detail(v, B, mstar, p, q):
    """Return (success, winning k list, pentagon-aligned winning k count)."""
    s = v[p] + v[q]
    thr = mstar * s
    hits = [k for k in range(1, s) if minres(v, k, s) * B > thr]
    return bool(hits), hits

def diagnose(v, B, mstar, p, q):
    """Obstruction taxonomy for a failing pair (p,q)."""
    N = v[p] + v[q]
    R = (mstar * N) // B + 1
    others = [t for t in range(4) if t not in (p, q)]
    eff = [v[p], v[others[0]], v[others[1]]]   # pair member + two others
    tags = []
    for t in others:
        y = v[t]
        g = gcd(y, N)
        if y % N == 0:                 # N divides y: runner stuck on grid N
            tags.append(f"stuck(v{t}={y})")
        elif g * 2 == N:
            tags.append(f"half(v{t}={y})")
        elif g * 3 == N:
            tags.append(f"third(v{t}={y})")
    if N % 5 == 0 and any(w % 5 == 0 for w in eff):
        tags.append("5blocked")
    G = sum(gcd(w, N) for w in eff)
    cover_margin = 6 * (R - 1) + G - N
    if not tags:
        if cover_margin >= 0:
            tags.append(f"covering({cover_margin})")
        else:
            tags.append("OTHER")
    small = "smallN" if N <= 12 else ""
    return N, R, tags, small

def main():
    points = []   # two-sided failing points
    for V in combinations(range(1, VMAX + 1), 4):
        g = 0
        for vi in V:
            g = gcd(g, vi)
        if g != 1:
            continue
        for B in range(6, BMAX + 1):
            mstar, args = argmax_grid(V, B)
            if mstar == 0 or not (5 * mstar < B):
                continue
            for x in args:
                rs = [(vi * x) % B for vi in V]
                L = [i for i in range(4) if rs[i] == mstar]
                R = [i for i in range(4) if rs[i] == B - mstar]
                if not (L and R) or len(L) + len(R) != 2:
                    continue
                points.append(dict(V=V, B=B, m=mstar, x=x, rs=rs,
                                   L=L, R=R, delta=B - 5 * mstar))

    print(f"two-sided exactly-2-binder failing points: {len(points)}")

    t1 = [0, 0]; t2 = [0, 0]; t3 = [0, 0]; t4 = [0, 0]
    taxonomy = Counter()
    t6 = Counter()
    t1_fail_dumps = []

    for pt in points:
        V, B, m = pt["V"], pt["B"], pt["m"]
        L, R = pt["L"], pt["R"]
        binders = L + R
        ints = [i for i in range(4) if i not in binders]
        Midx = max(range(4), key=lambda i: V[i])
        srt = sorted(range(4), key=lambda i: -V[i])

        pairs = {}
        for p in range(4):
            for q in range(p + 1, 4):
                ok, hits = rescue_detail(V, B, m, p, q)
                pairs[(p, q)] = (ok, hits)

        # T1: some pair containing the max-speed runner
        okM = any(pairs[tuple(sorted((Midx, o)))][0]
                  for o in range(4) if o != Midx)
        t1[0] += okM; t1[1] += 1
        if not okM and len(t1_fail_dumps) < 15:
            t1_fail_dumps.append(pt)

        # T2: the two largest speeds as the pair
        okTS = pairs[tuple(sorted((srt[0], srt[1])))][0]
        t2[0] += okTS; t2[1] += 1

        # T3: some mixed pair
        okMix = any(pairs[tuple(sorted((bd, it)))][0]
                    for bd in binders for it in ints)
        t3[0] += okMix; t3[1] += 1

        # T4: pentagon lemma — pairs with 5 | N, no effective speed 0 mod 5
        for (p, q), (ok, hits) in pairs.items():
            N = V[p] + V[q]
            if N % 5 == 0:
                others = [t for t in range(4) if t not in (p, q)]
                eff = [V[p], V[others[0]], V[others[1]]]
                if not any(w % 5 == 0 for w in eff):
                    t4[0] += ok; t4[1] += 1

        # T5: taxonomy of failing pairs; T6: winning-k pentagon alignment
        for (p, q), (ok, hits) in pairs.items():
            if not ok:
                N, R, tags, small = diagnose(V, B, m, p, q)
                for tg in tags:
                    taxonomy[tg.split("(")[0] + ("+" + small if small else "")] += 1
            elif Midx in (p, q):
                N = V[p] + V[q]
                pa = sum(1 for k in hits if (5 * k) % N in
                         (0, N // 5 if N % 5 == 0 else 1,
                          2 * N // 5 if N % 5 == 0 else 2))
                best = min(((5 * k) % N, N - (5 * k) % N) and
                           min((5 * k) % N, N - (5 * k) % N) for k in hits)
                t6["pentagon_exact" if best == 0 else
                   ("pentagon_near" if best <= 2 else "other")] += 1

    print("\n[T1] some pair containing the MAX-SPEED runner rescues: "
          f"{t1[0]}/{t1[1]}")
    print("[T2] the two-largest-speeds pair rescues: "
          f"{t2[0]}/{t2[1]}")
    print("[T3] some MIXED (binder+interior) pair rescues: "
          f"{t3[0]}/{t3[1]}")
    print("[T4] pentagon lemma (5|N, no effective speed 0 mod 5) rescues: "
          f"{t4[0]}/{t4[1]}")

    print("\n[T5] obstruction taxonomy over ALL failing pairs "
          "(a point contributes 6 pairs):")
    for k, c in taxonomy.most_common():
        print(f"   {k:<24} {c}")

    print("\n[T6] winning-k structure for max-speed pairs (pentagon alignment):")
    for k, c in t6.most_common():
        print(f"   {k:<24} {c}")

    if t1_fail_dumps:
        print("\n[T1] dumps of max-speed-pair rescue FAILURES:")
        for pt in t1_fail_dumps:
            V, B, m, x = pt["V"], pt["B"], pt["m"], pt["x"]
            print(f"   V={V} B={B} m={m} x={x} rs={pt['rs']} "
                  f"delta={pt['delta']}")
            for (p, q) in [(a, b) for a in range(4) for b in range(a + 1, 4)]:
                ok, hits = rescue_detail(V, B, m, p, q)
                N, R, tags, small = diagnose(V, B, m, p, q)
                print(f"      pair ({p},{q}) speeds ({V[p]},{V[q]}) N={N} "
                      f"R={R} ok={ok} tags={tags}{(' ' + small) if small else ''}")

if __name__ == "__main__":
    main()
