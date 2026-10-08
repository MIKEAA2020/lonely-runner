#!/usr/bin/env python3
"""
verify_top_pair_rescue_k3.py  —  LRC residue-language program.

Targeted example-generation / verification supporting the k=3 (n=4) TOP-PAIR
RESCUE proof attempt. This is NOT a census scan (per the reviewer directive):
it checks exactly the statements the proof attempt relies on, on a small base.

Statements under test (residue language, see worklog for vocabulary):

  E1  Covering law at every failing grid:
        B <= sum_i g_i (2*floor(m*/g_i) + 1),  g_i = gcd(v_i, B).
  E2  Binder pair-sum identity at two-sided argmax points:
        (v_u + v_w) x  ==  0 (mod B);  with h = gcd(x,B) = 1 also B | v_u+v_w.
  E3  Pair-sum-grid rescue: for pair (p,q), s = v_p+v_q, does some crossing time
        k/s satisfy   min_l ||v_l k||_s / s  >  m*/B   (cross-multiplied)?
      Stratified by shallowness delta = B - 5m*.
  E4  Specificity by pair role: interior-interior (top) pair vs binder pair vs
      mixed pairs (at two-sided, exactly-2-binder witnesses).
  E5  Binder / tie degeneracy counts.
  E6  Near-regular pentagon gap statistics at shallow two-sided witnesses
      (five gaps {m, g1, g2, g3, m}; sigma = g1 - g3 is the reflection deficit
      of the interior pair).
  E7  Dump of top-pair rescue failures with level-3 degeneracy analysis
      (level-3 problem: speeds {a, c, e} on grid N = v_i + v_j, target
      R = floor(m*N/B) + 1).

Key identity used throughout (pair-collapse): on grid s = v_p + v_q,
||v_p k||_s == ||v_q k||_s for all k, so the 4-runner minimum at a crossing
time of the pair equals the 3-effective-runner minimum.
"""

from math import gcd
from itertools import combinations
from collections import Counter, defaultdict

VMAX = 12          # speeds in [1, VMAX]
BMAX = 50          # grids in [6, BMAX]
MAX_DUMPS = 40
SHALLOW_MAX = 4    # delta <= this counts as shallow for E6

# ---------------------------------------------------------------- helpers

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

def rescue_hits(v, B, mstar, p, q):
    """All k in [1, s-1] with min_l ||v_l k||_s * B > mstar * s."""
    s = v[p] + v[q]
    thr = mstar * s
    return [k for k in range(1, s) if minres(v, k, s) * B > thr]

def level3_argmax(a, c, e, N):
    """Per-grid argmax of the 3-effective-runner problem; returns
    (best value, list of (k, residues))."""
    best, bk = -1, []
    for k in range(N):
        vals = []
        m = N
        for w in (a, c, e):
            y = (w * k) % N
            vals.append(y)
            yy = y if y <= N - y else N - y
            if yy < m:
                m = yy
        if m > best:
            best, bk = m, [(k, vals)]
        elif m == best:
            bk.append((k, vals))
    return best, bk

# ---------------------------------------------------------------- main

def main():
    n_v = 0
    n_grids = 0
    n_grids_m0 = 0
    n_grids_fail = 0
    delta_hist = Counter()
    e1_viol = []
    e2_points = 0
    e2_ident_viol = 0
    e2_h = Counter()
    e2_h1 = [0, 0]          # [h==1 points with B | v_u+v_w for all binder pairs, h==1 points]
    records = []

    for V in combinations(range(1, VMAX + 1), 4):
        g = 0
        for vi in V:
            g = gcd(g, vi)
        if g != 1:
            continue
        n_v += 1
        for B in range(6, BMAX + 1):
            n_grids += 1
            mstar, args = argmax_grid(V, B)
            if mstar == 0:
                n_grids_m0 += 1
                continue
            if not (5 * mstar < B):
                continue
            n_grids_fail += 1
            delta_hist[B - 5 * mstar] += 1

            # E1 covering law
            cov = 0
            for vi in V:
                gi = gcd(vi, B)
                cov += gi * (2 * (mstar // gi) + 1)
            if cov < B:
                e1_viol.append((V, B, mstar, cov))

            for x in args:
                rs = [(vi * x) % B for vi in V]
                cs = [y if y <= B - y else B - y for y in rs]
                L = [i for i in range(4) if rs[i] == mstar]
                R = [i for i in range(4) if rs[i] == B - mstar]
                two = bool(L) and bool(R)
                if two:
                    e2_points += 1
                    for u in L:
                        for w in R:
                            if (V[u] + V[w]) * x % B != 0:
                                e2_ident_viol += 1
                    h = gcd(x, B)
                    e2_h[h] += 1
                    if h == 1:
                        e2_h1[1] += 1
                        div_ok = all((V[u] + V[w]) % B == 0
                                     for u in L for w in R)
                        if div_ok:
                            e2_h1[0] += 1

                pr = {}
                for p in range(4):
                    for q in range(p + 1, 4):
                        pr[(p, q)] = bool(rescue_hits(V, B, mstar, p, q))

                order = sorted(range(4), key=lambda i: (-cs[i], i))
                top_pair = tuple(sorted((order[0], order[1])))
                tie_top = cs[order[1]] == cs[order[2]] if len(order) > 2 else False

                records.append(dict(
                    V=V, B=B, m=mstar, x=x, rs=rs, cs=cs, L=L, R=R,
                    two=two, n_bind=len(L) + len(R), delta=B - 5 * mstar,
                    pr=pr, top_pair=top_pair, top_ok=pr[top_pair],
                    any_ok=any(pr.values()), tie_top=tie_top))

    # ------------------------------------------------------------ aggregates

    def bucket(d):
        if d <= 4:
            return f"delta={d}"
        if d <= 8:
            return "delta 5-8"
        if d <= 16:
            return "delta 9-16"
        return "delta 17+"

    agg = defaultdict(Counter)
    for r in records:
        b = bucket(r["delta"])
        agg[b]["points"] += 1
        agg[b]["top_ok"] += r["top_ok"]
        agg[b]["any_ok"] += r["any_ok"]
        if r["two"]:
            agg[b]["two"] += 1
            agg[b]["two_top_ok"] += r["top_ok"]
            agg[b]["two_any_ok"] += r["any_ok"]
            if r["n_bind"] == 2:
                ints = [i for i in range(4) if i not in r["L"] and i not in r["R"]]
                ip = tuple(sorted((ints[0], ints[1])))
                agg[b]["n2"] += 1
                agg[b]["int_ok"] += r["pr"][ip]
                agg[b]["int_is_top"] += (ip == r["top_pair"])
                bp = any(r["pr"][tuple(sorted((u, w)))] for u in r["L"] for w in r["R"])
                mp = any(r["pr"][pq] for pq in r["pr"]
                         if (pq[0] in r["L"] + r["R"]) != (pq[1] in r["L"] + r["R"]))
                agg[b]["bind_ok"] += bp
                agg[b]["mixed_ok"] += mp

    # E5 degeneracy
    n_bind_hist = Counter(r["n_bind"] for r in records if r["two"])
    tie_hist = Counter(r["tie_top"] for r in records)
    two_hist = Counter(r["two"] for r in records)

    # E6 pentagon gaps (two-sided, exactly 2 binders, shallow)
    gapstats = []
    gap_anom = 0
    for r in records:
        if not (r["two"] and r["n_bind"] == 2 and r["delta"] <= SHALLOW_MAX):
            continue
        B, m = r["B"], r["m"]
        pos = sorted(r["rs"])
        if pos[0] != m or pos[3] != B - m:
            gap_anom += 1
            continue
        g1, g2, g3 = pos[1] - pos[0], pos[2] - pos[1], pos[3] - pos[2]
        gapstats.append(dict(V=r["V"], B=B, m=m, x=r["x"], pos=pos,
                             g1=g1, g2=g2, g3=g3, sigma=g1 - g3,
                             cs=r["cs"]))

    # E7 failures of the canonical top-pair rescue (two-sided points)
    fails = [r for r in records if r["two"] and not r["top_ok"]]

    # ------------------------------------------------------------ report

    print("=" * 72)
    print(f"k=3 TOP-PAIR RESCUE — targeted verification "
          f"(V<= {VMAX}, B in [6,{BMAX}], gcd(V)=1; NOT a census scan)")
    print("=" * 72)

    print("\n[E0] base counts")
    print(f"  speed sets                : {n_v}")
    print(f"  grids examined            : {n_grids}  (+{n_grids_m0} degenerate m*=0 skipped)")
    print(f"  failing grids (5m* < B)   : {n_grids_fail}")
    print(f"  failing witness points    : {len(records)}")
    print(f"  delta histogram (grid level): {dict(sorted(delta_hist.items()))}")

    print("\n[E1] covering law B <= sum g_i(2 floor(m*/g_i)+1) at failing grids")
    print(f"  violations: {len(e1_viol)}", e1_viol[:5])

    print("\n[E2] binder pair-sum identity at two-sided points")
    print(f"  two-sided points          : {e2_points}")
    print(f"  identity violations       : {e2_ident_viol}  (expect 0)")
    print(f"  h = gcd(x*,B) histogram   : {dict(sorted(e2_h.items()))}")
    print(f"  h==1: B | v_u+v_w holds   : {e2_h1[0]}/{e2_h1[1]}  (expect identity => all)")

    print("\n[E3/E4] pair-sum-grid rescue success by shallowness bucket")
    hdr = (f"  {'bucket':<12}{'pts':>6}{'two':>6}{'top%':>7}{'any%':>7}"
           f"{'|2bd':>6}{'int%':>7}{'int=top%':>9}{'bnd%':>7}{'mix%':>7}")
    print(hdr)
    order = ["delta=1", "delta=2", "delta=3", "delta=4",
             "delta 5-8", "delta 9-16", "delta 17+"]
    for b in order:
        c = agg[b]
        if c["points"] == 0:
            continue
        def pct(n, d):
            return f"{100.0*n/d:6.1f}" if d else "   -- "
        line = (f"  {b:<12}{c['points']:>6}{c['two']:>6}"
                f"{pct(c['top_ok'], c['points'])}{pct(c['any_ok'], c['points'])}"
                f"{c['n2']:>6}{pct(c['int_ok'], c['n2'])}"
                f"{pct(c['int_is_top'], c['n2'])}"
                f"{pct(c['bind_ok'], c['n2'])}{pct(c['mixed_ok'], c['n2'])}")
        print(line)
    print("  (top% = canonical top pair rescues; any% = some pair rescues;")
    print("   int% = interior pair at exactly-2-binder two-sided points;")
    print("   bnd% = binder pair; mix% = any mixed pair)")

    print("\n[E5] degeneracy")
    print(f"  two-sided among points        : {dict(two_hist)}")
    print(f"  binder count at two-sided pts : {dict(sorted(n_bind_hist.items()))}")
    print(f"  top-pair value tie            : {dict(tie_hist)}")

    print("\n[E6] near-regular pentagon gaps at shallow two-sided 2-binder points")
    print(f"  n = {len(gapstats)}  (positional anomalies skipped: {gap_anom})")
    if gapstats:
        import statistics as st
        g1s = [g["g1"] / g["m"] for g in gapstats]
        g2s = [g["g2"] / g["m"] for g in gapstats]
        g3s = [g["g3"] / g["m"] for g in gapstats]
        sig = [g["sigma"] for g in gapstats]
        big = [g for g in gapstats if sorted(g["cs"])[-1] >= 2 * g["m"] and
               sorted(g["cs"])[-2] >= 2 * g["m"]]
        print(f"  g1/m mean={st.mean(g1s):.3f}  g2/m mean={st.mean(g2s):.3f}  "
              f"g3/m mean={st.mean(g3s):.3f}")
        print(f"  sigma = g1-g3: mean={st.mean(sig):.2f}  "
              f"|sigma|<=2 in {sum(1 for s in sig if abs(s) <= 2)}/{len(sig)}")
        print(f"  both interior residuals >= 2m: {len(big)}/{len(gapstats)}")
        print("  examples (V, B, m, x, pos, gaps, sigma, residuals):")
        for g in gapstats[:12]:
            print(f"    V={g['V']} B={g['B']} m={g['m']} x={g['x']} "
                  f"pos={g['pos']} gaps=({g['g1']},{g['g2']},{g['g3']}) "
                  f"sigma={g['sigma']} cs={sorted(g['cs'])}")

    print(f"\n[E7] two-sided points where the canonical TOP-PAIR rescue FAILS: "
          f"{len(fails)}")
    shown = 0
    for r in fails:
        if shown >= MAX_DUMPS:
            break
        V, B, m, x = r["V"], r["B"], r["m"], r["x"]
        u = r["L"][0]
        w = r["R"][0]
        a, c = V[u], V[w]
        i, j = r["top_pair"]
        N = V[i] + V[j]
        e = V[i]
        R3 = (m * N) // B + 1
        m3, arg3 = level3_argmax(a, c, e, N)
        gcds = (gcd(a, N), gcd(c, N), gcd(e, N))
        g4 = gcd(gcd(a, c), gcd(e, N))
        g4 = gcd(g4, N)
        # level-3 two-sidedness at its own argmax
        l3_two = False
        for (k, vals) in arg3[:1]:
            l3_two = (any(vals[t] == m3 for t in range(3)) and
                      any(vals[t] == N - m3 for t in range(3)))
        ratio_ok = (m3 * B > m * N)   # level-3 beats the level-4 witness ratio?
        print(f"   FAIL V={V} B={B} m={m} x={x} delta={r['delta']} "
              f"rs={r['rs']} cs={r['cs']}")
        print(f"        binders a={a}(i={u}) c={c}(i={w}); top pair ({i},{j}) "
              f"speeds ({V[i]},{V[j]}) N={N} e={e}")
        print(f"        level-3: R={R3} m3*={m3} (beats witness ratio: {ratio_ok}) "
              f"gcds(a,c,e;N)={gcds} gcd_all={g4} l3_two_sided={l3_two}")
        shown += 1
    if len(fails) > MAX_DUMPS:
        print(f"   ... {len(fails) - MAX_DUMPS} more failures not shown")

    print("\n[SUMMARY]")
    tot = len(records)
    tot_two = two_hist.get(True, 0)
    tot_top = sum(1 for r in records if r["top_ok"])
    tot_top_two = sum(1 for r in records if r["two"] and r["top_ok"])
    tot_any = sum(1 for r in records if r["any_ok"])
    print(f"  all failing points          : top-pair rescue {tot_top}/{tot}; "
          f"any-pair {tot_any}/{tot}")
    print(f"  two-sided failing points    : top-pair rescue {tot_top_two}/{tot_two}")
    print(f"  shallow (delta<=4) points   : "
          f"{sum(1 for r in records if r['delta'] <= 4)}, "
          f"of which top-pair rescue ok: "
          f"{sum(1 for r in records if r['delta'] <= 4 and r['top_ok'])}")

if __name__ == "__main__":
    main()
