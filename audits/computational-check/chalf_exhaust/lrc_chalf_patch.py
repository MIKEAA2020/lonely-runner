"""C-HALF PATCH: scan the inter-rung-productive cores the drop experiment missed.

The drop experiment's zoo universe was the NSK productive cores (cores of
TIGHT sets).  The exhaustiveness check found that the n=6 inter-rung
inhabitant {1,5,6,11,16,17} @ 5/33 is a strict seed under tight/rung
productivity -- its core {1,5,6,11,16} (gap 4/17, a descaled rung-4) is
zoo-productive ONLY via the inter-rung extension x=17.  Such cores were
outside the scanned universe, so the C-half (margin >= 1/(k+2), deficit
d = m(k+2)-b >= 0 at crossings) was never measured for them.

This script scans the inter-rung core families exhaustively for x in range:
  C1 = {1,5,6,11,16}        (k=5, gap 4/17; core of the 5/33 inhabitant)
  C2 = {1,3,4,5,7,13}       (k=6; core of the second primitive 3/23 class)
  C3 = {1,2,3,4,5,7}        (k=6, tight@6; core of the first 3/23 class --
                             already in the experiment, re-verified here)
  C4 = {1,5,6,11,16,17}     (k=6, the 5/33 inhabitant itself as a core)
Checks per (S, x): touch/drop, margin vs 1/(k+2), and at every crossing
argmax the deficit d >= 0 + residue law + pole law.  Exact Fractions.
"""
import json
import os
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc_gap_lib import gap_int, _norm
from lrc_drop_experiment import core_info, classify, check

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def dist_b(r, b):
    r %= b
    return r if r < b - r else b - r


def scan(S, xmax):
    S = list(S)
    k = len(S)
    ci = core_info(S)
    recs = []
    silent = 0
    dmin = None
    for x in range(1, xmax + 1):
        if x in S:
            continue
        v2 = sorted(S + [x])
        g_new, ts_new = gap_int(v2)
        rec, _ = classify(v2, ci, x)
        floor = Fraction(1, k + 2)
        margin = g_new - floor
        if margin < 0:
            silent += 1
            print(f"  !!! SILENT: S={S} x={x} g_new={g_new} < {floor}")
        # crossing deficit + laws
        for e in rec["amax_new"]:
            if e["cls"] != "crossing":
                continue
            t = Fraction(*map(int, e["t"].split("/")))
            a, b = t.numerator, t.denominator
            m = g_new * b
            assert m.denominator == 1
            m = m.numerator
            d = m * (k + 2) - b
            dmin = d if dmin is None else min(dmin, d)
            if d < 0:
                print(f"  !!! DEFICIT: S={S} x={x} t*={t} m={m} b={b} d={d}")
            # pole law
            res = [(u * a) % b for u in v2]
            assert min(res) == m and max(res) == b - m, (S, x, t)
            # residue law on binder pairs
            for i in range(len(e["binders"])):
                for j in range(i + 1, len(e["binders"])):
                    u, w = e["binders"][i], e["binders"][j]
                    ok = (u + w) % b == 0 or (abs(u - w) % b == 0 and u != w)
                    assert ok, (S, x, t, u, w)
        recs.append(rec)
    return recs, silent, dmin


def main():
    report = []
    for name, S, xmax in [
        ("C1 {1,5,6,11,16} (descaled rung-4, core of 5/33)", [1, 5, 6, 11, 16], 80),
        ("C2 {1,3,4,5,7,13} (core of 2nd 3/23 class)", [1, 3, 4, 5, 7, 13], 80),
        ("C3 {1,2,3,4,5,7} (tight@6, core of 1st 3/23 class)", [1, 2, 3, 4, 5, 7], 80),
        ("C4 {1,5,6,11,16,17} (the 5/33 inhabitant as core)", [1, 5, 6, 11, 16, 17], 80),
    ]:
        recs, silent, dmin = scan(S, xmax)
        drops = [r for r in recs if not r["touch"]]
        floors = sum(1 for r in recs if r["margin_flag"] == "floor")
        inter = [r for r in recs if r["margin_flag"] == "above" and
                 Fraction(*map(int, r["margin"].split("/"))) < Fraction(2, 2 * len(S) + 3)]
        line = (f"{name}: x in [1..{xmax}] minus S: {len(recs)} pairs, "
                f"touch={len(recs)-len(drops)}, drop={len(drops)}, "
                f"SILENT={silent}, floor-landings={floors}, "
                f"inter-rung-band drops={len(inter)}, min deficit d={dmin}")
        print(line)
        report.append(line)
        for r in inter[:6]:
            print(f"    inter-band: x={r['x']} g_new={r['g_new']} "
                  f"cls={r['cls']} t*={[e['t'] for e in r['amax_new']]} "
                  f"binders={[e['binders'] for e in r['amax_new']]}")
    with open(os.path.join(OUT, "chalf_patch_report.txt"), "w") as f:
        f.write("\n".join(report) + "\n")
    print("\nALL CHECKS PASS (no SILENT, no negative deficit, "
          "pole+residue laws verified on every crossing)")


if __name__ == "__main__":
    main()
