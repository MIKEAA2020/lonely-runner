"""Case book: the complete residual case list for the deficit law at
n=4 (m=1..8) and n=5 (m=1..3) — every bridge config that SURVIVES the
orbit lemma (L2) must be killed by an off-grid pair-sum crossing (L3).
For each survivor we record: the residue config R, its circular gap
structure (NET FORM data), the killing pair/time/value, and tags.

This is the finite case load that a k=3 residue proof must discharge
beyond the orbit lemma.

Usage: python3 lrc_diag_casebook.py out.json
"""
import sys
import json
from collections import Counter
from fractions import Fraction
from math import gcd, comb
import itertools

import numpy as np

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_diag_gap import GapEngine, reduce_candidate
from lrc_diag_bridge import orbit_fail

CAP = 80


def circular_gaps(R, b, m):
    """Gaps between consecutive residues on the circle, with the
    pole-to-pole exterior gap 2m included (through 0)."""
    pts = sorted(R)
    gaps = [pts[i + 1] - pts[i] for i in range(len(pts) - 1)]
    gaps.append(b - pts[-1] + pts[0])
    return gaps


def casebook(eng):
    out = {}
    for n, m_list in [(4, range(1, 9)), (5, [1, 2, 3])]:
        cases = []
        orbit_pass = 0
        for m in m_list:
            lo, hi = m * (n + 1) + 1, min(n * (2 * m + 1) + 1, CAP)
            for b in range(lo, hi + 1):
                pool = list(range(m + 1, b - m))
                need = n - 2
                if need > len(pool):
                    continue
                units = [s for s in range(1, b) if gcd(s, b) == 1]
                target = Fraction(m, b)
                for interiors in itertools.combinations(pool, need):
                    R = [m] + list(interiors) + [b - m]
                    if len(orbit_fail(R, b, m)) == 0:
                        orbit_pass += 1
                        # canonical lift s=1 and the full unit family
                        realized = False
                        best_win = None
                        for s in units:
                            c = sorted((s * r) % b for r in R)
                            g = eng.gap(c)
                            if g == target:
                                realized = True
                                break
                        if realized:
                            cases.append(dict(m=m, b=b, R=R, realized=True))
                            continue
                        c1 = sorted(r % b for r in R)
                        g1, (mp, D, jp) = eng.gap(c1, with_argmax=True)
                        a2, b2, m2 = reduce_candidate(mp, D, jp)
                        # winning pair
                        pair = None
                        for u, w in itertools.combinations(c1, 2):
                            if u + w == D:
                                pair = (u, w)
                                break
                        gaps = circular_gaps(R, b, m)
                        interior_sum = b - 2 * m
                        avg_int = interior_sum / (n - 1)
                        cases.append(dict(
                            m=m, b=b, R=R, realized=False,
                            lift=c1, gap_lift=str(g1),
                            win_D=D, win_pair=pair, win_t=f"{jp}/{D}",
                            win_val=str(g1), win_b2=b2, win_m2=m2,
                            win_d=m2 * (n + 1) - b2,
                            gaps=gaps, max_gap=max(gaps),
                            avg_interior=round(avg_int, 3),
                            deficit_mag=b - m * (n + 1),
                        ))
        out[n] = dict(orbit_survivors=orbit_pass, cases=cases)
    return out


if __name__ == "__main__":
    eng = GapEngine(CAP)
    res = casebook(eng)
    with open(sys.argv[1], "w") as f:
        json.dump(res, f, indent=1)
    for n, r in res.items():
        cs = r["cases"]
        real = [c for c in cs if c["realized"]]
        mg = Counter(c.get("max_gap") for c in cs if not c["realized"])
        dm = Counter(c.get("deficit_mag") for c in cs if not c["realized"])
        print(f"n={n}: orbit survivors={r['orbit_survivors']}, "
              f"realized={len(real)}, killed_offgrid={len(cs)-len(real)}")
        print(f"   survivor max-gap hist: {sorted(mg.items())[:10]}")
        print(f"   survivor deficit-magnitude hist: {sorted(dm.items())[:12]}")
