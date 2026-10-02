"""Numerical verification of the deepseek/GLM Theorems 1 and 2 (item (c), R1).

Theorem 1:  min(||t||, ||2t||, ..., ||(n-1)t||, ||2nt||) <= 2/(2n+1)  all t.
Theorem 2:  equality forces t = +-2/(2n+1) mod 1; binding runners {1, 2n};
            every other runner sits at distance >= 3/(2n+1).

Machine check (exact, no floating point), n = 2..40:
  T1a  gap(family(n)) == 2/(2n+1)
  T2a  argmax set == {2/(2n+1), (2n-1)/(2n+1)}   (R1: uniqueness)
  T2b  at t = 2/(2n+1): binders == {1, 2n}
  T2c  at t = 2/(2n+1): every k in 2..n-1 has ||k t|| >= 3/(2n+1)
  T2d  the pair-distances of P = {0, t, 2t, ..., nt} at the argmax:
       gap between 0 and nt is delta/2, the other n gaps are exactly delta
       (the equality configuration of the proof).
"""
import sys
from fractions import Fraction

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap, dist_at, family

def main():
    ok = True
    for n in range(2, 41):
        V = family(n)
        d = 2 * n + 1
        delta = Fraction(2, d)
        g, ts = gap(V)
        line = f"n={n:2d} gap={g}"
        # T1a
        if g != delta:
            ok = False; line += f"  T1a FAIL (expected {delta})"; print(line); continue
        # T2a: unique argmax pair
        expect_ts = sorted([Fraction(2, d), Fraction(d - 2, d)])
        if ts != expect_ts:
            ok = False; line += f"  T2a FAIL argmaxes={ts}"
        # T2b/T2c at t = 2/(2n+1)
        t = Fraction(2, d)
        bind = sorted(u for u in V if dist_at(u, t) == delta)
        if bind != [1, 2 * n]:
            ok = False; line += f"  T2b FAIL binders={bind}"
        slack = min((dist_at(k, t) for k in range(2, n)), default=Fraction(1))
        if slack < Fraction(3, d):
            ok = False; line += f"  T2c FAIL slack={slack}"
        # T2d: equality configuration of P
        pts = [Fraction(k) * t for k in range(n + 1)]
        circ = sorted(set(pts))
        gaps_p = [b - a for a, b in zip(circ, circ[1:])] + [circ[0] + 1 - circ[-1]]
        gaps_p.sort()
        expect_gaps = sorted([delta / 2] + [delta] * n)
        if gaps_p != expect_gaps:
            ok = False; line += f"  T2d FAIL gaps={gaps_p}"
        if n <= 12 or n % 10 == 0 or not ok:
            print(line + ("  ok" if ok else ""))
    print("ALL CHECKS PASSED" if ok else "FAILURES PRESENT")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
