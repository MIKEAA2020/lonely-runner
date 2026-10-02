"""Rung-family check, extended: value, EXACT argmax set, binders, slack.

Claims (rung-family note, Theorems 1-3):
  (i)  gap({1,...,n-1, kn}) = k/(kn+1) for all n >= 2, k >= 1;
  (ii) k >= 2: the argmax set is EXACTLY {k/(kn+1), (kn+1-k)/(kn+1)}
       (unique up to reflection; NO gcd condition -- the no-wrap
       induction sigma(j+1) = sigma(j) + m replaces the sum argument);
  (iii) k = 1: the argmax set is exactly {a/(n+1) : gcd(a, n+1) = 1}
       (phi(n+1) maximizers, the units grid);
  (iv) binders at t = k/(kn+1) are exactly {1, kn};
  (v)  slack: every other speed sits at distance >= (k+1)/(kn+1) there.

Range: n = 2..16, k = 1..12 (exact rational arithmetic throughout).
"""
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap_int


def family_k(n, k):
    return list(range(1, n)) + [k * n]


def main():
    print("n   k   gap        = k/(kn+1)?   argmax set exact?   "
          "binders    slack>= (k+1)/(kn+1)?")
    nfail = 0
    for n in range(2, 17):
        for k in range(1, 13):
            V = family_k(n, k)
            g, ts = gap_int(V)
            D = k * n + 1
            target = Fraction(k, D)
            ok_val = (g == target)
            if k >= 2:
                expect = sorted([target, Fraction(D - k, D)])
                ok_amx = (ts == expect)
            else:
                units = sorted(Fraction(a, D) for a in range(1, D)
                               if gcd(a, D) == 1)
                ok_amx = (ts == units)
            t = target
            vals = {u: Fraction(min((u * t.numerator) % t.denominator,
                                    t.denominator -
                                    (u * t.numerator) % t.denominator),
                                t.denominator) for u in V}
            bind = sorted(u for u in V if vals[u] == target)
            ok_bind = (bind == sorted([1, k * n]))
            others = [vals[u] for u in V if vals[u] != target]
            slack = min(others) if others else Fraction(1)
            ok_slack = (slack >= Fraction(k + 1, D))
            ok = ok_val and ok_amx and ok_bind and ok_slack
            if not ok:
                nfail += 1
                print(f"FAIL n={n} k={k}: gap={g} target={target} "
                      f"argmax={ts} binders={bind} slack={slack}")
            elif n <= 6 or k == 2 or (n, k) in [(7, 3), (8, 3), (9, 4),
                                                (12, 5), (16, 12), (16, 2)]:
                print(f"{n:2d}  {k:2d}  {str(g):9s}  "
                      f"{'YES' if ok_val else 'NO'}          "
                      f"{'YES' if ok_amx else 'NO: ' + str(ts[:4])}          "
                      f"{bind}     "
                      f"{'YES' if ok_slack else 'NO'}")
    print(f"\nfailures: {nfail} "
          f"(range n=2..16, k=1..12 = {15*12} parameter pairs)")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
