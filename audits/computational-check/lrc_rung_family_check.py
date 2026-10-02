"""Rung-family check: does Theorem 1 generalize to k >= 2?

Claim (Theorem 1'): for n >= 2, k >= 2, gap({1,...,n-1,kn}) = k/(kn+1),
with binding runners {1, kn} at t = k/(kn+1).  The pigeonhole proof of
Theorem 1 goes through verbatim with 2n -> kn and delta = k/(kn+1):
  - classical step gives beta = ||nt|| <= 1/(n+1) < k/(kn+1) (k >= 2);
  - pair distances > delta except {0, nt};
  - 1 > beta + n*delta  =>  beta < 1/(kn+1)  =>  ||kn t|| <= k beta < delta.
The equality case (unique argmax) needs gcd(n(n+1)/2, kn+1) = 1 for the
sum argument, which can fail (e.g. n=k=3: gcd(6,10)=2) -- so we check
argmax uniqueness numerically per (n,k).
"""
import sys
from fractions import Fraction

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap, dist_at

def family_k(n, k):
    return list(range(1, n)) + [k * n]

def main():
    print("n   k   gap           k/(kn+1)      argmax unique at +-k/(kn+1)?")
    uniq_fail = []
    for n in range(2, 13):
        for k in range(2, 7):
            V = family_k(n, k)
            g, ts = gap(V)
            target = Fraction(k, k * n + 1)
            ok_val = (g == target)
            expect = sorted([target, Fraction(k * n + 1 - k, k * n + 1)])
            ok_uniq = (ts == expect)
            if not ok_val:
                print(f"VALUE FAIL n={n} k={k}: {g} != {target}")
            if not ok_uniq:
                uniq_fail.append((n, k, ts))
            if n <= 6 or k == 2 or not ok_uniq or not ok_val:
                print(f"{n:2d}  {k:2d}  {str(g):12s}  {str(target):12s} "
                      f" {'YES' if ok_uniq else 'NO: ' + str(ts)}")
    print(f"\nargmax-uniqueness failures (k >= 3): {uniq_fail if uniq_fail else 'none in range'}")
    # binding runners at the canonical time
    n, k = 7, 2
    t = Fraction(k, k * n + 1)
    V = family_k(n, k)
    bind = [u for u in V if dist_at(u, t) == Fraction(k, k * n + 1)]
    print(f"n=7 k=2 binders at t=k/(kn+1): {bind} (expect [1, {k*n}])")
    return 0

if __name__ == "__main__":
    sys.exit(main())
