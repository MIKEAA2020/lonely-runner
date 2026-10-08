#!/usr/bin/env python3
"""strand_vprime2_check.py -- verify Lemma V'' (the corrected orbit
invariant): the multiset of PAIRWISE-DIFFERENCE marginal sizes is an orbit
invariant, matched by the size permutation sigma.  p=19, orbit pair
S=(2,4,5,8) -> S2=(3,4,6,8) at delta=5 (the machine-verified bijection).
"""
import sys
from collections import Counter

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples
from scaling_closure_verify import solutions_m5

T, M = 8, 5
p = 19


def enum(p, fam, szts):
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    return {key[:5]: sols[key] for key in sorted(sols)}


def pairdiff_profile(lst):
    return tuple(sorted(len({(a[j] - a[i]) % p for a in lst})
                        for i in range(4) for j in range(i + 1, 4)))


def sigma_map(szt):
    sN, s2, s3, s4, s5 = szt
    return (s4, s3, sN, s5, s2)


def main():
    k, (on, off), cls, eps = profile_sizes(p, T)
    szts = size_tuples(M, on, off, p)
    recS = enum(p, (2, 4, 5, 8), szts)
    recS2 = enum(p, (3, 4, 6, 8), szts)
    mism = 0
    checked = 0
    for szt in szts:
        lst = recS.get(szt, [])
        if not lst:
            continue
        lst2 = recS2.get(sigma_map(szt), [])
        p1 = pairdiff_profile(lst)
        p2 = pairdiff_profile(lst2)
        checked += 1
        if p1 != p2 or len(lst) != len(lst2):
            mism += 1
            print('MISMATCH at szt %s -> %s: n %d vs %d, profiles %s vs %s'
                  % (szt, sigma_map(szt), len(lst), len(lst2), p1, p2))
    print('Lemma V\'\' pairwise-difference profile invariance: %d/%d '
          'sigma-matched size tuples match, %d mismatches'
          % (checked - mism, checked, mism))


if __name__ == '__main__':
    main()
