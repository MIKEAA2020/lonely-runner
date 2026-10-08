#!/usr/bin/env python3
"""strand_vprime_settle.py -- settle the Lemma V' question at mixed size
tuples: does the proof's bijection Sol(S,szt) -> Sol(S',szt.sigma) hold?
p=19, S=(2,4,5,8), delta=5, S'=(3,4,6,8).  For EVERY size tuple, construct
images of all solutions and test membership; compare per-key counts."""
import sys
from collections import Counter

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples
from scaling_closure_verify import solutions_m5

T, M = 8, 5
p = 19
S = (2, 4, 5, 8)
delta = 5
inv = pow(delta, -1, p)
S2 = (3, 4, 6, 8)


def enum(p, fam, szts):
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    return {key[:5]: sols[key] for key in sorted(sols)}


def sigma_map(szt):
    """Arc correspondence for delta=5, inv=4, S=(2,4,5,8) -> S2=(3,4,6,8):
    S slots (d2,d3,d4,d5)=(2,4,5,8), normalizer N (d=1), sizes
    (s_N, s_2, s_3, s_4, s_5).  Images: N -> d=4 (slot 1 of S2);
    d2=2 -> d=8 (slot 3 of S2); d3=4 -> raw lam 16 -> canon 3 (slot 0);
    d4=5=delta -> normalizer (slot N of S2); d5=8 -> raw lam 13 -> canon 6
    (slot 2).  So new sizes: (s_4 | s_3, s_N, s_5, s_2)."""
    sN, s2, s3, s4, s5 = szt
    return (s4, s3, sN, s5, s2)


def construct(a, szt):
    """Image of solution a under the proof recipe (sizes carried)."""
    sN, s2, s3, s4, s5 = szt
    a2, a3, a4, a5 = a
    lam = {2: (2 * inv) % p, 4: (4 * inv) % p, 8: (8 * inv) % p}
    # sizes of the source arcs, per arc identity
    sz = {'N': sN, 2: s2, 4: s3, 5: s4, 8: s5}
    img = {}
    # normalizer (d=1, start 0, size sN) -> d=inv*1=4, start 0, size sN
    # re-anchor later by -inv*a4
    img[4] = (0, sz['N'])
    # d2=2 -> lam=8 <= 9: no flip; start inv*a2, size s2
    img[8] = (inv * a2, sz[2])
    # d3=4 -> lam=16 > 9: flip to d=3; start inv*a3 + 16*(s3-1), size s3
    img[3] = (inv * a3 + lam[4] * (sz[4] - 1), sz[4])
    # d4=5=delta -> d=1, start inv*a4, size s4  (new normalizer)
    img[1] = (inv * a4, sz[5])
    # d5=8 -> lam=13 > 9: flip to d=6; start inv*a5 + 13*(s5-1), size s5
    img[6] = (inv * a5 + lam[8] * (sz[8] - 1), sz[8])
    # re-anchor: subtract the new normalizer's start
    anchor = img[1][0]
    # build the new (normalizer size, moving starts per S2 slot (3,4,6,8))
    newszt = (img[1][1], img[3][1], img[4][1], img[6][1], img[8][1])
    b = ((img[3][0] - anchor) % p, (img[4][0] - anchor) % p,
         (img[6][0] - anchor) % p, (img[8][0] - anchor) % p)
    return b, newszt


def covers(p, szt, fam, a):
    pts = set(range(szt[0]))
    for i, d in enumerate(fam):
        pts |= {(a[i] + t * d) % p for t in range(szt[i + 1])}
    return len(pts) == p


def main():
    k, (on, off), cls, eps = profile_sizes(p, T)
    szts = size_tuples(M, on, off, p)
    recS = enum(p, S, szts)
    recS2 = enum(p, S2, szts)
    tot_ok = tot_bad = 0
    print('per-key bijection test (S=%s -> S2=%s, delta=%d):' %
          (S, S2, delta))
    for szt in szts:
        lst = recS.get(szt, [])
        if not lst:
            continue
        newszt = sigma_map(szt)
        lst2 = recS2.get(newszt, [])
        set2 = set(lst2)
        ok = bad = 0
        badex = None
        for a in lst:
            b, nsz = construct(a, szt)
            if nsz != newszt:
                print('  SIZE MISMATCH at szt %s: constructed sizes %s, '
                      'expected %s' % (szt, nsz, newszt))
                bad += 1
                continue
            if b in set2:
                ok += 1
            else:
                bad += 1
                if badex is None:
                    badex = (a, b)
        tot_ok += ok
        tot_bad += bad
        if bad:
            print('  szt %s (n=%d) -> %s (n2=%d): ok=%d bad=%d  example %s'
                  % (szt, len(lst), newszt, len(lst2), ok, bad, badex))
    print('TOTAL: images in target %d, misses %d' % (tot_ok, tot_bad))
    # and the reverse direction count check per key
    mism = 0
    for szt in szts:
        n1 = len(recS.get(szt, []))
        n2 = len(recS2.get(sigma_map(szt), []))
        if n1 != n2:
            mism += 1
            if mism <= 8:
                print('  count mismatch: Sol(S,%s)=%d vs Sol(S2,%s)=%d'
                      % (szt, n1, sigma_map(szt), n2))
    print('per-key count mismatches under sigma: %d' % mism)
    # recompute the profile Counter comparison very explicitly
    def prof(szt, lst):
        return tuple(sorted(len({a[i] for a in lst}) for i in range(4)))
    c1 = Counter(prof(szt, lst) for szt, lst in recS.items() if lst)
    c2 = Counter(prof(szt, lst) for szt, lst in recS2.items() if lst)
    print('profile counters equal: %s' % (c1 == c2))
    print('  c1 - c2: %s' % dict(c1 - c2))
    print('  c2 - c1: %s' % dict(c2 - c1))


if __name__ == '__main__':
    main()
