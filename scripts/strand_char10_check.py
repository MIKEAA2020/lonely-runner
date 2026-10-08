#!/usr/bin/env python3
"""strand_char10_check.py -- verify the closure of the EXTENDED character
system: the 10 functionals {4 moving-start singletons} U {6 pairwise
differences (mod p)} on the start space.

Claim (Lemma V'''): under the renormalization action, the multiset of the
10 marginal sizes -- aggregated over all size tuples (sigma-free form) --
is an orbit invariant.  Mechanism (verified construction, delta=5, p=19):
image slots = (pair(a3,a4), singleton(a4), pair(a5,a4), pair(a2,a4));
image pairwise = (singleton(a3), pair(a3,a5), pair(a2,a3), singleton(a5),
singleton(a2), pair(a2,a5)) -- exactly the source's 10 functionals.
Checked on the full top orbits at p=17 and p=19.
"""
import sys
from collections import Counter

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples
from scaling_closure_verify import solutions_m5

T, M = 8, 5


def canon(x, p):
    x %= p
    return min(x, p - x)


def renorm_orbit(S, p):
    D = sorted(set([1] + list(S)))
    out = []
    for delta in D:
        inv = pow(delta, -1, p)
        out.append(tuple(sorted(canon(d * inv, p) for d in D if d != delta)))
    return sorted(set(out))


def enum(p, fam, szts):
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    return {key[:5]: sols[key] for key in sorted(sols)}


def char10_profile(lst, p):
    vals = [len({a[i] % p for a in lst}) for i in range(4)]
    vals += [len({(a[j] - a[i]) % p for a in lst})
             for i in range(4) for j in range(i + 1, 4)]
    return tuple(sorted(vals))


def main():
    for p, S in [(17, (2, 3, 4, 8)), (19, (2, 4, 5, 8))]:
        k, (on, off), cls, eps = profile_sizes(p, T)
        szts = size_tuples(M, on, off, p)
        orb = renorm_orbit(S, p)
        print('=== p=%d, orbit of %s: %s ===' % (p, S, orb))
        counters = []
        for s in orb:
            rec = enum(p, s, szts)
            c = Counter(char10_profile(lst, p)
                        for lst in rec.values() if lst)
            counters.append(c)
            print('  member %s: %d nonzero keys, total %d'
                  % (s, sum(c.values()), sum(len(v) for v in
                                             rec.values())))
        ok = all(c == counters[0] for c in counters[1:])
        print('  10-functional aggregate multiset invariant across orbit: '
              '%s' % ok)
        if not ok:
            for i, c in enumerate(counters[1:], 1):
                d1, d2 = counters[0] - c, c - counters[0]
                if d1 or d2:
                    print('    rep-only: %s' % dict(list(d1.items())[:6]))
                    print('    mem%d-only: %s' % (i, dict(list(d2.items())[:6])))
                    break


if __name__ == '__main__':
    main()
