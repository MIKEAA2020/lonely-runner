#!/usr/bin/env python3
"""
gt_pool_brute.py -- independent ground-truth check of pool_sat decisions.

For sampled distinct difference tuples, decide admissibility by PURE
enumeration (no cascade pruning, no arc-fit formula): for every ordered
size tuple with sum >= p and every start triple (a2,a3,a4), build the
five arc masks and check all p positions a5 of the fifth arc directly.

Convention (census_m5-compatible): a tuple is ADMISSIBLE iff some
placement has arcs 1..4 leaving remainder R with 1 <= |R| <= s5 and arc
5 CONTAINING R (some a5 with R & ~arc5 == 0).  R empty = bystander
(reported separately, not admissible).

Reads a pool_sat dump (tsv: "d2,d3,d4,d5<TAB>decision<TAB>steps"),
samples admissible + inadmissible tuples, verifies each, compares.
"""
import random
import sys

def rot(base, a, p, full):
    if a == 0:
        return base
    return ((base << a) | (base >> (p - a))) & full

def brute_admissible(p, T, m, tup):
    """tup = (d2..dm). Returns 1 admissible / 0 not / 2 bystander-only."""
    k = (p - p % T) // T
    eps = p % T
    d0 = (eps * eps - 1) // T
    M = k * eps + d0
    if 2 * M >= p - 1:
        on, off = 2 * k + 1, 2 * k + 2
    else:
        on, off = 2 * k, 2 * k + 1
    full = (1 << p) - 1
    D = (p - 1) // 2
    dd = (1,) + tuple(tup)          # d[1] = 1 pinned
    bystander = False
    for code in range(1 << m):
        s = tuple(off if (code >> i) & 1 else on for i in range(m))
        if sum(s) < p:
            continue
        bases = []
        for i in range(m):
            b = 0
            x = 0
            for j in range(s[i]):
                b |= 1 << x
                x += dd[i]
                if x >= p:
                    x -= p
            bases.append(b)
        m1 = bases[0]                 # a1 = 0
        b5 = bases[m - 1]
        # precompute rotations of arc5
        rots5 = [rot(b5, a, p, full) for a in range(p)]
        for a2 in range(p):
            m2 = rot(bases[1], a2, p, full)
            u12 = m1 | m2
            for a3 in range(p):
                u123 = u12 | rot(bases[2], a3, p, full)
                for a4 in range(p):
                    u = u123 | rot(bases[3], a4, p, full)
                    R = full & ~u
                    if R == 0:
                        bystander = True
                        continue
                    nr = R.bit_count()
                    if nr > s[m - 1]:
                        continue
                    for a5 in range(p):
                        if R & ~rots5[a5] == 0:
                            return 1
    return 2 if bystander else 0

def main():
    p = int(sys.argv[1]); T = int(sys.argv[2]); m = int(sys.argv[3])
    dump = sys.argv[4]
    n_sample = int(sys.argv[5]) if len(sys.argv) > 5 else 5
    dec = {}
    with open(dump) as f:
        for line in f:
            t, d, st = line.strip().split('\t')
            dec[tuple(int(x) for x in t.split(','))] = int(d)
    adm = [t for t, d in dec.items() if d == 1]
    inadm = [t for t, d in dec.items() if d == 0]
    rng = random.Random(20261007)
    sample = (rng.sample(adm, min(n_sample, len(adm))) +
              rng.sample(inadm, min(n_sample, len(inadm))))
    ok = True
    for tup in sample:
        got = brute_admissible(p, T, m, tup)
        want = dec[tup]
        tag = 'OK ' if got == want else 'MISMATCH'
        if got != want:
            ok = False
        print('  %s tuple %s: brute=%d pool_sat=%d' % (tag, tup, got, want))
        sys.stdout.flush()
    print('GT %s: %d/%d agree' % ('PASS' if ok else 'FAIL',
                                  len(sample) - (0 if ok else 1), len(sample)))

if __name__ == '__main__':
    main()
