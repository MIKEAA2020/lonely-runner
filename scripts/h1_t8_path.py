#!/usr/bin/env python3
"""
h1_t8_path.py -- the remaining top-class candidates at p=29 + the chain
benchmarks at p=19/p=29 (the H1 volume-law calibration).

The doubling path mod 29 (canonicalized): 2->4->8->13->3->6->12->5->10->
9->11->7->14->(1).  Runs one representative per consecutive quadruple
(the orbit law makes the 24 orders redundant for totals) + the chain
family (1,1,1,1) at p=19 and p=29 for the benchmark ratio.
"""
import sys
import time
from math import comb

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples
from h1_t8_prime import count_sols, inv_table, brute_sols

T, M = 8, 5
PATH = [2, 4, 8, 13, 3, 6, 12, 5, 10, 9, 11, 7, 14]


def fam_totals(p, fams):
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(M, on, off, p)
    inv = inv_table(p)
    res = {}
    for fam in fams:
        per = {}
        for szt in szts:
            per[szt] = count_sols(p, inv, tuple(szt) + tuple(fam))
        res[tuple(fam)] = per
    return res, szts


def main():
    # ---- p=29: the doubling-path quadruples ----
    quads = [tuple(PATH[i:i + 4]) for i in range(len(PATH) - 3)]
    # plus 3-of-a-path + neighbor candidates
    quads += [(2, 4, 8, 13), (4, 8, 13, 3), (2, 4, 13, 3), (2, 8, 13, 3)]
    print('=== p=29 doubling-path quadruples (one order per set) ===')
    per29, szts29 = fam_totals(29, quads)
    rows = []
    for fam, per in per29.items():
        tot = sum(per.values())
        bs = max(per.items(), key=lambda kv: kv[1])
        rows.append((tot, fam, bs))
    for tot, fam, bs in sorted(rows, reverse=True):
        print('  %s: total %6d, max/size %5d at %s (ov=%d)'
              % (fam, tot, bs[1], bs[0], sum(bs[0]) - 29))
    # chain benchmark at p=29
    print('  chain (1,1,1,1):')
    per, _ = fam_totals(29, [(1, 1, 1, 1)])
    for fam, pp in per.items():
        tot = sum(pp.values())
        bs = max(pp.items(), key=lambda kv: kv[1])
        print('    total %7d, max/size %6d at %s (ov=%d)'
              % (tot, bs[1], bs[0], sum(bs[0]) - 29))
        # GT on the top key
        key = tuple(bs[0]) + fam
        nb = brute_sols(29, key)
        print('    GT %s: brute %d vs cascade %d -> %s'
              % (key, nb, bs[1], 'MATCH' if nb == bs[1] else 'MISMATCH'))
    # GT one path quad top key
    tot, fam, bs = sorted(rows, reverse=True)[0]
    key = tuple(bs[0]) + fam
    nb = brute_sols(29, key)
    print('  GT %s: brute %d vs cascade %d -> %s'
          % (key, nb, bs[1], 'MATCH' if nb == bs[1] else 'MISMATCH'))

    # ---- p=19: chain benchmark ----
    print('=== p=19 chain (1,1,1,1) benchmark ===')
    per19, _ = fam_totals(19, [(1, 1, 1, 1)])
    for fam, pp in per19.items():
        tot = sum(pp.values())
        bs = max(pp.items(), key=lambda kv: kv[1])
        print('  total %7d, max/size %6d at %s (ov=%d)'
              % (tot, bs[1], bs[0], sum(bs[0]) - 19))
        key = tuple(bs[0]) + fam
        nb = brute_sols(19, key)
        print('  GT %s: brute %d vs cascade %d -> %s'
              % (key, nb, bs[1], 'MATCH' if nb == bs[1] else 'MISMATCH'))

    # ---- the summary fit ----
    print()
    print('=== the volume-law table (top distinct vs chain, per top size '
          'tuple) ===')
    print('p=17 k=2: distinct 876 @ ov=8  | chain 8880 @ ov=8   '
          'ratio %.4f' % (876 / 8880.0))
    print('p=19 k=2: distinct 108 @ ov=6  | chain (above)      ')
    print('p=29 k=3: distinct 186 @ ov=11 | chain (above)      ')


if __name__ == '__main__':
    main()
