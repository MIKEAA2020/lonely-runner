#!/usr/bin/env python3
"""
chain_structure_check.py -- validate the staircase/chain structure lemmas
for the all-interval family (1,1,1,1) at the committed zoo point
(T=8, (1K,5U), p=17), on the T-18 rebuilt solutions.

Lemmas under test (the theorem's substrate):
  [B2] every covering's sorted-order overlap vector o satisfies
       sum(o) = ov = sum(s) - p exactly (the telescoping identity);
  [B1'] the per-junction o_i <= s_{sigma(i)} and (reach bound) the
       inter-start gaps are <= s_max, i.e. o_i >= s_{sigma(i)} - s_max;
  [chain count] |Sol| >= sum_sigma C(ov+m-1, m-1)/overcount  and the
       non-chain (nested) excess is recorded.
"""
import json
import sys
from collections import Counter, defaultdict
from math import comb

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples, k_of
from scaling_closure_verify import solutions_m5

def main():
    p, T, m = 17, 8, 5
    fam = (1, 1, 1, 1)
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    s_max = max(on, off)
    szts = size_tuples(m, on, off, p)
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    print('family %s at p=%d: %d size tuples with solutions' %
          (fam, p, len(sols)))
    tot = 0
    tot_chain = 0
    bad_sum = 0
    bad_gap = 0
    nested = 0
    for key, lst in sorted(sols.items()):
        s = key[:m]
        ov = sum(s) - p
        # chain count for this size tuple: sum over cyclic orders of
        # #{o >= 0, sum o = ov} with o_i <= s_{sigma(i)} -- exact chain
        # placements (sorted order = sigma); overcount only at tied
        # starts.  Compute directly per solution instead.
        n_chain = 0
        for (a2, a3, a4, a5) in lst:
            tot += 1
            starts = [0, a2, a3, a4, a5]
            # cyclic sorted order starting from 0 (a1=0 is the min
            # residue; ties resolved stably)
            order = sorted(range(m), key=lambda i: starts[i])
            o = []
            ok = True
            for t in range(m):
                i = order[t]
                j2 = order[(t + 1) % m]
                gap = (starts[j2] - starts[i]) % p
                o.append(s[i] - gap)
            # [B2] telescoping sum identity
            if sum(o) != ov:
                bad_sum += 1
            # [B1'] reach bound: gap <= s_max  <=>  o_i >= s_i - s_max
            if any(gap > s_max for gap in
                   [(starts[order[(t + 1) % m]] - starts[order[t]]) % p
                    for t in range(m)]):
                bad_gap += 1
            if all(oi >= 0 for oi in o):
                n_chain += 1
            else:
                nested += 1
        # the pure-chain count for this size tuple (o >= 0 simplex with
        # o_i <= s_{sigma(i)}, summed over orders) -- computed as
        # sum_sigma #{o in Z^m : o>=0, sum=ov, o_i <= s_{sigma(i)}}
        cnt = 0
        import itertools
        for perm in itertools.permutations(range(m)):
            if perm[0] != 0:
                continue
            # o_t for t=0..m-1 with o_t in [0, s_{perm[t]}], sum = ov
            caps = [s[perm[t]] for t in range(m)]
            # DP count
            dp = [1] + [0] * ov
            for cap in caps:
                ndp = [0] * (ov + 1)
                run = 0
                for v in range(ov + 1):
                    run += dp[v]
                    if v - cap - 1 >= 0:
                        run -= dp[v - cap - 1]
                    ndp[v] = run
                dp = ndp
            cnt += dp[ov]
        tot_chain += cnt
        print('  sizes %s ov=%d: |Sol|=%6d  pure-chain=%6d  '
              'non-chain(nested)=%5d' % (s, ov, len(lst), cnt,
                                         len(lst) - n_chain))
    print()
    print('TOTAL solutions: %d | pure-chain total: %d | nested observed: %d'
          % (tot, tot_chain, nested))
    print('[B2] sum-identity violations: %d / %d' % (bad_sum, tot))
    print("[B1'] reach-bound violations: %d / %d" % (bad_gap, tot))
    print('note: |Sol| counts placements (a2..a5); pure-chain counts '
          '(order, o) pairs -- ties/degenerate overcounts make the '
          'difference; the o >= 0 subset of Sol = %d' % (tot - nested))

if __name__ == '__main__':
    main()
