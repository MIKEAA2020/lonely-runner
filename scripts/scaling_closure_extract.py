#!/usr/bin/env python3
"""
Scaling-closure lemma attempt -- calibration extraction (READ-ONLY).

Extracts from the COMMITTED outputs (out_higher_rung_step1.json,
out_higher_rung_step2.json) the three quantities the scaling-closure
lemma's proof needs as calibration:

  [POOL]   the T=8 (1K,5U) zoo pool counts at p=17/19 vs the counting law
           P = (D)_{T-4} (falling factorial of distinct canonical
           differences, D = (p-1)/2, n_d = T-4 free differences).
  [SOL]    the per-(family, size) census solution counts |Sol| for the
           committed T=7 and T=8-critical censuses vs the junction-volume
           law |Sol| <= (R(T) * k)^(m-2).
  [EXIT]   the committed mu-candidate counts (= |Sol| at the anchor fiber)
           from out_higher_rung_step2.json (cross-check of [SOL]).

No new census, no new primes, no cell-emptiness runs.  Read-only.
"""
import json
from math import comb

S1 = '/home/z/my-project/scripts/out_higher_rung_step1.json'
S2 = '/home/z/my-project/scripts/out_higher_rung_step2.json'

with open(S1) as f:
    step1 = json.load(f)
with open(S2) as f:
    step2 = json.load(f)

out = {}

# ----------------------------------------------------------------------
# [POOL] the zoo pool law P = (D)_{T-4}
# ----------------------------------------------------------------------
print('=== [POOL] T=8 (1K,5U) zoo pools vs the counting law ===')
pool_res = {}
for pp in sorted(step1['T8_m5_census'].keys()):
    dat = step1['T8_m5_census'][pp]
    p = int(pp[1:])
    D = (p - 1) // 2
    n_d = 4                      # T-4 free differences at T=8
    falling = 1
    for i in range(n_d):
        falling *= (D - 1 - i)   # values in [2, D] (excluding 1), ordered
    pool_res[pp] = {
        'p': p, 'k': dat['k'], 'eps': dat['eps'],
        'committed_families': dat['n_families'],
        'committed_distinct': dat['n_distinct_families'],
        'law_(D-1)_{T-4}': falling,
        'n_size_tuples_admissible': dat['n_families'] /
        max(1, dat['n_distinct_families']),
    }
    print('  %s: committed distinct families %d vs law (D-1)_{4} = %d  '
          '(k=%d, eps=%d)' % (pp, dat['n_distinct_families'], falling,
                              dat['k'], dat['eps']))
out['pool'] = pool_res

# ----------------------------------------------------------------------
# [SOL] junction-volume law |Sol| <= (R k)^(m-2) on the committed censuses
# ----------------------------------------------------------------------
print()
print('=== [SOL] junction-volume law on committed censuses ===')


def sol_stats(section, label, m_ap):
    """Per prime: distribution of per-(family,size) config counts n."""
    res = {}
    for pp, dat in sorted(step1[section].items(),
                          key=lambda kv: kv[1]['k']):
        k = dat['k']
        ns = [v['n'] for v in dat['famsize_detail'].values()]
        if not ns:
            continue
        ns.sort()
        nmax = ns[-1]
        med = ns[len(ns) // 2]
        # the junction law: |Sol| <= (R k)^(m-2) -> R_needed = nmax^(1/(m-2))/k
        R_needed = (nmax ** (1.0 / (m_ap - 2))) / k if k else None
        res[pp] = {'k': k, 'eps': dat['eps'],
                   'n_keys': len(ns), 'n_configs': sum(ns),
                   'n_max': nmax, 'n_median': med,
                   'R_needed_for_max': round(R_needed, 3) if R_needed else None}
        print('  %s %s: k=%2d eps=%d: %5d keys, max |Sol|=%6d, '
              'median=%5d, R(max)=%.2f'
              % (label, pp, k, dat['eps'], len(ns), nmax, med, R_needed))
    return res


out['sol_T7'] = sol_stats('T7_census', 'T=7', 4)
out['sol_T8c'] = sol_stats('T8critical_census', 'T=8c', 4)

# aggregate: the R that dominates all primes
for nm, sec in (('T=7', out['sol_T7']), ('T=8c', out['sol_T8c'])):
    Rs = [v['R_needed_for_max'] for v in sec.values()
          if v['R_needed_for_max']]
    if Rs:
        print('  -> %s: max over primes of R_needed = %.2f '
              '(junction law |Sol| <= (R k)^2 with this R)' % (nm, max(Rs)))

# ----------------------------------------------------------------------
# [EXIT] committed mu-candidate counts (anchor solution counts)
# ----------------------------------------------------------------------
print()
print('=== [EXIT] committed anchor |Sol| (mu-candidate counts) ===')
exit_res = {}
for key in ('exit_p17', 'exit_p29', 'exit_p43'):
    if key not in step2:
        continue
    d = {}
    for fam, v in step2[key].items():
        d[fam] = {'n_fibers': v['n_fibers'], 'best': v['best'],
                  'n_cand_mu': v['n_cand_mu']}
    exit_res[key] = d
    mu_counts = sorted(v['n_cand_mu'] for v in d.values())
    print('  %s: mu-candidates per family: %s (best exit max %d)'
          % (key, mu_counts, max(v['best'] for v in d.values())))
out['exit'] = exit_res

with open('/home/z/my-project/scripts/out_scaling_closure_extract.json',
          'w') as f:
    json.dump(out, f, indent=1)
print()
print('Wrote scripts/out_scaling_closure_extract.json')
