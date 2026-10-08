#!/usr/bin/env python3
"""t25b_reflection_addendum.py -- completes the T-25 session record:
(a) the reflection pair (j0 -> p-j0) two-fiber survival fraction at all
    three primes for the session's six committed families (Lemma reflect:
    survival fraction exactly 1 -- in automaton language: the deletion at
    the lambda = -1 step is empty);
(b) the p=29 homogenized reflection sub-shifts (capped out of the main
    run's first-8 window).
Committed scope only. Output: scripts/out_t25b_reflection.json.
"""
import json
import sys

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import ball_mask, k_of
from scaling_closure_verify import fiber_data_m5, sign_patterns_for
from t25_transfer_matrix_session import (build_sizedict, T,
                                          homogenized_pair_subshift, FAMS)

OUT = {'pairs': {}, 'subshifts': {}}

for p in (17, 19, 29):
    for label, fam in FAMS[p]:
        sizedict, _ = build_sizedict(p, fam)
        pat = sign_patterns_for(p, fam)[0]
        k = k_of(p, T)
        bm = ball_mask(p, k)
        U = [j for j in range(p) if not (bm >> j) & 1]
        fibers = []
        for j in U:
            sz, cc = fiber_data_m5(p, T, pat, j)
            if sz in sizedict:
                fibers.append((j, sz, cc))
        j0, sz0, c0 = fibers[0]
        S0 = sizedict[sz0]
        inv_j0 = pow(j0, -1, p)
        jr = p - j0                     # the reflection fiber
        if jr in [f[0] for f in fibers]:
            idx = [f[0] for f in fibers].index(jr)
            _, szr, ccr = fibers[idx]
            lam = (jr * inv_j0) % p     # = p-1
            tau = tuple((ccr[i] - lam * c0[i]) % p for i in range(4))
            Solr = sizedict[szr]
            surv = sum(1 for a in S0
                       if tuple((lam * a[i] + tau[i]) % p
                                for i in range(4)) in Solr)
            key = 'p%d_%s' % (p, label)
            OUT['pairs'][key] = {
                'j0': j0, 'j_reflect': jr, 'lambda': lam,
                'n_sol0': len(S0), 'n_sol_reflect': len(Solr),
                'survivors': surv,
                'survival_fraction': round(surv / float(len(S0)), 4),
                'sizes_match': list(sz0) == list(szr)}
            print('%s: reflection pair j0=%d -> j=%d (lam=-1): '
                  'survivors %d / %d = %.4f | sizes match: %s'
                  % (key, j0, jr, surv, len(S0),
                     surv / float(len(S0)), list(sz0) == list(szr)))
            if p == 29:
                h = homogenized_pair_subshift(p, lam, tau, Solr)
                h['j'] = jr
                OUT['subshifts'][key] = {
                    kk: h[kk] for kk in ('j', 'rho', 'n_inside_cycles',
                                         'inside_total',
                                         'asymptotic_fraction',
                                         'surv_T', 'sim_mismatches')}
                print('   p29 reflection subshift: rho=%d inside=%d(%d) '
                      'asymp=%.6f sim_bad=%d'
                      % (h['rho'], h['n_inside_cycles'], h['inside_total'],
                         h['asymptotic_fraction'], h['sim_mismatches']))
        else:
            print('p%d_%s: reflection fiber %d not admissible'
                  % (p, label, jr))

with open('/home/z/my-project/scripts/out_t25b_reflection.json', 'w') as f:
    json.dump(OUT, f, indent=1, default=str)
print('DONE; wrote scripts/out_t25b_reflection.json')
