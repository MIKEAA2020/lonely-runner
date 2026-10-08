#!/usr/bin/env python3
"""Retroactive provenance re-assertion for every number cited in the
type-mismatch paper.

Reads the two persisted run JSONs of the zonotope program
    scripts/out_lrc_zono_final.json      (zonotope-1)
    scripts/out_lrc_zono_rescan.json     (zonotope-2 re-check)
and re-asserts, by INDEPENDENT exact computation (the sigma-formula D of
certify_rung2_n5.py AND the V-norm lattice DFS D_exact of lrc_zono_lib,
plus lambda_exact), every witness value, rho value, lambda/delta value,
and Spearman correlation that the paper cites. Also re-computes the h*
strictness rows n = 3..10 from e-coefficients (via hstar_n9_n10.analyze).

Output: scripts/out_paper_provenance.json
        scripts/out_paper_provenance.txt   (human-readable ledger)
"""
import sys, json, math
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F

import certify_rung2_n5 as C
from lrc_zono_lib import make_forms, D_exact, lambda_exact

FIN = '/home/z/my-project/scripts/out_lrc_zono_final.json'
RES = '/home/z/my-project/scripts/out_lrc_zono_rescan.json'
OUTJ = '/home/z/my-project/scripts/out_paper_provenance.json'
OUTT = '/home/z/my-project/scripts/out_paper_provenance.txt'

ledger = []          # (item, claimed, recomputed, ok)


def add(item, claimed, recomputed):
    ok = (claimed == recomputed)
    ledger.append({'item': item, 'claimed': str(claimed),
                   'recomputed': str(recomputed), 'ok': ok})
    return ok


def D_two_ways(c, v):
    """D via sigma-formula AND via lib V-norm DFS; assert equal; return it."""
    a = C.D_at(c, v)
    b = D_exact(tuple(c), v, make_forms(v))
    assert a == b, ('two-way mismatch', v, c, a, b)
    return a


def main():
    fin = json.load(open(FIN))
    res = json.load(open(RES))

    # ---------------------------------------------------------------- lambda/delta
    for n in range(3, 11):
        for fam, v, lam_exp in [
                ('harmonic', list(range(1, n + 1)), F(1, n + 1)),
                ('rung2', list(range(1, n)) + [2 * n], F(2, 2 * n + 1))]:
            lam, _ = lambda_exact(v)
            add('lambda %s n=%d' % (fam, n), lam_exp, lam)

    # ---------------------------------------------------------------- witnesses
    # n=5, n=6 core witnesses (paper Table: the rho-refutation archive)
    WIT = [
        ((1, 2, 3, 4, 5), [F(1, 32), F(15, 32), F(3, 32), F(7, 8)], F(97, 288)),
        ((1, 2, 3, 4, 10), [F(1, 2), F(0), F(1, 2), F(1, 2)], F(7, 22)),
        ((1, 2, 3, 4, 5, 6), [F(0), F(5, 16), F(15, 16), F(3, 4), F(1, 2)],
         F(29, 80)),
        ((1, 2, 3, 4, 5, 7), [F(1, 16), F(9, 16), F(0), F(3, 8), F(0)],
         F(4, 11)),
    ]
    for v, w, val in WIT:
        add('witness D %s' % (v,), val, D_two_ways(w, v))

    # every Phase D (n=5) and Phase E (n=6) entry of the re-check
    for phase, label in ((res['D'], 'n5'), (res['E'], 'n6')):
        for e in phase:
            v = tuple(e['v'])
            if 'witness' in e:
                w = [F(s) for s in e['witness']]
                add('witness D %s (rescan %s m=%s)' % (v, label, e.get('m')),
                    F(e['witness_val']), D_two_ways(w, v))
            if 'lambda' in e:
                lam, _ = lambda_exact(v)
                add('lambda %s (rescan)' % (v,), F(e['lambda']), lam)
            if 'rho_lo_exact' in e and 'witness' in e:
                # rho_lo == witness value at the recorded argmax
                add('rho_lo %s (rescan)' % (v,), F(e['rho_lo_exact']),
                    D_two_ways(w, v))

    # n=3 sweep (final.json): rho at recorded argmax
    for e in fin['rho_n3_sweep']:
        v = tuple(e['v'])
        if 'argmax' in e and 'rho' in e:
            a = [F(s) for s in e['argmax']]
            add('rho_exact %s (n3 sweep)' % (v,), F(e['rho']),
                D_two_ways(a, v))
    # n=3 Phase B of rescan
    for e in res['B']:
        v = tuple(e['v'])
        if 'argmax' in e and 'rho_exact' in e:
            a = [F(s) for s in e['argmax']]
            add('rho_exact %s (rescan B)' % (v,), F(e['rho_exact']),
                D_two_ways(a, v))

    # n=4 certified rows: c* depth == delta == rho (the deepest-hole value)
    for e in res['C']:
        v = tuple(e['v'])
        lam, _ = lambda_exact(v)
        add('lambda %s (rescan C)' % (v,), F(e['lambda']), lam)
        cstar = tuple(F(1 - vj, 2) % 1 for vj in v[1:])
        d_cstar = D_two_ways(cstar, v)
        add('D(c*) %s (rescan C)' % (v,), F(e['delta']), d_cstar)

    # rho_n5_n6 enclosure lower ends (final.json)
    for e in fin['rho_n5_n6']:
        v = tuple(e['v'])
        w = [F(s) for s in e['witness']]
        add('witness D %s (final)' % (v,), F(e['witness_val']),
            D_two_ways(w, v))

    # ---------------------------------------------------------------- Spearman
    from scipy.stats import spearmanr
    for key, expect in (('correlation_n3', 0.32), ('correlation_n4', -0.31)):
        cc = fin[key]
        stored = cc.get('spearman',
                        cc.get('spearman_strictness_vs_delta_margin'))
        rho = spearmanr(cc['strictness'], cc['delta_margin']).statistic
        ok = bool(abs(rho - stored) < 1e-9 and
                  abs(round(rho, 2) - expect) < 0.005)
        ledger.append({'item': 'spearman %s' % key,
                       'claimed': '%.2f (stored %.6f)' % (expect, stored),
                       'recomputed': '%.6f' % rho, 'ok': ok})

    # strictness arrays independently recomputed from e-coefficients
    import hstar_n9_n10 as H
    n3_vs = [tuple(e['v']) for e in fin['hstar_n3_sweep']]
    n4_vs = [tuple([1, 2, 3, m]) for m in range(5, 13)]   # m = 5..12
    for key, vs in (('n3', n3_vs), ('n4', n4_vs)):
        assert len(vs) == len(fin['correlation_%s' % key]['strictness'])
        cc = fin['correlation_%s' % key]
        for i, v in enumerate(vs):
            row = H.analyze(list(v), 'sweep', do_dfs=False)
            ok = bool(abs(row['strictness_f'] - cc['strictness'][i]) < 1e-9)
            ledger.append({'item': 'strictness %s (v=%s)' % (key, v),
                           'claimed': '%.6f' % cc['strictness'][i],
                           'recomputed': '%.6f' % row['strictness_f'],
                           'ok': ok})
            lam, _ = lambda_exact(list(v))
            # delta margin (stored as thr - delta, positive = slack)
            n = len(v)
            dm = F(n - 1, 2 * (n + 1)) - (F(1, 2) - lam)
            ok2 = bool(abs(float(dm) - cc['delta_margin'][i]) < 1e-9)
            ledger.append({'item': 'delta_margin %s (v=%s)' % (key, v),
                           'claimed': '%.6f' % cc['delta_margin'][i],
                           'recomputed': '%.6f' % float(dm), 'ok': ok2})

    # h* strictness main table n=3..8 (both families) - recomputed and
    # compared with the stored table (strictness stored as decimal string)
    for e in fin['hstar_table']:
        v = tuple(e['v'])
        row = H.analyze(list(v), e.get('family', '?'), do_dfs=False)
        ok = bool(abs(row['strictness_f'] - float(e['strictness'])) < 1e-9)
        ledger.append({'item': 'strictness %s (%s)' % (v, e.get('family')),
                       'claimed': e['strictness'],
                       'recomputed': '%.6f' % row['strictness_f'],
                       'ok': ok})

    # ---------------------------------------------------------------- report
    nbad = sum(1 for e in ledger if not e['ok'])
    with open(OUTT, 'w') as f:
        f.write('PAPER PROVENANCE LEDGER  (%d items, %d failures)\n\n' %
                (len(ledger), nbad))
        for e in ledger:
            f.write('[%s] %-52s claimed=%-14s recomputed=%s\n' %
                    ('OK' if e['ok'] else 'FAIL', e['item'],
                     e['claimed'], e['recomputed']))
    json.dump({'items': ledger, 'failures': nbad}, open(OUTJ, 'w'), indent=1)
    print('items: %d, failures: %d' % (len(ledger), nbad))
    for e in ledger:
        if not e['ok']:
            print('FAIL:', e)
    print('written:', OUTJ, OUTT)


if __name__ == '__main__':
    main()
