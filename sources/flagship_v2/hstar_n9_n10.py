#!/usr/bin/env python3
"""h*-polynomial extension to n = 9, 10 (harmonic and rung-2 families),
plus a provenance re-computation of the n = 3..8 rows.

Every number in the paper's h*-table (Appendix A) traces to this run:
    output : scripts/out_hstar_n9_n10.json

Checks performed per instance (all exact unless stated):
  * e_k = sum_{|S|=k} gcd(v_{S^c})            (ehr_coeffs)
  * h* via Eulerian expansion == h* via the standard inversion from
    Ehr(0..d+1)                                (asserted equal)
  * Ehrhart reconstructed from h* at t = 0..d+2 (asserted equal)
  * DFS coset count at t = 1, 2 == Ehr(1), Ehr(2)  (independent count)
  * Ehr(-1) == (-1)^d * (interior count) sanity: signed reciprocity
    spot-asserted where an independent interior count exists
  * volume: e_{d} == sum v_j
  * lambda exact (breakpoint/crossover enumeration), asserted equal to
    1/(n+1) (harmonic) resp. 2/(2n+1) (rung-2)
  * real-rootedness: numpy companion-matrix roots AND (primary) mpmath
    polyroots at 40 dps on the exact integer coefficients
  * log-concavity h_k^2 >= h_{k-1} h_{k+1}, Newton inequalities
    (binomial-normalized), min ratio (the "strictness" column)
"""
import sys, json, math, time
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction
from lrc_zono_lib import (ehr_coeffs, hstar_from_ehr, hstar_by_inversion,
                          ehr_from_hstar, ehr_poly_eval, lambda_exact,
                          dfs_count)

OUT = '/home/z/my-project/scripts/out_hstar_n9_n10.json'


def min_ratio(hs):
    """min_k h_k^2 / (h_{k-1} h_{k+1})  (log-concavity strictness)."""
    best = None
    for k in range(1, len(hs) - 1):
        if hs[k - 1] and hs[k + 1]:
            r = Fraction(hs[k] * hs[k], hs[k - 1] * hs[k + 1])
            if best is None or r < best:
                best = r
    return best


def newton_check(hs):
    """Newton's inequalities: (h_k/C(d,k))^2 >= (h_{k-1}/C(d,k-1))(h_{k+1}/C(d,k+1))."""
    d = len(hs) - 1
    ok = True
    worst = None
    for k in range(1, len(hs) - 1):
        lhs = Fraction(hs[k] * hs[k] * math.comb(d, k - 1) * math.comb(d, k + 1))
        rhs = Fraction(hs[k - 1] * hs[k + 1] * math.comb(d, k) * math.comb(d, k))
        if lhs < rhs:
            ok = False
        r = lhs / rhs
        if worst is None or r < worst:
            worst = r
    return ok, worst


def roots_info(hs):
    """Real-rootedness: numpy first pass, mpmath (40 dps, exact int coeffs) as
    the primary verifier. Returns (roots, n_real_numpy, real_numpy,
    real_mpmath, roots_rounded)."""
    import numpy as np
    coeffs_desc = [int(x) for x in hs[::-1]]      # exact ints
    r = np.roots([float(x) for x in coeffs_desc])
    n_real_np = sum(1 for z in r
                    if abs(z.imag) <= 1e-6 * (1.0 + abs(z.real)))
    real_np = (n_real_np == len(r))
    mp_res = None
    try:
        import mpmath as mp
        mp.mp.dps = 45
        rr = mp.polyroots([mp.mpf(x) for x in coeffs_desc],
                          maxsteps=300, error=False)
        n_real_mp = sum(1 for z in rr
                        if abs(mp.im(z)) <= mp.mpf(10) ** (-15) * (1 + abs(mp.re(z))))
        mp_res = (n_real_mp == len(rr))
    except Exception as ex:                       # pragma: no cover
        mp_res = 'unavailable (%s)' % ex
    roots_rounded = sorted(round(float(z.real), 4) for z in r)
    return ([complex(z) for z in r], n_real_np, real_np, mp_res, roots_rounded)


def analyze(v, family, do_dfs=True):
    n = len(v)
    d = n - 1
    t0 = time.time()
    e = ehr_coeffs(v)
    h = hstar_from_ehr(e, d)
    # --- h* inversion cross-check
    Ehr_vals = [ehr_poly_eval(e, t) for t in range(d + 2)]
    h2 = hstar_by_inversion(Ehr_vals, d)
    assert h2 == h, ('h* inversion mismatch', v, h, h2)
    # --- reconstruction
    for t in range(d + 3):
        assert ehr_from_hstar(h, d, t) == ehr_poly_eval(e, t), ('recon', v, t)
    # --- independent DFS count
    dfs = {}
    if do_dfs:
        for t in (1, 2):
            tt = time.time()
            c = dfs_count(v, t)
            assert c == Ehr_vals[t], ('dfs', v, t, c, Ehr_vals[t])
            dfs['t%d' % t] = c
            dfs['sec_t%d' % t] = round(time.time() - tt, 2)
    # --- lambda
    lam, tstar = lambda_exact(v)
    if family in ('harmonic', 'rung2'):
        lam_exp = Fraction(1, n + 1) if family == 'harmonic' \
            else Fraction(2, 2 * n + 1)
        assert lam == lam_exp, ('lambda', v, lam, lam_exp)
    delta = Fraction(1, 2) - lam
    thr = Fraction(n - 1, 2 * (n + 1))
    # --- coefficient structure
    lc = all(h[k] * h[k] >= h[k - 1] * h[k + 1] for k in range(1, len(h) - 1))
    mr = min_ratio(h)
    nw, nw_worst = newton_check(h)
    vol_ok = (e[d] == sum(v))
    roots, n_real_np, real_np, real_mp, roots_r = roots_info(h)
    # --- signed reciprocity spot value
    ehr_m1 = ehr_poly_eval(e, -1)
    return {
        'family': family, 'v': v, 'n': n, 'd': d,
        'e_coeffs': e,
        'hstar': h, 'hstar_str': ' + '.join(
            ('%d t^%d' % (c, k)) if k else str(c)
            for k, c in enumerate(h) if c),
        'Ehr_0_to_d1': Ehr_vals,
        'Ehr_minus1': ehr_m1,
        'vol_sum_v': sum(v), 'vol_ok': vol_ok,
        'dfs': dfs,
        'lambda': str(lam), 'lambda_t': str(tstar),
        'delta': str(delta), 'thr': str(thr),
        'delta_margin': str(delta - thr),
        'degree': len(h) - 1,
        'real_rooted_numpy': real_np, 'n_real_numpy': n_real_np,
        'real_rooted_mpmath': real_mp,
        'roots_rounded': roots_r,
        'log_concave': lc, 'newton_ok': nw,
        'newton_worst_ratio': str(nw_worst),
        'strictness_min_ratio': str(mr),
        'strictness_f': float(mr),
        'strictness_times_n': float(mr) * n,
        'sec': round(time.time() - t0, 2),
    }


def families(n):
    return {'harmonic': list(range(1, n + 1)),
            'rung2': list(range(1, n)) + [2 * n]}


def main():
    rows = []
    print('%-9s %-2s %-58s %s' % ('family', 'n', 'h*(t)', 'real/logc/Newton strictness'))
    for n in range(3, 11):
        for fam, v in families(n).items():
            row = analyze(v, fam, do_dfs=True)
            rows.append(row)
            print('%-9s %-2d %-58s %s/%s/%s %s' % (
                fam, n, row['hstar_str'][:58],
                row['real_rooted_mpmath'], row['log_concave'],
                row['newton_ok'],
                ('%.4f  (n*S = %.1f)' % (row['strictness_f'],
                                         row['strictness_times_n']))))
    # strictness decay fit: S(n) ~ c/n over the harmonic family
    harm = [r for r in rows if r['family'] == 'harmonic']
    cvals = [r['strictness_times_n'] for r in harm]
    json.dump({'rows': rows,
               'strictness_decay': {
                   'harmonic_n_times_strictness': cvals,
                   'rung2_n_times_strictness':
                       [r['strictness_times_n'] for r in rows
                        if r['family'] == 'rung2']},
               'all_real_rooted': all(
                   (r['real_rooted_mpmath'] is True) for r in rows),
               'all_log_concave': all(r['log_concave'] for r in rows),
               'all_newton': all(r['newton_ok'] for r in rows)},
              open(OUT, 'w'), indent=1)
    print('\nall real-rooted (mpmath):', all(r['real_rooted_mpmath'] is True
                                            for r in rows))
    print('all log-concave:', all(r['log_concave'] for r in rows))
    print('all Newton:', all(r['newton_ok'] for r in rows))
    print('n*strictness (harmonic):', [round(x, 2) for x in cvals])
    print('written:', OUT)


if __name__ == '__main__':
    main()
