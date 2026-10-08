#!/usr/bin/env python3
"""
strand_rigidity_probe2.py -- T-22 second stage.

(a) DIRECT test of the Lemma V proof's per-size-assignment bijection claim:
    p=19, S=(2,4,5,8), delta=5, all-equal size tuple (5,5,5,5,5) (sigma
    trivial).  Construct the image of every solution under the proof's
    recipe (dilate by delta^-1, re-anchor at the new normalizer,
    re-parameterize reflected arcs) and test membership in
    Sol((3,4,6,8), (5,5,5,5,5)).  Localizes whether the committed proof
    step "the map is a bijection at each size assignment" holds.
(b) The gamma-envelope: max over admissible keys of |Pi_i|/ov (the
    marginal-to-overlap ratio), against the cascade threshold
    T/(T-6) = 4 at T=8.  Also the fullness fraction max |Pi_i|/p.
    Per-key records saved for the note.
(c) p=29 targeted marginals (raw + pairwise differences) for the two
    admissible committed sets.
"""
import json
import sys
from collections import defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples
from scaling_closure_verify import solutions_m5

T, M = 8, 5


def enum(p, fam, szts):
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    return {key[:5]: sols[key] for key in sorted(sols)}


def covers(p, szt, fam, a):
    """Direct covering test for start tuple a (moving arcs) + normalizer."""
    pts = set(range(szt[0]))
    for i, d in enumerate(fam):
        pts |= {(a[i] + t * d) % p for t in range(szt[i + 1])}
    return len(pts) == p


def main():
    OUT = {}

    # ---------- (a) the Lemma V per-assignment bijection test ----------
    p = 19
    S = (2, 4, 5, 8)
    delta = 5
    inv = pow(delta, -1, p)          # 4
    S2 = (3, 4, 6, 8)                # renorm_delta(S), verified earlier
    szt = (5, 5, 5, 5, 5)
    sols = enum(p, S, [szt])[szt]
    sols2 = enum(p, S2, [szt])[szt]
    set2 = set(sols2)
    print('=== (a) Lemma V per-assignment bijection test, p=%d ===' % p)
    print('  |Sol(%s)| = %d ; |Sol(%s)| = %d (all-5 szt, sigma trivial)'
          % (S, len(sols), S2, len(sols2)))
    # slots of S: d2=2, d3=4, d4=5(delta), d5=8; normalizer d=1 size 5.
    # after dilation by inv: normalizer -> d=inv? no: d=1*inv=4 (start 0);
    # d2=2 -> 2*inv=8 (start inv*a2); d3=4 -> 4*inv=16 -> canon 3,
    #   re-parameterize start inv*a3 + 16*(s-1);
    # d4=5=delta -> 1 (start inv*a4) = NEW NORMALIZER;
    # d5=8 -> 8*inv=32=13 -> canon 6, start inv*a5 + 13*(s-1).
    # re-anchor by -inv*a4.
    lam = {2: (2 * inv) % p, 4: (4 * inv) % p, 5: (5 * inv) % p,
           8: (8 * inv) % p}
    s = 5
    ok = bad = 0
    bad_examples = []
    for a in sols:
        a2, a3, a4, a5 = a
        img = {}
        # new normalizer: old delta arc, difference 1, start 0 after anchor
        # old normalizer (d=1): -> difference 4, start 0 - inv*a4
        img[4] = (0 - inv * a4) % p
        # old d2: difference 8 (no canon flip, 8 <= 9), start inv*a2 - inv*a4
        img[8] = (inv * a2 - inv * a4) % p
        # old d3: raw lam 16 -> canon 3 (flip), start inv*a3 + 16*(s-1) - inv*a4
        img[3] = (inv * a3 + lam[4] * (s - 1) - inv * a4) % p
        # old d5: raw lam 13 -> canon 6 (flip), start inv*a5 + 13*(s-1) - inv*a4
        img[6] = (inv * a5 + lam[8] * (s - 1) - inv * a4) % p
        b = (img[3], img[4], img[6], img[8])   # slots of S2 sorted (3,4,6,8)
        if b in set2:
            ok += 1
        else:
            bad += 1
            if len(bad_examples) < 3:
                bad_examples.append((a, b))
    print('  images landing in Sol(S2): %d / %d ; misses: %d'
          % (ok, len(sols), bad))
    for a, b in bad_examples:
        cov = covers(p, szt, S2, b)
        print('    a=%s -> b=%s ; direct covering test of b: %s'
              % (a, b, cov))
    # injectivity of the construction
    imgs = []
    for a in sols:
        a2, a3, a4, a5 = a
        b = (((inv * a3 + lam[4] * (s - 1) - inv * a4) % p,
              (0 - inv * a4) % p,
              (inv * a5 + lam[8] * (s - 1) - inv * a4) % p,
              (inv * a2 - inv * a4) % p))
        imgs.append(b)
    print('  construction injective: %s' % (len(set(imgs)) == len(imgs)))
    OUT['vtest'] = {'p': p, 'S': list(S), 'S2': list(S2), 'delta': delta,
                    'n_source': len(sols), 'n_target': len(sols2),
                    'images_in_target': ok, 'misses': bad,
                    'bad_examples': [list(a) + list(b)
                                     for a, b in bad_examples]}

    # ---------- (b) gamma envelopes + per-key records ----------
    print()
    print('=== (b) gamma-envelope vs the cascade threshold '
          'T/(T-6) = %g ===' % (T / (T - 6)))
    OUT['gamma'] = {}
    for p, sets_mode in [(17, 'full'), (19, 'full'), (29, 'targeted')]:
        k, (on, off), cls, eps = profile_sizes(p, T)
        szts = size_tuples(M, on, off, p)
        if sets_mode == 'full':
            from itertools import combinations
            D = (p - 1) // 2
            sets = [tuple(sorted(c)) for c in
                    combinations(range(2, D + 1), 4)]
        else:
            sets = []
            for x in range(3, (p - 1) // 2 + 1):
                if x not in (2, 4, 8):
                    sets.append(tuple(sorted((2, 4, 8, x))))
            path = [2, 4, 8, 13, 3, 6, 12, 5, 10, 9, 11, 7, 14]
            for i in range(len(path) - 3):
                sets.append(tuple(sorted(path[i:i + 4])))
            sets = sorted(set(sets))
        gmax = 0.0
        gkey = None
        fmax = 0.0
        perkey = []
        for s in sets:
            rec = enum(p, s, szts)
            for szt, lst in rec.items():
                if not lst:
                    continue
                ov = sum(szt) - p
                proj = [len({a[i] for a in lst}) for i in range(4)]
                mp = max(proj)
                g = mp / float(ov)
                if ov >= 2 and g > gmax:
                    gmax, gkey = g, (s, szt)
                fmax = max(fmax, mp / float(p))
                perkey.append({'set': list(s), 'szt': list(szt),
                               'ov': ov, 'n': len(lst), 'proj': proj})
        print('  p=%d (%d sets): gamma-env (ov>=2) = %.3f at %s ; '
              'fullness max |Pi|/p = %.3f ; admissible keys: %d'
              % (p, len(sets), gmax, gkey, fmax, len(perkey)))
        OUT['gamma']['p%d' % p] = {'gmax': gmax,
                                   'gkey': [list(gkey[0]), list(gkey[1])],
                                   'fmax': fmax, 'n_keys': len(perkey)}
        with open('/home/z/my-project/scripts/out_strand_gamma_p%d.json'
                  % p, 'w') as f:
            json.dump(perkey, f)

    # ---------- (c) p=29 targeted marginals ----------
    print()
    print('=== (c) p=29 targeted marginals ===')
    p = 29
    k, (on, off), cls, eps = profile_sizes(p, T)
    szts = size_tuples(M, on, off, p)
    OUT['p29_targeted'] = {}
    for s0 in [(2, 4, 7, 8), (2, 4, 8, 13)]:
        rec = enum(p, s0, szts)
        nz = [(szt, lst) for szt, lst in rec.items() if lst]
        if not nz:
            print('  %s: inadmissible' % (s0,))
            continue
        szt, lst = max(nz, key=lambda kv: len(kv[1]))
        ov = sum(szt) - p
        raw = [len({a[i] for a in lst}) for i in range(4)]
        pd = [len({(a[j] - a[i]) % p for a in lst})
              for i in range(4) for j in range(i + 1, 4)]
        print('  %s @ %s (ov=%d): n=%d' % (s0, szt, ov, len(lst)))
        print('    raw  marginals: %s  (ratios to ov: %s)'
              % (raw, [round(r / float(ov), 2) for r in raw]))
        print('    pairdiff marginals: %s' % pd)
        OUT['p29_targeted'][str(list(s0))] = {
            'szt': list(szt), 'ov': ov, 'n': len(lst), 'raw': raw,
            'raw_over_ov': [round(r / float(ov), 3) for r in raw],
            'pairdiff': pd}

    with open('/home/z/my-project/scripts/out_strand_rigidity2.json',
              'w') as f:
        json.dump(OUT, f, indent=1, default=int)
    print()
    print('wrote scripts/out_strand_rigidity2.json '
          '+ out_strand_gamma_p{17,19,29}.json')


if __name__ == '__main__':
    main()
