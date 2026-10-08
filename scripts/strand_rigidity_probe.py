#!/usr/bin/env python3
"""
strand_rigidity_probe.py -- T-22: the multi-character rigidity attack on
the strand census (H1g at T=8).

The reviewer's directive (post-T-21): attack the strand census, but with a
different tool -- the multi-character rigidity that was calibrated but not
developed.  In the committed record that phenomenon is the projection-
pinning / p-independent-counts law, calibrated at four sites and never
turned into a lemma:
  (i) T-7 zoo: solution sets ~k^2 but one-coordinate projections
      11,12,18,24,...,42,42,42 -- saturating at 42 while p grows;
 (ii) T-8c census: eleven families with placement counts 24/312 constant
      in p over seventeen primes (Conjecture SL's site);
(iii) the pinning-size lemma input: relative-offset sets <= 6 elements,
      counts 28/68 p-independent, verified p <= 113;
 (iv) T-8 zoo probes: coordinate projections [16,17,14,12] (top family)
      vs [9,9,7,8] (min class) at p=17.

This script develops the tool against the actual obligation:

  [A] FULL orbit-level rigidity census at p=17 and p=19: every
      renormalization orbit (Lemma V' representatives), all 32 size
      tuples -- per-character projection sizes |Pi_i(Sol)|, box product,
      joint density, volumes.  Exhaustive over the committed cells; no
      sampling, hence no intervals.
  [V'] Lemma V' machine verification: the multiset (over size tuples) of
      sorted projection profiles is identical across ALL members of each
      renormalization orbit (p=17: 35 sets; p=19: 70 sets).
  [B] TARGETED rigidity census at p=29 (k=3): the committed targeted
      classes -- every set containing the binary-cascade core {2,4,8},
      plus the doubling-path quadruples -- all size tuples, GT'd against
      the committed T-20 state where the ordered family was measured.
  [C] The k-trend of the pinning envelope: max |Pi_i| over admissible
      keys at p=17/19 (k=2) vs p=29 (k=3).  Saturating => the rigidity
      tool lives at the zoo cell; linear-in-p => it dies there.
  [D] Resonance classification of the measured sets: canonical ratio
      spectrum, resonant-pair count (lambda_bar <= 4), cascade flag --
      the mechanism correlate for the pinning lemma.
  [GT] Per-key solution counts cross-checked against the committed T-20
      census states (h1_t8_state.json, h1_t8_prime_state.json).
  [X] Chain-contrast control: the chain family (1,1,1,1) -- the proved
      full-dimensional class -- projected the same way, as the dichotomy
      control row.

No new cascade-depth measurements; no T=11; enumeration is confined to
committed cells and the committed p=29 targeted scope.
"""
import json
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples
from scaling_closure_verify import solutions_m5

T, M = 8, 5
OUT = {'primes': {}, 'meta': {}}


def canon(x, p):
    x %= p
    return min(x, p - x)


def renorm_orbit(S, p):
    """Renormalization orbit of difference set S (normalizer 1):
    {canon(d*delta^-1) : d in {1} u S, d != delta} over delta in {1} u S."""
    D = sorted(set([1] + list(S)))
    out = []
    for delta in D:
        inv = pow(delta, -1, p)
        S2 = tuple(sorted(canon(d * inv, p) for d in D if d != delta))
        out.append(S2)
    return out


def orbit_partition(sets, p):
    """Group sets into renormalization orbits; return list of member lists."""
    seen = {}
    orbits = []
    for s in sorted(sets):
        if s in seen:
            continue
        orb = sorted(set(renorm_orbit(s, p)))
        for m in orb:
            seen[m] = len(orbits)
        orbits.append(orb)
    return orbits


def ratio_data(S, p):
    """Resonance classification of D = {1} u S: canonical pairwise ratios,
    resonant-pair count (lambda_bar <= 4), min nontrivial ratio order."""
    D = sorted(set([1] + list(S)))
    lam = []
    for i in range(len(D)):
        for j in range(i + 1, len(D)):
            lam.append(canon(D[j] * pow(D[i], -1, p), p))
    lam = sorted(lam)
    res = sum(1 for x in lam if x <= 4)
    minord = None
    for i in range(len(D)):
        for j in range(len(D)):
            if i == j:
                continue
            r = (D[j] * pow(D[i], -1, p)) % p
            if r == 1:
                continue
            o = 1
            v = r
            while v != 1:
                v = (v * r) % p
                o += 1
            if minord is None or o < minord:
                minord = o
    return {'lam_spectrum': lam, 'resonant_pairs': res,
            'min_ratio_order': minord}


def probe_sets_p29(p):
    """The committed targeted classes at p=29: every set containing the
    binary-cascade core {2,4,8}, plus the doubling-path quadruples."""
    D = (p - 1) // 2
    sets = set()
    for x in range(3, D + 1):
        if x not in (2, 4, 8):
            sets.add(tuple(sorted((2, 4, 8, x))))
    path = [2, 4, 8, 13, 3, 6, 12, 5, 10, 9, 11, 7, 14]
    for i in range(len(path) - 3):
        sets.add(tuple(sorted(path[i:i + 4])))
    return sorted(sets)


def enumerate_set(p, fam, szts, label=''):
    """All size tuples for one ordered family: per-key n, projection sizes
    (raw starts), consecutive-relative marginals, box, density."""
    keys = [szt + fam for szt in szts]
    t0 = time.time()
    sols = solutions_m5(p, T, keys)
    rec = {}
    for key in sorted(sols):
        szt, lst = key[:5], sols[key]
        n = len(lst)
        proj = [len({a[i] for a in lst}) for i in range(4)]
        rel = [len({a[i] - a[i - 1] for a in lst}) for i in range(1, 4)]
        box = 1
        for q in proj:
            box *= q
        rec[str(szt)] = {'n': n, 'proj': proj, 'rel': rel,
                         'box': box, 'density': (n / box) if box else 0.0}
    return rec, time.time() - t0


def top_key(rec):
    """The size tuple with the largest n."""
    best, bn = None, -1
    for szt, r in rec.items():
        if r['n'] > bn:
            best, bn = szt, r['n']
    return best, bn


def main():
    # ---------------- Phase A + V' : p = 17, p = 19 (full) ----------------
    for p, stpath, stsec in [(17, '/home/z/my-project/scripts/h1_t8_state.json',
                              None),
                             (19, '/home/z/my-project/scripts/'
                              'h1_t8_prime_state.json', 'p19')]:
        k, (on, off), cls, eps = profile_sizes(p, T)
        szts = size_tuples(M, on, off, p)
        D = (p - 1) // 2
        allsets = [tuple(sorted(c)) for c in
                   _combinations(range(2, D + 1), 4)]
        orbits = orbit_partition(allsets, p)
        st = json.load(open(stpath))
        done = st['done'] if stsec is None else st[stsec]['done']

        print('=== p=%d (k=%d, eps=%d, %s, sizes %s/%s): %d sets, '
              '%d orbits, %d size tuples ==='
              % (p, k, eps, cls, on, off, len(allsets), len(orbits),
                 len(szts)))

        gt_hit = gt_miss = gt_skip = 0
        env_adm = 0            # pinning envelope over admissible keys
        env_key = None
        dens = []
        orbit_recs = []
        vprime_ok = True
        for oi, orb in enumerate(orbits):
            prof_counters = []
            tot = 0
            for s in orb:
                rec, dt = enumerate_set(p, s, szts)
                tot += sum(r['n'] for r in rec.values())
                # GT vs committed census (sorted-order family key)
                fk = str(list(s))
                if fk in done:
                    for szt, r in rec.items():
                        ck = str(tuple(ast_literal(szt)))
                        if ck in done[fk]:
                            if done[fk][ck] == r['n']:
                                gt_hit += 1
                            else:
                                gt_miss += 1
                        else:
                            gt_skip += 1
                else:
                    gt_skip += len(rec)
                # envelope + density over admissible keys
                for szt, r in rec.items():
                    if r['n'] > 0:
                        mp = max(r['proj'])
                        if mp > env_adm:
                            env_adm, env_key = mp, (s, szt)
                        if r['box']:
                            dens.append(r['density'])
                prof_counters.append(Counter(
                    tuple(sorted(r['proj'])) for r in rec.values()))
            # Lemma V': identical profile multiset across orbit members
            same = all(pc == prof_counters[0] for pc in prof_counters[1:])
            vprime_ok &= same
            rep = orb[0]
            orbit_recs.append({
                'rep': list(rep), 'members': [list(m) for m in orb],
                'total_volume': tot, 'vprime_profile_multiset_equal': same,
                'resonance': ratio_data(rep, p)})
        pos = [o for o in orbit_recs if o['total_volume'] > 0]
        print('  positive orbits: %d / %d; Lemma V\' profile multiset '
              'equal across members: %s' % (len(pos), len(orbit_recs),
                                             vprime_ok))
        print('  GT vs committed T-20 census: %d hit, %d miss, %d skip'
              % (gt_hit, gt_miss, gt_skip))
        print('  PINNING ENVELOPE (max |Pi_i| over admissible keys): %d at %s'
              % (env_adm, env_key if env_key else 'n/a'))
        print('  joint density n/box over admissible keys: median %.4f, '
              'max %.4f' % (sorted(dens)[len(dens) // 2], max(dens)))
        print('  top orbit: rep %s vol %d res_pairs %d lam %s'
              % (pos[0]['rep'], pos[0]['total_volume'],
                 pos[0]['resonance']['resonant_pairs'],
                 pos[0]['resonance']['lam_spectrum'][:8]))
        OUT['primes']['p%d' % p] = {
            'k': k, 'eps': eps, 'cls': cls, 'sizes': [on, off],
            'n_sets': len(allsets), 'n_orbits': len(orbits),
            'n_positive_orbits': len(pos), 'vprime_ok': vprime_ok,
            'gt_hit': gt_hit, 'gt_miss': gt_miss, 'gt_skip': gt_skip,
            'pinning_envelope': env_adm, 'envelope_key': [list(env_key[0]),
                                                          env_key[1]],
            'density_median': sorted(dens)[len(dens) // 2],
            'density_max': max(dens),
            'orbits': orbit_recs}

        # top-orbit marginal sets at the top size tuple (mechanism record)
        rec, _ = enumerate_set(p, tuple(pos[0]['rep']), szts)
        tk, bn = top_key(rec)
        fam = tuple(pos[0]['rep'])
        key0 = ast_literal(tk) + fam
        sols = solutions_m5(p, T, [key0])[key0]
        OUT['primes']['p%d' % p]['top_orbit_detail'] = {
            'rep': list(fam), 'top_szt': ast_literal(tk), 'n': bn,
            'proj': rec[tk]['proj'], 'rel': rec[tk]['rel'],
            'marginals': [sorted({a[i] for a in sols}) for i in range(4)]}

    # ---------------- Phase B : p = 29 (targeted, committed scope) --------
    p = 29
    k, (on, off), cls, eps = profile_sizes(p, T)
    szts = size_tuples(M, on, off, p)
    probe = probe_sets_p29(p)
    st = json.load(open('/home/z/my-project/scripts/h1_t8_prime_state.json'))
    done = st['p29']['done']
    print()
    print('=== p=%d (k=%d, eps=%d, %s, sizes %s/%s): %d targeted sets '
          '(committed classes: {2,4,8,x} + doubling path) ==='
          % (p, k, eps, cls, on, off, len(probe)))
    gt_hit = gt_miss = gt_skip = 0
    env_adm = 0
    env_key = None
    dens = []
    set_recs = []
    for s in probe:
        rec, dt = enumerate_set(p, s, szts)
        tot = sum(r['n'] for r in rec.values())
        fk = str(list(s))
        if fk in done:
            for szt, r in rec.items():
                ck = str(tuple(ast_literal(szt)))
                if ck in done[fk]:
                    if done[fk][ck] == r['n']:
                        gt_hit += 1
                    else:
                        gt_miss += 1
                else:
                    gt_skip += 1
        else:
            gt_skip += len(rec)
        for szt, r in rec.items():
            if r['n'] > 0:
                mp = max(r['proj'])
                if mp > env_adm:
                    env_adm, env_key = mp, (s, szt)
                if r['box']:
                    dens.append(r['density'])
        set_recs.append({'set': list(s), 'total_volume': tot,
                         'cascade': set(s) >= {2, 4, 8},
                         'resonance': ratio_data(s, p)})
    adm = [r for r in set_recs if r['total_volume'] > 0]
    print('  admissible among targeted: %d / %d' % (len(adm), len(probe)))
    print('  GT vs committed T-20 targeted census: %d hit, %d miss, %d skip'
          % (gt_hit, gt_miss, gt_skip))
    print('  PINNING ENVELOPE (max |Pi_i| over admissible keys): %d at %s'
          % (env_adm, env_key if env_key else 'n/a'))
    if dens:
        print('  joint density n/box: median %.4f, max %.4f'
              % (sorted(dens)[len(dens) // 2], max(dens)))
    print('  top targeted set: %s vol %d' %
          (max(adm, key=lambda r: r['total_volume'])['set'] if adm else None,
           max((r['total_volume'] for r in adm), default=0)))
    OUT['primes']['p29'] = {
        'k': k, 'eps': eps, 'cls': cls, 'sizes': [on, off],
        'n_probe_sets': len(probe), 'n_admissible': len(adm),
        'gt_hit': gt_hit, 'gt_miss': gt_miss, 'gt_skip': gt_skip,
        'pinning_envelope': env_adm,
        'envelope_key': [list(env_key[0]), env_key[1]] if env_key else None,
        'density_median': sorted(dens)[len(dens) // 2] if dens else None,
        'density_max': max(dens) if dens else None,
        'sets': set_recs}

    # ---------------- Phase X : the chain-contrast control ----------------
    print()
    print('=== chain-contrast control: family (1,1,1,1), top size tuple ===')
    chain = []
    for p in (17, 19, 29):
        k, (on, off), cls, eps = profile_sizes(p, T)
        szts = size_tuples(M, on, off, p)
        rec, _ = enumerate_set(p, (1, 1, 1, 1), szts)
        tk, bn = top_key(rec)
        ov = sum(ast_literal(tk)) - p
        chain.append({'p': p, 'k': k, 'top_szt': ast_literal(tk), 'ov': ov,
                      'n': bn, 'proj': rec[tk]['proj'],
                      'rel': rec[tk]['rel']})
        print('  p=%d k=%d szt=%s ov=%d: n=%d, projections %s, rel %s'
              % (p, k, tk, ov, bn, rec[tk]['proj'], rec[tk]['rel']))
    OUT['chain_control'] = chain

    with open('/home/z/my-project/scripts/out_strand_rigidity.json',
              'w') as f:
        json.dump(OUT, f, indent=1, default=int)
    print()
    print('wrote scripts/out_strand_rigidity.json')


def ast_literal(s):
    import ast
    return ast.literal_eval(s)


def _combinations(items, r):
    """itertools.combinations without the import noise in main path."""
    from itertools import combinations
    return combinations(items, r)


if __name__ == '__main__':
    main()
