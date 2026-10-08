#!/usr/bin/env python3
"""
Input E boundary census — the exact foot-count thresholds.

For each prime p (both classes), each critical s1 (I1 = [s1, p-1]), each
critical size s, each +-rep d:  min over a of feet(A(a,d,s)) where
feet = s - |A ∩ I1|.  Viable(f; s1,s) = {d : min_a feet <= f}.

Proposition E (extended trichotomy):  clean(f) ⟺ for ALL critical (s1,s):
Viable(f; s1,s) ⊆ {2, (p-1)/2}.

Questions answered (the reviewer's ask):
  Q1  one-foot (f=1) boundary at p ≡ 1 mod 6, all pairs
  Q2  two-foot (f=2) boundary at p ≡ 1 mod 6, all pairs   (note: 43)
  Q3  one-foot boundary at p ≡ 5 mod 6, all pairs         (note: 23)
  Q4  the CONSUMED pairs of Lemma A M4/M5 at p ≡ 1 mod 6:
      (s1,s) = (2k+2, 2k+1) [M4] and (2k+1, 2k+1) [M5], one-foot.
      -> decides whether Lemma A extends below p >= 43.
  Q5  Lemma B coverage at k=2 (p=17): the note's §2 proof consumes only
      the CONTAINED census (Lemma T) + T1/T2/T3; check the contained forms
      at (s1,s)=(6,6) are exactly {E,T} and |E∩T| >= 2.
"""
import json
import sys

import numpy as np

OUT = {}


def primes_6k(n):
    sieve = np.ones(n + 1, bool)
    sieve[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return [int(x) for x in np.nonzero(sieve)[0] if x % 6 in (1, 5)]


def critical_sizes(k, cls):
    return (2 * k + 1, 2 * k + 2) if cls == 5 else (2 * k, 2 * k + 1, 2 * k + 2)


def foot_table(p):
    """feetmin[(s1,s)][d] = min over a of feet; d over +-reps 1..(p-1)/2."""
    cls = p % 6
    k = (p - cls) // 6
    sizes = critical_sizes(k, cls)
    q = (p - 1) // 2
    a_grid = np.arange(p)[:, None]
    table = {}
    for s1 in sizes:
        for s in sizes:
            best = {}
            for d in range(1, q + 1):
                js = (d * np.arange(s)) % p
                arr = (a_grid + js[None, :]) % p
                tmax = int((arr >= s1).sum(axis=1).max())
                best[d] = s - tmax
            table[(s1, s)] = best
    return k, sizes, q, table


def viable_set(best, f):
    # +-distinct normalization: d1 = 1, so the OTHER APs have ||d|| >= 2;
    # Proposition E's viable set is over d in [2, q].
    return sorted(d for d, m in best.items() if m <= f and d >= 2)


def main():
    res = {}
    print('=== Input E boundary census ===')
    for p in primes_6k(419):
        if p < 7:
            continue
        k, sizes, q, table = foot_table(p)
        cls = p % 6
        row = {'k': k, 'cls': cls}
        # global cleanliness per f
        for f in (0, 1, 2):
            bad = {}
            for (s1, s), best in table.items():
                v = viable_set(best, f)
                if not set(v) <= {2, q}:
                    bad['%d_%d' % (s1, s)] = v
            row['clean_f%d' % f] = (len(bad) == 0)
            if bad:
                row['dirty_f%d' % f] = bad
        # consumed pairs (Lemma A M4/M5), one-foot
        if cls == 1 and k >= 1:
            cons = {}
            for (s1, s) in ((2 * k + 2, 2 * k + 1), (2 * k + 1, 2 * k + 1)):
                if (s1, s) in table:
                    v = viable_set(table[(s1, s)], 1)
                    cons['%d_%d' % (s1, s)] = {'viable': v,
                                               'clean': set(v) <= {2, q}}
            row['consumed_onefoot'] = cons
        res[p] = row
        flags = ''.join('F%d:%s' % (f, 'c' if row['clean_f%d' % f] else 'D')
                        for f in (0, 1, 2))
        extra = ''
        if 'consumed_onefoot' in row:
            extra = ' consumed:' + ','.join(
                '%s%s' % (kk, '+' if vv['clean'] else '-')
                for kk, vv in row['consumed_onefoot'].items())
        print('p=%3d (k=%2d, %d mod 6): %s%s' % (p, k, cls, flags, extra))
        sys.stdout.flush()

    # boundaries
    bnd = {}
    for cls in (1, 5):
        for f in (0, 1, 2):
            dirty = [p for p, r in res.items()
                     if r['cls'] == cls and not r['clean_f%d' % f]]
            clean = [p for p, r in res.items()
                     if r['cls'] == cls and r['clean_f%d' % f]]
            bnd['cls%d_f%d' % (cls, f)] = {
                'dirty_primes': dirty,
                'clean_from': min(clean) if clean else None}
    OUT['per_prime'] = {str(p): r for p, r in res.items()}
    OUT['boundaries'] = bnd

    # Q5: Lemma B at k=2 (p=17): contained size-(2k+2) forms at s1=2k+2
    p = 17
    k = 2
    s1 = s = 6
    forms = []
    for d in range(1, p):
        for a in range(p):
            A = {(a + d * j) % p for j in range(s)}
            if all(s1 <= x <= p - 1 for x in A):
                forms.append((a, d, tuple(sorted(A))))
    # dedupe by set
    sets = {}
    for a, d, A in forms:
        sets.setdefault(A, []).append((a, d))
    E = tuple(x for x in range(6, 17) if x % 2 == 0)
    O = tuple(x for x in range(6, 17) if x % 2 == 1)
    # T form: d = (p-1)/2 = 8: blocks of 3 with right endpoints differing 8
    Tb = tuple(sorted(list(range(6, 9)) + list(range(14, 17))))
    inter = len(set(E) & set(Tb))
    OUT['lemmaB_p17'] = {
        'contained_forms_at_(6,6)': {str(A): v for A, v in sets.items()},
        'n_forms': len(sets),
        'E_in': E in sets, 'T_in': Tb in sets,
        'E_T_intersection': inter,
        'O_two_consecutive': any(
            (y - x == 1) for x, y in zip(O, O[1:]))}
    print('LemmaB p=17: contained (6,6)-forms: %d sets; E∈%s T∈%s '
          '|E∩T|=%d' % (len(sets), E in sets, Tb in sets, inter))

    with open('out_foot_census.json', 'w') as f:
        json.dump(OUT, f, indent=1, default=str)
    print('written: out_foot_census.json')
    print('\nBOUNDARIES:')
    for kk, vv in bnd.items():
        print('  %s: dirty=%s clean_from=%s' %
              (kk, vv['dirty_primes'], vv['clean_from']))


if __name__ == '__main__':
    main()
