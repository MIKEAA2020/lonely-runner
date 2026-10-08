#!/usr/bin/env python3
"""
h1_t8_structure.py -- the H1-general attack at T=8, Stage 2 (T-20).

Probes the STRUCTURE of the distinct-family solution sets at the committed
zoo cell (T=8 (1K,5U), p=17):

  [V] renormalization-invariance verification: the 7 volume classes =
      the orbits of S -> {canon(d*delta^-1)} over delta in {1} U S
      (the covering problem is dilation-invariant; volume summed over
      size assignments is orbit-invariant).  Verified set-by-set.
  [E] essential-arc patterns per solution (which arcs are needed).
  [P] pairwise overlap matrices O_ij; the budget identity
      sum_{i<j} O_ij = ov + sum_x C(o(x),2) checked exactly.
  [S] the pair-overlap spacing caps: max_a O_ij vs s_j/b_min+1 form.
  [J] coordinate projection sizes.
  [G] the junction geometry: block decomposition of each arc in the
      natural coordinate; the interleaving pattern of the top class.

GT: solutions_m5 lists re-derived and cross-checked against the Stage-1
counts (independent counting path).
"""
import json
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples, ap_base, rot
from scaling_closure_verify import solutions_m5

P, T, M = 17, 8, 5


def canon(x, p):
    x %= p
    return min(x, p - x)


def renorm_orbit(S, p):
    """The renormalization orbit of difference set S (with normalizer 1):
    all {canon(d*delta^-1) : d in {1} U S, d != delta}, delta in {1} U S."""
    D = sorted(set([1] + list(S)))
    out = []
    for delta in D:
        inv = pow(delta, -1, p)
        S2 = tuple(sorted(canon(d * inv, p) for d in D if d != delta))
        out.append(S2)
    return out


def arcs_of(a, sz, fam, p):
    """The 5 arc point-sets: normalizer interval + 4 APs."""
    m = len(sz)
    ds = [1] + list(fam)
    starts = [0] + list(a)
    return [set((starts[i] + t * ds[i]) % p for t in range(sz[i]))
            for i in range(m)]


def blocks(A, p):
    """Maximal consecutive-integer runs of a point-set on Z_p."""
    if not A:
        return []
    S = sorted(A)
    runs = [[S[0]]]
    for x in S[1:]:
        if x == runs[-1][-1] + 1:
            runs[-1].append(x)
        else:
            runs.append([x])
    if len(runs) > 1 and runs[0][0] == 0 and runs[-1][-1] == p - 1:
        runs[-1].extend(runs.pop(0))
    return runs


def bmin(lam, L, p):
    for b in range(1, p):
        if canon(lam * b, p) <= L:
            return b
    return p


def main():
    p = P
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(M, on, off, p)
    with open('/home/z/my-project/scripts/h1_t8_state.json') as f:
        st = json.load(f)
    counts = st['done']
    import ast
    tot_by_set = defaultdict(list)
    for kk, per in counts.items():
        fam = tuple(ast.literal_eval(kk))
        tot_by_set[tuple(sorted(fam))].append(sum(per.values()))

    # ---------------- [V] renormalization orbits ----------------
    print('=== [V] renormalization orbits vs volume classes ===')
    vol_class = defaultdict(set)
    for s, vols in tot_by_set.items():
        vol_class[vols[0]].add(s)
    ok_all = True
    for vol, sets in sorted(vol_class.items(), key=lambda kv: -kv[0]):
        reps = set()
        seen = set()
        orbits = []
        for s in sorted(sets):
            if s in seen:
                continue
            orb = set(renorm_orbit(s, p))
            # orbit elements are 4-sets; keep those that are distinct-4
            # sets in [2, 8] (all are, by construction, unless collision)
            seen |= orb
            orbits.append(orb)
        same = (len(orbits) == 1 and sets == orbits[0])
        ok_all &= same
        print('  vol %6d: %d sets, %d orbit(s), orbit==class: %s  %s'
              % (vol, len(sets), len(orbits), same,
                 sorted(orbits[0]) if same else ''))
    print('  ALL classes are single renormalization orbits: %s' % ok_all)

    # ---------------- solution structure of the top class ----------------
    print()
    print('=== structure probes: top family (2,4,6,8) ===')
    fams_probe = [(2, 4, 6, 8), (2, 3, 4, 8), (2, 3, 5, 6), (2, 3, 4, 7)]
    labels = ['top-1 (binary cascade + spectator 6)',
              'top-1 (spectator 3)', 'min-class', 'mid-class']
    OUT = {'orbits_ok': ok_all, 'probes': {}}
    for fam, lab in zip(fams_probe, labels):
        keys = [szt + fam for szt in szts]
        t0 = time.time()
        sols = solutions_m5(p, T, keys)
        # GT vs stage-1 counts
        mism = 0
        for key, lst in sols.items():
            if len(lst) != counts[str(list(fam))][str(key[:5])]:
                mism += 1
        tot = sum(len(v) for v in sols.values())
        print('  --- %s %s: %d solutions / %d size tuples '
              '(%.1fs, GT mismatches vs Stage-1: %d) ---'
              % (fam, lab, tot, len(sols), time.time() - t0, mism))
        # concentrate on the top size tuple
        best = max(sols.items(), key=lambda kv: len(kv[1]))
        szt, lst = best
        sz = szt[:5]
        ov = sum(sz) - p
        print('      top size tuple %s (ov=%d): %d solutions'
              % (sz, ov, len(lst)))
        # [E] essential arcs
        pat_counter = Counter()
        for a in lst:
            As = arcs_of(a, sz, fam, p)
            ess = []
            for i in range(5):
                rest = set()
                for j in range(5):
                    if j != i:
                        rest |= As[j]
                ess.append(0 if len(rest) == p else 1)
            pat_counter[tuple(ess)] += 1
        print('      essential-arc patterns (1=needed): %s'
              % dict(pat_counter.most_common(8)))
        # [P] pairwise overlaps + budget identity
        bad_budget = 0
        O_sum = Counter()
        O_max = defaultdict(int)
        for a in lst:
            As = arcs_of(a, sz, fam, p)
            w = Counter()
            for A in As:
                for x in A:
                    w[x] += 1
            # budget identity: sum_{i<j} O_ij = ov + sum_x C(o(x),2)
            lhs = 0
            for i in range(5):
                for j in range(i + 1, 5):
                    Oij = len(As[i] & As[j])
                    lhs += Oij
                    O_sum[(i, j)] += Oij
                    O_max[(i, j)] = max(O_max[(i, j)], Oij)
            rhs = ov + sum(v * (v - 1) // 2 for v in
                           (w[x] - 1 for x in w))
            if lhs != rhs:
                bad_budget += 1
        print('      budget identity violations: %d / %d' %
              (bad_budget, len(lst)))
        print('      mean O_ij: %s' %
              {ij: round(O_sum[ij] / float(len(lst)), 2) for ij in O_sum})
        print('      max  O_ij: %s' % dict(O_max))
        # [S] spacing caps for each pair
        ds = [1] + list(fam)
        caps = {}
        for i in range(5):
            for j in range(5):
                if i == j:
                    continue
                lam = (ds[j] * pow(ds[i], -1, p)) % p
                caps[(i, j)] = (sz[j] // max(1, bmin(lam, sz[i], p)) + 1
                                if sz[j] <= sz[i] else
                                sz[j] // max(1, bmin(lam, sz[i], p)) + 1)
        cap_ok = all(O_max[(i, j)] <= caps[(i, j)]
                     for (i, j) in O_max if (i, j) in caps)
        print('      pair caps (i<j, oriented d_j/d_i): %s  all max<=cap: %s'
              % ({ij: caps.get(ij) for ij in sorted(O_max)}, cap_ok))
        # [J] projections
        proj = [len({a[i] for a in lst}) for i in range(4)]
        print('      coordinate projections |Pi_i|: %s' % proj)
        # [G] block structure of a few sample solutions
        print('      sample solutions (block decomposition):')
        for a in lst[:3]:
            As = arcs_of(a, sz, fam, p)
            print('        a=%s' % (a,))
            for i, A in enumerate(As):
                print('          arc %d (d=%d,s=%d): blocks %s'
                      % (i, ds[i], sz[i],
                         [(r[0], len(r)) for r in blocks(A, p)]))
        OUT['probes'][str(list(fam))] = {
            'label': lab, 'top_size': list(szt), 'ov': ov,
            'n_sols': len(lst), 'essential': {str(k_): v for k_, v in
                                              pat_counter.items()},
            'mean_O': {str(ij): O_sum[ij] / float(len(lst))
                       for ij in O_sum},
            'max_O': {str(ij): O_max[ij] for ij in O_max},
            'projections': proj, 'budget_bad': bad_budget,
            'total_gt_mismatch': mism,
        }
    with open('/home/z/my-project/scripts/out_h1_t8_structure.json', 'w') as f:
        json.dump(OUT, f, indent=1, default=int)
    print()
    print('wrote scripts/out_h1_t8_structure.json')


if __name__ == '__main__':
    main()
