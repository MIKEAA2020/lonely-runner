#!/usr/bin/env python3
"""
t26_multichar_rigidity_session.py -- T-26: the multi-character rigidity
DEVELOPED, as the arrangement-phase decomposition of the joint survivor
set (the object the T-22 calibration named: the sparsity is joint, not
marginal). Run as a proof attempt under the pre-committed criterion.

CHARTER (fixed before the work):
  [boundary vs Lemma M]  No transfer operators, no spectral quantities.
      T-25 proved the margin external to the lossless dynamics; the only
      in-program move left is a static decomposition of the per-fiber
      solution sets themselves. Lemma M is untouched and inherited.
  [compliance]           Not a fourth marginal/projection variation: no
      single-character law is sought. The object is the joint split
      X = arrangement (the 6 pairwise-difference characters, dim m-2)
      ⊕ phase (the diagonal direction). Marginals are computed only as
      GT anchors against the T-22 record.
  [pre-committed criterion]  The session succeeds iff it proves a
      uniform-in-k margin of |Sol| against the trivial box p^(m-1) --
      or lands the decomposition with the margin localized to one
      named factor (the honest expected outcome under the T-25
      mechanism). Outcome tree: (a) margin proved -> H1g discharged;
      (b) decomposition proved + margin localized -> the obligation
      compressed to its finest grain; (c) no law at all -> inconclusive,
      session not extended.
  [scope]                Committed primes/families only; no T=11; no new
      cells; targeted verification on committed objects (rule iii).

Session content:
  [0] Committed objects rebuilt + GT (committed top counts, family
      totals, covers() sample, T-22 marginal profiles).
  [1] LEMMA D (the phase-packing bound), proved and verified
      class-by-class: every arrangement class satisfies
          |C_u| <= s1 - span(H_u) + 1 <= s1 - h_u + 1 <= ov + 1,
      where H_u is the hole set of the moving union (phase-invariant).
  [2] The decomposition measured: N_arr (realized arrangement classes),
      class-size histograms, the arrangement/phase factorization of the
      committed margins 95.3x / 1,206.7x / 3,802.6x.
  [3] The arrangement-space ladder: p^3 -> |Arr_adm| (moving-moving
      overlap-admissible arrangements, the t-invariant constraints)
      -> N_arr (realized). The chain control (full arrangement space)
      vs the distinct families.
  [4] The orbit test (p=19, S=(2,4,5,8) -> S2=(3,4,6,8), delta=5):
      how the arrangement/phase split transforms under the
      renormalization (the C10 anchor-mixing).
  [5] The phase structure within classes (the run / pairing law).

Output: scripts/out_t26_session.json + stdout log.
"""
import json
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples, ap_base, rot
from scaling_closure_verify import solutions_m5

OUT = {}
T, M = 8, 5

FAMS = {
    17: [('chain', (1, 1, 1, 1)), ('top-distinct', (2, 3, 4, 8))],
    19: [('chain', (1, 1, 1, 1)), ('top-distinct', (2, 4, 5, 8))],
    29: [('committed-A', (2, 4, 7, 8)), ('committed-B', (3, 4, 8, 13)),
         ('chain', (1, 1, 1, 1))],
}

# committed GT anchors (T-22 / T-25 records)
GT_TOP = {  # top size-tuple counts
    (17, (2, 3, 4, 8)): 876, (19, (2, 4, 5, 8)): 108,
    (29, (2, 4, 7, 8)): 186, (29, (3, 4, 8, 13)): 186,
    (17, (1, 1, 1, 1)): 8880, (19, (1, 1, 1, 1)): 4680,
    (29, (1, 1, 1, 1)): 29760,
}
GT_TOTAL = {  # whole-family totals over all size tuples
    (17, (2, 3, 4, 8)): 9633, (19, (2, 4, 5, 8)): 801,
    (29, (2, 4, 7, 8)): 1311, (29, (3, 4, 8, 13)): 1311,
}
GT_MARG = {  # T-22 top-key 10-functional profiles: singletons / pairs
    (17, (2, 3, 4, 8)): ([12, 17, 17, 16], [17, 10, 17, 15, 14, 12]),
    (19, (2, 4, 5, 8)): ([4, 4, 9, 8], [4, 8, 11, 8, 12, 16]),
    (29, (2, 4, 7, 8)): ([9, 7, 10, 27], [7, 6, 17, 6, 15, 14]),
}


def covers(p, szt, fam, a):
    """Direct covering test (independent of solutions_m5)."""
    pts = set(range(szt[0]))
    for i, d in enumerate(fam):
        pts |= {(a[i] + t * d) % p for t in range(szt[i + 1])}
    return len(pts) == p


def build_sizedict(p, fam):
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(M, on, off, p)
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    sd = defaultdict(set)
    for key, lst in sols.items():
        sd[key[:5]].update(lst)
    return dict(sd), len(szts)


def moving_masks(p, szt, fam, a):
    full = (1 << p) - 1
    return [rot(ap_base(fam[i], szt[i + 1], p), a[i], p, full)
            for i in range(4)]


def arr_key(p, a):
    """Arrangement key: the coset of the diagonal (1,1,1,1)."""
    return ((a[0] - a[1]) % p, (a[0] - a[2]) % p, (a[0] - a[3]) % p)


def decompose(p, sols):
    classes = defaultdict(list)
    for a in sols:
        classes[arr_key(p, a)].append(a)
    for k in classes:
        classes[k].sort()
    return classes


def lemma_D_check(p, szt, fam, classes):
    """Class-by-class verification of |C_u| <= s1 - span(H_u) + 1
    <= s1 - h_u + 1 <= ov + 1 (H_u = holes of the moving union)."""
    s1 = szt[0]
    ov = sum(szt) - p
    rec = {'classes': len(classes), 'violations': 0,
           'sharp_equal': 0, 'max_class': 0, 'max_ov_bound': ov + 1,
           'h_min': None, 'h_max': None, 'span_min': None,
           'zero_h_classes': 0, 'zero_h_max': 0,
           'worst': None}
    for key, lst in classes.items():
        a0 = lst[0]
        U = 0
        for m in moving_masks(p, szt, fam, a0):
            U |= m
        holes = [i for i in range(p) if not (U >> i) & 1]
        h = len(holes)
        if h == 0:
            # the moving union alone covers: the covering bound is
            # vacuous here; record separately (bounded only by p).
            rec['zero_h_classes'] += 1
            rec['zero_h_max'] = max(rec['zero_h_max'], len(lst))
            rec['max_class'] = max(rec['max_class'], len(lst))
            rec['h_min'] = 0 if rec['h_min'] is None else min(rec['h_min'], 0)
            continue
        if h == 1:
            span = 0
        else:
            span = p - max((holes[(i + 1) % h] - holes[i]) % p
                           for i in range(h))
        b1 = s1 - span + 1
        b2 = s1 - h + 1
        if len(lst) > b1 or len(lst) > b2 or len(lst) > ov + 1:
            rec['violations'] += 1
            if rec['worst'] is None:
                rec['worst'] = {'key': list(key), 'size': len(lst),
                                'span_bound': b1, 'h_bound': b2}
        if len(lst) == b1:
            rec['sharp_equal'] += 1
        rec['max_class'] = max(rec['max_class'], len(lst))
        rec['h_min'] = h if rec['h_min'] is None else min(rec['h_min'], h)
        rec['h_max'] = h if rec['h_max'] is None else max(rec['h_max'], h)
        rec['span_min'] = span if rec['span_min'] is None \
            else min(rec['span_min'], span)
    return rec


def arr_admissible(p, szt, fam):
    """|Arr_adm|: gauge a2=0; count (x3,x4,x5) with all 6 moving-moving
    overlaps <= ov (the t-invariant constraints alone)."""
    full = (1 << p) - 1
    ov = sum(szt) - p
    m0 = ap_base(fam[0], szt[1], p)
    base3 = ap_base(fam[1], szt[2], p)
    base4 = ap_base(fam[2], szt[3], p)
    base5 = ap_base(fam[3], szt[4], p)
    cnt = 0
    for x3 in range(p):
        m3 = rot(base3, x3, p, full)
        if (m0 & m3).bit_count() > ov:
            continue
        for x4 in range(p):
            m4 = rot(base4, x4, p, full)
            if (m0 & m4).bit_count() > ov or (m3 & m4).bit_count() > ov:
                continue
            for x5 in range(p):
                m5 = rot(base5, x5, p, full)
                if (m0 & m5).bit_count() > ov:
                    continue
                if (m3 & m5).bit_count() > ov or (m4 & m5).bit_count() > ov:
                    continue
                cnt += 1
    return cnt


def phase_structure(p, classes):
    """Within-class phase sets: size histogram, gap histogram,
    pair-gap examples."""
    hist = Counter(len(lst) for lst in classes.values())
    gaps = Counter()
    examples = []
    for key, lst in classes.items():
        if len(lst) >= 2:
            base = lst[0][0]
            ts = sorted((a[0] - base) % p for a in lst)
            for i in range(len(ts) - 1):
                gaps[ts[i + 1] - ts[i]] += 1
            if len(examples) < 6:
                examples.append({'phases': ts,
                                 'arr': list(key)})
    return {'size_hist': {str(k): v for k, v in sorted(hist.items())},
            'gap_hist': {str(k): v for k, v in sorted(gaps.items())},
            'examples': examples}


def analyze_key(p, fam, szt, sols, deep=True):
    """Full decomposition analysis of one (family, size-tuple) key."""
    n = len(sols)
    ov = sum(szt) - p
    classes = decompose(p, sols)
    n_arr = len(classes)
    sizes = [len(v) for v in classes.values()]
    avg = n / n_arr if n_arr else 0.0
    rec = {
        'p': p, 'fam': list(fam), 'szt': list(szt), 'ov': ov, 'n_sol': n,
        'trivial_box': p ** 4, 'margin_over_trivial': round(p ** 4 / n, 1)
        if n else None,
        'n_arr': n_arr, 'arr_box': p ** 3,
        'arr_compression': round(p ** 3 / n_arr, 2) if n_arr else None,
        'avg_class': round(avg, 3), 'max_class': max(sizes) if sizes else 0,
        'phase_factor': round(p / avg, 2) if avg else None,
        'c1_ov_units': round(n / ov ** 4, 4) if ov else None,
        'narr_over_ov3': round(n_arr / ov ** 3, 4) if ov else None,
    }
    if deep and n:
        rec['lemma_D'] = lemma_D_check(p, szt, fam, classes)
        rec['phase_structure'] = phase_structure(p, classes)
        t0 = time.time()
        rec['arr_adm'] = arr_admissible(p, szt, fam)
        rec['arr_adm_sec'] = round(time.time() - t0, 1)
        rec['adm_to_realized'] = round(rec['arr_adm'] / n_arr, 2)
        # T-22 marginal GT anchors
        sing = [len({a[i] for a in sols}) for i in range(4)]
        pairs = [len({(a[j] - a[i]) % p for a in sols})
                 for i in range(4) for j in range(i + 1, 4)]
        rec['marginals'] = {'singletons': sing, 'pairs': pairs}
        box = 1
        for m in sing + pairs:
            box *= m
        rec['x_box_density'] = round(n / box, 5)
    return rec


def main():
    t_start = time.time()
    print('=== T-26: the arrangement-phase decomposition session ===')
    print('Charter: pre-committed criterion (uniform-in-k margin, or the')
    print('decomposition with the margin localized to one named factor).')
    print()
    OUT['gt'] = {}

    # build every committed family ONCE (solutions_m5 is the cost)
    STORED = {}
    for p in (17, 19, 29):
        for label, fam in FAMS[p]:
            t0 = time.time()
            sd, n_szts = build_sizedict(p, fam)
            STORED[(p, fam)] = (sd, n_szts)
            print('  built p=%d %-12s fam=%s: %d size tuples, %.1fs'
                  % (p, label, fam, n_szts, time.time() - t0))
            sys.stdout.flush()

    # ---------- [0] committed objects + GT ----------
    print()
    print('--- [0] committed objects rebuilt + GT ---')
    for p in (17, 19, 29):
        for label, fam in FAMS[p]:
            sd, n_szts = STORED[(p, fam)]
            total = sum(len(v) for v in sd.values())
            top_szt = max(sd, key=lambda s: len(sd[s])) if sd else None
            top_n = len(sd[top_szt]) if top_szt else 0
            gtt = GT_TOTAL.get((p, fam))
            gtp = GT_TOP.get((p, fam))
            ok_t = (gtt is None or total == gtt)
            ok_p = (gtp is None or top_n == gtp)
            OUT['gt']['%d_%s' % (p, label)] = {
                'fam': list(fam), 'n_szts': n_szts, 'total': total,
                'gt_total': gtt, 'total_ok': ok_t,
                'top_szt': list(top_szt) if top_szt else None,
                'top_n': top_n, 'gt_top': gtp, 'top_ok': ok_p}
            print('  p=%d %-12s fam=%s  total=%s (GT %s: %s)  top=%s @%s '
                  '(GT %s: %s)'
                  % (p, label, fam, total, gtt, 'OK' if ok_t else 'MISMATCH',
                     top_n, top_szt, gtp, 'OK' if ok_p else 'MISMATCH'))
            # covers() GT on a sample (first/last/mid solutions, top key)
            if top_szt and top_n:
                lst = sorted(sd[top_szt])
                sample = ([lst[0], lst[-1], lst[top_n // 2]]
                          if top_n >= 3 else lst)
                cov_ok = all(covers(p, top_szt, fam, a) for a in sample)
                OUT['gt']['%d_%s' % (p, label)]['covers_gt'] = cov_ok
                if not cov_ok:
                    print('    covers() GT FAILED for', sample)

    # ---------- [2] the decomposition on the committed top keys ----------
    print()
    print('--- [2] the decomposition: top keys of the committed families ---')
    OUT['top_keys'] = {}
    for p in (17, 19, 29):
        for label, fam in FAMS[p]:
            sd, _ = STORED[(p, fam)]
            if not sd:
                continue
            top_szt = max(sd, key=lambda s: len(sd[s]))
            rec = analyze_key(p, fam, top_szt, sorted(sd[top_szt]))
            OUT['top_keys']['%d_%s' % (p, label)] = rec
            print('  p=%d %-12s top szt=%s  |Sol|=%d  margin=%.1fx'
                  % (p, label, top_szt, rec['n_sol'],
                     rec['margin_over_trivial']))
            print('      N_arr=%d (p^3/N_arr=%.1fx)  avg class=%.3f  '
                  'max class=%d  phase factor=%.2f'
                  % (rec['n_arr'], rec['arr_compression'], rec['avg_class'],
                     rec['max_class'], rec['phase_factor']))
            print('      ov=%d  c1(ov-units)=%s  N_arr/ov^3=%s  '
                  'Lemma D: viol=%d sharp-eq=%d  h in [%s, %s]  '
                  'zero-h classes=%d (max %d)'
                  % (rec['ov'], rec['c1_ov_units'], rec['narr_over_ov3'],
                     rec['lemma_D']['violations'],
                     rec['lemma_D']['sharp_equal'],
                     rec['lemma_D']['h_min'], rec['lemma_D']['h_max'],
                     rec['lemma_D']['zero_h_classes'],
                     rec['lemma_D']['zero_h_max']))
            print('      Arr_adm=%d (adm/realized=%.2fx)  '
                  'class-size hist=%s'
                  % (rec['arr_adm'], rec['adm_to_realized'],
                     rec['phase_structure']['size_hist']))
            gm = GT_MARG.get((p, fam))
            if gm:
                m_ok = (rec['marginals']['singletons'] == gm[0]
                        and rec['marginals']['pairs'] == gm[1])
                OUT['top_keys']['%d_%s' % (p, label)]['marg_gt_ok'] = m_ok
                print('      T-22 marginal GT: %s (sing %s, pairs %s)'
                      % ('OK' if m_ok else 'MISMATCH',
                         rec['marginals']['singletons'],
                         rec['marginals']['pairs']))

    # ---------- [4] the orbit test (p=19) ----------
    print()
    print('--- [4] the orbit test: S=(2,4,5,8) -> S2=(3,4,6,8), d=5 ---')
    p = 19
    S, S2, delta = (2, 4, 5, 8), (3, 4, 6, 8), 5
    inv = pow(delta, -1, p)
    lam = {2: (2 * inv) % p, 4: (4 * inv) % p, 5: (5 * inv) % p,
           8: (8 * inv) % p}
    sdS, _ = STORED[(p, S)]
    sdS2, _ = build_sizedict(p, S2)   # S2 not in FAMS: build inline
    totS = sum(len(v) for v in sdS.values())
    totS2 = sum(len(v) for v in sdS2.values())
    topS = max(sdS, key=lambda s: len(sdS[s]))
    topS2 = max(sdS2, key=lambda s: len(sdS2[s]))
    recS = analyze_key(p, S, topS, sorted(sdS[topS]))
    recS2 = analyze_key(p, S2, topS2, sorted(sdS2[topS2]))
    OUT['orbit'] = {'p': p, 'S': list(S), 'S2': list(S2),
                    'total_S': totS, 'total_S2': totS2,
                    'top_S': recS, 'top_S2': recS2}
    print('  totals: %d vs %d (committed: 801 each)' % (totS, totS2))
    print('  top decompositions:')
    for tag, r in (('S ', recS), ('S2', recS2)):
        print('    %s top szt=%s |Sol|=%d N_arr=%d avg=%.3f max=%d'
              % (tag, r['szt'], r['n_sol'], r['n_arr'], r['avg_class'],
                 r['max_class']))
    # the image of the source diagonal direction under the renorm map
    s = topS[1]
    sols = sorted(sdS[topS])
    if sols:
        a = sols[0]

        def renorm(a):
            a2, a3, a4, a5 = a
            return ((inv * a3 + lam[4] * (s - 1) - inv * a4) % p,
                    (0 - inv * a4) % p,
                    (inv * a5 + lam[8] * (s - 1) - inv * a4) % p,
                    (inv * a2 - inv * a4) % p)
        b0 = renorm(a)
        # diagonal image: differences of renorm(a+t*1) - renorm(a)
        t = 1
        aT = tuple((x + t) % p for x in a)
        bT = renorm(aT)
        v_star = tuple((bT[i] - b0[i]) % p for i in range(4))
        # honest check: is v_star a scalar multiple of (1,1,1,1)?
        par = (v_star[0] != 0
               and len({(v_star[i] * pow(v_star[0], -1, p)) % p
                        for i in range(4)}) == 1)
        OUT['orbit']['diagonal_image'] = {
            'a0': list(a), 'b0': list(b0), 'v_star': list(v_star),
            'parallel_to_diagonal': par}
        print('  the source diagonal maps to direction %s -- parallel to '
              'the target diagonal: %s' % (v_star, par))
        # do source classes map into single target classes?
        clsS = decompose(p, sols)
        setS2 = set(sdS2[topS2]) if topS2 == topS else None
        if setS2 is not None and sols:
            mapped_single = 0
            tot = 0
            for key, lst in clsS.items():
                imgs = [renorm(x) for x in lst]
                keys2 = {arr_key(p, b) for b in imgs}
                tot += 1
                if len(keys2) == 1:
                    mapped_single += 1
            OUT['orbit']['classes_mapped_to_single_class'] = {
                'single': mapped_single, 'total': tot}
            print('  source classes mapping into a single target class: '
                  '%d / %d' % (mapped_single, tot))

    # ---------- [5] verdict data: the margin factorization ----------
    print()
    print('--- [5] the margin factorization (committed top keys) ---')
    print('  %-14s %8s %8s %10s %10s %10s' %
          ('key', '|Sol|', 'margin', 'arr-factor', 'phase-factor',
           'N_arr/ov^3'))
    OUT['factorization'] = []
    for p in (17, 19, 29):
        for label, fam in FAMS[p]:
            key = '%d_%s' % (p, label)
            r = OUT['top_keys'].get(key)
            if not r:
                continue
            row = [key, r['n_sol'], r['margin_over_trivial'],
                   r['arr_compression'], r['phase_factor'],
                   r['narr_over_ov3']]
            OUT['factorization'].append(row)
            print('  %-14s %8d %8.1f %10.2f %10.2f %10.4f'
                  % (key, row[1], row[2], row[3], row[4], row[5]))

    OUT['runtime_sec'] = round(time.time() - t_start, 1)
    print()
    print('runtime: %.1f s' % OUT['runtime_sec'])
    with open('/home/z/my-project/scripts/out_t26_session.json', 'w') as f:
        json.dump(OUT, f, indent=1, default=list)
    print('output: scripts/out_t26_session.json')


if __name__ == '__main__':
    main()
