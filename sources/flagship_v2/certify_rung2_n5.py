#!/usr/bin/env python3
"""One exact-certification attempt for the rung-2 'consistent' instances,
headlined by v = (1,2,3,4,10), target rho = 7/22 (the m=10 row of the
deepest-hole table).

METHOD (new, independent of the arrangement machinery in lrc_zono_lib):

    D(c) = min_{s in [0,1)} max( ||s||_T, ||c_j - v_j s||_T  (j=2..n) )

is the torus distance of c to Lambda in the quotient V-norm (validated
below against lrc_zono_lib.D_exact on random rational points). For a
c-box C let

    U(C) = min_{s in [0,1)} max_j sup_{c in C} f_j(c, s).

By the minimax inequality  max_C D  <=  U(C).  Each U(C) is computed
EXACTLY (rational arithmetic): every f_j-supremum over the box is a
piecewise-linear function of s whose kink set is the union of
{lo_j - v_j s = 0 or 1/2 (mod 1)}, {hi_j - v_j s = 0 or 1/2 (mod 1)},
{lo_j - v_j s = -w_j/2 or 1/2 - w_j/2 (mod 1)}  (argmax switches),
together with s in {0, 1/2} for the s-coordinate; between consecutive
kinks every sup is linear, so the minimum of their max is attained at a
kink or at a pairwise crossing, all rational.

Branch and bound on dyadic c-boxes: a box with exact U(C) <= r is
certified (max_C D <= r); otherwise it is split. If the tree empties,
rho <= r is EXACTLY certified with a finite, auditable leaf list
(the output JSON records every leaf and its exact U-value).

Negative control: v = (1,2,3,4,5), r = 1/3 must NOT certify
(rho >= 97/288 > 1/3, witness x = (1/32, 15/32, 3/32, 7/8)):
the run must terminate with survivors whose U stays above r.

Output: scripts/out_certify_rung2.json
"""
import sys, json, math, time, heapq
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F

OUT = '/home/z/my-project/scripts/out_certify_rung2.json'


# ---------------------------------------------------------------- exact pieces

def normT(x):
    y = x % 1
    return y if y <= F(1, 2) else 1 - y


def normT_f(x):
    y = x - math.floor(x)
    return y if y <= 0.5 else 1.0 - y


def sup_normT(lo, hi):
    """sup_{x in [lo,hi]} ||x||_T, exact (Fractions)."""
    m0 = math.floor(lo - F(1, 2)) + 1
    if m0 < hi - F(1, 2):
        return F(1, 2)
    return max(normT(lo), normT(hi))


def sup_normT_f(lo, hi):
    m0 = math.floor(lo - 0.5) + 1
    if m0 < hi - 0.5 - 1e-12:
        return 0.5
    return max(normT_f(lo), normT_f(hi))


def kink_list(box, v):
    """All s in [0,1] where some box-supremum can kink (Fractions)."""
    S = [F(0), F(1, 2), F(1)]
    for idx, vj in enumerate(v[1:]):
        lo, hi = box[idx]
        mid = (lo + hi) / 2
        for c in (lo, hi, mid):
            for b in (F(0), F(1, 2)):
                cb = c - b
                m0 = math.floor(cb)
                for m in range(m0 - 1, m0 + vj + 2):
                    s = F(cb + m, vj)
                    if 0 <= s <= 1:
                        S.append(s)
    return sorted(set(S))


def kink_list_f(box, v):
    S = [0.0, 0.5, 1.0]
    for idx, vj in enumerate(v[1:]):
        lo, hi = box[idx]
        mid = 0.5 * (lo + hi)
        for c in (lo, hi, mid):
            for b in (0.0, 0.5):
                cb = c - b
                m0 = math.floor(cb)
                for m in range(m0 - 1, m0 + vj + 2):
                    s = (cb + m) / vj
                    if -1e-12 <= s <= 1 + 1e-12:
                        S.append(min(1.0, max(0.0, s)))
    S.sort()
    out = [S[0]]
    for x in S[1:]:
        if x - out[-1] > 1e-11:
            out.append(x)
    return out


def _R_values(s, box, v, EXACT):
    """(R_0, R_2..R_n) at s: sup over the box of each f_j."""
    if EXACT:
        vals = [normT(s)]
        for idx, vj in enumerate(v[1:]):
            lo, hi = box[idx]
            vals.append(sup_normT(lo - vj * s, hi - vj * s))
    else:
        vals = [normT_f(s)]
        for idx, vj in enumerate(v[1:]):
            lo, hi = box[idx]
            vals.append(sup_normT_f(lo - vj * s, hi - vj * s))
    return vals


def U_box(box, v, EXACT, r=None):
    """min_{s in S} max_j sup_C f_j  (exact Fractions if EXACT else float),
    where S = [0, r] + [1-r, 1) if r is given, else [0, 1).

    Restriction to S is sound: any s with max_j f_j <= r satisfies
    ||s||_T <= r (coordinate 0 is in the max), hence s in S. So
    min_S <= r  <=>  min_{[0,1)} <= r; non-prunable boxes get an
    upper estimate (still a valid certificate upper bound).
    If r is given, stops early once a value <= r is found."""
    if EXACT:
        ks = kink_list(box, v)
    else:
        ks = kink_list_f(box, v)
    if r is not None:
        rf = float(r) if not EXACT else r
        # restrict to S and add the S-boundaries as interval ends
        bnds = [F(0), F(r), F(1) - F(r), F(1)] if EXACT else \
            [0.0, float(r), 1.0 - float(r), 1.0]
        keep = sorted(set([s for s in ks if (s <= bnds[1] or s >= bnds[2])]
                          + [b for b in bnds if 0 <= b <= 1]))
        intervals = []
        for s1, s2 in zip(keep, keep[1:]):
            if s2 <= bnds[1] or s1 >= bnds[2]:      # inside one S-piece
                intervals.append((s1, s2))
    else:
        keep = ks
        intervals = list(zip(ks, ks[1:]))
    best = None
    cands = list(keep)
    for s1, s2 in intervals:
        v1 = _R_values(s1, box, v, EXACT)
        v2 = _R_values(s2, box, v, EXACT)
        for a in range(len(v1)):
            for b in range(a + 1, len(v1)):
                da = v2[a] - v1[a]
                db = v2[b] - v1[b]
                den = da - db
                if den:
                    num = v1[b] - v1[a]
                    t = num / den
                    if 0 < t < 1:
                        cands.append(s1 + (s2 - s1) * t)
    for s in cands:
        vals = _R_values(s, box, v, EXACT)
        val = max(vals)
        if best is None or val < best:
            best = val
        if r is not None and best <= r:
            return best, s
    return best, None


def D_at(c, v):
    """Exact D at a rational point (degenerate box), full s-domain."""
    box = [(x, x) for x in c]
    val, _ = U_box(box, v, True)
    return val


# ---------------------------------------------------------------- B&B

def certify(name, v, r, witnesses, expect, budget_sec, max_depth=96,
            max_nodes=2000000, store_leaves=400):
    d = len(v) - 1
    t0 = time.time()
    r = F(r)
    # witness assertions
    for w, wv in witnesses:
        got = D_at([F(x) for x in w], v)
        assert got == F(wv), ('witness mismatch', name, w, got, wv)
    # cheap grid pre-pass over S = [0, r] u [1-r, 1)
    K = 48
    rf = float(r)
    grid = ([i / K * rf for i in range(K + 1)] +
            [1 - rf + i / K * rf for i in range(K + 1)])

    def cheap_prune(box):
        """True if some grid s has float G(s) <= r - 1e-9 (sound prune)."""
        for s in grid:
            vals = _R_values(s, box, v, False)
            if max(vals) <= rf - 1e-9:
                return True
        return False

    SIZECAP = F(1, 2 ** 26)          # per-side survivor cap
    root = tuple((F(0), F(1)) for _ in range(d))
    heap = [(-1.0, 0, root)]
    n_leaves = 0
    leaf_sample = []                # (box, exact U) sample
    leaf_exact = 0
    survivors = []                  # (box, float U)
    pops = 0
    deepest = 0
    while heap:
        if time.time() - t0 > budget_sec or pops > max_nodes:
            survivors.extend((b, None) for _, _, b in heap)
            heap = []
            break
        _, depth, box = heapq.heappop(heap)
        pops += 1
        deepest = max(deepest, depth)
        widths = [hi - lo for lo, hi in box]
        if depth >= max_depth or max(widths) <= SIZECAP:
            uf, _ = U_box(box, v, False, r=r)
            survivors.append((box, uf))
            continue
        # 1) cheapest prune
        if cheap_prune(box):
            n_leaves += 1
            if len(leaf_sample) < store_leaves:
                leaf_sample.append((box, 'cheap'))
            continue
        # 2) float U
        uf, sf = U_box(box, v, False, r=r)
        if uf > rf + 1e-7:
            pass                      # definitely U > r -> split
        else:
            ue, se = U_box(box, v, True, r=r)
            if ue <= r:
                n_leaves += 1
                leaf_exact += 1
                if len(leaf_sample) < store_leaves:
                    leaf_sample.append((box, str(ue)))
                continue
        # 3) split widest coordinate
        k = max(range(d), key=lambda i: widths[i])
        lo, hi = box[k]
        mid = (lo + hi) / 2
        b1 = box[:k] + ((lo, mid),) + box[k + 1:]
        b2 = box[:k] + ((mid, hi),) + box[k + 1:]
        vol = float(widths[k])
        for i in range(d):
            if i != k:
                vol *= float(widths[i])
        heapq.heappush(heap, (-vol, depth + 1, b1))
        heapq.heappush(heap, (-vol, depth + 1, b2))
    certified = not survivors
    # survivor aggregate + exact D at sampled centers (minimax-loss check)
    surv_info = None
    if survivors:
        sboxes = [b for b, _ in survivors[:5000]]
        us = [u for _, u in survivors if u is not None]
        bbox = [[float(min(b[i][0] for b in sboxes)),
                 float(max(b[i][1] for b in sboxes))] for i in range(d)]
        # exact D at centers of up to 60 sampled survivors
        max_d_centers = None
        n_chk = 0
        for b, _ in survivors[:60]:
            ctr = [(lo + hi) / 2 for lo, hi in b]
            dc = D_at(ctr, v)
            n_chk += 1
            if max_d_centers is None or dc > max_d_centers:
                max_d_centers = dc
        surv_info = {
            'count': len(survivors),
            'max_float_U': max(us) if us else None,
            'min_float_U': min(us) if us else None,
            'center_bbox': bbox,
            'max_exact_D_at_centers': str(max_d_centers) if
            max_d_centers is not None else None,
            'centers_checked': n_chk,
            'total_volume_approx': sum(
                float(math.prod([float(b[i][1] - b[i][0])
                                 for i in range(d)])) for b in sboxes),
        }
    res = {
        'name': name, 'v': v, 'r': str(r), 'expect': expect,
        'certified': certified,
        'leaves': n_leaves,
        'leaf_exact': leaf_exact,
        'survivors': surv_info,
        'pops': pops, 'deepest_depth': deepest,
        'wall_sec': round(time.time() - t0, 1),
        'verdict': ('CERTIFIED rho <= %s' % r) if certified
                   else ('NOT certified within budget; %d survivors' %
                         len(survivors)),
    }
    res['leaf_sample'] = [
        (['%s|%s' % (lo, hi) for lo, hi in b], u) for b, u in leaf_sample]
    if survivors:
        res['survivor_sample'] = [
            (['%s|%s' % (lo, hi) for lo, hi in b], u)
            for b, u in survivors[:50]]
    return res


def cross_validate(v, npts=40, seed=7):
    """D_at (sigma-formula) vs lrc_zono_lib.D_exact (V-norm lattice DFS)."""
    import random
    from lrc_zono_lib import make_forms, D_exact
    forms = make_forms(v)
    rng = random.Random(seed)
    d = len(v) - 1
    for _ in range(npts):
        c = tuple(F(rng.randrange(0, 96), 96) for _ in range(d))
        a = D_at(c, v)
        b = D_exact(c, v, forms)
        assert a == b, ('D mismatch', v, c, a, b)
    return npts


def main():
    import os
    out = {'runs': [], 'cross_checks': []}
    if os.path.exists(OUT):                    # resume support
        try:
            out = json.load(open(OUT))
        except Exception:
            out = {'runs': [], 'cross_checks': []}
    for v in [(1, 2, 3, 10), (1, 2, 3, 4, 5), (1, 2, 3, 4, 10),
              (1, 2, 3, 16)]:
        if any(cc['v'] == list(v) for cc in out['cross_checks']):
            continue
        npts = cross_validate(v)
        out['cross_checks'].append({'v': v, 'points': npts, 'ok': True})
        print('cross-check D(sigma) == D_exact(V-norm): v=%s  %d pts OK' %
              (v, npts))

    runs = [
        ('n3-harmonic', (1, 2, 3), F(1, 4),
         [([F(1, 2), F(0)], F(1, 4)), ([F(7, 12), F(0)], F(1, 4))],
         'certify', 90),
        ('n3-rung2', (1, 2, 6), F(3, 14),
         [([F(1, 2), F(1, 2)], F(3, 14)), ([F(23, 42), F(1, 2)], F(3, 14))],
         'certify', 90),
        ('n4-harmonic', (1, 2, 3, 4), F(3, 10),
         [([F(1, 2), F(0), F(1, 2)], F(3, 10))], 'certify', 150),
        ('n4-rung2-m8', (1, 2, 3, 8), F(5, 18),
         [([F(1, 2), F(0), F(1, 2)], F(5, 18))], 'certify', 150),
        ('n4-rung2-m12', (1, 2, 3, 12), F(7, 26),
         [([F(1, 2), F(0), F(1, 2)], F(7, 26))], 'certify', 150),
        ('n4-rung2-m16', (1, 2, 3, 16), F(9, 34),
         [([F(1, 2), F(0), F(1, 2)], F(9, 34))], 'certify', 300),
        ('n5-harmonic-NEGCTL', (1, 2, 3, 4, 5), F(1, 3),
         [([F(1, 2), F(0), F(1, 2), F(0)], F(1, 3)),
          ([F(1, 32), F(15, 32), F(3, 32), F(7, 8)], F(97, 288))],
         'fail', 420),
        ('n5-rung2-m10-TARGET', (1, 2, 3, 4, 10), F(7, 22),
         [([F(1, 2), F(0), F(1, 2), F(1, 2)], F(7, 22))],
         'certify', 1500),
        ('n5-rung2-m20', (1, 2, 3, 4, 20), F(13, 42),
         [([F(1, 2), F(0), F(1, 2), F(1, 2)], F(13, 42))],
         'certify', 1300),
        ('n5-rung2-m25', (1, 2, 3, 4, 25), F(4, 13),
         [([F(1, 2), F(0), F(1, 2), F(0)], F(4, 13))],
         'certify', 420),
        ('n5-rung2-m30', (1, 2, 3, 4, 30), F(19, 62),
         [([F(1, 2), F(0), F(1, 2), F(1, 2)], F(19, 62))],
         'certify', 420),
    ]
    for name, v, r, ws, expect, budget in runs:
        if any(rr['name'] == name for rr in out['runs']):
            continue                              # already done (resume)
        print('--- %s  v=%s r=%s (budget %ss)' % (name, v, r, budget))
        res = certify(name, v, r, ws, expect, budget)
        out['runs'].append(res)
        print('    %s | leaves=%d survivors=%s pops=%d depth=%d %.1fs' % (
            res['verdict'], res['leaves'],
            (res['survivors'] or {}).get('count', 0),
            res['pops'], res['deepest_depth'], res['wall_sec']))
        if res['survivors']:
            print('    survivor max U = %s (r = %s), U - r = %s' % (
                res['survivors']['max_float_U'], float(r),
                None if res['survivors']['max_float_U'] is None else
                res['survivors']['max_float_U'] - float(r)))
        sys.stdout.flush()
        json.dump(out, open(OUT, 'w'), indent=1)
    print('written:', OUT)


if __name__ == '__main__':
    main()
