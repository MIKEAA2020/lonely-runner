#!/usr/bin/env python3
"""The harmonic-maximum and boundary-row certification session (v4).

TARGETS
  (a) rho(1,2,3,4,5) exactly.  The archive's recorded lower bound is
      4183/12288 (the audit witness, superseding 607/1792 and 97/288);
      the deepest hole of the harmonic instance has never been certified
      from above.  The Lipschitz box certificate closes exactly at a
      nondegenerate PL argmax (the m = 10 mechanism: at the optimal
      sigma* the box on each binding side prunes with sup = r attained
      at the witness corner).
  (b) The two n=6 boundary rows: v = (1,2,3,4,5,9), witness
      (0,1/2,0,0,0), and v = (1,2,3,4,5,12), witness
      (1/16,11/16,1/8,7/16,7/8), each with D = 5/14 = thr_6 exactly.
      Question: rho >, =, or < 5/14?
  (c) The three n=5 boundary rows m = 6,7,8 (witness value exactly
      1/3 = thr_5): the same question at n=5, completing the boundary
      ledger.

METHOD (unchanged from the validated m = 10 second attempt,
scripts/certify_m10_v2.py; primitives in certify_rung2_n5.py):

    D(c) = min_{s in [0,1)} max( ||s||_T, ||c_j - v_j s||_T (j=2..n) )

Sound search screen on a 2x2048 s-grid over S = [0,r] u [1-r,1) (the
grid minimum is an upper bound for D, so no point above the threshold
can be missed); Lipschitz box certificate max_C D <= D(c0) + h; exact
rational minimax U(C); dyadic branch and bound; every sampled prune
re-confirmed exactly in phase `verify`.

NEW in this session: the harmonic search is seeded with the 30,000
stored centers of the negative-control survivor boxes (the deep-hole
clusters of the harmonic instance recorded in the m = 10 session), and
the argmax hunt runs batched uphill climbs + Nelder-Mead + exact
re-evaluation on everything near the frontier.

Phases: selftest hsearch hcert n6m9 n6m12 n5b verify analyze status
Output: scripts/out_certify_harmonic_v4.json (+ _survivors.json,
_state.json; checkpoint-resumed across ~460 s foreground slices).
Usage: python3 certify_harmonic_v4.py <phase> [<phase> ...]
"""
import sys, os, json, math, time, heapq
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F
import numpy as np

import certify_rung2_n5 as C          # validated primitives: U_box, D_at

OUT = '/home/z/my-project/scripts/out_certify_harmonic_v4.json'
OUT_SURV = '/home/z/my-project/scripts/out_certify_harmonic_v4_survivors.json'
STATE = '/home/z/my-project/scripts/out_certify_harmonic_v4_state.json'
NEGCTL_SURV = '/home/z/my-project/scripts/out_certify_m10_v2_survivors.json'

VH = (1, 2, 3, 4, 5)          # harmonic n=5  (target a)
V69 = (1, 2, 3, 4, 5, 9)      # n=6 boundary row m=9  (target b)
V612 = (1, 2, 3, 4, 5, 12)    # n=6 boundary row m=12 (target b)
V56 = (1, 2, 3, 4, 6)         # n=5 boundary rows (target c)
V57 = (1, 2, 3, 4, 7)
V58 = (1, 2, 3, 4, 8)
R6 = F(5, 14)                 # thr_6
R5 = F(1, 3)                  # thr_5
CHUNK = 512                   # rows per numpy batch
SLICE_SEC = 460
HEAP_CAP = 1500000            # large: the harmonic argmax band is wide
SURV_CAP = 40000

WITNESSES = [
    ('harm-97/288', VH, [F(1, 32), F(15, 32), F(3, 32), F(7, 8)], F(97, 288)),
    ('harm-607/1792', VH, [F(11, 256), F(163, 256), F(15, 256),
                           F(79, 256)], F(607, 1792)),
    ('harm-4183/12288', VH, [F(465, 4096), F(605, 1024), F(1019, 4096),
                             F(315, 4096)], F(4183, 12288)),
    ('n6m9-5/14', V69, [F(0), F(1, 2), F(0), F(0), F(0)], R6),
    ('n6m12-5/14', V612, [F(1, 16), F(11, 16), F(1, 8), F(7, 16),
                          F(7, 8)], R6),
    ('n6m6-29/80', (1, 2, 3, 4, 5, 6),
     [F(0), F(5, 16), F(15, 16), F(3, 4), F(1, 2)], F(29, 80)),
    ('n6m7-4/11', (1, 2, 3, 4, 5, 7),
     [F(1, 16), F(9, 16), F(0), F(3, 8), F(0)], F(4, 11)),
    ('n5m6-1/3', V56, [F(31, 32), F(17, 32), F(15, 16), F(1, 16)], R5),
    ('n5m7-1/3', V57, [F(5, 32), F(23, 32), F(5, 16), F(17, 32)], R5),
    ('n5m8-1/3', V58, [F(7, 8), F(7, 16), F(7, 8), F(3, 4)], R5),
    ('m10-7/22 (control)', (1, 2, 3, 4, 10),
     [F(1, 2), F(0), F(1, 2), F(1, 2)], F(7, 22)),
]


# ------------------------------------------------------------------ io helpers

def load_out():
    if os.path.exists(OUT):
        try:
            return json.load(open(OUT))
        except Exception:
            pass
    return {}


def dump_out(out):
    json.dump(out, open(OUT, 'w'), indent=1, default=str)


def load_state():
    if os.path.exists(STATE):
        try:
            return json.load(open(STATE))
        except Exception:
            pass
    return None


def dump_state(st):
    json.dump(st, open(STATE, 'w'), default=str)


def box_to_str(box):
    return ['%s|%s' % (lo, hi) for lo, hi in box]


def box_from_str(ss):
    out = []
    for s in ss:
        lo, hi = s.split('|')
        out.append((F(lo), F(hi)))
    return tuple(out)


# ------------------------------------------------------------------ numpy grids

_GRIDS = {}


def grids(rf):
    if rf not in _GRIDS:
        n = 8192
        gd = np.concatenate([np.linspace(0.0, rf, n),
                             np.linspace(1.0 - rf, 1.0, n)])
        g0d = np.abs(gd - np.round(gd))
        m = 2048
        gs = np.concatenate([np.linspace(0.0, rf, m),
                             np.linspace(1.0 - rf, 1.0, m)])
        g0s = np.abs(gs - np.round(gs))
        _GRIDS[rf] = (gd, g0d, gs, g0s)
    return _GRIDS[rf]


def G_point_batch(P, v, grid, g0):
    """min-grid UPPER bounds of D(c) for a row-chunked batch P (B x d).
    Sound: the grid is a subset of S, so min over grid >= D."""
    d = len(v) - 1
    out = np.empty(len(P))
    for i0 in range(0, len(P), CHUNK):
        Q = P[i0:i0 + CHUNK]
        M = np.broadcast_to(g0, (len(Q), len(grid))).astype(np.float64)
        for j in range(d):
            vj = float(v[j + 1])
            x = Q[:, j][:, None] - vj * grid[None, :]
            M = np.maximum(M, np.abs(x - np.round(x)))
        out[i0:i0 + len(Q)] = M.min(axis=1)
    return out


def G_box_sup(box, v, grid, g0):
    """max_j sup_{c in box} f_j(c, s) over a fixed s-grid (1-D array)."""
    M = g0.copy()
    for j, vj in enumerate(v[1:]):
        lo = float(box[j][0]); hi = float(box[j][1])
        a = lo - vj * grid
        b = hi - vj * grid
        na = np.abs(a - np.round(a))
        nb = np.abs(b - np.round(b))
        sup = np.maximum(na, nb)
        m0 = np.floor(a - 0.5) + 1.0        # half-integer strictly inside?
        sup = np.where(m0 < b - 0.5, 0.5, sup)
        M = np.maximum(M, sup)
    return M


def floatD_point(c, v):
    """Float D at a point via the validated kink evaluator (full s-domain).
    Error << 1e-9; used only to GATE exact confirmations."""
    return C.U_box([(x, x) for x in c], v, False, r=None)[0]


def argmin_sigma(c, v):
    """Exact argmin sigma and value of the runner-form max at a point c
    (full s-domain); for the paper's hand-check records."""
    box = [(x, x) for x in c]
    ks = C.kink_list(box, v)
    cands = list(ks)
    for s1, s2 in zip(ks, ks[1:]):
        v1 = C._R_values(s1, box, v, True)
        v2 = C._R_values(s2, box, v, True)
        for a in range(len(v1)):
            for b in range(a + 1, len(v1)):
                da = v2[a] - v1[a]
                db = v2[b] - v1[b]
                den = da - db
                if den:
                    t = (v1[b] - v1[a]) / den
                    if 0 < t < 1:
                        cands.append(s1 + (s2 - s1) * t)
    best, bs = None, None
    for s in cands:
        val = max(C._R_values(s, box, v, True))
        if best is None or val < best:
            best, bs = val, s
    return bs, best


# ------------------------------------------------------------------ phase: selftest

def phase_selftest(out):
    if 'selftest' in out:
        return out['selftest']
    res = {}
    t0 = time.time()
    for name, v, w, wv in WITNESSES:
        got = C.D_at(list(w), v)
        assert got == wv, ('witness mismatch', name, w, got, wv)
    res['witnesses_exact'] = len(WITNESSES)
    import random
    from lrc_zono_lib import make_forms, D_exact
    rng = random.Random(17)
    for v in (VH, V69, V612, V56):
        forms = make_forms(v)
        for _ in range(40):
            c = tuple(F(rng.randrange(0, 96), 96)
                      for _ in range(len(v) - 1))
            a = C.D_at(list(c), v)
            b = D_exact(c, v, forms)
            assert a == b, ('D mismatch', v, c, a, b)
    res['cross_checks'] = '160 points, four instances, exact agreement'
    adv = [
        (VH, F(1, 3), [F(11, 256), F(163, 256), F(15, 256), F(79, 256)]),
        (VH, F(1, 3), [F(465, 4096), F(605, 1024), F(1019, 4096),
                       F(315, 4096)]),
        ((1, 2, 3, 4, 5, 6), R6,
         [F(0), F(5, 16), F(15, 16), F(3, 4), F(1, 2)]),
        ((1, 2, 3, 4, 5, 7), R6,
         [F(1, 16), F(9, 16), F(0), F(3, 8), F(0)]),
    ]
    refusals = []
    for v, r, w in adv:
        h = F(1, 1024)
        box = [(x, x + F(1, 512)) for x in
               (F(math.floor(float(x) * 512), 512) for x in w)]
        c0 = [(lo + hi) / 2 for lo, hi in box]
        dc = C.D_at(c0, v)
        assert dc + h > r, ('ADVERSARIAL FAILURE: the Lipschitz prune '
                            'would kill a true deeper hole', v, r, w, dc)
        refusals.append({'v': list(v), 'r': str(r), 'center_D': str(dc)})
    res['adversarial_refusals'] = refusals
    # positive Lipschitz controls at each target threshold
    pruned = 0
    for v, r in [(VH, F(4183, 12288)), (V69, R6), (V612, R6), (V56, R5)]:
        ctr = [F(1, 4)] * (len(v) - 1)
        box = [(x, x + F(1, 512)) for x in ctr]
        c0 = [(lo + hi) / 2 for lo, hi in box]
        if C.D_at(c0, v) + F(1, 1024) <= r:
            pruned += 1
    assert pruned >= 3, ('positive Lipschitz control too weak', pruned)
    res['positive_controls_pruned'] = pruned
    res['sec'] = round(time.time() - t0, 1)
    out['selftest'] = res
    dump_out(out)
    print('[selftest] OK:', json.dumps(res)[:600])
    return res


# ------------------------------------------------------------------ phase: hsearch

def _exact_rows(cands, v, tag, res, budget_sec=600):
    """Exact D_at for candidate points; returns rows (fd, str, c_strs)
    sorted by depth, deepest first."""
    t0 = time.time()
    seen = set()
    rows = []
    for c in cands:
        if time.time() - t0 > budget_sec:
            break
        key = tuple(round(float(x), 9) for x in c)
        if key in seen:
            continue
        seen.add(key)
        cc = [x if isinstance(x, F) else F(x) for x in c]
        d = C.D_at(cc, v)
        rows.append((float(d), str(d), [str(x) for x in cc]))
    rows.sort(key=lambda t: -t[0])
    res['exact_' + tag] = {'n': len(rows),
                           'deepest': rows[0][1] if rows else None,
                           'deepest_c': rows[0][2] if rows else None}
    return rows


def phase_hsearch(out):
    """Harmonic argmax hunt (no fixed r: collect the frontier, then exact).
    Stage machine with checkpoint resume; returns None mid-slice."""
    if 'hsearch' in out:
        print('[hsearch] already done')
        return out['hsearch']
    v = VH
    rf = 0.35                       # screen domain S(0.35): sound for any
    gd, g0d, gs, g0s = grids(rf)    # candidate with D <= 0.35
    t0 = time.time()
    st = load_state() or {}
    if 'hsearch' not in st:
        st['hsearch'] = {'stage': 'seeds', 'chunk': 0, 'top': [],
                         'cand': [], 'best': None, 'extra': {}}
    S = st['hsearch']
    rng = np.random.default_rng(20261008)

    def save():
        S['top'] = S['top'][:4096]
        S['cand'] = S['cand'][:8192]
        if 'pairs' in S:
            S['pairs'] = S['pairs'][:2048]
        dump_state(st)
        print('[hsearch] slice save: stage=%s chunk=%d top=%d cand=%d'
              % (S['stage'], S.get('chunk', 0), len(S['top']),
                 len(S['cand'])))

    def screen(P):
        Dg = G_point_batch(P, v, gs, g0s)
        for i in np.argsort(-Dg)[:32]:
            S['top'].append((float(Dg[i]), [float(x) for x in P[i]]))
        hit = Dg > 0.336
        if hit.any() and len(S['cand']) < 8192:
            for i in np.nonzero(hit)[0][:256]:
                S['cand'].append([float(x) for x in P[i]])
        return Dg

    while True:
        stage = S['stage']
        if stage == 'seeds':
            boxes = json.load(open(NEGCTL_SURV))
            while S['chunk'] < len(boxes):
                lo, hi = S['chunk'], min(S['chunk'] + 1024, len(boxes))
                for sb in boxes[lo:hi]:
                    c = []
                    for s in sb:
                        a, b = s.split('|')
                        c.append((F(a) + F(b)) / 2)
                    d = C.D_at(c, v)
                    fd = float(d)
                    if S['best'] is None or fd > S['best'][0]:
                        S['best'] = (fd, [str(x) for x in c], str(d))
                    S['top'].append((fd, [float(x) for x in c]))
                S['chunk'] = hi
                if time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['top'].sort(key=lambda t: -t[0])
            S['top'] = S['top'][:2048]
            for nm, _, wc, wv in WITNESSES[:3]:
                S['top'].append((float(wv), [float(x) for x in wc]))
                fd = float(wv)
                if S['best'] is None or fd > S['best'][0]:
                    S['best'] = (fd, [str(x) for x in wc], str(wv))
            S['extra']['seed_boxes'] = len(boxes)
            S['extra']['seed_best_D'] = S['best'][2] if S['best'] else None
            S['stage'] = 'rand'; S['chunk'] = 0
        elif stage == 'rand':
            N = 1 << 21
            while S['chunk'] < N // 2048:
                P = rng.random((2048, 4))
                screen(P)
                S['chunk'] += 1
                if S['chunk'] % 64 == 0 and time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['stage'] = 'dyadic16'; S['chunk'] = 0
        elif stage == 'dyadic16':
            q = 16
            axis = np.arange(q) / q
            P = np.stack(np.meshgrid(*[axis] * 4, indexing='ij'),
                         -1).reshape(-1, 4)
            Dg = G_point_batch(P, v, gs, g0s)
            j = int(np.argmax(Dg))
            S['extra']['dyadic_q16_best_gridD'] = float(Dg[j])
            S['top'].append((float(Dg[j]), [float(x) for x in P[j]]))
            idx = np.argsort(-Dg)[:24]
            _exact_rows([[F(int(round(x * q)), q) for x in P[i]]
                         for i in idx], v, 'dyadic16', S['extra'], 120)
            S['stage'] = 'dyadic32'; S['chunk'] = 0
            if time.time() - t0 > SLICE_SEC - 60:
                save(); return None
        elif stage == 'dyadic32':
            q = 32
            axis = np.arange(q) / q
            P = np.stack(np.meshgrid(*[axis] * 4, indexing='ij'),
                         -1).reshape(-1, 4)
            for i0 in range(0, len(P), 2048):
                Dg = G_point_batch(P[i0:i0 + 2048], v, gs, g0s)
                for i in np.argsort(-Dg)[:8]:
                    S['top'].append((float(Dg[i]),
                                     [float(x) for x in P[i0 + i]]))
                hit = Dg > 0.336
                if hit.any() and len(S['cand']) < 8192:
                    for i in np.nonzero(hit)[0][:64]:
                        S['cand'].append(
                            [float(x) for x in P[i0 + i]])
                S['chunk'] += 1
                if S['chunk'] % 256 == 0 and time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['extra']['dyadic_q32_done'] = True
            S['stage'] = 'climbs'; S['chunk'] = 0
        elif stage == 'climbs':
            if 'X' not in S:
                S['top'].sort(key=lambda t: -t[0])
                seeds = [c for _, c in S['top'][:1536]]
                seeds += [[float(x) for x in rng.random(4)]
                          for _ in range(256)]
                X = np.array(seeds)
                S['X'] = X.tolist()
                S['DX'] = G_point_batch(X, v, gs, g0s).tolist()
                S['climb_evals'] = 0
            X = np.array(S['X'])
            DX = np.array(S['DX'])
            sig_idx = S['chunk']
            while sig_idx < 19:
                sig = 2.0 ** -(sig_idx + 4)
                for _round in range(8):
                    prop = []
                    for j in range(4):
                        for sgn in (+1.0, -1.0):
                            Y = X.copy()
                            Y[:, j] = (Y[:, j] + sgn * sig) % 1.0
                            prop.append(Y)
                    allprop = np.stack(prop)
                    vals = np.stack([G_point_batch(Y, v, gs, g0s)
                                     for Y in prop])
                    S['climb_evals'] += 8 * len(X)
                    bestv = vals.max(axis=0)
                    bestd = vals.argmax(axis=0)
                    acc = bestv > DX + 1e-12
                    if not acc.any():
                        break
                    idx = np.nonzero(acc)[0]
                    X[idx] = allprop[bestd[idx], idx]
                    DX[idx] = bestv[idx]
                if time.time() - t0 > SLICE_SEC:
                    S['X'] = X.tolist(); S['DX'] = DX.tolist()
                    S['chunk'] = sig_idx
                    save(); return None
                sig_idx += 1
                S['chunk'] = sig_idx
            S['X'] = X.tolist(); S['DX'] = DX.tolist()
            S['extra']['climb_starts'] = int(len(X))
            S['extra']['climb_evals'] = int(S['climb_evals'])
            S['extra']['climb_best_gridD'] = float(DX.max())
            S['stage'] = 'nm'; S['chunk'] = 0
        elif stage == 'nm':
            from scipy.optimize import minimize
            X = np.array(S['X']); DX = np.array(S['DX'])
            finals = X[np.argsort(-DX)[:128]]
            if 'polished' not in S:
                S['polished'] = []
            while S['chunk'] < 64:
                c0 = finals[S['chunk'] % len(finals)]

                def negD(c):
                    return -floatD_point(c, v)

                try:
                    rr = minimize(negD, np.array(c0, dtype=float),
                                  method='Nelder-Mead',
                                  options={'xatol': 1e-11, 'fatol': 1e-13,
                                           'maxiter': 800})
                    S['polished'].append(
                        [float(x) for x in np.mod(rr.x, 1.0)])
                except Exception:
                    pass
                S['chunk'] += 1
                if S['chunk'] % 4 == 0 and time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['extra']['nm_polished'] = len(S['polished'])
            S['extra']['nm_best_floatD'] = max(
                (floatD_point(c, v) for c in S['polished']), default=None)
            S['stage'] = 'refine'; S['chunk'] = 0
        elif stage == 'refine':
            if 'pairs' not in S:
                S['top'].sort(key=lambda t: -t[0])
                pairs = []
                for c in [c for _, c in S['top'][:300]]:
                    pairs.append([floatD_point(c, v), c])
                X = np.array(S['X']); DX = np.array(S['DX'])
                for c in X[np.argsort(-DX)[:128]]:
                    pairs.append([floatD_point(c, v), [float(x) for x in c]])
                for c in S['polished']:
                    pairs.append([floatD_point(c, v), c])
                for c in S['cand'][:4096]:
                    pairs.append([floatD_point(c, v), c])
                pairs.sort(key=lambda t: -t[0])
                S['pairs'] = pairs
            pairs = S['pairs']
            near = [c for fd, c in pairs if fd > 0.334][:1024]
            deep50 = [c for _, c in pairs[:50]]
            rows = _exact_rows(near + deep50, v, 'search', S['extra'], 600)
            S['extra']['n_refined'] = len(pairs)
            S['extra']['refined_best_floatD'] = (
                pairs[0][0] if pairs else None)
            for fd, ds, cs in rows:
                if S['best'] is None or fd > S['best'][0]:
                    S['best'] = (fd, cs, ds)
            if S['best'] is None or F(S['best'][2]) < F(4183, 12288):
                S['best'] = (float(F(4183, 12288)),
                             [str(x) for x in (F(465, 4096), F(605, 1024),
                                               F(1019, 4096),
                                               F(315, 4096))],
                             '4183/12288')
            S['extra']['best_exact'] = S['best'][2]
            S['extra']['best_c'] = S['best'][1]
            S['stage'] = 'done'
        else:
            break

    res = {
        'screen_grid': '2 x 2048 s-points on S(0.35) (sound upper bound)',
        'seed_boxes': S['extra'].get('seed_boxes'),
        'seed_best_D': S['extra'].get('seed_best_D'),
        'random_points': 1 << 21,
        'climb_starts': S['extra'].get('climb_starts'),
        'climb_evals': S['extra'].get('climb_evals'),
        'climb_best_gridD': S['extra'].get('climb_best_gridD'),
        'nm_polished': S['extra'].get('nm_polished'),
        'nm_best_floatD': S['extra'].get('nm_best_floatD'),
        'n_refined': S['extra'].get('n_refined'),
        'refined_best_floatD': S['extra'].get('refined_best_floatD'),
        'best_D': S['best'][2] if S['best'] else None,
        'best_c': S['best'][1] if S['best'] else None,
    }
    for k in ('exact_dyadic16', 'exact_search'):
        if k in S['extra']:
            res[k] = S['extra'][k]
    out['hsearch'] = res
    dump_out(out)
    st = load_state() or {}
    st.pop('hsearch', None)
    dump_state(st)
    print('[hsearch] COMPLETE:', json.dumps(res)[:700])
    return res


# -------------------------------------------------- phase: light search (fixed r)

def light_search(out, key, v, r, seed_witnesses):
    """Fixed-r instance search: random + dyadic + climbs + NM + exact,
    screened at r; any exact D > r is recorded as a refutation witness.
    Returns None while incomplete (re-invoke)."""
    if key in out:
        print('[%s] search already done' % key)
        return out[key]
    rf = float(r)
    gd, g0d, gs, g0s = grids(rf)
    d = len(v) - 1
    t0 = time.time()
    st = load_state() or {}
    if key + '_s' not in st:
        st[key + '_s'] = {'stage': 'rand', 'chunk': 0, 'top': [],
                          'cand': [], 'extra': {}}
    S = st[key + '_s']
    rng = np.random.default_rng(20261008 + d)

    def save():
        S['top'] = S['top'][:2048]
        S['cand'] = S['cand'][:4096]
        dump_state(st)
        print('[%s-s] slice save: stage=%s chunk=%d top=%d'
              % (key, S['stage'], S.get('chunk', 0), len(S['top'])))

    def screen(P):
        Dg = G_point_batch(P, v, gs, g0s)
        for i in np.argsort(-Dg)[:16]:
            S['top'].append((float(Dg[i]), [float(x) for x in P[i]]))
        hit = Dg > rf - 2e-3
        if hit.any() and len(S['cand']) < 4096:
            for i in np.nonzero(hit)[0][:128]:
                S['cand'].append([float(x) for x in P[i]])
        return Dg

    while True:
        stage = S['stage']
        if stage == 'rand':
            N = 1 << 20
            while S['chunk'] < N // 2048:
                screen(rng.random((2048, d)))
                S['chunk'] += 1
                if S['chunk'] % 64 == 0 and time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['stage'] = 'dyadic'; S['chunk'] = 0
        elif stage == 'dyadic':
            q = 16
            axis = np.arange(q) / q
            P = np.stack(np.meshgrid(*[axis] * d, indexing='ij'),
                         -1).reshape(-1, d)
            for i0 in range(0, len(P), 2048):
                Dg = G_point_batch(P[i0:i0 + 2048], v, gs, g0s)
                for i in np.argsort(-Dg)[:8]:
                    S['top'].append((float(Dg[i]),
                                     [float(x) for x in P[i0 + i]]))
                hit = Dg > rf - 2e-3
                if hit.any() and len(S['cand']) < 4096:
                    for i in np.nonzero(hit)[0][:64]:
                        S['cand'].append([float(x) for x in P[i0 + i]])
                S['chunk'] += 1
                if S['chunk'] % 256 == 0 and time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['stage'] = 'climbs'; S['chunk'] = 0
        elif stage == 'climbs':
            if 'X' not in S:
                S['top'].sort(key=lambda t: -t[0])
                seeds = [c for _, c in S['top'][:384]]
                seeds += [[float(x) for x in w] for w in seed_witnesses]
                seeds += [[float(x) for x in rng.random(d)]
                          for _ in range(128)]
                X = np.array(seeds)
                S['X'] = X.tolist()
                S['DX'] = G_point_batch(X, v, gs, g0s).tolist()
                S['climb_evals'] = 0
            X = np.array(S['X'])
            DX = np.array(S['DX'])
            sig_idx = S['chunk']
            while sig_idx < 15:
                sig = 2.0 ** -(sig_idx + 4)
                for _round in range(8):
                    prop = []
                    for j in range(d):
                        for sgn in (+1.0, -1.0):
                            Y = X.copy()
                            Y[:, j] = (Y[:, j] + sgn * sig) % 1.0
                            prop.append(Y)
                    allprop = np.stack(prop)
                    vals = np.stack([G_point_batch(Y, v, gs, g0s)
                                     for Y in prop])
                    S['climb_evals'] += 2 * d * len(X)
                    bestv = vals.max(axis=0)
                    bestd = vals.argmax(axis=0)
                    acc = bestv > DX + 1e-12
                    if not acc.any():
                        break
                    idx = np.nonzero(acc)[0]
                    X[idx] = allprop[bestd[idx], idx]
                    DX[idx] = bestv[idx]
                if time.time() - t0 > SLICE_SEC:
                    S['X'] = X.tolist(); S['DX'] = DX.tolist()
                    S['chunk'] = sig_idx
                    save(); return None
                sig_idx += 1
                S['chunk'] = sig_idx
            S['X'] = X.tolist(); S['DX'] = DX.tolist()
            S['extra']['climb_starts'] = int(len(X))
            S['extra']['climb_evals'] = int(S['climb_evals'])
            S['extra']['climb_best_gridD'] = float(DX.max())
            S['stage'] = 'nm'; S['chunk'] = 0
        elif stage == 'nm':
            from scipy.optimize import minimize
            X = np.array(S['X']); DX = np.array(S['DX'])
            finals = X[np.argsort(-DX)[:48]]
            if 'polished' not in S:
                S['polished'] = []
            while S['chunk'] < 32:
                c0 = finals[S['chunk'] % max(1, len(finals))]

                def negD(c):
                    return -floatD_point(c, v)

                try:
                    rr = minimize(negD, np.array(c0, dtype=float),
                                  method='Nelder-Mead',
                                  options={'xatol': 1e-11, 'fatol': 1e-13,
                                           'maxiter': 600})
                    S['polished'].append(
                        [float(x) for x in np.mod(rr.x, 1.0)])
                except Exception:
                    pass
                S['chunk'] += 1
                if S['chunk'] % 4 == 0 and time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['extra']['nm_polished'] = len(S['polished'])
            S['stage'] = 'refine'; S['chunk'] = 0
        elif stage == 'refine':
            if 'pairs' not in S:
                S['top'].sort(key=lambda t: -t[0])
                pairs = []
                for c in [c for _, c in S['top'][:200]]:
                    pairs.append([floatD_point(c, v), c])
                X = np.array(S['X']); DX = np.array(S['DX'])
                for c in X[np.argsort(-DX)[:64]]:
                    pairs.append([floatD_point(c, v), [float(x) for x in c]])
                for c in S['polished']:
                    pairs.append([floatD_point(c, v), c])
                for c in S['cand'][:2048]:
                    pairs.append([floatD_point(c, v), c])
                pairs.sort(key=lambda t: -t[0])
                S['pairs'] = pairs
            pairs = S['pairs']
            near = [c for fd, c in pairs if fd > rf - 1e-5][:512]
            deep30 = [c for _, c in pairs[:30]]
            rows, refuted = _exact_refute(near + deep30, v, r, 'search',
                                          S['extra'], 480)
            S['extra']['n_refined'] = len(pairs)
            S['extra']['refined_best_floatD'] = (
                pairs[0][0] if pairs else None)
            S['extra']['refuted'] = refuted
            S['stage'] = 'done'
        else:
            break

    res = {
        'screen_grid': '2 x 2048 s-points on S (sound upper bound)',
        'r': str(r), 'random_points': 1 << 20,
        'climb_starts': S['extra'].get('climb_starts'),
        'climb_evals': S['extra'].get('climb_evals'),
        'climb_best_gridD': S['extra'].get('climb_best_gridD'),
        'nm_polished': S['extra'].get('nm_polished'),
        'n_refined': S['extra'].get('n_refined'),
        'refined_best_floatD': S['extra'].get('refined_best_floatD'),
        'refuted': S['extra'].get('refuted'),
        'exact_search': S['extra'].get('exact_search'),
    }
    out[key] = res
    dump_out(out)
    st = load_state() or {}
    st.pop(key + '_s', None)
    dump_state(st)
    print('[%s-s] COMPLETE:' % key, json.dumps(res)[:500])
    return res


def _exact_refute(cands, v, r, tag, res, budget_sec=480):
    """Exact D_at for candidates; any D > r is a refutation witness."""
    t0 = time.time()
    seen = set()
    rows = []
    refuted = None
    for c in cands:
        if time.time() - t0 > budget_sec:
            break
        key = tuple(round(float(x), 9) for x in c)
        if key in seen:
            continue
        seen.add(key)
        cc = [x if isinstance(x, F) else F(x) for x in c]
        d = C.D_at(cc, v)
        rows.append((float(d), str(d), [str(x) for x in cc]))
        if d > F(r):
            refuted = {'c': [str(x) for x in cc], 'D': str(d)}
            break
    rows.sort(key=lambda t: -t[0])
    res['exact_' + tag] = {'n': len(rows),
                           'deepest': rows[0][1] if rows else None,
                           'deepest_c': rows[0][2] if rows else None}
    return rows, refuted


# ------------------------------------------------- sliding-sigma gap certificate

def gap_certificate(box, v, r):
    """Vertex-checked sliding-sigma interval certificate (sound, exact).

    For a fixed integer lift vector m = (m_2..m_n): a point c admits a
    sigma with every runner term <= r  iff  the interval
        [L(c), R(c)] = [max_j (c_j - m_j - r)/v_j, min_j (c_j - m_j + r)/v_j]
    intersects S = [0, r] u [1-r, 1] (the sigma-term constraint).  L is
    convex and R concave in c, so R - L is concave: its minimum over the
    box is attained at a vertex, as are the boundary conditions.  Hence
    checking the 2^d vertices certifies the WHOLE box -- with sigma
    allowed to SLIDE with c, which is exactly what boxes straddling a
    diagonal argmax band need (the single-sigma U-bound cannot do this).

    Returns the m-vector on success (exact Fractions), else None.
    Float pre-scan with early exit; every reported m is re-confirmed
    exactly with Fractions before returning."""
    d = len(v) - 1
    lo = [box[j][0] for j in range(d)]
    hi = [box[j][1] for j in range(d)]
    rfl = [float(x) for x in lo]
    rhl = [float(x) for x in hi]
    rf = float(r)
    vs = [float(v[j + 1]) for j in range(d)]
    # candidate lifts: need |c_j - v_j*sigma - m_j| <= r for some sigma in
    # [0,1) and c_j in [lo_j, hi_j]  =>  m_j in [lo_j - v_j - r, hi_j + r]
    msc = []
    for j in range(d):
        a = math.floor(lo[j] - v[j + 1] - r)
        b = math.floor(hi[j] + r)
        msc.append(range(a, b + 1))
    nverts = 1 << d
    cand = []
    # ---- float pre-scan over m-combos with early exit at first vertex
    for m in _iter_combos(msc):
        ok = True
        maxL = -1e18
        minR = 1e18
        for t in range(nverts):
            L = -1e18
            R = 1e18
            for j in range(d):
                cj = rfl[j] if not (t >> j) & 1 else rhl[j]
                a = (cj - m[j] - rf) / vs[j]
                if a > L:
                    L = a
                b = (cj - m[j] + rf) / vs[j]
                if b < R:
                    R = b
            if L > R + 1e-12:
                ok = False
                break
            if L > maxL:
                maxL = L
            if R < minR:
                minR = R
        if ok and ((maxL <= rf + 1e-12 and minR >= -1e-12) or
                   (maxL <= 1.0 + 1e-12 and minR >= 1.0 - rf - 1e-12)):
            cand.append(m)
            if len(cand) >= 8:
                break
    if not cand:
        return None
    # ---- exact confirmation with Fractions
    for m in cand:
        ok = True
        maxL = None
        minR = None
        for t in range(nverts):
            L = None
            R = None
            for j in range(d):
                cj = lo[j] if not (t >> j) & 1 else hi[j]
                a = (cj - m[j] - r) / v[j + 1]
                b = (cj - m[j] + r) / v[j + 1]
                if L is None or a > L:
                    L = a
                if R is None or b < R:
                    R = b
            if L > R:
                ok = False
                break
            if maxL is None or L > maxL:
                maxL = L
            if minR is None or R < minR:
                minR = R
        if ok and ((maxL <= r and minR >= 0) or
                   (maxL <= 1 and minR >= 1 - r)):
            return list(m)
    return None


def _iter_combos(msc):
    """Iterate over the product of integer ranges (itertools-free)."""
    if not msc:
        yield ()
        return
    idx = [0] * len(msc)
    lens = [len(rng) for rng in msc]
    if any(L == 0 for L in lens):
        return
    while True:
        yield tuple(rng[i] for rng, i in zip(msc, idx))
        p = len(msc) - 1
        while p >= 0:
            idx[p] += 1
            if idx[p] < lens[p]:
                break
            idx[p] = 0
            p -= 1
        if p < 0:
            return


# --------------------------------------------- triple-identity band certificate

def _affine_range(box, coeffs, const):
    """Exact range of an affine form over a box (interval arithmetic)."""
    lo = const
    hi = const
    for j, b in coeffs.items():
        if b > 0:
            lo += b * box[j][0]
            hi += b * box[j][1]
        elif b < 0:
            lo += b * box[j][1]
            hi += b * box[j][0]
    return lo, hi


def _identity_triple(box, v, r, triple, EXACT):
    """One triple of runners; returns True if the identity certificate
    holds (exactly, if EXACT else float-gated).  Conditions:
    (i) each pair-crossing sigma_i(c) in [0,1) over the box;
    (ii) branch validity: |form_i| in [0, 1/2], sign-consistent;
    (iii) every non-pair runner term and the sigma-term <= r at each dip;
    (iv) the weighted identity constant S <= (sum w) * r, where the
    c-coefficients of sum w_i r_i cancel exactly."""
    x, y, z = triple
    vx, vy, vz = v[x + 1], v[y + 1], v[z + 1]
    pairs = [((x, y), vz * (vx + vy)),
             ((x, z), vy * (vx + vz)),
             ((y, z), vx * (vy + vz))]
    conv = (lambda t: t) if EXACT else (lambda t: float(t))
    half = conv(F(1, 2)) if EXACT else 0.5
    rf = conv(r)
    options = []
    for (a, b), w in pairs:
        va = conv(v[a + 1])
        vb = conv(v[b + 1])
        V = va + vb
        lo_a, hi_a = conv(box[a][0]), conv(box[a][1])
        lo_b, hi_b = conv(box[b][0]), conv(box[b][1])
        P = lo_a + lo_b
        Q = hi_a + hi_b
        opts = []
        # M = m_a + m_b: sigma = (c_a + c_b - M)/V in [0, 1] for all c
        # (sigma = 1 is sigma = 0 on the torus: the closed end is allowed)
        Mlo = Q - V                     # M >= Mlo
        Mhi = P                         # M <= Mhi
        Mi = math.floor(Mlo)
        while Mi <= Mhi:
            M = Mi
            if (P - M) >= 0 and (Q - M) <= V:
                # base range of (vb*c_a - va*c_b)
                bvals = [vb * s_a - va * s_b for s_a in (lo_a, hi_a)
                         for s_b in (lo_b, hi_b)]
                base_lo, base_hi = min(bvals), max(bvals)
                # m_a: form = (base + va*M - V*m_a)/V in [-1/2, 1/2]
                # -> V*m_a in [base_lo + va*M - V/2, base_hi + va*M + V/2]
                mlo = math.floor((base_lo + va * M - V * half) / V)
                mhi = math.floor((base_hi + va * M + V * half) / V)
                for ma in range(mlo, mhi + 1):
                    f_lo = base_lo + va * M - V * ma
                    f_hi = base_hi + va * M - V * ma
                    if f_lo >= 0 and f_hi <= V * half:
                        s = 1
                    elif f_hi <= 0 and f_lo >= -V * half:
                        s = -1
                    else:
                        continue
                    opts.append((M, ma, s))
            Mi += 1
        if not opts:
            return False
        options.append((opts, (a, b), w, va, vb, V))
    # combine one option per pair
    import itertools as _it
    for combo in _it.product(*[o[0] for o in options]):
        # affine form of S = sum w_i * s_i * form_i; require c-coeffs vanish
        coeffs = {}
        S = 0
        ok = True
        for (M, ma, s), (opts_, (a, b), w, va, vb, V) in zip(combo, options):
            mb = M - ma
            # form_i = (vb*c_a - va*c_b + va*mb - vb*ma)/V
            for j, cf in ((a, vb / V * s * w), (b, -va / V * s * w)):
                coeffs[j] = coeffs.get(j, 0) + cf
            S += (va * mb - vb * ma) / V * s * w
        for j, cf in coeffs.items():
            if abs(cf) > (0 if EXACT else 1e-9):
                ok = False
                break
        if not ok:
            continue
        tol = 0.0 if EXACT else 1e-9
        if S > (options[0][2] + options[1][2] + options[2][2]) * rf + tol:
            continue
        # (iii) non-binding checks per pair
        all_ok = True
        for (M, ma, s), (opts_, (a, b), w, va, vb, V) in zip(combo, options):
            # sigma-term: sup |sigma| <= r  over the box
            sig_coeffs = {a: conv(1) / V, b: conv(1) / V}
            sig_const = -conv(M) / V
            slo, shi = _affine_range_f(box, sig_coeffs, sig_const, EXACT)
            sup_sig = C.sup_normT(slo, shi) if EXACT else \
                _sup_normT_f(slo, shi)
            if sup_sig > (r if EXACT else rf):
                all_ok = False
                break
            # other runners j (c-coordinates) at this dip
            for j in range(len(v) - 1):
                if j in (a, b):
                    continue
                vj = conv(v[j + 1])
                cf = {j: conv(1), a: -vj / V, b: -vj / V}
                cst = vj * conv(M) / V
                lo2, hi2 = _affine_range_f(box, cf, cst, EXACT)
                sup_j = C.sup_normT(lo2, hi2) if EXACT else \
                    _sup_normT_f(lo2, hi2)
                if sup_j > (r if EXACT else rf):
                    all_ok = False
                    break
            if not all_ok:
                break
        if all_ok:
            return True
    return False


def _affine_range_f(box, coeffs, const, EXACT):
    if EXACT:
        return _affine_range(box, coeffs, const)
    lo = const
    hi = const
    for j, b in coeffs.items():
        if b > 0:
            lo += b * float(box[j][0])
            hi += b * float(box[j][1])
        elif b < 0:
            lo += b * float(box[j][1])
            hi += b * float(box[j][0])
    return lo, hi


def _sup_normT_f(lo, hi):
    import math as _m
    m0 = _m.floor(lo - 0.5) + 1
    if m0 < hi - 0.5:
        return 0.5
    return max(abs(lo - round(lo)), abs(hi - round(hi)))


def identity_certificate(box, v, r):
    """Triple-identity band certificate (sound).  Float-gated, then
    exact re-confirmation.  Returns True if max_C D <= r is certified."""
    import itertools as _it
    d = len(v) - 1
    if d < 3:
        return False
    for triple in _it.combinations(range(d), 3):
        try:
            if _identity_triple(box, v, r, triple, False):
                if _identity_triple(box, v, r, triple, True):
                    return True
        except Exception:
            continue
    return False


# ------------------------------------------------------------------ phase: B&B

def certify_v2(name, v, r, witnesses, budget_sec, out,
               max_nodes=6000000, max_depth=120, sizecap=F(1, 2 ** 30),
               lattice=None):
    """Sliced dyadic (or lattice-aligned) branch and bound; certificate:
    tree empties with every box pruned by the sound cascade (cheap grid
    sup / grid-Lipschitz / exact minimax U).  Deterministic pop order
    (-volume, depth, box).  When `lattice` = q is given, split points
    snap to the i/q grid of the instance's dip lattice, so boxes align
    with the rational slice edges of the argmax structure (dyadic
    splitting alone can strand boxes straddling non-dyadic edges)."""
    d = len(v) - 1
    r = F(r)
    rf = float(r)
    gd, g0d, gs, g0s = grids(rf)
    t0 = time.time()
    for w, wv in witnesses:
        got = C.D_at([F(x) for x in w], v)
        assert got == F(wv), ('witness mismatch', name, w, got, wv)
    K = 48
    grid_c = np.array([i / K * rf for i in range(K + 1)] +
                      [1 - rf + i / K * rf for i in range(K + 1)])
    # merge in the 470-division grid (the harmonic dip lattice: the three
    # tied dips of the argmax family sit at sigma = integer/235, the
    # c_2 band edges at integer/470); points kept only inside S.
    q = 470
    gq = np.array([i / q for i in range(q + 1)])
    gq = gq[(gq <= rf + 1e-12) | (gq >= 1 - rf - 1e-12)]
    grid_c = np.unique(np.concatenate([grid_c, gq]))
    g0_c = np.abs(grid_c - np.round(grid_c))

    st = load_state()
    if st and st.get('bb', {}).get('name') == name:
        B = st['bb']
        cnt = B['cnt']
        cnt.setdefault('gap', 0)
        cnt.setdefault('identity', 0)
        survivors = [box_from_str(s) for s in B['survivors']]
        heap = [(-e[0], e[1], box_from_str(e[2])) for e in B['heap']]
        heapq.heapify(heap)
        deepest = B['deepest']
        elapsed = B['elapsed']
        prune_sample = B['prune_sample']
        budget_hit = B['budget_hit']
        expiry_i = B['expiry_i']
        stage = B['stage']
    else:
        cnt = {'cheap': 0, 'gridlips': 0, 'gap': 0, 'identity': 0,
               'exactU': 0, 'splits': 0, 'pops': 0}
        survivors = []
        heap = [(-1.0, 0, tuple((F(0), F(1)) for _ in range(d)))]
        deepest = 0
        elapsed = 0.0
        prune_sample = []
        budget_hit = False
        expiry_i = 0
        stage = 'tree'
    t0adj = time.time() - elapsed

    def save():
        st2 = load_state() or {}
        hlist = [[-vol, dep, box_to_str(b)] for vol, dep, b in heap]
        if len(hlist) > HEAP_CAP:
            hlist.sort(key=lambda e: -e[0])          # largest volume first
            tail = len(hlist) - HEAP_CAP
            hlist = hlist[:HEAP_CAP]
            cnt['surv_overflow'] = cnt.get('surv_overflow', 0) + tail
        slist = [box_to_str(b) for b in survivors]
        if len(slist) > SURV_CAP:
            cnt['surv_overflow'] = (cnt.get('surv_overflow', 0)
                                    + len(slist) - SURV_CAP)
            slist = slist[:SURV_CAP]
        st2['bb'] = {
            'name': name, 'v': list(v), 'r': str(r),
            'budget_sec': budget_sec, 'cnt': cnt,
            'survivors': slist,
            'heap': hlist,
            'deepest': deepest, 'elapsed': time.time() - t0adj,
            'prune_sample': prune_sample, 'budget_hit': budget_hit,
            'expiry_i': expiry_i, 'stage': stage}
        dump_state(st2)
        print('[%s] slice save: stage=%s pops=%d heap=%d surv=%d(+%d) '
              'elapsed=%.0fs' % (name, stage, cnt['pops'], len(hlist),
                                 len(slist), cnt.get('surv_overflow', 0),
                                 time.time() - t0adj))

    def maybe_sample(box, which):
        if len(prune_sample) < 32 or cnt[which] % 64 == 0:
            prune_sample.append([which, box_to_str(box)])

    def test_and_maybe_prune(box):
        """Sound cascade.  True -> certified max_C D <= r."""
        G = G_box_sup(box, v, grid_c, g0_c)
        if G.min() <= rf - 1e-9:
            cnt['cheap'] += 1
            maybe_sample(box, 'cheap')
            return True
        c0 = [(lo + hi) / 2 for lo, hi in box]
        h = max(hi - lo for lo, hi in box) / 2
        Gc = G_point_batch(np.array([[float(x) for x in c0]]),
                           v, gd, g0d)[0]
        if float(Gc) + float(h) <= rf - 1e-9:
            cnt['gridlips'] += 1
            maybe_sample(box, 'gridlips')
            return True
        mgap = gap_certificate(box, v, r)
        if mgap is not None:
            cnt['gap'] += 1
            if len(prune_sample) < 32 or cnt['gap'] % 64 == 0:
                prune_sample.append(['gap', box_to_str(box), mgap])
            return True
        if identity_certificate(box, v, r):
            cnt['identity'] += 1
            if len(prune_sample) < 32 or cnt['identity'] % 64 == 0:
                prune_sample.append(['identity', box_to_str(box)])
            return True
        uf, _ = C.U_box(box, v, False, r=r)
        if uf is not None and uf <= rf + 1e-7:
            ue, _ = C.U_box(box, v, True, r=r)
            if ue is not None and ue <= r:
                cnt['exactU'] += 1
                return True
        return False

    def split_at(lo, hi):
        if lattice:
            q = lattice
            i_lo = math.floor(lo * q) + 1
            i_hi = math.ceil(hi * q) - 1
            if i_lo <= i_hi:
                i = min(max(round(float((lo + hi) / 2) * q), i_lo), i_hi)
                return F(i, q)
        return (lo + hi) / 2

    while stage == 'tree':
        if not heap:
            break
        if (time.time() - t0adj) > budget_sec or cnt['pops'] > max_nodes:
            budget_hit = True
            stage = 'expiry'
            break
        if time.time() - t0 > SLICE_SEC:
            save(); return None
        _, depth, box = heapq.heappop(heap)
        cnt['pops'] += 1
        deepest = max(deepest, depth)
        widths = [hi - lo for lo, hi in box]
        if test_and_maybe_prune(box):
            continue
        if max(widths) <= sizecap or depth >= max_depth:
            if len(survivors) < SURV_CAP:
                survivors.append(box)
            else:
                cnt['surv_overflow'] = cnt.get('surv_overflow', 0) + 1
            continue
        k = max(range(d), key=lambda i: widths[i])
        lo, hi = box[k]
        mid = split_at(lo, hi)
        b1 = box[:k] + ((lo, mid),) + box[k + 1:]
        b2 = box[:k] + ((mid, hi),) + box[k + 1:]
        vol = float(widths[k])
        for i in range(d):
            if i != k:
                vol *= float(widths[i])
        cnt['splits'] += 1
        heapq.heappush(heap, (-vol, depth + 1, b1))
        heapq.heappush(heap, (-vol, depth + 1, b2))

    if stage == 'tree':                 # heap emptied naturally
        stage = 'done'
    elif stage == 'expiry':
        # test expired heap boxes before calling them survivors
        while expiry_i < len(heap):
            if time.time() - t0 > SLICE_SEC:
                save(); return None
            _, _, box = heap[expiry_i]
            if (time.time() - t0adj) < budget_sec + 300 and \
                    test_and_maybe_prune(box):
                pass
            elif len(survivors) < SURV_CAP:
                survivors.append(box)
            else:
                cnt['surv_overflow'] = cnt.get('surv_overflow', 0) + 1
            expiry_i += 1
        heap = []
        stage = 'done'

    certified = not survivors and cnt.get('surv_overflow', 0) == 0
    res = {
        'name': name, 'v': list(v), 'r': str(r), 'certified': certified,
        'budget_hit': budget_hit,
        'survivors': len(survivors) + cnt.get('surv_overflow', 0),
        'counts': cnt, 'deepest_depth': deepest,
        'wall_sec': round(time.time() - t0adj, 1),
        'verdict': ('CERTIFIED rho <= %s' % r) if certified else
                   ('NOT certified; %d survivors' %
                    (len(survivors) + cnt.get('surv_overflow', 0))),
        'survivor_boxes': [box_to_str(b) for b in survivors[:3000]],
        'n_survivor_boxes_stored': len(survivors),
    }
    out[name] = res
    dump_out(out)
    if survivors:
        json.dump([box_to_str(b) for b in survivors], open(OUT_SURV, 'w'))
    st2 = load_state() or {}
    st2.pop('bb', None)
    st2['bb_done_' + name] = {'certified': certified,
                              'prune_sample': prune_sample}
    dump_state(st2)
    print('[%s] %s | pops=%d cheap=%d gridlips=%d exactU=%d surv=%d %.1fs'
          % (name, res['verdict'], cnt['pops'], cnt['cheap'],
             cnt['gridlips'], cnt['exactU'], len(survivors),
             res['wall_sec']))
    return res


# ------------------------------------------------------------------ phase drivers

def phase_hcert(out):
    if 'hcert' in out:
        return out['hcert']
    hs = out.get('hsearch') or {}
    best_D = F(hs.get('best_D') or '4183/12288')
    best_c = [F(x) for x in (hs.get('best_c') or
                             ['465/4096', '605/1024', '1019/4096',
                              '315/4096'])]
    if best_D < F(4183, 12288):
        best_D = F(4183, 12288)
        best_c = [F(465, 4096), F(605, 1024), F(1019, 4096), F(315, 4096)]
    # The balance-family witness discovered by the structural analysis of
    # the search summit: three tied dips (speed pairs (3,5), (4,5), (3,4))
    # with the weighted identity 32 r_A + 27 r_B + 35 r_C = 32 force the
    # tie value 16/47; witness c = (230, 92, 231, 175)/235 (k = 35 of the
    # family c(k) = (2k+160, 3k-13, 4k+91, 5k)/235, k = 20..35).
    FAM_D = F(16, 47)
    FAM_C = [F(230, 235), F(92, 235), F(231, 235), F(175, 235)]
    if FAM_D > best_D:
        best_D, best_c = FAM_D, FAM_C
        out.setdefault('hsearch', {})['family_witness'] = {
            'D': '16/47', 'c': [str(x) for x in FAM_C],
            'family': 'c(k) = (2k+160, 3k-13, 4k+91, 5k)/235, k=20..35',
            'sigmas': ['(k-31)/235', '(k-16)/235', '(k+179)/235'],
            'identity': '32 r_(3,5) + 27 r_(4,5) + 35 r_(3,4) = 32',
            'note': ('the three dips tie only at 16/47; witness verified '
                     'exactly (D_at == 16/47)')}
        dump_out(out)
    return certify_v2('hcert', VH, best_D, [(best_c, best_D)], 21600, out,
                      lattice=470)


def phase_n6m9(out):
    if 'n6m9' not in out:
        r = light_search(out, 'n6m9', V69, R6,
                         [[F(0), F(1, 2), F(0), F(0), F(0)]])
        if r is None:
            return None
    if 'n6m9cert' in out:
        return out['n6m9cert']
    refuted = (out['n6m9'] or {}).get('refuted')
    if refuted:
        out['n6m9cert'] = {'name': 'n6m9cert', 'v': list(V69),
                           'r': str(R6), 'certified': False,
                           'refuted_by_search': refuted,
                           'verdict': 'REFUTED: deeper witness D=%s at c=%s'
                                      % (refuted['D'], refuted['c'])}
        dump_out(out)
        print('[n6m9cert]', out['n6m9cert']['verdict'])
        return out['n6m9cert']
    return certify_v2('n6m9cert', V69, R6,
                      [([F(0), F(1, 2), F(0), F(0), F(0)], R6)],
                      7200, out)


def phase_n6m12(out):
    if 'n6m12' not in out:
        r = light_search(out, 'n6m12', V612, R6,
                         [[F(1, 16), F(11, 16), F(1, 8), F(7, 16),
                           F(7, 8)]])
        if r is None:
            return None
    if 'n6m12cert' in out:
        return out['n6m12cert']
    refuted = (out['n6m12'] or {}).get('refuted')
    if refuted:
        out['n6m12cert'] = {'name': 'n6m12cert', 'v': list(V612),
                            'r': str(R6), 'certified': False,
                            'refuted_by_search': refuted,
                            'verdict': 'REFUTED: deeper witness D=%s at c=%s'
                                       % (refuted['D'], refuted['c'])}
        dump_out(out)
        print('[n6m12cert]', out['n6m12cert']['verdict'])
        return out['n6m12cert']
    return certify_v2('n6m12cert', V612, R6,
                      [([F(1, 16), F(11, 16), F(1, 8), F(7, 16),
                         F(7, 8)], R6)],
                      7200, out)


def phase_n5b(out):
    rows = [('n5m6', V56, [F(31, 32), F(17, 32), F(15, 16), F(1, 16)]),
            ('n5m7', V57, [F(5, 32), F(23, 32), F(5, 16), F(17, 32)]),
            ('n5m8', V58, [F(7, 8), F(7, 16), F(7, 8), F(3, 4)])]
    for key, v, w in rows:
        if key not in out:
            r = light_search(out, key, v, R5, [w])
            if r is None:
                return None
        if key + 'cert' not in out:
            refuted = (out[key] or {}).get('refuted')
            if refuted:
                out[key + 'cert'] = {
                    'name': key + 'cert', 'v': list(v), 'r': str(R5),
                    'certified': False, 'refuted_by_search': refuted,
                    'verdict': 'REFUTED: deeper witness D=%s at c=%s'
                               % (refuted['D'], refuted['c'])}
                dump_out(out)
                print('[%scert]' % key, out[key + 'cert']['verdict'])
                continue
            r = certify_v2(key + 'cert', v, R5, [(w, R5)], 14400, out)
            if r is None:
                return None
    return True


# ------------------------------------------------------------------ phase: verify

def phase_verify(out):
    """Exact re-confirmation of sampled stage-0/1 prunes of every
    completed B&B run + exact survivor-center analysis of open runs."""
    t0 = time.time()
    st = load_state() or {}
    res = {}
    for key in [k for k in st if k.startswith('bb_done_')]:
        name = key[8:]
        if name not in out:
            continue
        rec = st[key]
        v = tuple(out[name]['v'])
        r = F(out[name]['r'])
        n_ok = 0
        n_bad = 0
        for entry in rec['prune_sample']:
            which = entry[0]
            sbox = entry[1]
            box = box_from_str(sbox)
            if which == 'cheap':
                # re-confirm against the SAME merged grid the prune used
                rf = float(r)
                K = 48
                grid_c = np.array([i / K * rf for i in range(K + 1)] +
                                  [1 - rf + i / K * rf
                                   for i in range(K + 1)])
                q = 470
                gq = np.array([i / q for i in range(q + 1)])
                gq = gq[(gq <= rf + 1e-12) | (gq >= 1 - rf - 1e-12)]
                grid_c = np.unique(np.concatenate([grid_c, gq]))
                G = G_box_sup(box, v, grid_c,
                              np.abs(grid_c - np.round(grid_c)))
                j = int(np.argmin(G))
                s = _grid_sigma(grid_c[j], r)
                vals = [C.normT(s)] + [
                    C.sup_normT(box[i][0] - vj * s, box[i][1] - vj * s)
                    for i, vj in enumerate(v[1:])]
                if max(vals) <= r:
                    n_ok += 1
                else:
                    n_bad += 1
            else:
                # gridlips / gap / identity samples
                c0 = [(lo + hi) / 2 for lo, hi in box]
                h = max(hi - lo for lo, hi in box) / 2
                if which == 'gap':
                    mg = gap_certificate(box, v, r)
                    if mg is not None:
                        n_ok += 1
                    else:
                        n_bad += 1
                elif which == 'identity':
                    if identity_certificate(box, v, r):
                        n_ok += 1
                    else:
                        n_bad += 1
                elif C.D_at(c0, v) + h <= r:
                    n_ok += 1
                else:
                    n_bad += 1
        res[name] = {'sampled_prunes': len(rec['prune_sample']),
                     'exact_confirmed': n_ok, 'exact_failed': n_bad}
        assert n_bad == 0, ('VERIFY FAILURE', name)
    # survivor deep analysis for every open run
    for name in [k for k in out if isinstance(out.get(k), dict)
                 and out[k].get('certified') is False
                 and 'survivor_boxes' in out[k]
                 and not out[k].get('refuted_by_search')]:
        v = tuple(out[name]['v'])
        r = F(out[name]['r'])
        boxes = out[name]['survivor_boxes']
        rows = []
        tA = time.time()
        boxes = sorted(boxes, key=lambda sb: -math.prod(
            [float(F(s.split('|')[1]) - F(s.split('|')[0])) for s in sb]))
        for i, sb in enumerate(boxes):
            if i > 2000 and time.time() - tA > 300:
                break
            if i > 400 and time.time() - t0 > SLICE_SEC:
                break
            box = box_from_str(sb)
            c0 = [(lo + hi) / 2 for lo, hi in box]
            row = {'box': sb, 'center_D': str(C.D_at(c0, v))}
            if i < 100:
                ue, _ = C.U_box(box, v, True, r=r)
                row['exact_U'] = str(ue)
            rows.append(row)
        out[name]['survivor_rows'] = rows
        depths = [F(x['center_D']) for x in rows]
        out[name]['max_center_D'] = str(max(depths)) if depths else None
        us = [F(x['exact_U']) for x in rows if 'exact_U' in x]
        out[name]['enclosure_upper'] = (
            str(max([r] + us)) if us else str(r))
        res[name + '_survivor_analysis'] = {
            'boxes_analyzed': len(rows),
            'max_center_D': out[name]['max_center_D'],
            'enclosure_upper': out[name]['enclosure_upper']}
    out['verify'] = res
    dump_out(out)
    print('[verify]', json.dumps(res)[:800])
    return res


def _grid_sigma(x, r):
    """Rational reconstruction of a merged-grid float point in S."""
    q = 470
    i = round(x * q)
    if abs(x * q - i) < 1e-9 and 0 <= i <= q:
        return F(i, q)
    rf = float(r)
    if rf > 0 and x <= rf + 1e-12:
        i = round(x / rf * 48)
        return F(i, 48) * F(r)
    i = round((x - (1 - rf)) / rf * 48)
    return (1 - F(r)) + F(i, 48) * F(r)


# ------------------------------------------------------------------ phase: analyze

def phase_analyze(out):
    res = {}
    hs = out.get('hsearch', {})
    res['harmonic_search'] = {
        'best_D': hs.get('best_D'), 'best_c': hs.get('best_c'),
        'seed_best_D': hs.get('seed_best_D'),
        'climb_evals': hs.get('climb_evals'),
        'n_refined': hs.get('n_refined')}
    hc = out.get('hcert', {})
    if hc:
        if hc.get('certified'):
            wv = F(hc['r'])
            fw = (out.get('hsearch') or {}).get('family_witness') or {}
            if F(fw.get('D') or '0') == wv:
                wc = [F(x) for x in fw['c']]
            else:
                wc = [F(x) for x in (hs.get('best_c') or [])]
            got = C.D_at(wc, VH) if len(wc) == 4 else None
            fam = fw.get('family')
            res['harmonic'] = {
                'value': str(wv), 'certified': True,
                'witness_reattains': (got == wv),
                'witness_c': [str(x) for x in wc] if len(wc) == 4 else None,
                'witness_family': fam,
                'identity': fw.get('identity'),
                'counts': hc.get('counts'), 'wall_sec': hc.get('wall_sec'),
                'statement': ('rho(1,2,3,4,5) = %s exactly: the '
                              'branch-and-bound certified rho <= %s and '
                              'the balance-family witness attains it.' %
                              (wv, wv))}
        else:
            res['harmonic'] = {
                'certified': False,
                'lower_bound': hs.get('best_D'),
                'max_center_D': hc.get('max_center_D'),
                'enclosure_upper': hc.get('enclosure_upper'),
                'survivors': hc.get('survivors')}
    for key, vname, r in [('n6m9cert', 'rho(1,2,3,4,5,9)', R6),
                          ('n6m12cert', 'rho(1,2,3,4,5,12)', R6),
                          ('n5m6cert', 'rho(1,2,3,4,6)', R5),
                          ('n5m7cert', 'rho(1,2,3,4,7)', R5),
                          ('n5m8cert', 'rho(1,2,3,4,8)', R5)]:
        rec = out.get(key)
        if not rec:
            continue
        if rec.get('refuted_by_search'):
            res[key] = {'v': rec.get('v'), 'refuted': True,
                        'witness': rec['refuted_by_search']}
        elif rec.get('certified'):
            res[key] = {
                'v': rec.get('v'), 'value': str(r), 'certified': True,
                'counts': rec.get('counts'),
                'wall_sec': rec.get('wall_sec'),
                'statement': ('%s = %s exactly (boundary row resolved as '
                              'equality: the tree emptied at r = %s with '
                              'the recorded witness attaining it).' %
                              (vname, r, r))}
        else:
            res[key] = {'v': rec.get('v'), 'certified': False,
                        'max_center_D': rec.get('max_center_D'),
                        'enclosure_upper': rec.get('enclosure_upper')}
    res['verify'] = out.get('verify')
    # sigma records for the paper's hand checks
    sig = {}
    for nm, v, w, wv in WITNESSES:
        s, val = argmin_sigma(list(w), v)
        assert val == wv, ('argmin sigma mismatch', nm, val, wv)
        sig[nm] = {'sigma': str(s), 'D': str(val)}
    for key in ('n6m9cert', 'n6m12cert', 'n5m6cert', 'n5m7cert', 'n5m8cert'):
        rec = out.get(key)
        if rec and rec.get('refuted_by_search'):
            c = [F(x) for x in rec['refuted_by_search']['c']]
            s, val = argmin_sigma(c, tuple(rec['v']))
            sig[key + '-refutation'] = {'sigma': str(s), 'D': str(val)}
    if hs.get('best_c'):
        c = [F(x) for x in hs['best_c']]
        s, val = argmin_sigma(c, VH)
        sig['harmonic-argmax'] = {'sigma': str(s), 'D': str(val)}
    fw = (out.get('hsearch') or {}).get('family_witness') or {}
    if fw.get('c'):
        c = [F(x) for x in fw['c']]
        s, val = argmin_sigma(c, VH)
        sig['harmonic-family-witness'] = {'sigma': str(s), 'D': str(val)}
        # the three tied dips of the family witness (hand-check record)
        dips = {}
        for sname, sguess in [('A', F(4, 235)), ('B', F(19, 235)),
                               ('C', F(214, 235))]:
            terms = [C.normT(sguess)] + [
                C.normT(c[j] - v * sguess)
                for j, v in enumerate((2, 3, 4, 5))]
            dips[sname] = {'sigma': str(sguess),
                           'terms_over_235': [str(t * 235)
                                              for t in terms]}
        sig['harmonic-family-dips'] = dips
    res['sigma_records'] = sig
    out['analysis'] = res
    dump_out(out)
    print('[analyze]', json.dumps(res, default=str)[:1200])
    return res


def phase_status(out):
    st = load_state()
    print('phases done:', [k for k in out])
    if st:
        for k, v in st.items():
            if k.endswith('_s') or k == 'hsearch':
                print(k, ': stage=%s chunk=%s top=%d'
                      % (v['stage'], v.get('chunk'), len(v.get('top', []))))
            elif k.startswith('bb_done_'):
                print(k, ': certified=%s sample=%d'
                      % (v['certified'], len(v['prune_sample'])))
            elif k == 'bb':
                print('bb state: name=%s stage=%s pops=%d heap=%d surv=%d '
                      'elapsed=%.0f' % (v['name'], v['stage'],
                                        v['cnt']['pops'], len(v['heap']),
                                        len(v['survivors']), v['elapsed']))


PHASES = {'selftest': phase_selftest, 'hsearch': phase_hsearch,
          'hcert': phase_hcert, 'n6m9': phase_n6m9,
          'n6m12': phase_n6m12, 'n5b': phase_n5b,
          'verify': phase_verify, 'analyze': phase_analyze,
          'status': phase_status}

if __name__ == '__main__':
    args = sys.argv[1:] or ['all']
    for a in args:
        if a == 'all':
            for p in ['selftest', 'hsearch', 'hcert', 'n6m9', 'n6m12',
                      'n5b', 'verify', 'analyze']:
                out = load_out()
                r = PHASES[p](out)
                if r is None:
                    print('--- phase %s needs more slices; re-invoke' % p)
                    break
        else:
            out = load_out()
            PHASES[a](out)
