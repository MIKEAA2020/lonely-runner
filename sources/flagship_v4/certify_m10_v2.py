#!/usr/bin/env python3
"""The second, audit-driven certification attempt for v = (1,2,3,4,10),
target rho = 7/22 (the m=10 row of the n=5 deepest-hole ladder).

WHY A SECOND ATTEMPT.  The audit (audits/flaghship.txt, opus points on the
harmonic hole bound and on m=10) criticized the first attempt
(out_certify_rung2.json) twice:  (a) the search evidence behind the
"consistent" label was thin -- 288 Nelder-Mead starts, 7.4 s -- while 40
random local searches sufficed to beat the published harmonic lower bound;
(b) the "quantized minimax loss" reading of the 21,153 unprocessed residual
boxes was unproved: residual volume below 6e-7 does not rule out a deeper
hole.

THIS SCRIPT answers both.

SEARCH UPGRADE (phase `search`).  4.2M uniform random points, the full
{k/16}^4 and {k/32}^4 dyadic grids, a 500k-point torsion neighborhood, and
2,048 batched random uphill climbs (19 sigma levels, 8 directions), all
screened on a fixed 2x2048-point s-grid over S = [0,r] u [1-r,1) (sound:
the grid minimum is an upper bound for D, so no point with D > r can be
missed), then Nelder-Mead polish of the top summits and EXACT rational
re-evaluation of every candidate near or above 7/22.  Any exact D > 7/22
refutes the row and is recorded as a witness.

THE LIPSCHITZ PRUNE (phases `negctl`, `n4`, `target`).  D is 1-Lipschitz in
the sup norm: for fixed s every term ||c_j - v_j s||_T moves by at most
|c_j - c'_j| while the ||s||_T term does not move, so G(c,s) <= G(c',s) +
||c-c'||_inf for every s; minimizing over s gives
|D(c) - D(c')| <= ||c-c'||_inf.  Hence for a dyadic c-box C with center c0
and half-width h = max_j (hi_j - lo_j)/2:

    max_{c in C} D(c) <= D(c0) + h,

and C is CERTIFIED (max_C D <= r) whenever D(c0) + h <= r.  At the residual
width of the first attempt (1/512 per side, h = 1/1024) this closes the
observed gap: 1135/3584 + 1/1024 = 0.317661... < 7/22 = 0.318181... .

Branch and bound cascade per box (every step sound):
  0. cheap grid prune (the first attempt's 97-point S-grid): some grid s
     has sup_{c in C} max_j f_j(c,s) <= r - 1e-9;
  1. grid-Lipschitz prune: min over a 2x8192-point float s-grid of
     G(c0, s) is an upper bound for D(c0) (grid points are reals in S);
     if bound + h <= r - 1e-9 the box is certified;
  2. float U(C) <= r + 1e-7 -> exact rational U(C) <= r (the first
     attempt's minimax box bound);
  3. split the widest coordinate.
Stages 0 and 1 are float-implemented with >= 5 orders of magnitude of
numerical slack; a 1-in-64 sample of their prunes (plus the first 32 of
each) is re-confirmed exactly in phase `verify`.  Size-capped and
budget-expired boxes are TESTED before being called survivors (the first
attempt declared 21,153 heap boxes survivors untested -- exactly what the
audit flagged).  Full survivor lists are persisted.

Controls: the harmonic negative control v = (1,2,3,4,5), r = 1/3 must NOT
certify; boxes around the three known refuting witnesses (607/1792,
4183/12288, 6/19) must REFUSE the Lipschitz prune; the previously certified
n=4 instances must re-certify under the new cascade.

SLICING.  The sandbox kills background processes, so every phase runs in
~460 s foreground slices with a checkpoint state
(scripts/out_certify_m10_v2_state.json); re-invoking the same phase resumes
exactly where the previous slice stopped.  The branch-and-bound pop
sequence is slice-independent (heap priority is fully determined by
(-volume, depth, box)), so a completed certificate is deterministic.

Output: scripts/out_certify_m10_v2.json (+ out_certify_m10_v2_survivors.json)
Usage: python3 certify_m10_v2.py selftest|search|negctl|n4|target|verify|analyze|status
"""
import sys, os, json, math, time, heapq
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F
import numpy as np

import certify_rung2_n5 as C          # validated primitives: U_box, D_at

OUT = '/home/z/my-project/scripts/out_certify_m10_v2.json'
OUT_SURV = '/home/z/my-project/scripts/out_certify_m10_v2_survivors.json'
STATE = '/home/z/my-project/scripts/out_certify_m10_v2_state.json'

V10 = (1, 2, 3, 4, 10)
R10 = F(7, 22)
CHUNK = 512                            # rows per numpy batch (3 GB RAM box)
SLICE_SEC = 460
HEAP_CAP = 80000                       # persisted heap entries
SURV_CAP = 30000                       # persisted survivor boxes


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

_GRIDS = {}                            # rf -> (dense grid, g0, search grid, g0)


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


# ------------------------------------------------------------------ phase: selftest

def phase_selftest(out):
    if 'selftest' in out:
        return out['selftest']
    res = {}
    t0 = time.time()
    W = [
        (V10, [F(1, 2), F(0), F(1, 2), F(1, 2)], F(7, 22)),
        ((1, 2, 3, 4, 5), [F(11, 256), F(163, 256), F(15, 256), F(79, 256)],
         F(607, 1792)),
        ((1, 2, 3, 4, 5), [F(465, 4096), F(605, 1024), F(1019, 4096),
                           F(315, 4096)], F(4183, 12288)),
        ((1, 2, 3, 4, 5), [F(1, 32), F(15, 32), F(3, 32), F(7, 8)],
         F(97, 288)),
        ((1, 2, 3, 4, 15), [F(7, 8), F(5, 16), F(3, 4), F(9, 16)],
         F(6, 19)),
    ]
    for v, w, wv in W:
        got = C.D_at(list(w), v)
        assert got == wv, ('witness mismatch', v, w, got, wv)
    res['witnesses_exact'] = len(W)
    import random
    from lrc_zono_lib import make_forms, D_exact
    rng = random.Random(11)
    for v in (V10, (1, 2, 3, 4, 5)):
        forms = make_forms(v)
        for _ in range(40):
            c = tuple(F(rng.randrange(0, 96), 96)
                      for _ in range(len(v) - 1))
            a = C.D_at(list(c), v)
            b = D_exact(c, v, forms)
            assert a == b, ('D mismatch', v, c, a, b)
    res['cross_checks'] = '80 points, two instances, exact agreement'
    adv = [
        ((1, 2, 3, 4, 5), F(1, 3),
         [F(11, 256), F(163, 256), F(15, 256), F(79, 256)]),
        ((1, 2, 3, 4, 5), F(1, 3),
         [F(465, 4096), F(605, 1024), F(1019, 4096), F(315, 4096)]),
        ((1, 2, 3, 4, 15), F(5, 16),
         [F(7, 8), F(5, 16), F(3, 4), F(9, 16)]),
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
    pruned = 0
    for ctr in ([F(1, 4)] * 4, [F(1, 8), F(3, 8), F(1, 8), F(3, 8)],
                [F(3, 4), F(1, 16), F(3, 4), F(1, 16)]):
        box = [(x, x + F(1, 512)) for x in ctr]
        c0 = [(lo + hi) / 2 for lo, hi in box]
        if C.D_at(c0, V10) + F(1, 1024) <= R10:
            pruned += 1
    assert pruned >= 2, ('positive Lipschitz control too weak', pruned)
    res['positive_controls_pruned'] = pruned
    res['sec'] = round(time.time() - t0, 1)
    out['selftest'] = res
    dump_out(out)
    print('[selftest] OK:', json.dumps(res)[:500])
    return res


# ------------------------------------------------------------------ phase: search

def _exact_check(cands, v, tag, res, budget_sec=900):
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
        if v == V10 and d > R10:
            refuted = {'c': [str(x) for x in cc], 'D': str(d)}
            break
    rows.sort(key=lambda t: -t[0])
    res['exact_checks_' + tag] = {
        'n': len(rows), 'deepest': rows[0][1] if rows else None,
        'deepest_c': rows[0][2] if rows else None, 'refuted': refuted}
    return rows, refuted


def phase_search(out):
    """Stage-machine search with checkpoint resume.  Returns None while
    incomplete (re-invoke to continue)."""
    if 'search' in out:
        print('[search] already done')
        return out['search']
    v = V10
    rf = float(R10)
    gd, g0d, gs, g0s = grids(rf)
    t0 = time.time()
    st = load_state()
    if not st:
        st = {}
    if 'search' not in st:
        st['search'] = {'stage': 'rand', 'chunk': 0, 'top': [],
                        'cand': [], 'best_rand': -1.0, 'extra': {}}
    S = st['search']
    rng = np.random.default_rng(20261008)

    def save():
        S['top'] = S['top'][:4096]
        S['cand'] = S['cand'][:8192]
        dump_state(st)
        print('[search] slice save: stage=%s chunk=%d elapsed=%.0fs'
              % (S['stage'], S.get('chunk', 0), time.time() - t0))

    def screen(P):
        Dg = G_point_batch(P, v, gs, g0s)
        for i in np.argsort(-Dg)[:32]:
            S['top'].append((float(Dg[i]), [float(x) for x in P[i]]))
        hit = Dg > rf - 2e-3
        if hit.any() and len(S['cand']) < 8192:
            for i in np.nonzero(hit)[0][:256]:
                S['cand'].append([float(x) for x in P[i]])
        return Dg

    while True:
        stage = S['stage']
        if stage == 'rand':
            N = 1 << 22
            while S['chunk'] < N // 2048:
                P = rng.random((2048, 4))
                Dg = screen(P)
                i = int(np.argmax(Dg))
                if Dg[i] > S['best_rand']:
                    S['best_rand'] = float(Dg[i])
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
            _exact_check([[F(int(round(x * q)), q) for x in P[i]]
                          for i in idx], v, 'dyadic16', S['extra'], 120)
            S['stage'] = 'dyadic32'; S['chunk'] = 0
            if time.time() - t0 > SLICE_SEC - 60:
                save(); return None
        elif stage == 'dyadic32':
            q = 32
            axis = np.arange(q) / q
            P = np.stack(np.meshgrid(*[axis] * 4, indexing='ij'),
                         -1).reshape(-1, 4)
            Dg = G_point_batch(P, v, gs, g0s)
            j = int(np.argmax(Dg))
            S['extra']['dyadic_q32_best_gridD'] = float(Dg[j])
            for i in np.argsort(-Dg)[:64]:
                S['top'].append((float(Dg[i]), [float(x) for x in P[i]]))
            S['stage'] = 'torsion'; S['chunk'] = 0
            if time.time() - t0 > SLICE_SEC - 180:
                save(); return None
        elif stage == 'torsion':
            cstar = np.array([0.5, 0.0, 0.5, 0.5])
            while S['chunk'] < 500000 // 2048:
                P = (cstar + rng.uniform(-0.02, 0.02, (2048, 4))) % 1.0
                screen(P)
                S['chunk'] += 1
                if S['chunk'] % 16 == 0 and time.time() - t0 > SLICE_SEC:
                    save(); return None
            S['stage'] = 'oldcenters'; S['chunk'] = 0
        elif stage == 'oldcenters':
            try:
                old = json.load(open(
                    '/home/z/my-project/scripts/out_certify_rung2.json'))
                run = [r for r in old['runs']
                       if r['name'] == 'n5-rung2-m10-TARGET'][0]
                cents = []
                for sbox, _ in run.get('survivor_sample', []):
                    c = []
                    for s in sbox:
                        lo, hi = s.split('|')
                        c.append((F(lo) + F(hi)) / 2)
                    cents.append(c)
                _exact_check(cents, v, 'old_survivor_centers',
                             S['extra'], 300)
                S['extra']['old_survivor_centers_checked'] = len(cents)
            except Exception as e:
                S['extra']['old_survivor_centers_error'] = str(e)
            S['stage'] = 'climbs'; S['chunk'] = 0
        elif stage == 'climbs':
            if 'X' not in S:
                S['top'].sort(key=lambda t: -t[0])
                seeds = [c for _, c in S['top'][:1536]]
                seeds += [[float(x) for x in rng.random(4)]
                          for _ in range(256)]
                cstar = np.array([0.5, 0.0, 0.5, 0.5])
                seeds += [[float(x) for x in
                           ((cstar + rng.uniform(-0.03, 0.03, 4)) % 1.0)]
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
                c0 = finals[S['chunk']]

                def negD(c):
                    return -floatD_point(c, v)

                try:
                    rr = minimize(negD, np.array(c0, dtype=float),
                                  method='Nelder-Mead',
                                  options={'xatol': 1e-11, 'fatol': 1e-13,
                                           'maxiter': 800})
                    S['polished'].append([float(x) for x in np.mod(rr.x, 1.0)])
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
                for c in [c for _, c in S['top'][:200]]:
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
            near = [c for fd, c in pairs if fd > rf - 1e-5][:1024]
            deep50 = [c for _, c in pairs[:50]]
            S['extra']['n_refined'] = len(pairs)
            S['extra']['refined_best_floatD'] = pairs[0][0] if pairs else None
            rows, refuted = _exact_check(near + deep50, v, 'search',
                                         S['extra'], 900)
            S['extra']['refuted'] = refuted
            S['stage'] = 'done'
        else:
            break

    res = {'screen_grid': '2 x 2048 s-points on S (sound upper bound)',
           'r': str(R10), 'random_points': 1 << 22,
           'torsion_neighborhood_points': 500000,
           'climb_starts': S['extra'].get('climb_starts'),
           'climb_evals': S['extra'].get('climb_evals'),
           'climb_best_gridD': S['extra'].get('climb_best_gridD'),
           'nm_polished': S['extra'].get('nm_polished'),
           'nm_best_floatD': S['extra'].get('nm_best_floatD'),
           'n_refined': S['extra'].get('n_refined'),
           'refined_best_floatD': S['extra'].get('refined_best_floatD'),
           'refuted': S['extra'].get('refuted'),
           'old_survivor_centers_checked':
               S['extra'].get('old_survivor_centers_checked'),
           'dyadic_q16_best_gridD':
               S['extra'].get('dyadic_q16_best_gridD'),
           'dyadic_q32_best_gridD':
               S['extra'].get('dyadic_q32_best_gridD')}
    for k in ('exact_checks_dyadic16', 'exact_checks_old_survivor_centers',
              'exact_checks_search'):
        if k in S['extra']:
            res[k] = S['extra'][k]
    res['random_best_gridD'] = S['best_rand']
    res['sec'] = 'multi-slice (see state history)'
    out['search'] = res
    dump_out(out)
    st.pop('search', None)
    dump_state(st)
    print('[search] COMPLETE:', json.dumps(res)[:700])
    return res


# ------------------------------------------------------------------ phase: B&B

def certify_v2(name, v, r, witnesses, budget_sec, out,
               max_nodes=4000000, max_depth=110, sizecap=F(1, 2 ** 30)):
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
    g0_c = np.abs(grid_c - np.round(grid_c))

    st = load_state()
    if st and st.get('bb', {}).get('name') == name:
        B = st['bb']
        cnt = B['cnt']
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
        cnt = {'cheap': 0, 'gridlips': 0, 'exactU': 0, 'splits': 0,
               'pops': 0}
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
        uf, _ = C.U_box(box, v, False, r=r)
        if uf is not None and uf <= rf + 1e-7:
            ue, _ = C.U_box(box, v, True, r=r)
            if ue is not None and ue <= r:
                cnt['exactU'] += 1
                return True
        return False

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
        mid = (lo + hi) / 2
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


def phase_negctl(out):
    if 'negctl' in out:
        return out['negctl']
    return certify_v2(
        'negctl', (1, 2, 3, 4, 5), F(1, 3),
        [([F(1, 2), F(0), F(1, 2), F(0)], F(1, 3)),
         ([F(11, 256), F(163, 256), F(15, 256), F(79, 256)], F(607, 1792)),
         ([F(465, 4096), F(605, 1024), F(1019, 4096), F(315, 4096)],
          F(4183, 12288))],
        300, out)


def phase_n4(out):
    if 'n4harm' not in out:
        r = certify_v2('n4harm', (1, 2, 3, 4), F(3, 10),
                       [([F(1, 2), F(0), F(1, 2)], F(3, 10))], 300, out)
        if r is None:
            return None
    if 'n4m16' not in out:
        r = certify_v2('n4m16', (1, 2, 3, 16), F(9, 34),
                       [([F(1, 2), F(0), F(1, 2)], F(9, 34))], 300, out)
        if r is None:
            return None
    ok = out['n4harm']['certified'] and out['n4m16']['certified']
    print('[n4] re-certified both:', ok)
    return ok


def phase_target(out):
    if 'target' in out:
        return out['target']
    return certify_v2(
        'target', V10, R10,
        [([F(1, 2), F(0), F(1, 2), F(1, 2)], F(7, 22))],
        5400, out)


# ------------------------------------------------------------------ phase: verify

def phase_verify(out):
    """Exact re-confirmation of the sampled stage-0/1 prunes of every
    completed B&B run + exact survivor-center analysis."""
    t0 = time.time()
    st = load_state() or {}
    res = {}
    for key in [k for k in st if k.startswith('bb_done_')]:
        name = key[8:]
        rec = st[key]
        v = tuple(out[name]['v'])
        r = F(out[name]['r'])
        n_ok = 0
        n_bad = 0
        for which, sbox in rec['prune_sample']:
            box = box_from_str(sbox)
            if which == 'cheap':
                G = G_box_sup(box, v, *grids(float(r))[2:])
                j = int(np.argmin(G))
                s = F(grids(float(r))[2][j])
                vals = [C.normT(s)] + [
                    C.sup_normT(box[i][0] - vj * s, box[i][1] - vj * s)
                    for i, vj in enumerate(v[1:])]
                if max(vals) <= r:
                    n_ok += 1
                else:
                    n_bad += 1
            else:
                c0 = [(lo + hi) / 2 for lo, hi in box]
                h = max(hi - lo for lo, hi in box) / 2
                if C.D_at(c0, v) + h <= r:
                    n_ok += 1
                else:
                    n_bad += 1
        res[name] = {'sampled_prunes': len(rec['prune_sample']),
                     'exact_confirmed': n_ok, 'exact_failed': n_bad}
        assert n_bad == 0, ('VERIFY FAILURE', name)
    # survivor deep analysis for the target (if open)
    if 'target' in out and not out['target']['certified']:
        v = V10
        boxes = json.load(open(OUT_SURV)) if os.path.exists(OUT_SURV) else []
        rows = []
        tA = time.time()
        boxes = sorted(boxes, key=lambda sb: -math.prod(
            [float(F(s.split('|')[1]) - F(s.split('|')[0])) for s in sb]))
        for i, sb in enumerate(boxes):
            if i > 2000 and time.time() - tA > 300:
                break
            if i > 500 and time.time() - t0 > SLICE_SEC:
                break
            box = box_from_str(sb)
            c0 = [(lo + hi) / 2 for lo, hi in box]
            row = {'box': sb, 'center_D': str(C.D_at(c0, v))}
            if i < 100:
                ue, _ = C.U_box(box, v, True, r=R10)
                row['exact_U'] = str(ue)
            rows.append(row)
        out['target']['survivor_rows'] = rows
        depths = [F(x['center_D']) for x in rows]
        out['target']['max_center_D'] = str(max(depths))
        us = [F(x['exact_U']) for x in rows if 'exact_U' in x]
        out['target']['enclosure_upper'] = str(max([R10] + us))
    out['verify'] = res
    dump_out(out)
    print('[verify]', json.dumps(res))
    return res


# ------------------------------------------------------------------ phase: analyze

def phase_analyze(out):
    res = {'r': str(R10), 'v': list(V10)}
    tgt = out.get('target', {})
    res['certified'] = tgt.get('certified')
    if tgt.get('certified'):
        wv = C.D_at([F(1, 2), F(0), F(1, 2), F(1, 2)], V10)
        assert wv == R10
        res['value_equals_r'] = True
        res['counts'] = tgt.get('counts')
        res['wall_sec'] = tgt.get('wall_sec')
        res['statement'] = ('rho(1,2,3,4,10) = 7/22 exactly: the '
                            'branch-and-bound certified rho <= 7/22 and '
                            'the torsion center (1/2,0,1/2,1/2) attains '
                            '7/22.')
    else:
        res['max_center_D'] = tgt.get('max_center_D')
        res['enclosure_upper'] = tgt.get('enclosure_upper')
    sr = out.get('search', {})
    res['search'] = {k: sr.get(k) for k in
                     ('random_points', 'climb_starts', 'climb_evals',
                      'nm_polished', 'refuted', 'n_refined')}
    neg = out.get('negctl', {})
    if neg:
        res['negctl_refused'] = not neg.get('certified', True)
        assert res['negctl_refused'], 'NEGATIVE CONTROL CERTIFIED -- BUG'
    else:
        res['negctl_refused'] = None
    res['n4_recertified'] = (out.get('n4harm', {}).get('certified') and
                             out.get('n4m16', {}).get('certified'))
    res['verify'] = out.get('verify')
    out['analysis'] = res
    dump_out(out)
    print('[analyze]', json.dumps(res, default=str)[:900])
    return res


def phase_status(out):
    st = load_state()
    print('phases done:', [k for k in out])
    if st:
        for k, v in st.items():
            if k == 'search':
                print('search state: stage=%s chunk=%s top=%d cand=%d'
                      % (v['stage'], v.get('chunk'), len(v.get('top', [])),
                         len(v.get('cand', []))))
            elif k.startswith('bb_done_'):
                print(k, ': certified=%s sample=%d'
                      % (v['certified'], len(v['prune_sample'])))
            elif k == 'bb':
                print('bb state: name=%s stage=%s pops=%d heap=%d surv=%d '
                      'elapsed=%.0f' % (v['name'], v['stage'],
                                        v['cnt']['pops'], len(v['heap']),
                                        len(v['survivors']), v['elapsed']))


PHASES = {'selftest': phase_selftest, 'search': phase_search,
          'negctl': phase_negctl, 'n4': phase_n4, 'target': phase_target,
          'verify': phase_verify, 'analyze': phase_analyze,
          'status': phase_status}

if __name__ == '__main__':
    args = sys.argv[1:] or ['all']
    for a in args:
        if a == 'all':
            for p in ['selftest', 'search', 'negctl', 'n4', 'target',
                      'verify', 'analyze']:
                out = load_out()
                r = PHASES[p](out)
                if r is None:
                    print('--- phase %s needs more slices; re-invoke' % p)
                    break
        else:
            out = load_out()
            PHASES[a](out)
