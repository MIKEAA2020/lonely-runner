#!/usr/bin/env python3
"""Close out the exact covering radii of the refutation instances and the
binding-structure (weight lemma) data, with persisted JSON evidence.

Method notes (established in-session):
- The pair-form characters chi_{ab}(y) = b*y_a - a*y_b (mod 1) are lattice-
  invariant circle maps on the quotient torus; they satisfy the linear
  identity  c*chi_{ab} + a*chi_{bc} = b*chi_{ac}  for any triple a<b<c.
- At a triple-binding argmax (binding pairs (a,b),(b,c),(a,c) with common
  value rho and sign pattern (+,+,-)), the identity forces
  2(ab+ac+bc)*rho = integer, i.e. rho = k / (2(ab+ac+bc)).
- Harmonic n=5 triple {3,4,5}: 2*47 = 94, claimed rho = 32/94 = 16/47.
  m=12 triple {2,5,12}: 2*94 = 188, claimed rho = 68/188 = 17/47.

Phases: grid47, bb_hunt, tiny_exact, structure, certify, m9_hunt, fresh
Output: /home/z/my-project/scripts/out_rho_closeout.json (incremental)
"""
import sys, json, math, time, heapq
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F
from itertools import combinations, product

import lrc_zono_lib as L
import certify_rung2_n5 as BB

OUT = '/home/z/my-project/scripts/out_rho_closeout.json'


def load_out():
    try:
        return json.load(open(OUT))
    except Exception:
        return {}


def save_out(d):
    json.dump(d, open(OUT, 'w'), indent=1)


# ---------------------------------------------------------------- z-set

def tight_Z0(v, forms, M):
    """Lattice points within V-distance M of the unit box (tight filter)."""
    d = len(v) - 1
    out = []
    rng = range(-2, 4)
    for z in product(rng, repeat=d):
        ok = True
        for coeffs, denom in forms:
            lo = hi = 0
            for k, cc in enumerate(coeffs):
                if cc:
                    t = cc * z[k]
                    interval = (min(t - cc, t), max(t - cc, t))
                    lo += interval[0]
                    hi += interval[1]
            # min over x in [0,1]^d of |form(z-x)| <= M required
            if lo > 0 or hi < 0:
                if min(abs(F(lo, denom)), abs(F(hi, denom))) > M:
                    ok = False
                    break
            # else 0 in [lo,hi] -> min 0, fine
        if ok:
            out.append(z)
    return out


# ---------------------------------------------------------------- float grid scan

def grid_scan(v, B, topk=400, Mextra=0.05):
    """Float D over the B^d grid on [0,1)^d; exact D at the top-k points.
    Returns dict with float max, exact max over top-k, argmax (exact)."""
    import numpy as np
    forms = L.make_forms(v)
    d = len(v) - 1
    rho_est = 0.5
    Z0 = tight_Z0(v, forms, rho_est + Mextra)
    print('  tight Z0 size:', len(Z0))
    axes = [np.arange(B) / B for _ in range(d)]
    fdata = []
    for coeffs, denom in forms:
        supp = [k for k in range(d) if coeffs[k] != 0]
        fdata.append((supp, [coeffs[k] / denom for k in supp]))
    best = np.full((B,) * d, np.inf)
    for z in Z0:
        acc = None
        for supp, cf in fdata:
            term = None
            for k, a in zip(supp, cf):
                sh = [1] * d
                sh[k] = B
                comp = (axes[k].reshape(sh) - z[k]) * a
                term = comp if term is None else term + comp
            term = np.abs(term)
            acc = term if acc is None else np.maximum(acc, term)
        best = np.minimum(best, acc)
        del acc
    fmax = float(best.max())
    print('  float grid max:', fmax)
    flat = best.flatten()
    kk = min(topk, flat.size)
    idxs = np.argpartition(flat, -kk)[-kk:]
    cands = []
    for fi in idxs:
        idx = np.unravel_index(fi, best.shape)
        x = tuple(F(int(i), B) for i in idx)
        cands.append((float(flat[fi]), x))
    # dedup near-duplicates, keep distinct points
    cands.sort(key=lambda t: -t[0])
    seen = []
    picked = []
    for fval, x in cands:
        xf = [float(t) for t in x]
        if any(max(abs(a - b) for a, b in zip(xf, s)) < 1.0 / B for s in seen):
            continue
        seen.append(xf)
        picked.append(x)
        if len(picked) >= 80:
            break
    best_exact = F(-1)
    best_x = None
    t0 = time.time()
    for x in picked:
        val = L.D_exact(x, v, forms)
        if val > best_exact:
            best_exact, best_x = val, x
    print('  exact max over top-%d distinct: %s at %s (%.1fs)' %
          (len(picked), best_exact, best_x, time.time() - t0))
    return {'B': B, 'float_max': fmax, 'n_z': len(Z0),
            'exact_max': str(best_exact), 'argmax': [str(t) for t in best_x],
            'top_points': [[str(t) for t in x] for x in picked[:20]]}


# ---------------------------------------------------------------- tiny local exact

def tiny_exact(v, center, hw, radius=0.45):
    """Local exact arrangement in a tiny box center +- hw (Fractions)."""
    forms = L.make_forms(v)
    d = len(v) - 1
    box = tuple((F(center[k]) - F(hw), F(center[k]) + F(hw)) for k in range(d))
    M = L.M_corner(v, forms)
    Z0 = L.enumerate_Z0(v, forms, M)
    t0 = time.time()
    val, arg, nh, nv = L._local_exact(v, forms, Z0, box, radius=radius)
    print('  tiny box %s..%s -> %s at %s (%d hyps %d cands %.1fs)' %
          (box[0], box[-1], val, arg, nh, nv, time.time() - t0))
    return {'box': [str(b) for b in box], 'val': str(val),
            'arg': [str(a) for a in arg], 'n_hyps': nh, 'n_cands': nv,
            'sec': round(time.time() - t0, 1)}


# ---------------------------------------------------------------- structure

def binding_structure(v, y):
    forms = L.make_forms(v)
    d = len(v) - 1
    n = len(v)
    rho = L.D_exact(y, v, forms)
    rho2 = BB.D_at(list(y), v)
    assert rho == rho2, ('D mismatch', rho, rho2)
    near = []
    base = [int(round(float(yk))) for yk in y]
    # adaptive per-coordinate z-range: |c_j| <= (1+v_j)*rho is the axis extent
    ranges = []
    for k in range(d):
        ext = int((1 + v[k + 1]) * float(rho)) + 1
        ranges.append(range(-ext, ext + 1))
    for off in product(*ranges):
        z = tuple(base[k] + off[k] for k in range(d))
        val = L.norm_V(tuple(yk - zk for yk, zk in zip(y, z)), forms)
        if val == rho:
            near.append(z)
    bindings = []
    for z in near:
        yz = tuple(yk - zk for yk, zk in zip(y, z))
        fi = 0
        for j in range(2, n + 1):
            val = L.form_value(forms[fi], yz)
            if abs(val) == rho:
                bindings.append({'z': list(z), 'pair': [1, v[j - 1]],
                                 'value': str(val)})
            fi += 1
        for i in range(2, n + 1):
            for j in range(i + 1, n + 1):
                val = L.form_value(forms[fi], yz)
                if abs(val) == rho:
                    bindings.append({'z': list(z), 'pair': [v[i - 1], v[j - 1]],
                                     'value': str(val)})
                fi += 1
    chis = []
    for i in range(1, n):
        for j in range(i + 1, n):
            a, b = v[i], v[j]
            chi = F(b) * y[i - 1] - F(a) * y[j - 1]
            chis.append({'pair': [a, b], 'chi_mod1': str(chi % 1)})
    box = [(x, x) for x in y]
    ks = BB.kink_list(box, v)
    best_s, best_val = None, None
    for s in ks:
        val = max(BB.normT(s), max(BB.normT(y[j] - v[j + 1] * s) for j in range(d)))
        if best_val is None or val < best_val:
            best_val, best_s = val, s
    runners = [{'speed': 1, 'res': str(BB.normT(best_s))}]
    for j in range(d):
        runners.append({'speed': v[j + 1], 'res': str(BB.normT(y[j] - v[j + 1] * best_s))})
    return {'v': list(v), 'y': [str(x) for x in y], 'rho': str(rho),
            'sigma_star': str(best_s), 'nearest_z': [list(z) for z in near],
            'bindings': bindings, 'chi_mod1': chis, 'runners': runners}


# ---------------------------------------------------------------- B&B hunt

def bb_survivors(v, r, budget_sec, sizecap_bits=18, max_depth=60):
    d = len(v) - 1
    t0 = time.time()
    r = F(r)
    K = 48
    rf = float(r)
    grid = ([i / K * rf for i in range(K + 1)] +
            [1 - rf + i / K * rf for i in range(K + 1)])

    def cheap_prune(box):
        for s in grid:
            vals = BB._R_values(s, box, v, False)
            if max(vals) <= rf - 1e-9:
                return True
        return False

    SIZECAP = F(1, 2 ** sizecap_bits)
    root = tuple((F(0), F(1)) for _ in range(d))
    heap = [(-1.0, 0, root)]
    survivors = []
    pops = 0
    while heap:
        if time.time() - t0 > budget_sec or pops > 600000:
            survivors.extend(b for _, _, b in heap)
            break
        _, depth, box = heapq.heappop(heap)
        pops += 1
        widths = [hi - lo for lo, hi in box]
        if depth >= max_depth or max(widths) <= SIZECAP:
            survivors.append(box)
            continue
        if cheap_prune(box):
            continue
        uf, sf = BB.U_box(box, v, False, r=r)
        if uf <= rf + 1e-7:
            ue, se = BB.U_box(box, v, True, r=r)
            if ue <= r:
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
        heapq.heappush(heap, (-vol, depth + 1, b1))
        heapq.heappush(heap, (-vol, depth + 1, b2))
    return survivors, pops, round(time.time() - t0, 1)


# ---------------------------------------------------------------- main

if __name__ == '__main__':
    phase = sys.argv[1]
    out = load_out()
    if phase == 'grid47':
        v = tuple(int(x) for x in sys.argv[2].split(','))
        B = int(sys.argv[3]) if len(sys.argv) > 3 else 47
        print('grid scan v=%s B=%d' % (v, B))
        res = grid_scan(v, B)
        out['grid47_%s_B%d' % (sys.argv[2], B)] = res
    elif phase == 'tiny_exact':
        v = tuple(int(x) for x in sys.argv[2].split(','))
        center = [F(s) for s in sys.argv[3].split(',')]
        hw = F(sys.argv[4])
        res = tiny_exact(v, center, hw)
        out['tiny_%s_%s' % (sys.argv[2], sys.argv[3].replace('/', '_'))] = res
    elif phase == 'structure':
        v = tuple(int(x) for x in sys.argv[2].split(','))
        y = tuple(F(s) for s in sys.argv[3].split(','))
        out['structure_%s_%s' % (sys.argv[2], sys.argv[3].replace('/', '_'))] = \
            binding_structure(v, y)
        print(json.dumps(out['structure_%s_%s' % (sys.argv[2], sys.argv[3].replace('/', '_'))], indent=1))
    elif phase == 'bb_hunt':
        v = tuple(int(x) for x in sys.argv[2].split(','))
        r = F(sys.argv[3])
        budget = float(sys.argv[4]) if len(sys.argv) > 4 else 600
        eps = F(1, 2 ** 22)
        surv, pops, sec = bb_survivors(v, r - eps, budget)
        bbox = None
        if surv:
            bbox = [[str(min(b[k][0] for b in surv)),
                     str(max(b[k][1] for b in surv))] for k in range(len(v) - 1)]
        print('survivors:', len(surv), 'pops', pops, 'sec', sec, 'bbox', bbox)
        out['bb_hunt_%s' % sys.argv[2]] = {
            'v': list(v), 'r': str(r - eps), 'n_survivors': len(surv),
            'pops': pops, 'sec': sec, 'survivor_bbox': bbox,
            'survivor_boxes_sample': [[str(x) for x in b] for b in surv[:30]],
        }
    elif phase == 'certify':
        name = sys.argv[2]
        v = tuple(int(x) for x in sys.argv[3].split(','))
        r = F(sys.argv[4])
        witness = [F(s) for s in sys.argv[5].split(',')]
        wval = F(sys.argv[6])
        budget = float(sys.argv[7]) if len(sys.argv) > 7 else 2400
        got = BB.D_at(witness, v)
        assert got == wval, ('witness mismatch', got, wval)
        res = BB.certify(name, v, r, [(witness, wval)], 'certify', budget)
        rec = {k: val for k, val in res.items()
               if k not in ('leaf_sample', 'survivor_sample')}
        if res.get('survivors'):
            rec['survivors'] = res['survivors']
        out['certify_' + name] = rec
        print(name, res['verdict'], '| leaves', res['leaves'],
              '| pops', res['pops'], '|', res['wall_sec'], 's')
    save_out(out)
    print('saved', OUT)
