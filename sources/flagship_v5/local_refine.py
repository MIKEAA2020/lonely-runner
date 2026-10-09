#!/usr/bin/env python3
"""Local float-grid refinement around a seed witness, then exact top-K.
Usage: local_refine.py v_csv seed_csv [halfwidth] [K]
Prints float max, float argmax, exact max over top distinct points."""
import sys, time
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F
import numpy as np
import lrc_zono_lib as L


def tight_zloc(v, forms, box, Mcut):
    """Lattice points z with min over box of some |form(z-x)| <= Mcut,
    enumerated over a generous integer range (axis extents included)."""
    d = len(v) - 1
    out = []
    rngs = []
    for k in range(d):
        ext = int((1 + v[k + 1]) * Mcut) + 2
        lo = int(np.floor(min(box[k]))) - ext
        hi = int(np.ceil(max(box[k]))) + ext
        rngs.append(range(lo, hi + 1))
    from itertools import product
    for z in product(*rngs):
        ok = True
        for coeffs, denom in forms:
            lo = hi = 0.0
            for k, cc in enumerate(coeffs):
                if cc:
                    t = cc * z[k]
                    lo += min(t - cc * box[k][1], t - cc * box[k][0])
                    hi += max(t - cc * box[k][1], t - cc * box[k][0])
            if lo > 0 or hi < 0:
                if min(abs(lo / denom), abs(hi / denom)) > Mcut:
                    ok = False
                    break
        if ok:
            out.append(z)
    return out


def main():
    v = tuple(int(x) for x in sys.argv[1].split(','))
    seed = [F(s) for s in sys.argv[2].split(',')]
    hw = F(sys.argv[3]) if len(sys.argv) > 3 else F(3, 250)
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 64
    forms = L.make_forms(v)
    d = len(v) - 1
    box = tuple((float(s - hw), float(s + hw)) for s in seed)
    zloc = tight_zloc(v, forms, box, 0.5)
    print('zloc:', len(zloc))
    axes = [np.linspace(box[k][0], box[k][1], K) for k in range(d)]
    fdata = []
    for coeffs, denom in forms:
        supp = [k for k in range(d) if coeffs[k] != 0]
        fdata.append((supp, [coeffs[k] / denom for k in supp]))
    best = np.full((K,) * d, np.inf)
    t0 = time.time()
    for z in zloc:
        acc = None
        for supp, cf in fdata:
            term = None
            for k, a in zip(supp, cf):
                sh = [1] * d
                sh[k] = K
                comp = (axes[k].reshape(sh) - z[k]) * a
                term = comp if term is None else term + comp
            term = np.abs(term)
            acc = term if acc is None else np.maximum(acc, term)
        best = np.minimum(best, acc)
    print('grid %.1fs max %s' % (time.time() - t0, best.max()))
    idx = np.unravel_index(np.argmax(best), best.shape)
    fam = [float(axes[k][idx[k]]) for k in range(d)]
    print('float argmax', fam)
    flat = best.flatten()
    kk = min(60, flat.size)
    idxs = np.argpartition(flat, -kk)[-kk:]
    bex, bx = F(-1), None
    for fi in idxs:
        i = np.unravel_index(fi, best.shape)
        x = tuple(F(int(round(float(axes[k][i[k]]) * 10 ** 7)), 10 ** 7) for k in range(d))
        val = L.D_exact(x, v, forms)
        if val > bex:
            bex, bx = val, x
    print('exact max over top-%d (1e7 approx): %s at %s' % (kk, bex, bx))


if __name__ == '__main__':
    main()
