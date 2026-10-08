"""Stronger witness search for the n=5 consistent candidates (m=10,20,25,30)
and n=4 m=16: multi-start Nelder-Mead on float D (256 random + 16 torsion +
16 grid seeds), exact re-evaluation of refined optima.
A fail needs only ONE point with exact D > delta; if none is found, the
'consistent' label is upgraded to 'consistent (strong multi-start search)'."""
import sys, json, time
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
import lrc_zono_lib as L
from lrc_zono_lib import lambda_exact, make_forms, cstar, D_exact, M_corner, enumerate_Z0

def float_D_factory(v, forms, Z0):
    d = len(v) - 1
    C = np.array([[coeffs[k] / denom for k in range(d)] for coeffs, denom in forms],
                 dtype=float)                      # (nforms, d)
    Zarr = np.array(Z0, dtype=float)               # (nz, d)
    def Df(x):
        x = np.mod(np.asarray(x, dtype=float), 1.0)
        vals = np.abs((x - Zarr) @ C.T)            # (nz, nforms)
        return float(np.min(np.max(vals, axis=1)))
    return Df

def strong_search(v, n_random=256, n_grid=16):
    d = len(v) - 1
    forms = make_forms(v)
    M = M_corner(v, forms)
    Z0 = enumerate_Z0(v, forms, M)
    lam, _ = lambda_exact(v)
    delta = F(1, 2) - lam
    Df = float_D_factory(v, forms, Z0)
    rng = np.random.default_rng(12345 + len(v) * 1000 + v[-1])
    starts = [tuple(rng.uniform(0, 1, d)) for _ in range(n_random)]
    for bits in range(2 ** d):
        starts.append(tuple(0.5 if (bits >> k) & 1 else 0.0 for k in range(d)))
    # grid seeds (coarse K=16 to bound cost)
    K = 16
    g = L._grid_D_numpy(v, forms, Z0, K)
    flat = g.flatten()
    for fi in np.argpartition(flat, -n_grid)[-n_grid:]:
        idx = np.unravel_index(fi, g.shape)
        starts.append(tuple(int(i) / K for i in idx))
    refined = []
    t0 = time.time()
    for s in starts:
        r = minimize(lambda x: -Df(x), np.array(s), method="Nelder-Mead",
                     options={"maxiter": 500, "xatol": 1e-9, "fatol": 1e-12})
        refined.append(tuple(np.mod(r.x, 1.0)))
    # dedupe on 2^-14 grid, exact re-eval of the distinct top candidates
    uniq = {}
    for x in refined:
        key = tuple(int(round(c * 2 ** 14)) for c in x)
        uniq[key] = x
    vals = []
    for key, x in uniq.items():
        fval = Df(x)
        vals.append((fval, x))
    vals.sort(key=lambda p: -p[0])
    best_val = F(-1)
    best_x = None
    for fval, x in vals[:30]:
        xr = tuple(F(int(round(c * 2 ** 21)), 2 ** 21) for c in x)
        vv = D_exact(xr, v, forms)
        if vv > best_val:
            best_val, best_x = vv, xr
    return {"v": list(v), "d": d, "delta": str(delta),
            "best_val": str(best_val), "best_x": [str(t) for t in best_x],
            "beats_delta": best_val > delta,
            "margin": str(best_val - delta), "n_starts": len(starts),
            "n_uniq": len(uniq), "sec": round(time.time() - t0, 1)}

if __name__ == "__main__":
    out = []
    for v in [(1, 2, 3, 4, 10), (1, 2, 3, 4, 20), (1, 2, 3, 4, 25),
              (1, 2, 3, 4, 30), (1, 2, 3, 16)]:
        rec = strong_search(v)
        out.append(rec)
        print(f"v={v}: delta={rec['delta']} best={rec['best_val']} "
              f"beats_delta={rec['beats_delta']} margin={rec['margin']} "
              f"starts={rec['n_starts']} uniq={rec['n_uniq']} [{rec['sec']}s]",
              flush=True)
    with open("/home/z/my-project/scripts/out_lrc_zono_rescan_strong.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("STRONG SEARCH DONE")
