"""n=4 sweep (1,2,3,m) m=5..12: h* + rho enclosures with exact witnesses."""
import sys, json, time
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
import numpy as np
import lrc_zono_lib as L
from lrc_zono_run2 import hstar_analysis   # module-level func (safe import)

def rho_enclosure(v, K, K_loc=25, rounds=3):
    forms = L.make_forms(v)
    d = len(v) - 1
    M = L.M_corner(v, forms)
    Z0 = L.enumerate_Z0(v, forms, M)
    grid = L._grid_D_numpy(v, forms, Z0, K)
    fmax = float(grid.max()); h = 1.0 / K
    rho_hi = fmax + h / 2 + 1e-6
    tor = []
    for bits in range(2 ** d):
        x = tuple(F(1, 2) if (bits >> k) & 1 else F(0) for k in range(d))
        tor.append(L.D_exact(x, v, forms))
    tor_max = max(t for t in tor)
    flat = grid.flatten()
    idxs = np.argpartition(flat, -50)[-50:]
    exact_pts = []
    for fi in idxs:
        idx = np.unravel_index(fi, grid.shape)
        x = tuple(F(int(i), K) for i in idx)
        exact_pts.append((x, L.D_exact(x, v, forms)))
    best_x, best_val = max(exact_pts, key=lambda p: p[1])
    center = [float(c) for c in best_x]
    rng = 1.5 * h
    for r in range(rounds):
        axes = [np.linspace(center[k] - rng, center[k] + rng, K_loc) for k in range(d)]
        meshes = np.meshgrid(*axes, indexing="ij")
        pts = np.stack([m.ravel() for m in meshes], axis=1)
        vals = np.full(len(pts), np.inf)
        zloc = [z for z in Z0 if max(abs(z[k] - center[k]) for k in range(d)) < 1.5 * (rng + rho_hi)]
        for z in zloc:
            acc = None
            for coeffs, denom in forms:
                s = None
                for k, a in enumerate(coeffs):
                    if a:
                        t = (pts[:, k] - z[k]) * (a / denom)
                        s = t if s is None else s + t
                s = np.abs(s)
                acc = s if acc is None else np.maximum(acc, s)
            vals = np.minimum(vals, acc)
        i = int(np.argmax(vals))
        center = pts[i].tolist()
        rng = rng * 2 / (K_loc - 1) * 1.2
    xr = tuple(F(round(c * 2 ** 21), 2 ** 21) for c in center)
    vr = L.D_exact(xr, v, forms)
    if vr > best_val: best_x, best_val = xr, vr
    rho_lo = max(best_val, tor_max)
    return {"rho_hi": rho_hi, "rho_lo_exact": rho_lo, "rho_lo_f": float(rho_lo),
            "tor_max": str(tor_max), "witness": best_x,
            "witness_val": str(best_val), "fmax": fmax}

sweep4 = []; sweep4_h = []
for m in range(5, 13):
    v = (1, 2, 3, m)
    t0 = time.time()
    rec_h = hstar_analysis(v)
    sweep4_h.append(rec_h)
    lam, _ = L.lambda_exact(v)
    delta = F(1, 2) - lam
    thr = F(3, 10)
    res = rho_enclosure(v, K=40)
    rec = {"v": list(v), "m": m, "n": 4, "lambda": str(lam), "delta": str(delta),
           "thr": str(thr), "rho_lo_exact": str(res["rho_lo_exact"]),
           "rho_f": res["rho_lo_f"], "rho_hi": res["rho_hi"],
           "deepest_hole": res["rho_lo_exact"] == delta,
           "rho_leq_thr": res["rho_lo_exact"] <= thr,
           "rho_minus_thr_f": float(res["rho_lo_exact"] - thr),
           "witness": [str(a) for a in res["witness"]],
           "strictness_f": rec_h["strictness_f"],
           "delta_margin_f": rec_h["delta_margin_f"],
           "sec": round(time.time() - t0, 1)}
    sweep4.append(rec)
    print(f"(1,2,3,{m:2d}): lam={lam} delta={delta} rho_lo={res['rho_lo_exact']} "
          f"({res['rho_lo_f']:.6f}) hi={res['rho_hi']:.6f} deepest={rec['deepest_hole']} "
          f"rho<=thr={rec['rho_leq_thr']} strict={rec_h['strictness_f']:.3f} [{rec['sec']}s]", flush=True)

with open("/home/z/my-project/scripts/out_lrc_zono2b.json", "w") as f:
    json.dump({"rho_n4_sweep_enclosure": sweep4, "hstar_n4_sweep": sweep4_h},
              f, indent=1, default=str)
print("SWEEP DONE", flush=True)
