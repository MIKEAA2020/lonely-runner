"""Part 3: rho enclosures for d=4,5 (n=5,6 core instances) with exact witnesses."""
import sys, json, time
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
import numpy as np
import lrc_zono_lib as L

def rho_enclosure(v, K, K_loc=25, rounds=3, verbose=True):
    """Two-sided rho for d in {4,5}:
    global float grid (K) -> certified upper fmax + h/2;
    exact witnesses: torsion + top grid points; then local multi-scale
    float refinement around the best witness (exact re-evaluated)."""
    forms = L.make_forms(v)
    d = len(v) - 1
    M = L.M_corner(v, forms)
    Z0 = L.enumerate_Z0(v, forms, M)
    grid = L._grid_D_numpy(v, forms, Z0, K)
    fmax = float(grid.max())
    h = 1.0 / K
    rho_hi = fmax + h / 2 + 1e-6

    # exact torsion values
    tor = []
    for bits in range(2 ** d):
        x = tuple(F(1, 2) if (bits >> k) & 1 else F(0) for k in range(d))
        tor.append((x, L.D_exact(x, v, forms)))
    tor_max = max(v_ for _, v_ in tor)

    # exact at top-50 grid points
    flat = grid.flatten()
    idxs = np.argpartition(flat, -50)[-50:]
    exact_pts = []
    for fi in idxs:
        idx = np.unravel_index(fi, grid.shape)
        x = tuple(F(int(i), K) for i in idx)
        exact_pts.append((x, L.D_exact(x, v, forms)))
    best_x, best_val = max(exact_pts, key=lambda p: p[1])

    # local multi-scale refinement (float grid, restricted z set)
    center = [float(c) for c in best_x]
    rng = 1.5 * h
    for r in range(rounds):
        axes = [np.linspace(center[k] - rng, center[k] + rng, K_loc)
                for k in range(d)]
        meshes = np.meshgrid(*axes, indexing="ij")
        pts = np.stack([m.ravel() for m in meshes], axis=1)
        vals = np.full(len(pts), np.inf)
        zloc = [z for z in Z0
                if max(abs(z[k] - center[k]) for k in range(d)) < 1.5 * (rng + rho_hi)]
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
        best_f = float(vals[i])
    # exact at refined center
    xr = tuple(F(round(c * 2 ** 21), 2 ** 21) for c in center)
    vr = L.D_exact(xr, v, forms)
    if vr > best_val:
        best_x, best_val = xr, vr
    rho_lo = max(best_val, tor_max)
    return {
        "v": list(v), "d": d, "K": K,
        "rho_hi": rho_hi, "rho_lo_exact": rho_lo,
        "rho_lo_f": float(rho_lo),
        "enclosure_width": rho_hi - float(rho_lo),
        "tor_max": tor_max, "tor_max_f": float(tor_max),
        "witness": best_x, "witness_val": best_val,
        "witness_val_f": float(best_val),
        "fmax": fmax, "Z0": len(Z0),
    }

out = []
for v, K in [((1, 2, 3, 4, 5), 32), ((1, 2, 3, 4, 10), 32),
             ((1, 2, 3, 4, 5, 6), 16), ((1, 2, 3, 4, 5, 12), 16)]:
    t0 = time.time()
    lam, _ = L.lambda_exact(v)
    delta = F(1, 2) - lam
    thr = F(len(v) - 1, 2 * (len(v) + 1))
    res = rho_enclosure(v, K)
    res.update({"lambda": str(lam), "delta": str(delta), "thr": str(thr),
                "delta_f": float(delta), "thr_f": float(thr),
                "deepest_hole_possible": res["rho_lo_exact"] == delta,
                "rho_gt_thr": res["rho_lo_exact"] > thr,
                "rho_gt_delta": res["rho_lo_exact"] > delta,
                "sec": round(time.time() - t0, 1)})
    out.append(res)
    print(f"v={v}:", flush=True)
    print(f"  lam={lam} delta={delta} thr={thr}", flush=True)
    print(f"  rho_lo(exact)={res['rho_lo_exact']} ({res['rho_lo_f']:.6f}) "
          f"rho_hi={res['rho_hi']:.6f} width={res['enclosure_width']:.6f}", flush=True)
    print(f"  witness={res['witness']} tor_max={res['tor_max_f']:.6f}", flush=True)
    print(f"  rho>thr: {res['rho_gt_thr']}  rho>delta: {res['rho_gt_delta']} "
          f"deepest-hole: {res['deepest_hole_possible']}  [{res['sec']}s]", flush=True)

with open("/home/z/my-project/scripts/out_lrc_zono3.json", "w") as f:
    json.dump(out, f, indent=1, default=str)
print("PART 3 DONE", flush=True)
