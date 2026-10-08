"""n=4 sweep (1,2,3,m): h* + exact rho (d=3 grid+local). Standalone."""
import sys, json, time, math
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
from lrc_zono_lib import (ehr_coeffs, hstar_from_ehr, hstar_by_inversion,
                           ehr_from_hstar, ehr_poly_eval, dfs_count,
                           lambda_exact, make_forms, rho_semid2)

def spearman(xs, ys):
    n = len(xs)
    if n < 3: return None
    def ranks(a):
        order = sorted(range(n), key=lambda i: a[i])
        r = [0.0]*n; i = 0
        while i < n:
            j = i
            while j+1 < n and a[order[j+1]] == a[order[i]]: j += 1
            avg = (i+j)/2 + 1
            for k in range(i, j+1): r[order[k]] = avg
            i = j+1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx)/n, sum(ry)/n
    num = sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    den = math.sqrt(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))
    return num/den if den else None

def hstar_analysis(v):
    n = len(v); d = n-1
    e = ehr_coeffs(v)
    hs = hstar_from_ehr(e, d)
    Ehr_vals = [ehr_poly_eval(e, t) for t in range(d+2)]
    assert hs == hstar_by_inversion(Ehr_vals, d)
    assert [ehr_from_hstar(hs, d, t) for t in range(d+2)] == Ehr_vals
    for t in (1, 2):
        assert dfs_count(v, t) == Ehr_vals[t]
    deg = len(hs)-1
    import sympy as sp
    z = sp.symbols("z")
    p = sp.Poly(list(reversed(hs)), z)
    try:
        rr = p.real_roots(); n_real = len(rr)
        roots_str = [str(r.evalf(6)) for r in rr]
    except Exception as ex:
        n_real, roots_str = None, [f"err"]
    lc_ok = True; ratios = []
    for k in range(1, deg):
        if hs[k-1]*hs[k+1] > 0:
            r = hs[k]**2/(hs[k-1]*hs[k+1]); ratios.append(r)
            if r < 1: lc_ok = False
        elif hs[k] == 0 and (hs[k-1] > 0 or hs[k+1] > 0):
            lc_ok = False; ratios.append(None)
        else: ratios.append(None)
    newton_ok = True; newton_ratios = []
    from math import comb
    for k in range(1, deg):
        if hs[k-1]*hs[k+1] > 0:
            r = hs[k]**2*comb(d,k-1)*comb(d,k+1)/(hs[k-1]*hs[k+1]*comb(d,k)**2)
            newton_ratios.append(r)
            if r < 1: newton_ok = False
        else: newton_ratios.append(None)
    strictness = min([r for r in ratios if r is not None], default=None)
    lam, tstar = lambda_exact(v)
    delta = F(1,2) - lam
    thr = F(n-1, 2*(n+1))
    return {"e_coeffs": e, "hstar_str": [str(x) for x in hs], "degree": deg,
            "roots": roots_str, "n_real_roots": n_real,
            "real_rooted": (n_real == deg) if n_real is not None else None,
            "log_concave": lc_ok, "newton_ok": newton_ok,
            "strictness": str(strictness) if strictness is not None else None,
            "strictness_f": float(strictness) if strictness is not None else None,
            "lambda": str(lam), "delta": str(delta), "thr": str(thr),
            "delta_margin_f": float(thr - delta)}

if __name__ == "__main__":
    mlo, mhi = int(sys.argv[1]), int(sys.argv[2])
else:
    mlo, mhi = 0, 0
sweep4 = []; sweep4_h = []
for m in (range(mlo, mhi+1) if mlo else []):
    v = (1, 2, 3, m)
    t0 = time.time()
    rec_h = hstar_analysis(v)
    sweep4_h.append(rec_h)
    forms = make_forms(v)
    lam, _ = lambda_exact(v)
    delta = F(1,2) - lam
    thr = F(3, 10)
    res = rho_semid2(v, forms, K=40)
    rec = {"v": list(v), "m": m, "n": 4,
           "lambda": str(lam), "delta": str(delta), "thr": str(thr),
           "rho_lo_exact": str(res["rho_lo_exact"]), "rho_f": float(res["rho_lo_exact"]),
           "rho_exact": str(res["rho_exact"]),
           "away_certified": res["away_certified"],
           "rho_hi_grid": res["rho_hi_grid"],
           "deepest_hole": res["rho_lo_exact"] == delta,
           "rho_leq_thr": res["rho_lo_exact"] <= thr,
           "rho_minus_thr_f": float(res["rho_lo_exact"] - thr),
           "argmax": [str(a) for a in res["argmax"]] if res["argmax"] else None,
           "n_boxes": res["n_boxes"], "sec": round(time.time()-t0, 1),
           "strictness_f": rec_h["strictness_f"],
           "delta_margin_f": rec_h["delta_margin_f"]}
    sweep4.append(rec)
    print(f"(1,2,3,{m:2d}): lam={lam} delta={delta} rho={res['rho_lo_exact']} "
          f"exact={res['rho_exact'] is not None} away={res['away_certified']} "
          f"deepest={rec['deepest_hole']} rho<=thr={rec['rho_leq_thr']} "
          f"strict={rec_h['strictness_f']:.3f} [{rec['sec']}s]", flush=True)

with open(f"/home/z/my-project/scripts/out_lrc_zono2_{mlo}_{mhi}.json", "w") as f:
    json.dump({"rho_n4_sweep": sweep4, "hstar_n4_sweep": sweep4_h}, f, indent=1, default=str)
print("CHUNK DONE", flush=True)
