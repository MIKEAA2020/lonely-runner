"""Driver: LR-zonotope Ehrhart h* / covering-radius program (audit 5-step plan).

Steps:
 1. h* for n=3..8, families {1..n} (harmonic) and {1..n-1, 2n} (rung-2),
    with DFS cross-validation; real-rootedness (exact, sympy Sturm),
    log-concavity, Newton-normalized ratios, strictness.
 2. lambda/delta/threshold margins.
 3. rho: n=3 exact (d=2 global arrangement) incl. sweep (1,2,m);
    n=4 exact (d=3 grid+local) incl. sweep (1,2,3,m);
    n=5,6 (d=4,5) grid+local or enclosure.
 4. deepest-hole test everywhere.
 5. correlation (Spearman) between h* strictness and margins.
Outputs: scripts/out_lrc_zono.json + download/lrc_zonotope_ehrhart_report.md
"""
import sys, json, time, math
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
from itertools import combinations
from lrc_zono_lib import (ehr_coeffs, hstar_from_ehr, hstar_by_inversion,
                          ehr_from_hstar, ehr_poly_eval, dfs_count,
                          lambda_exact, make_forms, cstar, D_exact,
                          norm_V, M_corner, enumerate_Z0,
                          rho_exact_d2, rho_semid2, D_grid_float)

T0 = time.time()
OUT = {}

def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    def ranks(a):
        order = sorted(range(n), key=lambda i: a[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and a[order[j + 1]] == a[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else None

# ---------------------------------------------------------------- h* analysis

def hstar_analysis(v):
    n = len(v)
    d = n - 1
    e = ehr_coeffs(v)
    hs = hstar_from_ehr(e, d)
    Ehr_vals = [ehr_poly_eval(e, t) for t in range(d + 2)]
    hs2 = hstar_by_inversion(Ehr_vals, d)
    assert hs == hs2, f"h* mismatch {v}"
    assert [ehr_from_hstar(hs, d, t) for t in range(d + 2)] == Ehr_vals
    for t in (1, 2):
        assert dfs_count(v, t) == Ehr_vals[t], f"DFS mismatch {v} t={t}"
    deg = len(hs) - 1
    # roots
    import sympy as sp
    z = sp.symbols("z")
    p = sp.Poly(list(reversed(hs)), z)   # sympy wants descending
    try:
        rr = p.real_roots()
        n_real = len(rr)
        roots_str = [str(r.evalf(6)) for r in rr]
    except Exception as ex:
        n_real, roots_str = None, [f"err:{ex}"]
    real_rooted = (n_real == deg) if n_real is not None else None
    # log-concavity
    lc_ok = True
    ratios = []
    for k in range(1, deg):
        if hs[k - 1] * hs[k + 1] > 0:
            r = hs[k] ** 2 / (hs[k - 1] * hs[k + 1])
            ratios.append(r)
            if r < 1:
                lc_ok = False
        elif hs[k] == 0 and (hs[k - 1] > 0 or hs[k + 1] > 0):
            lc_ok = False  # internal zero
            ratios.append(None)
        else:
            ratios.append(None)
    newton_ok = True
    newton_ratios = []
    for k in range(1, deg):
        if hs[k - 1] * hs[k + 1] > 0:
            from math import comb
            num = hs[k] ** 2 * comb(d, k - 1) * comb(d, k + 1)
            den = hs[k - 1] * hs[k + 1] * comb(d, k) ** 2
            r = num / den
            newton_ratios.append(r)
            if r < 1:
                newton_ok = False
        else:
            newton_ratios.append(None)
    unimodal = all(hs[i] <= hs[i + 1] or all(hs[j] >= hs[j + 1] for j in range(i, deg))
                   for i in range(deg)) if deg >= 2 else True
    strictness = min([r for r in ratios if r is not None], default=None)
    lam, tstar = lambda_exact(v)
    delta = F(1, 2) - lam
    thr = F(n - 1, 2 * (n + 1))
    return {
        "v": list(v), "n": n, "d": d,
        "e_coeffs": e, "Ehr_0_to_d+1": Ehr_vals,
        "hstar": hs, "hstar_str": [str(x) for x in hs], "degree": deg,
        "roots": roots_str, "n_real_roots": n_real,
        "real_rooted": real_rooted, "log_concave": lc_ok,
        "lc_ratios": [str(r) if r is not None else None for r in ratios],
        "newton_ok": newton_ok,
        "newton_ratios": [str(r) if r is not None else None for r in newton_ratios],
        "unimodal": unimodal,
        "strictness": str(strictness) if strictness is not None else None,
        "strictness_f": float(strictness) if strictness is not None else None,
        "lambda": str(lam), "lambda_f": float(lam), "tstar": str(tstar),
        "delta": str(delta), "delta_f": float(delta),
        "thr": str(thr), "thr_f": float(thr),
        "delta_margin": str(thr - delta), "delta_margin_f": float(thr - delta),
        "lambda_margin": str(lam - F(1, n + 1)),
        "lambda_margin_f": float(lam - F(1, n + 1)),
        "Ehr_minus1_interior": ehr_poly_eval(e, -1),
        "vol_sum_v": sum(v),
    }

print("=== Step 1-2: h* tables (both families, n=3..8) ===", flush=True)
families = []
for n in range(3, 9):
    families.append(("harmonic", tuple(range(1, n + 1))))
    families.append(("rung2", tuple(list(range(1, n)) + [2 * n])))
hstar_table = []
for fam, v in families:
    rec = hstar_analysis(v)
    rec["family"] = fam
    hstar_table.append(rec)
    print(f"{fam:9s} n={rec['n']}: h*={rec['hstar_str']} deg={rec['degree']} "
          f"real_rooted={rec['real_rooted']} LC={rec['log_concave']} "
          f"Newton={rec['newton_ok']} strict={rec['strictness']} "
          f"lam={rec['lambda']} margin={rec['delta_margin']}", flush=True)
OUT["hstar_table"] = hstar_table

print("\n=== Step 3-4: rho n=3 (d=2, exact) both families + sweep (1,2,m) ===",
      flush=True)
rho_d2 = []
sweep3 = []
def rho_rec_d2(v, fam, tag):
    forms = make_forms(v)
    lam, _ = lambda_exact(v)
    delta = F(1, 2) - lam
    thr = F(len(v) - 1, 2 * (len(v) + 1))
    t0 = time.time()
    rho, arg, nhyp, ncand = rho_exact_d2(v, forms)
    rec = {
        "v": list(v), "family": fam, "tag": tag, "n": len(v),
        "lambda": str(lam), "delta": str(delta), "thr": str(thr),
        "rho": str(rho), "rho_f": float(rho),
        "argmax": [str(a) for a in arg],
        "deepest_hole": rho == delta,
        "rho_minus_delta": str(rho - delta), "rho_minus_delta_f": float(rho - delta),
        "rho_leq_thr": rho <= thr, "rho_minus_thr": str(rho - thr),
        "rho_minus_thr_f": float(rho - thr),
        "n_hyp": nhyp, "n_cand": ncand, "sec": round(time.time() - t0, 1),
    }
    return rec

for fam in ("harmonic", "rung2"):
    n = 3
    v = tuple(range(1, n + 1)) if fam == "harmonic" else (1, 2, 6)
    rec = rho_rec_d2(v, fam, "core")
    rho_d2.append(rec)
    print(f"{fam:9s} n=3: rho={rec['rho']} delta={rec['delta']} deepest={rec['deepest_hole']} "
          f"rho-thr={rec['rho_minus_thr']} [{rec['sec']}s]", flush=True)
for m in range(3, 15):
    v = (1, 2, m)
    rec = rho_rec_d2(v, "sweep(1,2,m)", f"m={m}")
    sweep3.append(rec)
    print(f"(1,2,{m:2d}): rho={rec['rho']:>8s} delta={rec['delta']:>8s} "
          f"deepest={rec['deepest_hole']} rho<=thr={rec['rho_leq_thr']} "
          f"rho-thr={rec['rho_minus_thr_f']:+.6f} [{rec['sec']}s]", flush=True)
OUT["rho_n3_core"] = rho_d2
OUT["rho_n3_sweep"] = sweep3

print("\n=== h* for the n=3 sweep (for correlation) ===", flush=True)
sweep3_h = []
for m in range(3, 15):
    rec = hstar_analysis((1, 2, m))
    rec["family"] = "sweep(1,2,m)"
    sweep3_h.append(rec)
OUT["hstar_n3_sweep"] = sweep3_h

print("\n=== Step 5a: correlation at n=3 (12 instances) ===", flush=True)
xs = [r["strictness_f"] for r in sweep3_h]
ys = [r["delta_margin_f"] for r in sweep3_h]
zs = [r["rho_f"] if "rho_f" in r else None for r in sweep3_h]
rho3 = [r["rho_f"] for r in sweep3]
sp_delta = spearman(xs, ys)
sp_rho_thr = spearman(xs, [1.0 if r["rho_leq_thr"] else 0.0 for r in sweep3])
print(f"  strictness values: {[round(x,3) for x in xs]}")
print(f"  delta margins:     {[round(y,4) for y in ys]}")
print(f"  Spearman(strictness, delta-margin) = {sp_delta}")
OUT["correlation_n3"] = {
    "strictness": xs, "delta_margin": ys,
    "spearman_strictness_vs_delta_margin": sp_delta,
    "rho_values": rho3,
}

print(f"\n[total {time.time()-T0:.0f}s] saving partial JSON", flush=True)
with open("/home/z/my-project/scripts/out_lrc_zono.json", "w") as f:
    json.dump(OUT, f, indent=1, default=str)
print("PART 1 DONE")
