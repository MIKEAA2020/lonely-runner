"""zonotope-2: audit re-check of the "deepest-hole iff n|m" law.

The audit flagged (all confirmed against the persisted data):
 1. The report stated the law "iff n|m for n=3,4,5 (rung-2 m=2n included, and
    harmonic m=n)" while its own n=5 harmonic entry (m=5, 5|5) FAILS exactly
    (witness 97/288 > 1/3 = delta) -- internally inconsistent as stated.
 2. The n=4 "core exact" values (3/10, 5/18) were hardcoded in the assembler
    without a persisted run.
 3. The positive entries n=4 m=8/12 and n=5 m=10 are enclosure-level
    "consistent" labels (rho_lo == delta), NOT certified holds.
 4. There was no n=5 m-scan beyond {5, 10} and no n=6 scan beyond {6, 12}.

Phases (selectable via argv):
  A <none>          exact re-assertions: refutation witnesses, D(c*)=delta, n=3 core
  B <mlo> <mhi>     n=3 sweep extension, exact rho (rho_exact_d2)
  C <mlo> <mhi>     n=4 (1,2,3,m) rho_semid2 sweep (away-certified exact, PERSISTED)
  D <mlo> <mhi>     n=5 (1,2,3,4,m) scan: enclosures + exact witnesses + exact
                    strict-local-max cone test at c*
  E <mlist>         n=6 (1,..,5,m) enclosures + exact witnesses + cone test
  digest            merge partials -> out_lrc_zono_rescan.json + console digest
"""
import sys, json, time, math
import numpy as np
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
import lrc_zono_lib as L
from lrc_zono_lib import (lambda_exact, make_forms, form_value, norm_V, cstar,
                           D_exact, rho_exact_d2, rho_semid2, M_corner,
                           enumerate_Z0, _grid_D_numpy)

OUT = "/home/z/my-project/scripts/out_lrc_zono_rescan"

# ---------------------------------------------------------------- cone test (exact FM)

def _fm_feasible(rows, d):
    """Exact Fourier-Motzkin feasibility of {a_i . w <= b_i}.
    rows: list of (coeffs tuple of ints, rhs Fraction)."""
    cur = {}
    def norm(c, b):
        c = tuple(int(x) for x in c)
        g = 0
        for x in c:
            g = math.gcd(g, abs(x))
        if g > 1:
            c = tuple(x // g for x in c)
            b = b / g
        return c, b
    def push(dct, c, b):
        if all(x == 0 for x in c):
            return b >= 0          # 0 <= b: contradiction iff b < 0
        if c not in dct or b < dct[c]:
            dct[c] = b
        return True
    for c, b in rows:
        c, b = norm(c, b)
        if not push(cur, c, b):
            return False
    for v in range(d):
        P, N, Z = [], [], []
        for c, b in cur.items():
            a = c[v]
            if a > 0:
                P.append((c, b))
            elif a < 0:
                N.append((c, b))
            else:
                Z.append((c, b))
        nxt = {}
        ok = True
        for c, b in Z:
            if not push(nxt, c, b):
                ok = False
                break
        if not ok:
            return False
        for ci, bi in P:
            avi = ci[v]
            for cj, bj in N:
                avj = cj[v]
                mi, mj = -avj, avi          # both positive
                cnew = tuple(mi * ci[k] + mj * cj[k] for k in range(d))
                bnew = mi * bi + mj * bj
                cnew, bnew = norm(cnew, bnew)
                if not push(nxt, cnew, bnew):
                    return False
        cur = nxt
        if not cur:
            break
    return True

def cone_trivial(rows, d):
    """{w : a_i . w <= 0 for all i} == {0} ?  (exact).
    Returns (trivial, witness_dir_or_None). witness_dir: float list or None."""
    base = [(tuple(int(x) for x in a), F(0)) for a in rows]
    for k in range(d):
        for s in (1, -1):
            e = [0] * d
            e[k] = -s                     # (-s) w_k <= -1  <=>  s w_k >= 1
            if _fm_feasible(base + [(tuple(e), F(-1))], d):
                w = _cone_witness(rows, d, k, s)
                return False, w
    return True, None

def _cone_witness(rows, d, k, s):
    try:
        from scipy.optimize import linprog
        import numpy as np
        A = np.array([[float(x) for x in a] for a in rows] + [[0.0 if j != k else -float(s)
                                                               for j in range(d)]])
        b = np.array([0.0] * len(rows) + [-1.0])
        r = linprog(c=np.zeros(d), A_ub=A, b_ub=b, bounds=[(None, None)] * d,
                    method="highs")
        if r.x is None:
            return None
        w = list(r.x)
        m = max(abs(t) for t in w) or 1.0
        w = [t / m for t in w]
        return w
    except Exception:
        return None

# ---------------------------------------------------------------- contacts + local max at c*

def contacts_cstar(v, forms, cs, delta):
    """All z in Z^d with ||c*-z||_V == delta, with active forms (value == delta).
    (D(c*) = delta is the min, so <= delta means == delta.)"""
    d = len(v) - 1
    beta = [1 + v[j + 1] for j in range(d)]
    ranges = []
    for kk in range(d):
        a = math.ceil(cs[kk] - delta * beta[kk]) - 1
        b = math.floor(cs[kk] + delta * beta[kk]) + 1
        ranges.append(range(a, b + 1))
    out = []
    import itertools as it
    for z in it.product(*ranges):
        diff = tuple(cs[kk] - z[kk] for kk in range(d))
        mx = F(0)
        vals = []
        for coeffs, denom in forms:
            val = form_value((coeffs, denom), diff)
            vals.append(val)
            a = abs(val)
            if a > mx:
                mx = a
        if mx == delta:
            act = []
            for fi, (coeffs, denom) in enumerate(forms):
                if abs(vals[fi]) == delta:
                    sgn = 1 if vals[fi] > 0 else -1
                    act.append((sgn, coeffs))
            out.append((z, act))
        elif mx < delta:
            raise AssertionError("contact with norm < delta (D(c*) != delta?)")
    return out

def local_max_certificate(v, forms, cs, delta):
    """Exact strict-local-max test at c* via the ascent cone.
    Returns dict: contacts count, trivial cone?, extracted dir, D(c*+eps*w)."""
    d = len(v) - 1
    conts = contacts_cstar(v, forms, cs, delta)
    rows = []
    for z, act in conts:
        for sgn, coeffs in act:
            rows.append(tuple(sgn * c for c in coeffs))
    trivial, w = cone_trivial(rows, d)
    res = {"n_contacts": len(conts), "n_rows": len(rows),
           "cstar_strict_local_max_exact": bool(trivial)}
    if not trivial and w is not None:
        probes = {}
        for eps in (F(1, 64), F(1, 256)):
            x = tuple((cs[k] + eps * F(round(w[k] * 4096), 4096)) % 1 for k in range(d))
            probes[str(eps)] = str(D_exact(x, v, forms))
        res["probe_D_cstar_plus_eps_w"] = probes
        res["delta"] = str(delta)
    return res

# ---------------------------------------------------------------- enclosure (from run3, copied)

def rho_enclosure(v, K, K_loc=25, rounds=3):
    forms = L.make_forms(v)
    d = len(v) - 1
    M = L.M_corner(v, forms)
    Z0 = L.enumerate_Z0(v, forms, M)
    grid = L._grid_D_numpy(v, forms, Z0, K)
    fmax = float(grid.max())
    h = 1.0 / K
    rho_hi = fmax + h / 2 + 1e-6
    tor = []
    for bits in range(2 ** d):
        x = tuple(F(1, 2) if (bits >> k) & 1 else F(0) for k in range(d))
        tor.append((x, L.D_exact(x, v, forms)))
    tor_max = max(v_ for _, v_ in tor)
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
    xr = tuple(F(round(c * 2 ** 21), 2 ** 21) for c in center)
    vr = L.D_exact(xr, v, forms)
    if vr > best_val:
        best_x, best_val = xr, vr
    rho_lo = max(best_val, tor_max)
    return {"rho_hi": rho_hi, "rho_lo_exact": rho_lo, "fmax": fmax,
            "tor_max": str(tor_max), "witness": best_x, "witness_val": str(best_val),
            "Z0": len(Z0), "forms": forms}

# ---------------------------------------------------------------- phases

def base_instance_record(v):
    n = len(v)
    lam, tstar = lambda_exact(v)
    delta = F(1, 2) - lam
    thr = F(n - 1, 2 * (n + 1))
    forms = make_forms(v)
    cs = cstar(v)
    dval = D_exact(cs, v, forms)
    assert dval == delta, f"D(c*)={dval} != delta={delta} for {v}"
    return {"v": list(v), "n": n, "d": n - 1, "lambda": str(lam), "delta": str(delta),
            "thr": str(thr), "D_cstar": str(dval), "asserts": "D(c*)==delta OK"}

def phase_A():
    import numpy as np  # noqa
    out = []
    # n=3 core exact (provenance re-establishment)
    for v, want in [((1, 2, 3), F(1, 4)), ((1, 2, 6), F(3, 14))]:
        forms = make_forms(v)
        rho, arg, nh, nc = rho_exact_d2(v, forms)
        rec = base_instance_record(v)
        assert rho == want, f"{v}: rho={rho} != persisted {want}"
        rec.update({"rho_exact": str(rho), "argmax": [str(a) for a in arg],
                    "deepest_hole": rho == rec_delta(rec), "note": "re-asserted exact"})
        out.append(rec)
    # n=5/6 refutation witnesses (exact)
    wit = [
        ((1, 2, 3, 4, 5), (F(1, 32), F(15, 32), F(3, 32), F(7, 8)), F(97, 288)),
        ((1, 2, 3, 4, 10), (F(1, 2), F(0), F(1, 2), F(1, 2)), F(7, 22)),
        ((1, 2, 3, 4, 5, 6), (F(0), F(5, 16), F(15, 16), F(3, 4), F(1, 2)), F(29, 80)),
        ((1, 2, 3, 4, 5, 12), (F(1, 16), F(11, 16), F(1, 8), F(7, 16), F(7, 8)), F(5, 14)),
    ]
    for v, x, want in wit:
        forms = make_forms(v)
        val = D_exact(x, v, forms)
        rec = base_instance_record(v)
        assert val == want, f"{v}: witness value {val} != recorded {want}"
        delta = F(1, 2) - lambda_exact(v)[0]
        rec.update({"witness": [str(t) for t in x], "witness_val": str(val),
                    "witness_beats_delta": val > delta,
                    "witness_beats_thr": val > rec_thr(rec)})
        out.append(rec)
    with open(f"{OUT}_A.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    for r in out:
        print(f"v={r['v']}: {r.get('note', '')} rho={r.get('rho_exact', r.get('witness_val'))}"
              f" delta={r['delta']}", flush=True)
    print("PHASE A DONE")

def rec_delta(rec):
    return F(*map(int, rec["delta"].split("/"))) if "/" in rec["delta"] else F(int(rec["delta"]))

def rec_thr(rec):
    return F(*map(int, rec["thr"].split("/"))) if "/" in rec["thr"] else F(int(rec["thr"]))

def phase_B(mlo, mhi):
    out = []
    for m in range(mlo, mhi + 1):
        v = (1, 2, m)
        t0 = time.time()
        forms = make_forms(v)
        lam, _ = lambda_exact(v)
        delta = F(1, 2) - lam
        thr = F(1, 4)
        rho, arg, nh, nc = rho_exact_d2(v, forms)
        rec = {"v": list(v), "m": m, "n": 3, "lambda": str(lam), "delta": str(delta),
               "rho_exact": str(rho), "deepest_hole": rho == delta,
               "rho_leq_thr": rho <= thr, "n_div_m": (m % 3 == 0),
               "argmax": [str(a) for a in arg], "n_hyp": nh, "sec": round(time.time() - t0, 1)}
        out.append(rec)
        print(f"(1,2,{m:2d}): rho={rho} delta={delta} deepest={rec['deepest_hole']} "
              f"3|m={rec['n_div_m']} agree={rec['deepest_hole'] == rec['n_div_m']} "
              f"[{rec['sec']}s]", flush=True)
    with open(f"{OUT}_B_{mlo}_{mhi}.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("PHASE B DONE")

def _C_save(rec, mode):
    fn = f"{OUT}_C_{mode}.json"
    try:
        data = json.load(open(fn))
    except Exception:
        data = []
    data.append(rec)
    with open(fn, "w") as f:
        json.dump(data, f, indent=1, default=str)

def phase_C(mlist, K=48):
    """n=4 (1,2,3,m): hold-candidates (m % 4 == 0) get the full away-certified
    rho_semid2 treatment; fail-entries get the cheap exact-witness assertion
    (enclosure + local refinement, exact re-evaluated) -- a fail needs only a
    witness > delta, which is exact and rigorous."""
    import numpy as np
    for m in mlist:
        v = (1, 2, 3, m)
        t0 = time.time()
        forms = make_forms(v)
        lam, _ = lambda_exact(v)
        delta = F(1, 2) - lam
        thr = F(3, 10)
        cs = cstar(v)
        assert D_exact(cs, v, forms) == delta
        if m % 4 == 0:
            res = rho_semid2(v, forms, K=K)
            rho_lo = res["rho_lo_exact"]
            # independent cross-check: certified exact value must sit under any grid bound
            if res["rho_exact"] is not None:
                M = M_corner(v, forms)
                Z0 = enumerate_Z0(v, forms, M)
                g72 = _grid_D_numpy(v, forms, Z0, 72)
                hi72 = float(g72.max()) + 1.0 / 144 + 1e-6
                assert float(res["rho_exact"]) <= hi72 + 1e-9, \
                    f"{v}: exact {res['rho_exact']} > grid72 bound {hi72}"
            rec = {"v": list(v), "m": m, "n": 4, "lambda": str(lam), "delta": str(delta),
                   "thr": str(thr), "mode": "semid2-exact",
                   "rho_lo_exact": str(rho_lo),
                   "rho_exact": str(res["rho_exact"]) if res["rho_exact"] is not None else None,
                   "away_certified": bool(res["away_certified"]),
                   "rho_hi_grid": res["rho_hi_grid"], "grid_K": K,
                   "deepest_hole_certified": (res["rho_exact"] is not None
                                              and res["rho_exact"] == delta),
                   "deepest_hole_refuted": rho_lo > delta,
                   "rho_leq_thr": (res["rho_exact"] is not None
                                   and res["rho_exact"] <= thr),
                   "rho_leq_thr_lo": rho_lo <= thr,
                   "n_div_m": (m % 4 == 0),
                   "argmax": [str(a) for a in res["argmax"]] if res["argmax"] else None,
                   "n_boxes": res["n_boxes"], "Z0": res["Z0_size"],
                   "sec": round(time.time() - t0, 1)}
            _C_save(rec, "exact")
            print(f"(1,2,3,{m:2d}) [exact]: delta={delta} rho_lo={rho_lo} "
                  f"exact={rec['rho_exact']} away={rec['away_certified']} "
                  f"deepest_cert={rec['deepest_hole_certified']} "
                  f"refuted={rec['deepest_hole_refuted']} 4|m=True [{rec['sec']}s]",
                  flush=True)
        else:
            enc = rho_enclosure(v, K=40)
            rho_lo = enc["rho_lo_exact"]
            forms_e = enc.pop("forms")
            rec = {"v": list(v), "m": m, "n": 4, "lambda": str(lam), "delta": str(delta),
                   "thr": str(thr), "mode": "witness-fail",
                   "rho_lo_exact": str(rho_lo), "rho_hi": enc["rho_hi"],
                   "witness": [str(t) for t in enc["witness"]],
                   "witness_val": str(enc["witness_val"]),
                   "tor_max": str(enc["tor_max"]),
                   "deepest_hole_refuted": rho_lo > delta,
                   "rho_leq_thr_lo": rho_lo <= thr,
                   "rho_gt_thr": rho_lo > thr,
                   "n_div_m": False,
                   "sec": round(time.time() - t0, 1)}
            _C_save(rec, "fail")
            flag = "" if rho_lo > delta else "  <<< WARNING: no witness beats delta"
            print(f"(1,2,3,{m:2d}) [witness]: delta={delta} rho_lo={rho_lo} "
                  f"wit={rec['witness_val']} refuted={rec['deepest_hole_refuted']} "
                  f"4|m=False [{rec['sec']}s]{flag}", flush=True)
    print("PHASE C DONE")

def phase_D(mlo, mhi, K=32):
    out = []
    for m in range(mlo, mhi + 1):
        v = (1, 2, 3, 4, m)
        t0 = time.time()
        lam, _ = lambda_exact(v)
        delta = F(1, 2) - lam
        thr = F(4, 12)
        res = rho_enclosure(v, K)
        rho_lo = res["rho_lo_exact"]
        forms = res.pop("forms")
        cs = cstar(v)
        lmc = local_max_certificate(v, forms, cs, delta)
        rec = {"v": list(v), "m": m, "n": 5, "lambda": str(lam), "delta": str(delta),
               "thr": str(thr),
               "rho_lo_exact": str(rho_lo), "rho_hi": res["rho_hi"],
               "fmax": res["fmax"], "tor_max": res["tor_max"],
               "witness": [str(t) for t in res["witness"]],
               "witness_val": res["witness_val"],
               "deepest_hole_refuted": rho_lo > delta,
               "deepest_hole_consistent": (rho_lo == delta
                                           and res["fmax"] <= float(delta) + 1e-7),
               "rho_gt_thr": rho_lo > thr,
               "rho_geq_thr": rho_lo >= thr,
               "n_div_m": (m % 5 == 0),
               "local_max": lmc,
               "Z0": res["Z0"], "sec": round(time.time() - t0, 1)}
        out.append(rec)
        print(f"(1,2,3,4,{m:2d}): delta={delta} rho_lo={rho_lo} "
              f"refuted={rec['deepest_hole_refuted']} cons={rec['deepest_hole_consistent']} "
              f"5|m={rec['n_div_m']} c*localmax={lmc['cstar_strict_local_max_exact']} "
              f"[{rec['sec']}s]", flush=True)
    with open(f"{OUT}_D_{mlo}_{mhi}.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("PHASE D DONE")

def _E_save(rec):
    fn = f"{OUT}_E.json"
    try:
        data = json.load(open(fn))
    except Exception:
        data = []
    data.append(rec)
    with open(fn, "w") as f:
        json.dump(data, f, indent=1, default=str)

def phase_E(mlist, K=16):
    out = []
    for m in mlist:
        v = (1, 2, 3, 4, 5, m)
        t0 = time.time()
        lam, _ = lambda_exact(v)
        delta = F(1, 2) - lam
        thr = F(5, 14)
        res = rho_enclosure(v, K)
        rho_lo = res["rho_lo_exact"]
        forms = res.pop("forms")
        cs = cstar(v)
        lmc = local_max_certificate(v, forms, cs, delta)
        rec = {"v": list(v), "m": m, "n": 6, "lambda": str(lam), "delta": str(delta),
               "thr": str(thr),
               "rho_lo_exact": str(rho_lo), "rho_hi": res["rho_hi"],
               "fmax": res["fmax"], "tor_max": res["tor_max"],
               "witness": [str(t) for t in res["witness"]],
               "witness_val": res["witness_val"],
               "deepest_hole_refuted": rho_lo > delta,
               "deepest_hole_consistent": (rho_lo == delta
                                           and res["fmax"] <= float(delta) + 1e-7),
               "rho_gt_thr": rho_lo > thr,
               "n_div_m": (m % 6 == 0),
               "local_max": lmc,
               "Z0": res["Z0"], "sec": round(time.time() - t0, 1)}
        out.append(rec)
        _E_save(rec)
        print(f"(1..5,{m:2d}): delta={delta} rho_lo={rho_lo} "
              f"refuted={rec['deepest_hole_refuted']} cons={rec['deepest_hole_consistent']} "
              f"6|m={rec['n_div_m']} c*localmax={lmc['cstar_strict_local_max_exact']} "
              f"[{rec['sec']}s]", flush=True)
    print("PHASE E DONE")

def phase_digest():
    import glob
    allp = {}
    for fn in sorted(glob.glob(f"{OUT}_A.json")):
        allp["A"] = json.load(open(fn))
    for fn in sorted(glob.glob(f"{OUT}_B_*.json")):
        allp.setdefault("B", []).extend(json.load(open(fn)))
    for mode in ("exact", "fail"):
        try:
            allp.setdefault("C", []).extend(json.load(open(f"{OUT}_C_{mode}.json")))
        except Exception:
            pass
    for fn in sorted(glob.glob(f"{OUT}_D_*.json")):
        allp.setdefault("D", []).extend(json.load(open(fn)))
    try:
        allp["E"] = json.load(open(f"{OUT}_E.json"))
    except Exception:
        pass
    with open(f"{OUT}.json", "w") as f:
        json.dump(allp, f, indent=1, default=str)
    print("== DIGEST ==")
    def _status(r):
        if r.get("deepest_hole_certified"):
            return "CERTIFIED-HOLD"
        if r.get("deepest_hole_refuted"):
            return "REFUTED(exact witness)"
        if r.get("deepest_hole_consistent"):
            return "consistent(uncertified)"
        if "deepest_hole" in r:
            return "EXACT-HOLD" if r["deepest_hole"] else "EXACT-FAIL"
        return "?"
    for tag in ("B", "C", "D", "E"):
        for r in sorted(allp.get(tag, []), key=lambda r: (r.get("n", 0), r.get("m", 0))):
            extra = ""
            if r.get("rho_gt_thr"):
                extra = "  [rho>thr: COVERING REFUTED]"
            print(f"{tag} n={r['n']} m={r.get('m'):>2}: {_status(r):24s} n|m={str(r.get('n_div_m')):5s}{extra}")
    print("saved ->", f"{OUT}.json")

if __name__ == "__main__":
    import numpy as np  # noqa: F401  (lib needs it)
    ph = sys.argv[1] if len(sys.argv) > 1 else "A"
    if ph == "A":
        phase_A()
    elif ph == "B":
        phase_B(int(sys.argv[2]), int(sys.argv[3]))
    elif ph == "C":
        phase_C([int(x) for x in sys.argv[2].split(",")])
    elif ph == "D":
        phase_D(int(sys.argv[2]), int(sys.argv[3]))
    elif ph == "E":
        phase_E([int(x) for x in sys.argv[2].split(",")])
    elif ph == "digest":
        phase_digest()
