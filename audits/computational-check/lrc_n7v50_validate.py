#!/usr/bin/env python3
"""
lrc_n7v50_validate.py -- validation layers for the sharded n=7, V=50 sweep
(16 shards, C(50,7) = 99,884,400 vectors), and census extraction.

Layers (resumable; run as:  python3 lrc_n7v50_validate.py A <lo> <hi> [shard idxs]
                                 python3 lrc_n7v50_validate.py B|C|D|E|F):
  A  per-shard full re-verification of every dump record against ALL polytope
     constraints (numpy, exact integers), plus structural checks (v increasing,
     box bounds, feasible flag == gap*(n+1) >= 1, witness == floor(t v)),
     rank-geometry checks (first/last record == unrank(lo)/unrank(hi-1)),
     and census extraction (tight / rung-2 (2/15) / inter-rung (1/8, 2/15)).
     Marker file per shard: out/n7_V50_s<i>_verify.json.
  B  aggregate: totals == C(50,7), spectrum == sum of shard spectra, tight
     census == predicted 14 (AP d<=7, A d<=4, B d<=3), inter-rung census ==
     the 4 known vectors, rung-2 list, V=38-vs-V=50 reconciliation.
  C  independent exact reference solver (pure Python, imported from
     lrc_gap_family_check.py) on every sample.txt record (~1000) + every
     tight / inter-rung / rung-2 vector.
  D  HiGHS MILP cross-solve (scipy) on tight + inter-rung + random samples.
  E  Theorem 6 anchors: (1,2,3,4,5,6,14) and multiples have gap exactly
     2/15 in the dump records.
  F  sample.txt <-> dump.bin consistency per shard.
"""
import glob
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

SCR = Path("/home/z/my-project/scripts")
OUT = SCR / "out"
sys.path.insert(0, str(SCR))
from lrc_gap_family_check import gap_exact_pure  # noqa: E402

N, V, SHARDS = 7, 50, 16
TOTAL = math.comb(V, N)
REC_DTYPE = np.dtype([
    ("v", "u1", 8), ("k", "u1", 8),
    ("ta", "<u2"), ("tb", "<u2"), ("num", "<u2"), ("den", "<u2"),
    ("feas", "u1"), ("pad", "u1", 7)])


def shard_range(idx):
    base, rem = divmod(TOTAL, SHARDS)
    lo = idx * base + min(idx, rem)
    hi = (idx + 1) * base + min(idx + 1, rem)
    return lo, hi


def unrank(r, n, V):
    v, x = [], 1
    for pos in range(n):
        while True:
            cnt = math.comb(V - x, n - 1 - pos) if V - x >= n - 1 - pos else 0
            if r < cnt:
                v.append(x)
                x += 1
                break
            r -= cnt
            x += 1
    return v


def verify_point_py(vlist, klist):
    n = len(vlist)
    for i in range(n):
        if not (0 <= klist[i] <= vlist[i] - 1):
            return False
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if (n + 1) * (vlist[j] * klist[i] - vlist[i] * klist[j]) \
                    > n * vlist[i] - vlist[j]:
                return False
    return True


# ------------------------------------------------------------------ A
def layer_A(indices):
    t0 = time.time()
    out = {}
    for idx in indices:
        marker = OUT / f"n7_V50_s{idx}_verify.json"
        if marker.exists():
            out[idx] = json.loads(marker.read_text())
            continue
        lo, hi = shard_range(idx)
        cnt = hi - lo
        path = OUT / f"n7_V50_s{idx}_dump.bin"
        fsize = path.stat().st_size
        assert fsize == cnt * 32, (idx, fsize, cnt * 32)
        tight_l, rung2_l, inter_l = [], [], []
        nfeas = nviol = 0
        with open(path, "rb") as f:
            first = last = None
            pos = 0
            while True:
                buf = f.read(500_000 * 32)
                if not buf:
                    break
                recs = np.frombuffer(buf, dtype=REC_DTYPE)
                m = len(recs)
                v = recs["v"][:, :N].astype(np.int64)
                k = recs["k"][:, :N].astype(np.int64)
                num = recs["num"].astype(np.int64)
                den = recs["den"].astype(np.int64)
                feas = recs["feas"].astype(np.int64)
                ta = recs["ta"].astype(np.int64)
                tb = recs["tb"].astype(np.int64)
                # structural
                assert (np.diff(v, axis=1) > 0).all()
                assert (v >= 1).all() and (v <= V).all()
                assert (k >= 0).all() and (k <= v - 1).all()
                # feasible flag consistency
                assert (feas == (num * (N + 1) >= den)).all()
                # all polytope constraints on the stored witness
                bad = np.zeros(m, dtype=bool)
                for i in range(N):
                    ki = k[:, i]
                    vi = v[:, i]
                    for j in range(N):
                        if i == j:
                            continue
                        D = v[:, j] * ki - vi * k[:, j]
                        bad |= (N + 1) * D > N * vi - v[:, j]
                nviol += int(bad.sum())
                # witness == floor(t v)
                assert ((k * tb[:, None]) <= (ta[:, None] * v)).all()
                assert ((ta[:, None] * v) < ((k + 1) * tb[:, None])).all()
                nfeas += int(feas.sum())
                # census
                tsel = num * (N + 1) == den
                rsel = num * 15 == 2 * den
                isel = (num * (N + 1) > den) & (num * 15 < 2 * den)
                for sel, acc in ((tsel, tight_l), (rsel, rung2_l),
                                 (isel, inter_l)):
                    for r in recs[sel]:
                        acc.append([int(x) for x in r["v"][:N]])
                if pos == 0:
                    first = [int(x) for x in recs["v"][0][:N]]
                last = [int(x) for x in recs["v"][-1][:N]]
                pos += m
        assert pos == cnt, (idx, pos, cnt)
        # rank geometry: dump order == lexicographic rank interval
        assert first == unrank(lo, N, V), (idx, first, unrank(lo, N, V))
        assert last == unrank(hi - 1, N, V), (idx, last)
        res = {"shard": idx, "count": cnt, "feasible": nfeas,
               "constraint_violations": nviol, "tight": tight_l,
               "rung2_2_15": rung2_l, "inter_rung": inter_l}
        marker.write_text(json.dumps(res))
        out[idx] = res
        print(f"  shard {idx}: {cnt} records, feasible={nfeas}, "
              f"violations={nviol}, tight={len(tight_l)}, "
              f"2/15={len(rung2_l)}, inter={len(inter_l)}", flush=True)
    print(f"layer A PASS for shards {sorted(out)} "
          f"({round(time.time()-t0,1)}s)")
    return out


# ------------------------------------------------------------------ B
def layer_B():
    t0 = time.time()
    A = {}
    for idx in range(SHARDS):
        A[idx] = json.loads((OUT / f"n7_V50_s{idx}_verify.json").read_text())
    tot = sum(a["count"] for a in A.values())
    feas = sum(a["feasible"] for a in A.values())
    viol = sum(a["constraint_violations"] for a in A.values())
    assert tot == TOTAL == 99_884_400
    assert feas == tot and viol == 0
    tight = [tuple(v) for a in A.values() for v in a["tight"]]
    rung2 = [tuple(v) for a in A.values() for v in a["rung2_2_15"]]
    inter = [tuple(v) for a in A.values() for v in a["inter_rung"]]
    # predicted tight census
    pred = []
    for d in range(1, 8):
        pred.append(tuple(d * i for i in range(1, 8)))
    for d in range(1, 5):
        pred.append(tuple(d * x for x in (1, 2, 3, 4, 5, 7, 12)))
    for d in range(1, 4):
        pred.append(tuple(d * x for x in (1, 4, 5, 6, 7, 11, 13)))
    assert sorted(tight) == sorted(pred), (sorted(tight), sorted(pred))
    # inter-rung: exactly the four known vectors
    known = sorted([(1, 2, 3, 4, 5, 7, 18), (2, 4, 6, 8, 10, 14, 36),
                    (1, 3, 4, 5, 7, 13, 18), (2, 6, 8, 10, 14, 26, 36)])
    assert sorted(inter) == known, sorted(inter)
    # rung-2 must contain the Theorem 6 family multiples d = 1..3
    fam = [(1, 2, 3, 4, 5, 6, 14), (2, 4, 6, 8, 10, 12, 28),
           (3, 6, 9, 12, 15, 18, 42)]
    for fvec in fam:
        assert fvec in rung2, fvec
    # spectrum aggregate vs shard spectra
    spec = {}
    for fp in sorted(glob.glob(str(OUT / "n7_V50_s*_spectrum.txt"))):
        for line in open(fp):
            num, den, c = line.split()
            key = (int(num), int(den))
            spec[key] = spec.get(key, 0) + int(c)
    assert sum(spec.values()) == TOTAL
    assert spec.get((1, 8)) == len(tight) == 14
    assert spec.get((3, 23)) == len(inter) == 4
    assert spec.get((2, 15)) == len(rung2)
    # V=38 -> V=50 reconciliation of the bottom spectrum
    v38 = {"1/8": 10, "3/23": 4, "2/15": 21, "5/37": 1, "3/22": 1,
           "4/29": 5, "5/36": 2}
    v50 = {"1/8": spec.get((1, 8), 0), "3/23": spec.get((3, 23), 0),
           "2/15": spec.get((2, 15), 0), "5/37": spec.get((5, 37), 0),
           "3/22": spec.get((3, 22), 0), "4/29": spec.get((4, 29), 0),
           "5/36": spec.get((5, 36), 0)}
    print("  bottom spectrum V=38 -> V=50:", v38, "->", v50)
    print(f"  rung-2 (2/15) vectors ({len(rung2)}):")
    for r in sorted(rung2):
        print("   ", r)
    res = {"layer": "B", "status": "PASS", "total": tot, "feasible": feas,
           "violations": viol, "tight": len(tight), "inter_rung": len(inter),
           "rung2_2_15": len(rung2), "v38_bottom": v38, "v50_bottom": v50,
           "rung2_list": [list(r) for r in sorted(rung2)],
           "seconds": round(time.time() - t0, 1)}
    print("layer B PASS:", json.dumps({k: v for k, v in res.items()
                                       if k not in ("rung2_list",)}))
    return res


# ------------------------------------------------------------------ C
def layer_C():
    t0 = time.time()
    vecs = []
    # every sample record
    nsamp = 0
    for fp in sorted(glob.glob(str(OUT / "n7_V50_s*_sample.txt"))):
        for line in open(fp):
            d = {}
            for tok in line.split():
                key, _, val = tok.partition("=")
                d[key] = val
            vecs.append(([int(x) for x in d["v"].strip("()").split(",")],
                         [int(x) for x in d["k"].strip("()").split(",")],
                         tuple(int(x) for x in d["t"].split("/")),
                         tuple(int(x) for x in d["gap"].split("/")),
                         int(d["feasible"])))
            nsamp += 1
    # tight / inter / rung-2 vectors from layer A markers
    A = {}
    for idx in range(SHARDS):
        A[idx] = json.loads((OUT / f"n7_V50_s{idx}_verify.json").read_text())
    census = [tuple(v) for a in A.values() for key in
              ("tight", "rung2_2_15", "inter_rung") for v in a[key]]
    checked = 0
    for item in vecs:
        vlist, klist, tfrac, gapfrac, feas = item
        bn, bd = gap_exact_pure(vlist)
        assert bn * gapfrac[1] == gapfrac[0] * bd, (vlist, bn, bd, gapfrac)
        assert feas == 1
        assert verify_point_py(vlist, klist)
        checked += 1
    for vlist in census:
        bn, bd = gap_exact_pure(vlist)
        assert bn * 8 >= bd, vlist            # feasible (LRC holds on it)
        checked += 1
    res = {"layer": "C", "status": "PASS", "samples_rechecked": nsamp,
           "census_vectors_rechecked": len(census),
           "total_reference_solves": checked,
           "seconds": round(time.time() - t0, 1)}
    print("layer C PASS:", json.dumps(res))
    return res


# ------------------------------------------------------------------ D
def layer_D(n_milp=170):
    from scipy.optimize import milp, LinearConstraint, Bounds
    import random
    t0 = time.time()
    random.seed(20261002)
    A = {}
    for idx in range(SHARDS):
        A[idx] = json.loads((OUT / f"n7_V50_s{idx}_verify.json").read_text())
    census = [list(v) for a in A.values() for key in
              ("tight", "rung2_2_15", "inter_rung") for v in a[key]]
    samples = []
    for fp in sorted(glob.glob(str(OUT / "n7_V50_s*_sample.txt"))):
        for line in open(fp):
            for tok in line.split():
                if tok.startswith("v="):
                    samples.append([int(x) for x in
                                    tok[2:].strip("()").split(",")])
    picks = census + random.sample(samples, min(n_milp - len(census),
                                                len(samples)))
    rows = []
    for vlist in picks:
        n = len(vlist)
        Amat, ub = [], []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                row = [0] * n
                row[i] = 8 * vlist[j]
                row[j] = -8 * vlist[i]
                Amat.append(row)
                ub.append(7 * vlist[i] - vlist[j])
        con = LinearConstraint(np.array(Amat), -np.inf, np.array(ub))
        bnd = Bounds([0] * n, [x - 1 for x in vlist])
        r = milp(c=np.zeros(n), constraints=con, bounds=bnd,
                 integrality=np.ones(n))
        assert r.status == 0, (vlist, r.status, r.message)
        ksol = np.round(r.x).astype(int)
        assert verify_point_py(vlist, list(ksol))
        rows.append(vlist)
    res = {"layer": "D", "status": "PASS", "milp_solves": len(rows),
           "seconds": round(time.time() - t0, 1)}
    print("layer D PASS:", json.dumps(res))
    return res


# ------------------------------------------------------------------ E
def layer_E():
    t0 = time.time()
    A = {}
    for idx in range(SHARDS):
        A[idx] = json.loads((OUT / f"n7_V50_s{idx}_verify.json").read_text())
    rung2 = [tuple(v) for a in A.values() for v in a["rung2_2_15"]]
    for d in (1, 2, 3):
        fam = tuple(d * x for x in (1, 2, 3, 4, 5, 6, 14))
        assert fam in rung2
        bn, bd = gap_exact_pure(list(fam))
        assert bn * 15 == 2 * bd
    res = {"layer": "E", "status": "PASS",
           "theorem6_anchors": ["(1,2,3,4,5,6,14)", "(2,...,28)",
                                "(3,...,42)"],
           "seconds": round(time.time() - t0, 1)}
    print("layer E PASS:", json.dumps(res))
    return res


# ------------------------------------------------------------------ F
def layer_F():
    t0 = time.time()
    checked = 0
    for idx in range(SHARDS):
        lo, hi = shard_range(idx)
        samples = {}
        for line in open(OUT / f"n7_V50_s{idx}_sample.txt"):
            d = {}
            for tok in line.split():
                key, _, val = tok.partition("=")
                d[key] = val
            vlist = [int(x) for x in d["v"].strip("()").split(",")]
            samples[tuple(vlist)] = d
        with open(OUT / f"n7_V50_s{idx}_dump.bin", "rb") as f:
            pos = 0
            while True:
                buf = f.read(99_733 * 32)
                if not buf:
                    break
                recs = np.frombuffer(buf, dtype=REC_DTYPE)
                for r in recs[:1]:
                    vlist = tuple(int(x) for x in r["v"][:N])
                    if vlist in samples:
                        s = samples[vlist]
                        assert s["t"] == f"{int(r['ta'])}/{int(r['tb'])}"
                        assert s["gap"] == f"{int(r['num'])}/{int(r['den'])}"
                        assert s["feasible"] == str(int(r["feas"]))
                        checked += 1
                pos += len(recs)
    res = {"layer": "F", "status": "PASS", "sample_dump_matches": checked,
           "seconds": round(time.time() - t0, 1)}
    print("layer F PASS:", json.dumps(res))
    return res


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "ABCDEF"
    if which == "A":
        if len(sys.argv) >= 4:
            idxs = list(range(int(sys.argv[2]), int(sys.argv[3])))
        else:
            idxs = list(range(SHARDS))
        layer_A(idxs)
        return
    results = []
    for ch in which:
        results.append(globals()[f"layer_{ch}"]())
    with open(OUT / "n7v50_validation.jsonl", "a") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")


if __name__ == "__main__":
    main()
