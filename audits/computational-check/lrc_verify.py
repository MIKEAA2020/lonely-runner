#!/usr/bin/env python3
"""
Cross-validation and analysis driver for the LRC ILP computational check.

Layers of validation:
  L1  Full numpy re-verification of every witness k in every dump record.
  L2  Independent exact reference implementation of gap(v) (pure Python,
      integer arithmetic) on random samples, compared record-by-record.
  L3  Literal ILP enumeration (C brute mode over the entire k-box) for all
      v with v_n <= 12, compared with the run-mode verdicts.
  L4  HiGHS MILP (scipy.optimize.milp) on random samples + doc rows.
  L5  Verification of the source audit document's hand tables and claims,
      including the overlap-improvement bound (doc next-step #2).
  L6  Structural checks: scaling invariance, bounded-speed k=0 validity,
      LP-relaxation explicit point k_i = (v_i - v_1)/((n+1) v_1).

Outputs: out/stats.json, out/analysis.txt
"""
import json
import hashlib
import math
import random
import subprocess
from fractions import Fraction
from itertools import combinations
from collections import Counter

import numpy as np

OUT = "/home/z/my-project/scripts/out"
BIN = "/home/z/my-project/scripts/lrc_ilp_check"
V = 50
SEED = 20261002

DTYPE = np.dtype([
    ("v", "u1", (8,)),
    ("k", "u1", (8,)),
    ("ta", "<u2"), ("tb", "<u2"),
    ("num", "<u2"), ("den", "<u2"),
    ("feasible", "u1"), ("pad", "u1", (7,)),
])
assert DTYPE.itemsize == 32


# ---------------------------------------------------------------- basics
def load(n):
    a = np.fromfile(f"{OUT}/n{n}_V50_dump.bin", dtype=DTYPE)
    return a


def verify_point(v, k):
    """Exact integer check of k in P(v)."""
    n = len(v)
    for i in range(n):
        if k[i] < 0 or k[i] > v[i] - 1:
            return False
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if (n + 1) * (v[j] * k[i] - v[i] * k[j]) > n * v[i] - v[j]:
                return False
    return True


def tight_pairs(v, k):
    """Pairs (i,j) where the doc's condition holds with exact equality:
    v_i(1+k_j) - v_j k_i == (v_i + v_j)/(n+1)."""
    n = len(v)
    out = []
    for i in range(n):
        for j in range(n):
            if (n + 1) * (v[i] * (1 + k[j]) - v[j] * k[i]) == v[i] + v[j]:
                out.append((i + 1, j + 1))
    return out


def min_dist(a, b, v):
    """min_i ||a*v_i/b|| as Fraction."""
    m = min(min((a * x) % b, b - (a * x) % b) for x in v)
    return Fraction(m, b)


def gap_ref(v):
    """Independent exact reference: max over candidate t of min_i ||t v_i||."""
    n = len(v)
    cands = set()
    for x in v:
        for a in range(1, 2 * x, 2):
            cands.add((a, 2 * x))
    for x, y in combinations(v, 2):
        for a in range(0, x + y + 1):
            cands.add((a, x + y))
        for a in range(0, y - x + 1):
            cands.add((a, y - x))
    bn, bd = 0, 1
    for a, b in cands:
        m = min(min((a * x) % b, b - (a * x) % b) for x in v)
        if m * bd > bn * b:
            bn, bd = m, b
    return Fraction(bn, bd)


def milp_feasible(v):
    """HiGHS MILP feasibility solve; returns witness tuple or None."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = len(v)
    rows, ubs = [], []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            row = np.zeros(n)
            row[i] = (n + 1) * v[j]
            row[j] = -(n + 1) * v[i]
            rows.append(row)
            ubs.append(n * v[i] - v[j])
    A = np.array(rows, dtype=float) if rows else None
    cons = LinearConstraint(A, -np.inf, np.array(ubs, float)) if rows else None
    res = milp(c=np.zeros(n),
               constraints=cons,
               integrality=np.ones(n),
               bounds=Bounds(np.zeros(n), np.array(v) - 1.0))
    if not res.success or res.x is None:
        return None
    k = tuple(int(round(x)) for x in res.x)
    assert verify_point(v, k), "HiGHS returned an infeasible point!"
    return k


# ---------------------------------------------------------------- L1
def l1_full_recheck(n, a):
    v = a["v"][:, :n].astype(np.int64)
    k = a["k"][:, :n].astype(np.int64)
    N = len(a)
    ok = np.ones(N, dtype=bool)
    ok &= (k >= 0).all(axis=1)
    ok &= (k <= v - 1).all(axis=1)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            lhs = (n + 1) * (v[:, j] * k[:, i] - v[:, i] * k[:, j])
            rhs = n * v[:, i] - v[:, j]
            ok &= lhs <= rhs
    bad = int((~ok).sum())
    feas_ok = bool((a["feasible"].astype(np.int64) *
                    (a["num"].astype(np.int64) * (n + 1)
                     >= a["den"].astype(np.int64))).all())
    return bad, feas_ok


# ---------------------------------------------------------------- L3
def l3_brute_crosscheck(n, vmax=12):
    """Literal k-box enumeration vs run-mode verdicts for v_n <= vmax."""
    lines, vs = [], []
    for comb in combinations(range(1, vmax + 1), n):
        vs.append(comb)
        lines.append(" ".join(map(str, comb)))
    proc = subprocess.run([BIN, "brute", str(n)], input="\n".join(lines),
                          capture_output=True, text=True)
    res = {}
    for line in proc.stdout.splitlines():
        left, right = line.split(") feasible_integer_points=")
        vv = tuple(int(x) for x in left.split("(")[1].split(","))
        cnt = int(right.split(" ")[0])
        first = None
        if "first_k=(" in line:
            first = tuple(int(x) for x in
                          line.split("first_k=(")[1].split(")")[0].split(","))
        res[vv] = (cnt, first)
    dump = {tuple(r["v"][:n]): r for r in load(n)}
    agree = sum(1 for vv in vs
                if (res[vv][0] > 0) == bool(dump[vv]["feasible"]))
    return dict(n=n, vmax=vmax, vectors=len(vs), agree=agree,
                mismatches=[list(vv) for vv in vs
                            if (res[vv][0] > 0) != bool(dump[vv]["feasible"])],
                total_solutions=sum(c for c, _ in res.values()))


# ---------------------------------------------------------------- L5 doc checks
DOC_TABLE1 = [  # (v, claimed k, claimed lonely time)
    ((1, 2), (0, 0), Fraction(1, 3)),
    ((2, 3), (0, 0), Fraction(1, 6)),
    ((1, 2, 3), (0, 0, 0), Fraction(1, 4)),
    ((1, 3, 4), (0, 0, 1), Fraction(3, 7)),
]
DOC_TABLE2 = [  # (v, claimed k, claimed tight annotation)
    ((1, 2), (0, 0), "(1,2) equality"),
    ((1, 3), (0, 1), "(1,2) equality"),
    ((1, 2, 3), (0, 0, 0), "(1,3) equality"),
    ((1, 3, 5), (0, 0, 1), "(1,3) equality"),
    ((1, 2, 3, 5), (0, 0, 0, 1), None),
    ((1, 2, 5, 10), (0, 0, 1, 2), None),
    ((1, 3, 7, 15), (0, 0, 1, 3), None),
    ((1, 4, 9, 20), (0, 1, 2, 6), "(4,3) near-tight"),
    ((1, 5, 9, 15), (0, 1, 2, 3), None),
]


def l5_doc_checks(dumps):
    out = {"table1": [], "table2": [], "errata": []}
    for v, k, t in DOC_TABLE1:
        n = len(v)
        kw = verify_point(v, k)
        tw = min_dist(t.numerator, t.denominator, v) >= Fraction(1, n + 1)
        row = dict(v=list(v), k=list(k), t=f"{t.numerator}/{t.denominator}",
                   k_valid=bool(kw), t_valid=bool(tw),
                   gap=str(gap_ref(v)))
        out["table1"].append(row)
        if not kw:
            out["errata"].append(
                f"Doc Table 1 row v={v}: claimed witness k={k} is NOT feasible in "
                f"P(v) (it violates constraint (i,j)=(3,2): 3*1 - 4*0 = 3 > 9/4). "
                f"The claimed lonely time t={t} itself IS valid; a correct witness "
                f"is (0,1,1), and indeed gap(v) = {row['gap']}.")
    for v, k, ann in DOC_TABLE2:
        n = len(v)
        kw = verify_point(v, k)
        tp = tight_pairs(v, k)
        row = dict(v=list(v), k=list(k), k_valid=bool(kw),
                   actual_tight_pairs=[list(p) for p in tp],
                   claimed=ann)
        out["table2"].append(row)
        if ann and "equality" in ann:
            ci, cj = int(ann[1]), int(ann[3])
            claimed = (ci, cj) in tp
            if not claimed:
                out["errata"].append(
                    f"Doc Table 2 row v={v}: claimed tight constraint {ann}, but the "
                    f"actual equality pair(s) for k={k} are {tp}.")
    return out


def l5_claim_checks(dumps):
    res = {}
    # bounded-speed: v_n <= n v_1  => k=0 valid and feasible
    for n in N_LIST:
        a = dumps[n]
        v = a["v"][:, :n].astype(np.int64)
        bs = v[:, n - 1] <= n * v[:, 0]
        k0_ok = np.ones(len(a), dtype=bool)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                k0_ok &= 0 <= n * v[:, i] - v[:, j]   # (n+1)(0-0) <= n v_i - v_j
        agree = bool((k0_ok[bs]).all())
        res[f"bounded_speed_n{n}"] = dict(
            count=int(bs.sum()),
            all_feasible=bool(a["feasible"][bs].all()),
            k0_valid_iff_condition=bool((k0_ok == bs).all()))
    # overlap improvement (doc claim (b) / next-step #2)
    for n in N_LIST:
        a = dumps[n]
        v = a["v"][:, :n].astype(np.int64)
        num = a["num"].astype(np.int64)
        den = a["den"].astype(np.int64)
        vn = v[:, n - 1]
        ov_num, ov_den = vn, 2 * (n * vn - (n - 1))
        viol = (num * ov_den < den * ov_num).sum()
        plain_viol = (num * 2 * n < den).sum()
        cert = ((ov_num * (n + 1) >= ov_den)).sum()
        res[f"overlap_n{n}"] = dict(
            gap_ge_overlap_violations=int(viol),
            gap_ge_1over2n_violations=int(plain_viol),
            overlap_certifies_conjecture_count=int(cert))
    # doc: (1,3,5,7) "gap >= 1/4"
    res["gap_1357"] = str(gap_ref((1, 3, 5, 7)))
    # doc: (1,3,5) t=1/4 k=(0,0,1) equality at (1,3)?
    res["tight_pairs_135"] = [list(p) for p in tight_pairs((1, 3, 5), (0, 0, 1))]
    return res


# ---------------------------------------------------------------- L6
def l6_structural(dumps):
    res = {}
    rng = random.Random(SEED)
    for n in N_LIST:
        a = dumps[n]
        v = a["v"][:, :n].astype(np.int64)
        # scaling invariance gap(dv) == gap(v)
        gaps = {}
        for r, vv in zip(a, v):
            gaps[tuple(vv)] = (int(r["num"]), int(r["den"]))
        bad = 0
        checked = 0
        for vv, (num, den) in gaps.items():
            for d in range(2, V // vv[-1] + 1):
                scaled = tuple(d * x for x in vv)
                if scaled in gaps:
                    checked += 1
                    if gaps[scaled][0] * den != gaps[scaled][1] * num:
                        bad += 1
        # LP explicit point k_i = (v_i - v_1)/((n+1) v_1) -- sampled exact check
        lp_ok = 0
        sample_idx = rng.sample(range(len(a)), min(500, len(a)))
        for idx in sample_idx:
            vv = tuple(int(x) for x in v[idx])
            if n == 1:
                klp = [Fraction(0)]
            else:
                klp = [Fraction(vv[i] - vv[0], (n + 1) * vv[0]) for i in range(n)]
            good = all(0 <= klp[i] <= vv[i] - 1 for i in range(n))
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    lhs = Fraction((n + 1)) * (vv[j] * klp[i] - vv[i] * klp[j])
                    if lhs > n * vv[i] - vv[j]:
                        good = False
            lp_ok += good
        res[f"struct_n{n}"] = dict(
            scaling_pairs_checked=checked, scaling_violations=bad,
            lp_point_samples=len(sample_idx), lp_point_valid=lp_ok)
    return res


# ---------------------------------------------------------------- L2/L4 samples
def l2_reference_samples(dumps):
    rng = random.Random(SEED)
    res = {}
    per_n = {1: 50, 2: 1225, 3: 2000, 4: 600, 5: 600}
    for n in N_LIST:
        a = dumps[n]
        recs = {tuple(int(x) for x in r["v"][:n]): r for r in a}
        pool = list(recs.keys())
        idx = rng.sample(range(len(pool)), min(per_n[n], len(pool)))
        bad = 0
        tight_bad = 0
        for i in idx:
            vv = pool[i]
            g = gap_ref(vv)
            r = recs[vv]
            gd = Fraction(int(r["num"]), int(r["den"]))
            if g != gd:
                bad += 1
            if (g >= Fraction(1, n + 1)) != bool(r["feasible"]):
                tight_bad += 1
        res[f"ref_n{n}"] = dict(samples=len(idx), gap_mismatches=bad,
                                feasibility_mismatches=tight_bad)
    return res


def l4_milp_samples(dumps):
    rng = random.Random(SEED)
    res = {}
    per_n = {1: 30, 2: 100, 3: 100, 4: 100, 5: 100}
    extra = [(1, 2), (1, 3), (2, 3), (1, 2, 3), (1, 3, 4), (1, 3, 5),
             (1, 2, 3, 4), (1, 3, 4, 7), (1, 2, 3, 4, 5), (1, 3, 4, 5, 9),
             (1, 2, 5, 10), (1, 3, 7, 15), (1, 4, 9, 20), (1, 5, 9, 15)]
    for n in N_LIST:
        a = dumps[n]
        recs = {tuple(int(x) for x in r["v"][:n]): r for r in a}
        pool = [k for k in recs if len(k) == n]
        idx = rng.sample(range(len(pool)), min(per_n[n], len(pool)))
        cases = [pool[i] for i in idx] + [e for e in extra if len(e) == n]
        bad = 0
        for vv in cases:
            k = milp_feasible(vv)
            feas = k is not None
            if feas != bool(recs[vv]["feasible"]):
                bad += 1
        res[f"milp_n{n}"] = dict(solved=len(cases), mismatches=bad)
    return res


# ---------------------------------------------------------------- stats
def stats(dumps):
    res = {}
    for n in N_LIST:
        a = dumps[n]
        num = a["num"].astype(np.int64)
        den = a["den"].astype(np.int64)
        tightm = num * (n + 1) == den
        # reduced gap distribution
        cnt = Counter()
        for nn, dd in zip(num.tolist(), den.tolist()):
            g = math.gcd(nn, dd)
            cnt[(nn // g, dd // g)] += 1
        dist = sorted(cnt.items(), key=lambda kv: -kv[1])[:12]
        # second-tightest (smallest gap strictly above 1/(n+1))
        above = num * (n + 1) > den
        if above.any():
            an, ad = num[above], den[above]
            best = None
            for nn, dd, r in zip(an.tolist(), ad.tolist(), a[above]):
                val = Fraction(nn, dd)
                if best is None or val < best[0]:
                    best = (val, [tuple(int(x) for x in r["v"][:n])])
                elif val == best[0] and len(best[1]) < 8:
                    best[1].append(tuple(int(x) for x in r["v"][:n]))
            second = dict(value=str(best[0]), examples=[list(x) for x in best[1]])
        else:
            second = None
        res[f"stats_n{n}"] = dict(
            total=int(len(a)),
            feasible=int(a["feasible"].sum()),
            tight=int(tightm.sum()),
            tight_fraction=str(Fraction(1, n + 1)),
            gap_distribution=[[f"{k[0]}/{k[1]}", v] for k, v in dist],
            second_tightest=second)
    return res


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


N_LIST = [1, 2, 3, 4, 5]


def main():
    dumps = {n: load(n) for n in N_LIST}
    report = {}

    print("L1: full numpy re-verification of every witness ...")
    l1 = {}
    for n in N_LIST:
        bad, feas_ok = l1_full_recheck(n, dumps[n])
        l1[f"n{n}"] = dict(witness_violations=bad, feasible_flag_consistent=feas_ok)
        print(f"  n={n}: violations={bad} flag_ok={feas_ok}")
    report["L1_full_witness_recheck"] = l1

    print("L2: independent Fraction reference on samples ...")
    report["L2_reference"] = l2_reference_samples(dumps)
    print(" ", json.dumps(report["L2_reference"]))

    print("L3: literal brute-force k-box enumeration (v_n<=12) ...")
    report["L3_brute"] = [l3_brute_crosscheck(n) for n in [2, 3, 4, 5]]
    for r in report["L3_brute"]:
        print(f"  n={r['n']}: vectors={r['vectors']} agree={r['agree']} "
              f"solutions={r['total_solutions']}")

    print("L4: HiGHS MILP on samples ...")
    report["L4_milp"] = l4_milp_samples(dumps)
    print(" ", json.dumps(report["L4_milp"]))

    print("L5: source-audit document claim checks ...")
    report["L5_doc"] = l5_doc_checks(dumps)
    report["L5_claims"] = l5_claim_checks(dumps)

    print("L6: structural checks ...")
    report["L6_structural"] = l6_structural(dumps)

    print("Stats ...")
    report["stats"] = stats(dumps)

    report["checksums"] = {
        f"n{n}_dump": sha256(f"{OUT}/n{n}_V50_dump.bin") for n in N_LIST}
    report["seed"] = SEED

    with open(f"{OUT}/stats.json", "w") as f:
        json.dump(report, f, indent=1, default=str)
    print("\nWrote", f"{OUT}/stats.json")
    # headline
    tot = sum(report["stats"][f"stats_n{n}"]["total"] for n in N_LIST)
    feas = sum(report["stats"][f"stats_n{n}"]["feasible"] for n in N_LIST)
    print(f"HEADLINE: total={tot} feasible={feas} "
          f"counterexamples={tot - feas}")


if __name__ == "__main__":
    main()
