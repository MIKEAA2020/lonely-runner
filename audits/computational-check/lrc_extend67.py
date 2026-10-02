#!/usr/bin/env python3
"""
lrc_extend67.py -- Extension of the computational check to n=6 (V=50) and n=7 (V=38).

Reads the dumps produced by the committed solver (lrc_ilp_check.c, unchanged),
re-verifies every witness against all polytope constraints (layer L1), compares
against an independent pure-Python reference solver on samples (L2), runs the
literal k-box brute enumeration for small v_n (L3), cross-solves with HiGHS
MILP (L4), and computes the census statistics: tight sets, second-tightest gap,
gap distribution, rung-2 minimizer structure (modular witnesses at denominator
2n+1), and targeted lookups for the theory-critical families.

Output: out/stats67.json
"""
import itertools
import json
import math
import random
import subprocess
import sys
from fractions import Fraction

import numpy as np

OUT = "/home/z/my-project/scripts/out"
BIN = "/home/z/my-project/scripts/lrc_ilp_check"
SEED = 20261002

# ---- record dtype: 32-byte packed Rec from lrc_ilp_check.c -----------------
DT = np.dtype({
    "names": ["v", "k", "ta", "tb", "num", "den", "feasible", "pad"],
    "formats": ["(8,)u1", "(8,)u1", "<u2", "<u2", "<u2", "<u2", "u1", "(7,)u1"],
    "offsets": [0, 8, 16, 18, 20, 22, 24, 25],
    "itemsize": 32,
})

CHUNK = 2_000_000


def gap_value(num, den):
    g = math.gcd(int(num), int(den))
    return Fraction(int(num) // g, int(den) // g)


def iter_dump(path, n):
    """Yield (chunk_of_records) reading the dump incrementally."""
    with open(path, "rb") as f:
        while True:
            recs = np.fromfile(f, dtype=DT, count=CHUNK)
            if recs.size == 0:
                break
            yield recs


def verify_chunk(n, recs):
    """L1: verify every witness k against ALL constraints of P(v). Exact int64."""
    v = recs["v"][:, :n].astype(np.int64)
    k = recs["k"][:, :n].astype(np.int64)
    bad_box = np.zeros(recs.size, dtype=bool)
    bad_box |= (k < 0).any(axis=1)
    bad_box |= (k > v - 1).any(axis=1)
    bad_pair = np.zeros(recs.size, dtype=bool)
    for i in range(n):
        vi = v[:, i]
        ki = k[:, i]
        for j in range(n):
            if i == j:
                continue
            lhs = (n + 1) * (v[:, j] * ki - vi * k[:, j])
            rhs = n * vi - v[:, j]
            bad_pair |= lhs > rhs
    num = recs["num"].astype(np.int64)
    den = recs["den"].astype(np.int64)
    feas_flag = recs["feasible"].astype(np.int64)
    feas_calc = (num * (n + 1) >= den).astype(np.int64)
    flag_mismatch = feas_flag != feas_calc
    ok = ~(bad_box | bad_pair)
    return ok, flag_mismatch


def analyze(n, V, prefix, interesting_vectors):
    """Full L1 + census scan of one dump."""
    total = 0
    feasible = 0
    n_viol = 0
    n_flagmis = 0
    key_counts = {}          # reduced gap key -> count
    examples = {}            # key -> list of (v tuple, k tuple, t frac str, gap)
    lookups = {}             # encoded v -> record dict for interesting vectors
    want = set()
    for vec in interesting_vectors:
        want.add(encode_v(vec))
    second_run = None

    for recs in iter_dump(f"{OUT}/{prefix}_dump.bin", n):
        ok, flagmis = verify_chunk(n, recs)
        n_viol += int((~ok).sum())
        n_flagmis += int(flagmis.sum())
        total += recs.size
        feasible += int(recs["feasible"].sum())

        num = recs["num"].astype(np.int64)
        den = recs["den"].astype(np.int64)
        g = np.gcd(num, den)
        rk = (num // g) * 256 + (den // g)  # den < 256 in range
        keys, counts = np.unique(rk, return_counts=True)
        for key, cnt in zip(keys.tolist(), counts.tolist()):
            key_counts[key] = key_counts.get(key, 0) + cnt

        # examples for interesting keys: append across chunks (cap 8)
        for key in keys.tolist():
            ex = examples.setdefault(key, [])
            if len(ex) < 8:
                mask = rk == key
                idx = np.nonzero(mask)[0][: 8 - len(ex)]
                for r in idx:
                    ex.append(record_info(n, recs[r]))

        # targeted lookups
        enc = encode_rows(recs, n)
        m = np.isin(enc, list(want))
        for r in np.nonzero(m)[0]:
            e = int(enc[r])
            if e not in lookups:
                lookups[e] = record_info(n, recs[r])

    # gap values as Fractions
    def key_to_frac(key):
        return Fraction(key // 256, key % 256)

    delta = Fraction(1, n + 1)
    values = sorted((key_to_frac(k), c) for k, c in key_counts.items())
    # second-tightest: min value strictly above delta with count>0
    above = [(val, c) for val, c in values if val > delta]
    second_val = above[0][0] if above else None
    second_count = above[0][1] if above else 0
    second_examples = examples.get(
        (second_val.numerator * 256 + second_val.denominator) if second_val else -1, []
    )

    dist_top = [
        (str(v), c) for v, c in sorted(values, key=lambda x: -x[1])[:15]
    ]
    values_str = [(str(v), c) for v, c in values]

    return {
        "n": n, "V": V, "total": total, "feasible": feasible,
        "infeasible": total - feasible,
        "witness_violations": n_viol,
        "feasible_flag_mismatches": n_flagmis,
        "tight_count": sum(c for v, c in values if v == delta),
        "second_tightest": {
            "value": str(second_val),
            "count": second_count,
            "examples": second_examples,
        },
        "gap_distribution_top15": dist_top,
        "gap_values_all": values_str,
        "lookups": {str(tuple(sorted(decode_v(e)))): r
                    for e, r in lookups.items()},
    }


def encode_v(vec):
    e = 0
    for x in vec:
        e = e * 64 + x
    return e


def decode_v(e):
    out = []
    while e:
        out.append(e % 64)
        e //= 64
    return list(reversed(out))


def encode_rows(recs, n):
    v = recs["v"][:, :n].astype(np.int64)
    e = np.zeros(recs.size, dtype=np.int64)
    for i in range(n):
        e = e * 64 + v[:, i]
    return e


def record_info(n, rec):
    return {
        "v": [int(x) for x in rec["v"][:n]],
        "k": [int(x) for x in rec["k"][:n]],
        "t": f"{int(rec['ta'])}/{int(rec['tb'])}",
        "gap": f"{int(rec['num'])}/{int(rec['den'])}",
    }


# ---------------- L2: independent reference solver --------------------------
def ref_gap(n, v):
    """Independent exact reference: same breakpoint theorem, fresh code path."""
    best = Fraction(0, 1)
    best_t = Fraction(0, 1)
    cands = set()
    for vi in v:
        for a in range(1, 2 * vi, 2):
            cands.add(Fraction(a, 2 * vi))
    for i in range(n):
        for j in range(i + 1, n):
            s = v[i] + v[j]
            for a in range(0, s + 1):
                cands.add(Fraction(a, s))
            d = v[j] - v[i]
            for a in range(0, d + 1):
                cands.add(Fraction(a, d))
    for t in cands:
        b = t.denominator
        a = t.numerator
        m = min(Fraction(min(a * vi % b, b - a * vi % b), b) for vi in v)
        if m > best:
            best = m
            best_t = t
    return best, best_t


def l2_reference(n, V, prefix, samples):
    random.seed(SEED + n)
    mism_gap = 0
    mism_feas = 0
    checked = 0
    recs_all = []
    # collect sampled records by reservoir over the dump
    for recs in iter_dump(f"{OUT}/{prefix}_dump.bin", n):
        recs_all.append(recs)
    alld = np.concatenate(recs_all) if recs_all else np.array([])
    idx = random.sample(range(alld.size), min(samples, alld.size))
    for r in idx:
        rec = alld[r]
        v = [int(x) for x in rec["v"][:n]]
        g, t = ref_gap(n, v)
        gd = gap_value(rec["num"], rec["den"])
        feas = int(rec["feasible"])
        if g != gd:
            mism_gap += 1
        if feas != (1 if g >= Fraction(1, n + 1) else 0):
            mism_feas += 1
        checked += 1
    return {"samples": checked, "gap_mismatches": mism_gap,
            "feasibility_mismatches": mism_feas}


# ---------------- L3: literal brute enumeration -----------------------------
def l3_brute(n, vmax):
    vecs = list(itertools.combinations(range(1, vmax + 1), n))
    inp = "\n".join(" ".join(str(x) for x in v) for v in vecs) + "\n"
    p = subprocess.run([BIN, "brute", str(n)], input=inp,
                       capture_output=True, text=True, timeout=560)
    lines = p.stdout.strip().splitlines()
    total_solutions = 0
    with_sol = 0
    agree = 0
    firsts = {}
    for line, v in zip(lines, vecs):
        # v=(...) feasible_integer_points=N first_k=(...)
        try:
            cnt = int(line.split("feasible_integer_points=")[1].split()[0])
            fk = line.split("first_k=(")[1].rstrip(")").strip()
            firsts[v] = (cnt, fk)
        except Exception:
            firsts[v] = (None, None)
        total_solutions += cnt if cnt else 0
        if cnt and cnt > 0:
            with_sol += 1
        if cnt is not None and cnt > 0:
            agree += 1
    return {"vmax": vmax, "vectors": len(vecs),
            "with_solution": with_sol, "agree": agree,
            "total_solutions": total_solutions, "firsts": firsts}


def l3_compare(n, firsts):
    """Brute-side agreement: every brute first_k verified exactly against P(v);
    dump side says every vector is feasible, so agreement == all counts > 0."""
    bad = []
    for v, (cnt, fk) in firsts.items():
        if cnt is None or cnt <= 0:
            bad.append(list(v))
            continue
        k = [int(x) for x in fk.split(",")]
        if not verify_point_py(n, list(v), k):
            bad.append(list(v))
    return {"vectors_checked": len(firsts),
            "empty_or_invalid": len(bad), "examples": bad[:10]}


def verify_point_py(n, v, k):
    for i in range(n):
        if k[i] < 0 or k[i] > v[i] - 1:
            return False
        for j in range(n):
            if i == j:
                continue
            if (n + 1) * (v[j] * k[i] - v[i] * k[j]) > n * v[i] - v[j]:
                return False
    return True


# ---------------- L4: HiGHS MILP cross-solve ---------------------------------
def milp_solve(n, v):
    from scipy.optimize import milp, LinearConstraint, Bounds
    rows, lbs, ubs = [], [], []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            row = np.zeros(n)
            row[i] = (n + 1) * v[j]
            row[j] = -(n + 1) * v[i]
            rows.append(row)
            lbs.append(-np.inf)
            ubs.append(n * v[i] - v[j])
    A = np.array(rows, dtype=float)
    con = LinearConstraint(A, lbs, ubs)
    lb = np.zeros(n)
    ub = np.array([v[i] - 1 for i in range(n)], dtype=float)
    res = milp(c=np.zeros(n), constraints=con, bounds=Bounds(lb, ub),
               integrality=np.ones(n))
    return res


def l4_milp(n, V, prefix, n_random, extra_vectors):
    random.seed(SEED + 100 * n)
    recs_all = []
    for recs in iter_dump(f"{OUT}/{prefix}_dump.bin", n):
        recs_all.append(recs)
    alld = np.concatenate(recs_all)
    idx = random.sample(range(alld.size), min(n_random, alld.size))
    vecs = [tuple(int(x) for x in alld[r]["v"][:n]) for r in idx]
    vecs += [tuple(v) for v in extra_vectors]
    solved = 0
    mism = 0
    bad_points = 0
    for v in vecs:
        res = milp_solve(n, list(v))
        ok = res.success
        if ok:
            k = np.round(res.x).astype(int)
            # verify exactly
            good = all(0 <= k[i] <= v[i] - 1 for i in range(n))
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    if (n + 1) * (v[j] * k[i] - v[i] * k[j]) > n * v[i] - v[j]:
                        good = False
            if not good:
                bad_points += 1
        # dump says feasible (all are, in practice) -> MILP must succeed
        if not ok:
            mism += 1
        solved += 1
    return {"solved": solved, "milp_failures": mism,
            "invalid_returned_points": bad_points}


# ---------------- modular witness machinery (theory) -------------------------
def modular_witness(n, v, k_rung):
    """Exists c in [1, kn] with c*v_i mod (kn+1) in [k, kn+1-k] for all i?
    Returns (True, c) or (False, None)."""
    q = k_rung * n + 1
    for c in range(1, q):
        if all(k_rung <= (c * vi) % q <= q - k_rung for vi in v):
            return True, c
    return False, None


# ============================================================================
def main():
    stats = {}

    interesting6 = [
        (1, 2, 3, 4, 5, 6),        # AP tight
        (1, 3, 4, 5, 6, 11),       # old sporadic family, predicted NOT tight
        (1, 2, 3, 4, 5, 12),       # rung-2 family {1..n-1, 2n}
        (2, 4, 6, 8, 10, 24),      # scaled rung-2 family
        (1, 3, 4, 5, 6, 7),        # census comparison
        (1, 2, 3, 4, 5, 7),
    ]
    interesting7 = [
        (1, 2, 3, 4, 5, 6, 7),        # AP tight
        (1, 2, 3, 4, 5, 7, 12),       # new family A
        (1, 4, 5, 6, 7, 11, 13),      # new family B
        (1, 3, 4, 5, 6, 7, 13),       # old sporadic family, predicted NOT tight
        (1, 2, 3, 4, 5, 6, 14),       # rung-2 family {1..n-1, 2n}
        (2, 4, 6, 8, 10, 12, 28),     # scaled rung-2
    ]

    print("analyzing n=6 ...", flush=True)
    stats["n6"] = analyze(6, 50, "n6_V50", interesting6)
    print("analyzing n=7 ...", flush=True)
    stats["n7"] = analyze(7, 38, "n7_V38", interesting7)

    print("L2 reference n=6 ...", flush=True)
    stats["l2_n6"] = l2_reference(6, 50, "n6_V50", 250)
    print("L2 reference n=7 ...", flush=True)
    stats["l2_n7"] = l2_reference(7, 38, "n7_V38", 250)

    print("L3 brute n=6 (v_n<=12) ...", flush=True)
    b6 = l3_brute(6, 12)
    stats["l3_n6"] = {kk: v for kk, v in b6.items() if kk != "firsts"}
    stats["l3_n6_compare"] = l3_compare(6, b6["firsts"])
    print("L3 brute n=7 (v_n<=10) ...", flush=True)
    b7 = l3_brute(7, 10)
    stats["l3_n7"] = {kk: v for kk, v in b7.items() if kk != "firsts"}
    stats["l3_n7_compare"] = l3_compare(7, b7["firsts"])

    print("L4 MILP n=6 ...", flush=True)
    stats["l4_n6"] = l4_milp(6, 50, "n6_V50", 40, interesting6)
    print("L4 MILP n=7 ...", flush=True)
    stats["l4_n7"] = l4_milp(7, 38, "n7_V38", 40, interesting7)

    # modular witness analysis for the second rung on rung-2 examples
    for key in ("n6", "n7"):
        n = stats[key]["n"]
        st = stats[key]
        examples = st["second_tightest"]["examples"]
        wit = []
        for ex in examples:
            ok, c = modular_witness(n, ex["v"], 2)
            wit.append({"v": ex["v"], "t": ex["t"], "witness_c": c})
        st["second_tightest"]["modular_witness_c"] = wit

    with open(f"{OUT}/stats67.json", "w") as f:
        json.dump(stats, f, indent=1)

    # console summary
    for key in ("n6", "n7"):
        st = stats[key]
        print(f"[{key}] total={st['total']} feasible={st['feasible']} "
              f"infeasible={st['infeasible']} tight={st['tight_count']} "
              f"viol={st['witness_violations']} "
              f"second_tightest={st['second_tightest']['value']} "
              f"(x{st['second_tightest']['count']})")
    print(json.dumps({k: v for k, v in stats.items()
                      if k.startswith(("l2", "l3", "l4"))}, indent=1))
    print("lookups n6:", json.dumps(stats["n6"]["lookups"], indent=1))
    print("lookups n7:", json.dumps(stats["n7"]["lookups"], indent=1))


if __name__ == "__main__":
    main()
