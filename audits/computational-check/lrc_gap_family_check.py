#!/usr/bin/env python3
"""
lrc_gap_family_check.py

Machine verification of THEOREM 6:  gap(1, 2, ..., n-1, 2n) = 2/(2n+1)  for
every n >= 2, together with the two ingredients of its proof.

Notation (matches audits/ladder analysis.txt): ||x|| = distance to the nearest
integer, gap(v) = max_t min_i ||t v_i||, S_n = {1, ..., n-1, 2n} (n speeds),
delta = 2/(2n+1).

Proof ingredients being machine-checked here:

  Step 1 (argmax localization).  g(t) = min_{v in S_n} ||vt|| is continuous
  piecewise linear; its maximum is attained at a sawtooth peak t=odd/(2v), at
  a crossing t=j/(v_i+v_j) or t=j/(v_j-v_i), or at t=0.  Hence evaluating g
  on the full candidate set computes the exact gap.

  Step 2 (residue kill lemma).  For every reduced denominator m' in
  [2, 3n-1] and every j' coprime to m':
      m' <= n-1            ->  some v in S_n has  v j' = 0   (mod m')
      m' | 2n              ->  v = 2n gives            0
      n+1 <= m' <= 2n-1    ->  some v gives  v j' = +-1 (mod m')
      2n+1 <= m' <= 3n-1   ->  some v gives  v j' = +-1 or +-2 (mod m')
  (min_{v in S_n} ||v j'/m'|| <= 0/0/1/0/2 respectively).

Layers:
  L1  exact family-gap check on the full candidate set (peaks + sum-crossings
      + diff-crossings + t=0), exhaustive n = 2..200, then samples
      {250, 400, 700, 1000} (diff-crossings dropped there: they are
      same-sign crossings and provably never maxima; they are fully covered
      by L2's denominator range anyway).  Certifies: g(t) <= 2/(2n+1) at
      every candidate, and g(2/(2n+1)) = 2/(2n+1) exactly.
  L2  residue kill lemma, exhaustive n = 2..250, samples {400, 800}.
  L3  cross-check against the COMMITTED solver binary lrc_ilp_check.v1.orig
      (untouched build of the committed source): family representatives
      (1,2,...,n-1,2n), n = 2..7, appear in its sweeps with gap exactly
      2/(2n+1); and per-vector gap agreement with this script's independent
      Python implementation on ALL vectors of run 3 6, run 4 8, run 5 10.
  L4  pure-Python (no numpy) exact reference for the family, n = 2..40.
  L5  dense-grid sanity scan (2,000,000 points), n = 2..12, with the
      2n-Lipschitz correction.

All arithmetic is exact (integers).  Any failure raises AssertionError.
Results appended as JSON lines to out/thm6_family_check.jsonl.

Usage:  python3 lrc_gap_family_check.py L1 L2 ...   (default: all layers)
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

SCR = Path("/home/z/my-project/scripts")
OUT = SCR / "out"
OUT.mkdir(exist_ok=True)
JSONL = OUT / "thm6_family_check.jsonl"


def speeds(n):
    """S_n = {1, ..., n-1, 2n} as a plain list."""
    return list(range(1, n)) + [2 * n]


# ---------------------------------------------------------------- candidates
def candidate_arrays(n, peaks=True, sums=True, diffs=True):
    """Candidate points t = A/B (int64 arrays, unreduced fractions):
    sawtooth peaks odd/(2v); crossings j/(v_i+v_j), j/(v_j-v_i); and t=0."""
    S = speeds(n)
    parts_a, parts_b = [], []
    if peaks:
        for v in S:
            parts_a.append(np.arange(1, 2 * v, 2, dtype=np.int64))
            parts_b.append(np.full(v, 2 * v, dtype=np.int64))
    V = np.array(S, dtype=np.int64)
    i, j = np.triu_indices(len(S), k=1)
    if sums:
        for m in np.unique(V[i] + V[j]):
            m = int(m)
            parts_a.append(np.arange(0, m + 1, dtype=np.int64))
            parts_b.append(np.full(m + 1, m, dtype=np.int64))
    if diffs:
        for d in np.unique(V[j] - V[i]):
            d = int(d)
            parts_a.append(np.arange(0, d + 1, dtype=np.int64))
            parts_b.append(np.full(d + 1, d, dtype=np.int64))
    parts_a.append(np.zeros(1, dtype=np.int64))
    parts_b.append(np.ones(1, dtype=np.int64))
    return np.concatenate(parts_a), np.concatenate(parts_b)


def eval_candidates(n, A, B, chunk_entries=2_000_000):
    """M[c] = min_{v in S_n} dist(A[c]*v mod B[c], 0); also count violations
    of M*(2n+1) > 2*B (i.e. g(t) > delta at a candidate point)."""
    V = np.array(speeds(n), dtype=np.int64)
    q = 2 * n + 1
    M = np.empty(len(A), dtype=np.int64)
    bad = 0
    chunk = max(2000, chunk_entries // max(len(V), 1))
    for s in range(0, len(A), chunk):
        a = A[s:s + chunk]
        b = B[s:s + chunk]
        R = (a[:, None] * V[None, :]) % b[:, None]
        D = np.minimum(R, b[:, None] - R)
        m = D.min(axis=1)
        M[s:s + chunk] = m
        viol = m * q > 2 * b
        if viol.any():
            idx = np.nonzero(viol)[0]
            bad += len(idx)
            for c in idx[:5]:
                print(f"  VIOLATION n={n} t={a[c]}/{b[c]} "
                      f"value={m[c]}/{b[c]} > 2/{q}", flush=True)
    return M, bad


def witness_value(n):
    """min_{v in S_n} ||v * 2/(2n+1)|| as an exact integer numerator
    (denominator 2n+1).  Theorem 6 lower bound: must equal 2."""
    V = np.array(speeds(n), dtype=np.int64)
    q = 2 * n + 1
    R = (2 * V) % q
    D = np.minimum(R, q - R)
    return int(D.min())


# ---------------------------------------------------------------- L1
def layer_L1():
    t0 = time.time()
    exhaustive = range(2, 201)
    samples = [250, 400, 700, 1000]
    stats = {}
    maximizer_report = {}
    for n in exhaustive:
        A, B = candidate_arrays(n, peaks=True, sums=True, diffs=True)
        M, bad = eval_candidates(n, A, B)
        assert bad == 0, f"L1: candidate above delta at n={n}"
        w = witness_value(n)
        assert w == 2, f"L1: witness numerator {w} != 2 at n={n}"
        q = 2 * n + 1
        eq = M * q == 2 * B          # candidates attaining delta exactly
        if n <= 40:
            # report the maximizer set (soft data; assert +-2/(2n+1) present)
            tset = sorted({(int(a), int(b)) for a, b in zip(A[eq], B[eq])})
            maximizer_report[n] = [
                f"{a}/{b}" for a, b in tset if 0 < a < b]
            assert (2, q) in tset or (q - 2, q) in tset, (n, tset)
        stats[n] = len(A)
    sample_stats = {}
    for n in samples:
        A, B = candidate_arrays(n, peaks=True, sums=True, diffs=False)
        M, bad = eval_candidates(n, A, B)
        assert bad == 0, f"L1(sample): candidate above delta at n={n}"
        assert witness_value(n) == 2, f"L1(sample): witness at n={n}"
        sample_stats[n] = len(A)
    res = {"layer": "L1", "status": "PASS",
           "exhaustive_n": "2..200", "candidates_max": max(stats.values()),
           "samples": {str(k): v for k, v in sample_stats.items()},
           "maximizers_n<=40": {str(k): v for k, v in maximizer_report.items()},
           "seconds": round(time.time() - t0, 1)}
    print("L1 PASS:", json.dumps({k: v for k, v in res.items()
                                  if k != "maximizers_n<=40"}))
    print("  maximizers (n<=10):",
          {k: v for k, v in list(maximizer_report.items())[:9]})
    return res


# ---------------------------------------------------------------- L2
def layer_L2():
    t0 = time.time()
    checked = 0
    units_total = 0
    for n in list(range(2, 251)) + [400, 800]:
        V = np.array(speeds(n), dtype=np.int64)
        for m in range(2, 3 * n):            # m' in [2, 3n-1]
            if m <= n - 1:
                bound = 0
            elif (2 * n) % m == 0:
                bound = 0                     # m' | 2n (includes m'=n, 2n)
            elif m <= 2 * n - 1:
                bound = 1                     # n+1 <= m' <= 2n-1
            else:
                bound = 2                     # 2n+1 <= m' <= 3n-1
            J = np.arange(1, m, dtype=np.int64)
            U = J[np.gcd(J, m) == 1]
            if len(U) == 0:
                continue
            R = (V[None, :] * U[:, None]) % m
            D = np.minimum(R, m - R)
            M = D.min(axis=1)
            assert (M <= bound).all(), \
                f"L2: kill lemma fails n={n} m'={m} maxM={M.max()} bound={bound}"
            checked += 1
            units_total += len(U)
    res = {"layer": "L2", "status": "PASS",
           "exhaustive_n": "2..250 + samples 400,800",
           "denominators_checked": checked, "units_checked": units_total,
           "seconds": round(time.time() - t0, 1)}
    print("L2 PASS:", json.dumps(res))
    return res


# ---------------------------------------------------------------- L3
REC_DTYPE = np.dtype([
    ("v", "u1", 8), ("k", "u1", 8),
    ("ta", "<u2"), ("tb", "<u2"), ("num", "<u2"), ("den", "<u2"),
    ("feas", "u1"), ("pad", "u1", 7)])


def gap_exact_pure(vlist):
    """Pure-Python exact gap via the candidate theorem (no numpy).
    Returns (num, den), the maximum of min_i ||a v_i / b|| over candidates,
    unreduced."""
    n = len(vlist)
    best_n, best_d = 0, 1
    # peaks
    for v in vlist:
        b = 2 * v
        for a in range(1, b, 2):
            m = min(min((a * x) % b, b - (a * x) % b) for x in vlist)
            if m * best_d > best_n * b:
                best_n, best_d = m, b
    # crossings (sums and diffs)
    for i in range(n):
        for j in range(i + 1, n):
            for b in (vlist[i] + vlist[j], vlist[j] - vlist[i]):
                for a in range(0, b + 1):
                    m = min(min((a * x) % b, b - (a * x) % b) for x in vlist)
                    if m * best_d > best_n * b:
                        best_n, best_d = m, b
    return best_n, best_d


def layer_L3():
    t0 = time.time()
    binary = str(SCR / "lrc_ilp_check.v1.orig")
    family_hits = {}
    agree = 0
    for k in range(2, 8):
        prefix = str(OUT / f"thm6_k{k}")
        r = subprocess.run([binary, "run", str(k), str(2 * k), prefix],
                           check=True, capture_output=True, text=True)
        recs = np.fromfile(f"{prefix}_dump.bin", dtype=REC_DTYPE)
        target = speeds(k)
        sel = np.all(recs["v"][:, :k] == np.array(target, dtype=np.uint8),
                     axis=1)
        assert sel.sum() == 1, f"L3: family vector not found n={k}"
        rec = recs[sel][0]
        num, den = int(rec["num"]), int(rec["den"])
        assert num * (2 * k + 1) == 2 * den, \
            f"L3: committed solver disagrees at n={k}: {num}/{den}"
        assert int(rec["feas"]) == 1
        family_hits[k] = f"{num}/{den}"
    # per-vector agreement with the committed solver on complete small sweeps
    for k, V in [(3, 6), (4, 8), (5, 10)]:
        prefix = str(OUT / f"thm6_agree_{k}_{V}")
        subprocess.run([binary, "run", str(k), str(V), prefix],
                       check=True, capture_output=True, text=True)
        recs = np.fromfile(f"{prefix}_dump.bin", dtype=REC_DTYPE)
        for rec in recs:
            vlist = [int(x) for x in rec["v"][:k]]
            bn, bd = gap_exact_pure(vlist)
            assert bn * int(rec["den"]) == int(rec["num"]) * bd, \
                f"L3: gap disagreement v={vlist}: mine {bn}/{bd} " \
                f"vs solver {int(rec['num'])}/{int(rec['den'])}"
            agree += 1
    res = {"layer": "L3", "status": "PASS",
           "committed_solver_family_gaps": {str(k): v for k, v
                                            in family_hits.items()},
           "per_vector_agreement_sweeps": "run 3 6, run 4 8, run 5 10",
           "per_vector_agreements": agree,
           "seconds": round(time.time() - t0, 1)}
    print("L3 PASS:", json.dumps(res))
    return res


# ---------------------------------------------------------------- L4
def layer_L4():
    t0 = time.time()
    for n in range(2, 41):
        bn, bd = gap_exact_pure(speeds(n))
        assert bn * (2 * n + 1) == 2 * bd, \
            f"L4: n={n} pure reference gap {bn}/{bd} != 2/{2*n+1}"
    res = {"layer": "L4", "status": "PASS", "exhaustive_n": "2..40",
           "seconds": round(time.time() - t0, 1)}
    print("L4 PASS:", json.dumps(res))
    return res


# ---------------------------------------------------------------- L5
def layer_L5():
    t0 = time.time()
    N = 2_000_000
    for n in range(2, 13):
        V = np.array(speeds(n), dtype=np.float64)
        q = 2 * n + 1
        best = 0.0
        for s in range(0, N, 400_000):
            T = np.arange(s, min(s + 400_000, N)) / N
            R = (T[:, None] * V[None, :]) % 1.0
            D = np.minimum(R, 1.0 - R)
            best = max(best, float(D.min(axis=1).max()))
        # g is 2n-Lipschitz: true max <= grid max + 2n/N; also grid <= true.
        assert best <= 2.0 / q + 1e-9, (n, best)
        assert best >= 2.0 / q - 2.0 * n / N - 1e-9, (n, best)
    res = {"layer": "L5", "status": "PASS", "exhaustive_n": "2..12",
           "grid_points": N, "seconds": round(time.time() - t0, 1)}
    print("L5 PASS:", json.dumps(res))
    return res


LAYERS = {"L1": layer_L1, "L2": layer_L2, "L3": layer_L3,
          "L4": layer_L4, "L5": layer_L5}


def main():
    which = sys.argv[1:] or list(LAYERS)
    results = []
    for name in which:
        print(f"=== {name} ===", flush=True)
        results.append(LAYERS[name]())
    with open(JSONL, "a") as f:
        for r in results:
            f.write(json.dumps(r) + "\n")
    print("ALL REQUESTED LAYERS PASS:", ",".join(which))


if __name__ == "__main__":
    main()
