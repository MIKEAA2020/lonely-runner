"""Test 2 (exhaustive): deficit-law scan over ALL n-subsets of [1..M].

For every set V: exact gap via the pair-sum engine; deficit d = m*(n+1)-b
at the (reduced) argmax; pole-law / covering-wall / binder-count stats;
critical census (gap == 1/(n+1)); deficit histogram; violation hunt
(gap < 1/(n+1) ⟺ d < 0 — any hit is either a bug or a counterexample and
is re-verified with the exact reference implementation).

Usage: python3 lrc_diag_scan.py n M out.json
"""
import sys
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd
import time

import numpy as np

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_diag_gap import GapEngine, reduce_candidate, gap_int


def scan(n, M):
    eng = GapEngine(M)
    bound = Fraction(1, n + 1)
    total = 0
    d_hist = Counter()
    near_gaps = Counter()          # gap values in the near-critical band
    violations = []
    critical_norm = Counter()      # gcd-normalized sets at gap == bound
    near_crit_sets = []            # d <= 1 specimens (normalized, capped)
    pole_fail = 0
    cover_fail = 0
    binder_counts = Counter()
    min_gap = None
    min_sets = []
    t0 = time.time()
    for combo in combinations(range(1, M + 1), n):
        V = list(combo)
        g, (m, D, j) = eng.gap(V, with_argmax=True)
        a, b, mr = reduce_candidate(m, D, j)
        d = mr * (n + 1) - b
        total += 1
        d_hist[d] += 1
        if min_gap is None or g < min_gap:
            min_gap, min_sets = g, [tuple(V)]
        elif g == min_gap:
            if len(min_sets) < 400:
                min_sets.append(tuple(V))
        if g < bound or d < 0:
            # re-verify with the independent reference before reporting
            gref, _ = gap_int(V)
            violations.append((V, str(g), str(gref), d, (a, b, mr)))
        if g == bound:
            gv = gcd(*V) if n > 1 else V[0]
            critical_norm[tuple(x // gv for x in V)] += 1
            if len(critical_norm) <= 400:
                pass
        if d <= 1 and len(near_crit_sets) < 300:
            gv = gcd(*V)
            near_crit_sets.append((tuple(x // gv for x in V), str(g), d))
        if g <= Fraction(5, 4 * (n + 1)):
            near_gaps[str(g)] += 1
        # pole law + covering wall + binder census at this argmax
        res = np.array([(v * a) % b for v in V], dtype=np.int64)
        dist = np.minimum(res, b - res)
        mmin = int(dist.min())
        if mmin != mr:
            pole_fail += 1  # sanity: reduced m must equal min dist
        if not (mr in dist.tolist() or (b - mr) in dist.tolist()):
            # residues must attain a pole (mr or b-mr); mr==b-mr counts once
            if 2 * mr != b:
                pole_fail += 1
        if b > n * (2 * mr + 1) + 1:
            cover_fail += 1
        binder_counts[int((dist == mmin).sum())] += 1
    stats = dict(
        n=n, M=M, total_sets=total,
        min_gap=str(min_gap),
        min_gap_sets=len(min_sets),
        min_gap_examples=[list(s) for s in min_sets[:20]],
        bound=str(bound),
        violations=violations,
        violation_count=len(violations),
        critical_normalized_count=len(critical_norm),
        critical_normalized_examples=[list(k) for k in list(critical_norm)[:40]],
        critical_total_instances=sum(critical_norm.values()),
        d_histogram={str(k): v for k, v in sorted(d_hist.items())},
        d_min=min(d_hist) if d_hist else None,
        near_critical_examples=[[list(s), g, int(dd)] for s, g, dd in near_crit_sets[:40]],
        near_gap_band={k: v for k, v in sorted(near_gaps.items(), key=lambda kv: Fraction(kv[0]))[:25]},
        pole_law_failures=pole_fail,
        covering_wall_failures=cover_fail,
        binder_count_hist={str(k): v for k, v in sorted(binder_counts.items())},
        seconds=round(time.time() - t0, 1),
    )
    return stats


if __name__ == "__main__":
    n = int(sys.argv[1])
    M = int(sys.argv[2])
    outp = sys.argv[3]
    stats = scan(n, M)
    with open(outp, "w") as f:
        json.dump(stats, f, indent=1)
    print(f"n={n} M={M}: {stats['total_sets']} sets in {stats['seconds']}s; "
          f"min_gap={stats['min_gap']} (bound {stats['bound']}); "
          f"violations={stats['violation_count']}; "
          f"d_min={stats['d_min']}; "
          f"critical(normalized)={stats['critical_normalized_count']}; "
          f"pole_fail={stats['pole_law_failures']} "
          f"cover_fail={stats['covering_wall_failures']}")
