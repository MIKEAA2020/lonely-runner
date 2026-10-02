"""Fast exact gap engine + argmax residue configs for the two-diagnostics run.

Pair-sum candidate theorem (proved in the audit): for a set V of distinct
positive speeds, the argmax of f(t) = min_i ||v_i t|| on [0,1) is attained
at some t = j/(u+w) with u,w in V.  (Peaks odd/(2v) and diff-crossings
j/(w-u) are never argmaxes except when the time coincides with a pair-sum
time; validated against lrc_gap_lib.gap_int, which uses the full candidate
set of Lemma 6.1.)

GapEngine(M) precomputes the dist tensor T[j, D, v] = min(v*j mod D,
D - v*j mod D) for all pair sums D <= 2M-1, so per-set evaluation is a
pure gather + reductions.

Deficit bookkeeping (reviewer's frame): at an argmax t* = a/b (reduced)
with m = b * gap(V), the deficit is d = m*(n+1) - b; the deficit law
(d >= 0) is LRC in residue form.
"""
import sys
from fractions import Fraction
from math import gcd
from itertools import combinations

import numpy as np

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap_int  # exact reference implementation


class GapEngine:
    def __init__(self, M):
        self.M = int(M)
        Dmax = 2 * self.M + 1
        # T[j, D, v] for 0 <= j < D <= Dmax, 1 <= v <= M
        self.T = np.full((Dmax + 1, Dmax + 1, self.M + 1), -1, dtype=np.int16)
        for D in range(2, Dmax + 1):
            j = np.arange(D, dtype=np.int64)
            v = np.arange(1, self.M + 1, dtype=np.int64)
            R = (j[:, None] * v[None, :]) % D
            self.T[0:D, D, 1:] = np.minimum(R, D - R)

    def gap(self, speeds, with_argmax=False):
        """Exact gap (Fraction). If with_argmax, also (m, D, j) of one
        maximizing pair-sum candidate: t = j/D, m = min residue dist."""
        V = sorted(set(int(x) for x in speeds))
        n = len(V)
        if n == 0:
            raise ValueError("empty speed set")
        if n == 1:
            # argmax value 1/2 at t = odd/(2v)
            if with_argmax:
                return Fraction(1, 2), (1, 2 * V[0], 1)
            return Fraction(1, 2)
        Va = np.array(V, dtype=np.int64)
        iu, iw = np.triu_indices(n, k=1)
        D = Va[iu] + Va[iw]
        csum = np.cumsum(D)
        total = int(csum[-1])
        starts = np.concatenate(([0], csum[:-1]))
        pair_idx = np.repeat(np.arange(len(D), dtype=np.int64), D)
        jflat = np.arange(total, dtype=np.int64) - starts[pair_idx]
        Df = D[pair_idx]
        # gather precomputed distances: shape (total, n)
        dist = self.T[jflat[:, None], Df[:, None], Va[None, :]]
        mvec = dist.min(axis=1)
        vals = mvec / Df  # float shortlist, then exact resolution
        j0 = int(np.argmax(vals))
        cand = np.nonzero(vals >= vals[j0] - 1e-9)[0]
        best = None
        for jj in cand:
            m, d = int(mvec[jj]), int(Df[jj])
            if best is None or m * best[1] > best[0] * d:
                best = (m, d, int(jflat[jj]))
        gapv = Fraction(best[0], best[1])
        if with_argmax:
            return gapv, best
        return gapv


def reduce_candidate(m, D, j):
    """Reduce t = j/D to lowest terms; returns (a, b, m_reduced) with the
    invariant g = gcd(j, D) dividing m (proved: g | j and g | D force
    g | (v*j mod D) for every speed v)."""
    g = gcd(j, D)
    a, b = j // g, D // g
    assert m % g == 0, (m, D, j)
    return a, b, m // g


def argmax_configs(speeds):
    """All argmax residue configs via the reference implementation.
    Returns list of dicts: t, a, b, m, residues, binders, d, cover_ok."""
    V = sorted(set(int(x) for x in speeds))
    n = len(V)
    g, ts = gap_int(V)
    out = []
    for t in ts:
        a, b = t.numerator, t.denominator
        res = {v: (v * a) % b for v in V}
        dists = {v: min(res[v], b - res[v]) for v in V}
        m = min(dists.values())
        assert Fraction(m, b) == g, (V, t, m, b, g)
        binders = [v for v in V if dists[v] == m]
        d = m * (n + 1) - b
        cover_ok = b <= n * (2 * m + 1) + 1
        # pole law: residues must attain both m and b - m (generic 2-binder)
        poles = (m in dists.values()) and ((b - m) in dists.values())
        out.append(dict(t=t, a=a, b=b, m=m, residues=res, dists=dists,
                        binders=binders, d=d, cover_ok=cover_ok,
                        pole_law=poles, gap=g, n=n))
    return out


def _random_sets(rng, count, nmax=8, Mmax=40):
    for _ in range(count):
        n = int(rng.integers(2, nmax + 1))
        M = int(rng.integers(n, Mmax + 1))  # ensure M >= n (distinct values exist)
        while True:
            S = set(int(x) for x in rng.integers(1, M + 1, size=n))
            if len(S) == n:
                break
        yield sorted(S)


def validate(count=2500, seed=0, verbose=True):
    """Cross-validate GapEngine against lrc_gap_lib.gap_int; also verify
    the pair-sum property of every reference argmax time."""
    rng = np.random.default_rng(seed)
    eng = GapEngine(60)
    checked = 0
    for V in _random_sets(rng, count):
        g_fast = eng.gap(V)
        g_ref, ts_ref = gap_int(V)
        assert g_fast == g_ref, ("GAP MISMATCH", V, g_fast, g_ref)
        for t in ts_ref:
            ok = any((u + w) % t.denominator == 0
                     for u, w in combinations(V, 2))
            assert ok, ("PAIR-SUM PROPERTY FAILS", V, t, ts_ref)
        # argmax extraction agrees with reference on value and reduction
        g2, (m, D, j) = eng.gap(V, with_argmax=True)
        assert g2 == g_ref
        a, b, mr = reduce_candidate(m, D, j)
        assert Fraction(mr, b) == g_ref, (V, m, D, j)
        checked += 1
    # known values (committed census data cross-checks)
    known = [
        ([1], Fraction(1, 2)),
        ([1, 2], Fraction(1, 3)),
        ([1, 3], Fraction(1, 2)),
        ([1, 2, 3], Fraction(1, 4)),
        (list(range(1, 9)), Fraction(1, 9)),
        ([1, 2, 3, 5], Fraction(1, 4)),      # hand-worked in the audit
        ([1, 2, 3, 4, 5, 7, 18], Fraction(3, 23)),  # aug=7 census vector
        ([1, 5, 6, 11, 16, 17], Fraction(5, 33)),   # N=6 exception (Task 7)
    ]
    for V, expect in known:
        got = eng.gap(V)
        assert got == expect, ("KNOWN-VALUE FAIL", V, got, expect)
        gref, _ = gap_int(V)
        assert gref == expect, ("REFERENCE FAIL", V, gref, expect)
    if verbose:
        print(f"validate: {checked} random sets OK; "
              f"{len(known)} known census values OK; "
              "pair-sum candidate theorem holds on all argmaxes")
    return True


if __name__ == "__main__":
    validate()
