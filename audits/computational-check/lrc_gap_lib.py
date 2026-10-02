"""Exact gap machinery for the lonely-runner zoo (Python reference).

Two independent implementations:
  gap_int(v)   -- integer arithmetic, mirrors lrc_ilp_check.c solve()
  gap_frac(v)  -- Fraction arithmetic, independent formulation
Both return (gap, argmaxes) with gap a Fraction and argmaxes the sorted
list of ALL candidate times attaining the gap.

Candidate theorem (repo, Lemma 6.1 / solver header): the max of
g(t) = min_i ||t v_i|| over [0,1] is attained among
  peaks      t = odd/(2 v_i),
  crossings  t(v_i+v_j) in Z  or  t(v_j-v_i) in Z,
plus t in {0,1}; valleys t = m/v_i give g = 0.
"""
from fractions import Fraction
from math import gcd, floor


def _candidates_int(v):
    """Yield (a, b) pairs with t = a/b, b > 0, covering all candidates."""
    out = []
    for x in v:                       # peaks
        for a in range(1, 2 * x, 2):
            out.append((a, 2 * x))
    n = len(v)
    for i in range(n):                # crossings
        for j in range(i + 1, n):
            s = v[i] + v[j]
            for a in range(0, s + 1):
                out.append((a, s))
            d = v[j] - v[i]
            for a in range(0, d + 1):
                out.append((a, d))
    out.append((0, 1))
    out.append((1, 1))
    return out


def _dist_mod(a, b):
    """dist(a, b) = min(a mod b, b - a mod b) for b > 0, a >= 0."""
    r = a % b
    return r if r < b - r else b - r


def gap_int(v):
    """(gap, argmaxes). Integer arithmetic; dedup argmaxes by value."""
    best_n, best_d = 0, 1
    amax = set()
    for a, b in _candidates_int(v):
        m = min(_dist_mod(a * x, b) for x in v)
        lhs, rhs = m * best_d, best_n * b
        if lhs > rhs:
            best_n, best_d = m, b
            amax = {(a, b)}
        elif lhs == rhs:
            amax.add((a, b))
    gap = Fraction(best_n, best_d)
    ts = sorted({Fraction(a, b) for a, b in amax})
    return gap, ts


def _norm(x):
    """||x|| for Fraction x: distance to nearest integer."""
    f = x - floor(x)
    if f > Fraction(1, 2):
        f = 1 - f
    return f


def gap_frac(v):
    """(gap, argmaxes). Fraction arithmetic reference."""
    cands = {Fraction(a, b) for a, b in _candidates_int(v)}
    best = Fraction(0)
    ts = []
    for t in sorted(cands):
        m = min(_norm(x * t) for x in v)
        if m > best:
            best, ts = m, [t]
        elif m == best:
            ts.append(t)
    return best, ts


def gap(v):
    """Consensus gap: both implementations must agree."""
    g1, t1 = gap_int(v)
    g2, t2 = gap_frac(v)
    assert g1 == g2, (v, g1, g2)
    assert t1 == t2, (v, t1, t2)
    return g1, t1


def binders(v, t):
    """Speeds u in v with ||u t|| = gap(v) at argmax t (must be one)."""
    g, _ = gap_int(v)
    # t may be any Fraction; normalize to denominator dividing lcm grid
    return [u for u in v if _norm(u * t) == g]


def dist_at(u, t):
    """||u t|| for Fraction t."""
    return _norm(u * t)


def family(n):
    """The Theorem 6 family {1,...,n-1, 2n}."""
    return list(range(1, n)) + [2 * n]
