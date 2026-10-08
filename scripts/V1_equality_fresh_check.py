#!/usr/bin/env python3
"""
V1_equality_fresh_check.py — Non-circular verification of the Equality Theorem
(M(V) = max over pair-sum grids), for the reviewer's directive 2.

Independence: M_unrestricted is computed ONLY from the folklore breakpoint set
  D = {2 v_p} U {v_p + v_q} U {|v_p - v_q|}   (peaks + curve intersections),
which is exhaustive for vertices of the PL lower envelope g = min_p ||v_p t||.
M_pairsum is computed ONLY on pair-sum lattices k/(v_u+v_w).

Checks per set V (gcd 1, distinct speeds):
  C1  M_pairsum == M_unrestricted  (exact integer arithmetic)
  C2  some unrestricted ARGMAX time lies on a pair-sum lattice (numerator check)
  C3  all speeds odd  <=>  M == 1/2   (the pole/parity case)
  C4  if M < 1/2: at the sum-grid argmax, binding is two-sided (residues m and
      N-m both present among binders)  [spot-check on a sample]
  C5  count argmaxes that are difference-only (folklore needs v_i - v_j) even
      though a sum argmax also exists  [structural sharpening stat]

Corpus: all 4-sets v<=12; random sample of 5-sets v<=20; targeted edge sets.
"""
from math import gcd
from itertools import combinations
from fractions import Fraction
import random

def dist_int(x, d):
    r = x % d
    return r if r <= d - r else d - r

def candidates(V):
    """Folklore breakpoint denominators for the full n-runner problem."""
    ds = set()
    for v in V:
        ds.add(2 * v)
    for (a, b) in combinations(V, 2):
        ds.add(a + b)
        if a != b:
            ds.add(abs(a - b))
    return sorted(ds)

def g_at(V, num, den):
    return min(dist_int(v * num, den) for v in V)

def M_unrestricted(V):
    """Exact M via folklore breakpoints. Returns (Fraction, argmax list)."""
    best = Fraction(0)
    args = []
    for d in candidates(V):
        for m in range(d):
            val = g_at(V, m, d)
            fr = Fraction(val, d)
            if fr > best:
                best, args = fr, [(m, d)]
            elif fr == best and (m, d) not in args:
                args.append((m, d))
    return best, args

def M_pairsum(V):
    """M via pair-sum lattices only. Returns (Fraction, (pair, k) list)."""
    best = Fraction(0)
    args = []
    for (p, q) in combinations(range(len(V)), 2):
        N = V[p] + V[q]
        for k in range(1, N):
            val = g_at(V, k, N)
            fr = Fraction(val, N)
            if fr > best:
                best, args = fr, [((p, q), k)]
            elif fr == best and ((p, q), k) not in args:
                args.append(((p, q), k))
    return best, args

def on_sum_lattice(m, d, V):
    """Is t=m/d (as given) equal to some k/N with N=v_u+v_w?  (m/d = k/N iff
    N*m/d integer... careful: reduce m/d first.)"""
    from math import gcd as G
    h = G(m, d)
    num, den = m // h, d // h
    for (a, b) in combinations(V, 2):
        N = a + b
        if (N * num) % den == 0:
            return True, (a + b)
    return False, None

def binding_sides(V, num, den):
    """At t=num/den: residues of binders. Two-sided iff some binder residue is
    m and another is den-m (m = min residual)."""
    res = [ (v * num) % den for v in V ]
    dists = [ min(r, den - r) for r in res ]
    m = min(dists)
    if m == 0 or 2 * m == den:
        return ('degenerate', m)
    sides = set(r for r, dd in zip(res, dists) if dd == m)
    # side +1: r == m ; side -1: r == den - m
    has_p = any(r == m for r in sides)
    has_m = any(r == den - m for r in sides)
    return ('two-sided' if (has_p and has_m) else 'one-sided'), m

def run_corpus(sets, label):
    n_sets = len(sets)
    c1_fail, c2_fail, c3_fail, c4_fail = [], [], [], []
    diff_only_count = 0
    two_sided_ok = 0
    for V in sets:
        Mu, args_u = M_unrestricted(V)
        Mp, args_p = M_pairsum(V)
        if Mu != Mp:
            c1_fail.append((V, Mu, Mp))
            continue
        # C2: some unrestricted argmax on a pair-sum lattice
        ok2 = any(on_sum_lattice(m, d, V)[0] for (m, d) in args_u)
        if not ok2:
            c2_fail.append((V, args_u[:4], Mp))
        # C5: difference-only argmaxes while sum argmax exists
        onsum = [a for a in args_u if on_sum_lattice(*a, V)[0]]
        if len(onsum) < len(args_u):
            diff_only_count += 1
        # C3: parity
        all_odd = all(v % 2 == 1 for v in V)
        if all_odd != (Mu == Fraction(1, 2)):
            c3_fail.append((V, Mu))
        # C4: two-sided binding at a sum-grid argmax (only when M < 1/2)
        if Mu < Fraction(1, 2):
            for ((p, q), k) in args_p:
                N = V[p] + V[q]
                side, m = binding_sides(V, k, N)
                if side == 'two-sided':
                    two_sided_ok += 1
                    break
            else:
                c4_fail.append((V, args_p[:2]))
    print(f"[{label}] sets={n_sets}")
    print(f"  C1 equality M_pairsum == M_unrestricted : "
          f"{n_sets - len(c1_fail)}/{n_sets}"
          + (f"  FAILS: {c1_fail[:3]}" if c1_fail else ""))
    print(f"  C2 some argmax on pair-sum lattice      : "
          f"{n_sets - len(c2_fail)}/{n_sets}"
          + (f"  FAILS: {c2_fail[:2]}" if c2_fail else ""))
    print(f"  C3 all-odd <=> M == 1/2                 : "
          f"{n_sets - len(c3_fail)}/{n_sets}"
          + (f"  FAILS: {c3_fail[:3]}" if c3_fail else ""))
    print(f"  C4 two-sided binding at sum-grid argmax : "
          f"{two_sided_ok}/{sum(1 for V in sets if M_unrestricted(V)[0] < Fraction(1,2))}"
          + (f"  FAILS: {c4_fail[:2]}" if c4_fail else ""))
    print(f"  C5 sets with difference-only argmaxes too: {diff_only_count}"
          f"  (folklore's |v_i-v_j| denominators appear but never needed)")
    return not (c1_fail or c2_fail or c3_fail or c4_fail)

def main():
    random.seed(20261003)
    # n=4 corpus, v <= 12
    from math import gcd as _g
    sets4 = [V for V in combinations(range(1, 13), 4) if _g(*V) == 1]
    ok4 = run_corpus(sets4, "n=4, v<=12, gcd 1")
    # n=5 random sample, v <= 20
    pool5 = list(combinations(range(1, 21), 5))
    sample5 = random.sample(pool5, 400)
    sets5 = [V for V in sample5 if _g(*V) == 1]
    ok5 = run_corpus(sets5, "n=5, v<=20, random 400")
    # targeted edge sets
    edge = [(1, 2), (1, 3), (2, 3), (1, 2, 3), (1, 3, 5, 7), (1, 3, 5, 7, 9),
            (1, 2, 3, 4), (1, 3, 4, 5), (1, 3, 4, 7), (2, 6, 8, 10, 11),
            (1, 2, 3, 16), (1, 6, 11, 16), (3, 5, 7, 11, 13)]
    oke = run_corpus(edge, "targeted edge sets (incl. n=2, all-odd, tight)")
    print()
    print("OVERALL:", "ALL CHECKS PASS" if (ok4 and ok5 and oke) else "FAILURES PRESENT")

if __name__ == "__main__":
    main()
