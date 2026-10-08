#!/usr/bin/env python3
"""
E3_recovery_diagnostic.py — Reviewer Step 3: can the classical n=4 theorem be
RECOVERED through the reduction, i.e. does the pair-sum grid capture classical
3-runner content, or does (TAU) rest on grid-specific values?

For every 4-speed set V (gcd 1, v <= 16) and every pair (p,q), N = v_p+v_q,
effective triple eff = (v_p, v_y, v_z):

  - grid max m3 = max_k min_w ||w k||_N   (the (TAU) certificate value)
  - UNRESTRICTED 3-runner optimum mu(eff) = max_t min_w ||w t||, computed
    EXACTLY by breakpoint enumeration (candidates t = m/d, d in
    {2a,2b,2c,a+b,a+c,b+c,|a-b|,|a-c|,|b-c|}) — the max of a PL function
    with no flat pieces is attained at a breakpoint.
  - grid-visible optimum mu_vis: max over UNRESTRICTED-OPTIMAL breakpoints
    that lie on the pair-sum grid (t = m/d = k/N  iff  d | N*m).

Pair classification (the recovery tiers):
  classical : mu_vis == mu   -> the grid captures the classical 3-runner
              optimum (>= 1/4 by the classical k=3 theorem); the (TAU)
              certificate for this pair is classical content + divisibility.
  semi      : certifies (TAU), mu_vis < mu, but the grid argmax k IS an
              unrestricted breakpoint (some 2wk = 0 or (wi+-wj)k = 0 mod N)
              with value in [1/5, 1/4): breakpoint structure known, value
              bound is grid-specific.
  gridonly  : certifies, mu_vis < mu, argmax not a breakpoint: fully
              grid-specific certificate.
  fail      : does not certify (5*m3 < N).
  degen     : effective triple stuck (0 mod N) or colliding mod N (flagged).

Set classification: classical if SOME pair is classical; else semi; else
gridonly; else fail.

Sanity assertions (bug-catchers):
  A1: m3/N <= mu for every pair (grid is a subset of times).
  A2: mu >= 1/4 for every non-colliding, non-stuck triple (classical k=3).
"""
from math import gcd
from itertools import combinations
from collections import Counter

VMAX = 16

def dist_d(x, d):
    r = x % d
    return r if r <= d - r else d - r

def unrestricted_opt(eff):
    """Exact max_t min_w ||w t||. Returns (mu_num, mu_den, [(d, m, val), ...])."""
    a, b, c = eff
    dens = sorted({d for d in (2*a, 2*b, 2*c, a+b, a+c, b+c,
                               abs(a-b), abs(a-c), abs(b-c)) if d > 0})
    best_num, best_den = 0, 1
    args = []
    for d in dens:
        for m in range(d):
            val = min(dist_d(w * m, d) for w in (a, b, c))
            # value = val / d ; compare with best_num/best_den
            if val * best_den > best_num * d:
                best_num, best_den = val, d
                args = [(d, m, val)]
            elif val * best_den == best_num * d:
                args.append((d, m, val))
    return best_num, best_den, args

def grid_max(N, eff):
    best, argk = 0, 0
    half = N // 2
    for k in range(1, half + 1):
        m = N
        for w in eff:
            y = (w * k) % N
            d = y if y <= N - y else N - y
            if d < m:
                m = d
        if m > best:
            best, argk = m, k
    return best, argk

def is_unrestricted_breakpoint(k, N, eff):
    ws = eff
    for w in ws:
        if (2 * w * k) % N == 0:
            return True
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            if ((ws[i] + ws[j]) * k) % N == 0 or (abs(ws[i] - ws[j]) * k) % N == 0:
                return True
    return False

def breakpoint_grid_max(N, eff):
    """max over GRID k that are unrestricted breakpoints, of min residual."""
    best = 0
    for k in range(1, N):
        if is_unrestricted_breakpoint(k, N, eff):
            m = N
            for w in eff:
                y = (w * k) % N
                d = y if y <= N - y else N - y
                if d < m:
                    m = d
            if m > best:
                best = m
    return best

def main():
    total = 0
    set_class = Counter()
    pair_class = Counter()
    nonclassical_examples = []
    mu_violations = []   # A2 violations
    grid_gt_mu = []      # A1 violations
    for V in combinations(range(1, VMAX + 1), 4):
        g = 0
        for x in V:
            g = gcd(g, x)
        if g != 1:
            continue
        total += 1
        per_pair = []
        setbest = None
        for (p, q) in combinations(range(4), 2):
            N = V[p] + V[q]
            others = [V[i] for i in range(4) if i not in (p, q)]
            eff = (V[p], others[0], others[1])
            m3, argk = grid_max(N, eff)
            mu_num, mu_den, args = unrestricted_opt(eff)
            # A1: m3/N <= mu
            if m3 * mu_den > mu_num * N:
                grid_gt_mu.append((V, (p, q), m3, N, mu_num, mu_den))
            # A2: mu >= 1/4 for genuine triples
            stuck = any(w % N == 0 for w in eff)
            collide = ((eff[0] % N) == (eff[1] % N) or (eff[0] % N) == (eff[2] % N)
                       or (eff[1] % N) == (eff[2] % N))
            if not stuck and not collide and mu_num * 4 < mu_den:
                mu_violations.append((V, eff, mu_num, mu_den))
            # grid-visible optimum: is some UNRESTRICTED-OPTIMAL breakpoint a
            # grid point?  t = m/d is on the grid (1/N)Z  iff  d | N*m.
            classical = any((N * m) % d == 0 for (d, m, val) in args)
            cert = 5 * m3 >= N
            bp_max = breakpoint_grid_max(N, eff)
            if stuck or collide:
                cls = 'degen-stuck' if stuck else 'degen-collide'
            elif classical:
                cls = 'classical'
            elif cert and 5 * bp_max >= N:
                cls = 'semi'
            elif cert:
                cls = 'gridonly'
            else:
                cls = 'fail'
            pair_class[cls] += 1
            per_pair.append((V, p, q, N, eff, m3, argk, cert, cls,
                             (mu_num, mu_den)))
        certs = [t for t in per_pair if t[7]]
        if any(t[8] == 'classical' for t in per_pair):
            set_class['classical'] += 1
        elif any(t[8] == 'semi' for t in per_pair if t[7]):
            set_class['semi'] += 1
        elif certs:
            set_class['gridonly'] += 1
            if len(nonclassical_examples) < 10:
                nonclassical_examples.append(per_pair)
        else:
            set_class['fail'] += 1

    print(f"4-speed sets (gcd 1, v <= {VMAX}): {total}")
    print(f"SET classification: {dict(set_class)}")
    print(f"PAIR classification (6 pairs each): {dict(pair_class)}")
    print(f"A1 violations (grid > unrestricted): {len(grid_gt_mu)}"
          + (f" e.g. {grid_gt_mu[:3]}" if grid_gt_mu else ""))
    print(f"A2 violations (mu < 1/4 genuine triple): {len(mu_violations)}"
          + (f" e.g. {mu_violations[:3]}" if mu_violations else ""))
    print("\ngridonly SET examples (certificate is fully grid-specific):")
    for pp in nonclassical_examples[:6]:
        V = pp[0][0]
        print(f"  V = {V}")
        for (Vv, p, q, N, eff, m3, argk, cert, cls, mu) in pp:
            print(f"    pair ({Vv[p]},{Vv[q]}) N={N} eff={eff} m3={m3} "
                  f"cert={cert} cls={cls} mu={mu[0]}/{mu[1]} argk={argk}")

if __name__ == "__main__":
    main()
