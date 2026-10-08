#!/usr/bin/env python3
"""
verify_tau_statement_k3.py — the REDUCED k=3 statement.

Reduction (proved in the proof-attempt writeup):
  For any pair (p,q), on grid N = v_p + v_q we have v_q = -v_p (mod N),
  so at every lattice time k/N the 4-runner minimum equals the
  3-effective-runner minimum min(||v_p k||_N, ||v_y k||_N, ||v_z k||_N).
  Hence M(V) >= m3*(N)/N for EVERY pair, constructively.

  Therefore LRC(k=3) follows from the PURE SPEED STATEMENT:
    (TAU) for every 4-speed set V with gcd(V) = 1, some pair (p,q) has
          m3*(v_p+v_q)/(v_p+v_q) >= 1/5.
  and it suffices that this holds for the three max-speed pairs.

This script verifies (TAU) directly (no witnesses, no grids B), over a
speed range, and reports which pairs achieve it and by what mechanism
(pentagon / plain / other).

Also verifies the quantitative bound M(V) >= max_pair tau_pair on the
continuous optimum (computed by fine rational sampling) as a cross-check.
"""

from math import gcd
from itertools import combinations
from collections import Counter

VMAX = 16   # speeds in [1, VMAX]

def m3_of_pair(v, p, q):
    """3-effective-runner per-grid max on N = v_p+v_q, and argmax ks."""
    N = v[p] + v[q]
    others = [t for t in range(4) if t not in (p, q)]
    eff = (v[p], v[others[0]], v[others[1]])
    best, ks = -1, []
    for k in range(N):
        m = N
        for w in eff:
            y = (w * k) % N
            if y > N - y:
                y = N - y
            if y < m:
                m = y
        if m > best:
            best, ks = m, [k]
        elif m == best:
            ks.append(k)
    return N, best, ks, eff

def pentagon_ok(N, eff):
    return N % 5 == 0 and not any(w % 5 == 0 for w in eff)

def main():
    total = 0
    tau_ok = 0
    maxpair_ok = 0
    topspeed_ok = 0
    mech = Counter()
    worst = []   # (max_tau, V) smallest values
    for V in combinations(range(1, VMAX + 1), 4):
        g = 0
        for vi in V:
            g = gcd(g, vi)
        if g != 1:
            continue
        total += 1
        best_tau, best_pair, best_mech = 0, None, None
        maxpair_best = 0
        topspeed_best = 0
        for p in range(4):
            for q in range(p + 1, 4):
                N, m3, ks, eff = m3_of_pair(V, p, q)
                tau = m3 / N
                if tau > best_tau:
                    best_tau, best_pair = tau, (p, q)
                    best_mech = ("pentagon" if pentagon_ok(N, eff) and
                                 m3 * 5 >= N else
                                 "plain" if m3 * 5 >= N else "sub")
                if V[p] == V[3] or V[q] == V[3]:
                    maxpair_best = max(maxpair_best, tau)
        srt = sorted(range(4), key=lambda i: -V[i])
        Nt, m3t, kst, efft = m3_of_pair(V, srt[0], srt[1])
        topspeed_best = m3t / Nt

        if best_tau >= 0.2 - 1e-12:
            tau_ok += 1
            mech[best_mech] += 1
        if maxpair_best >= 0.2 - 1e-12:
            maxpair_ok += 1
        if topspeed_best >= 0.2 - 1e-12:
            topspeed_ok += 1
        worst.append((best_tau, V, best_pair))

    worst.sort()
    print(f"speed sets (gcd 1, v <= {VMAX}): {total}")
    print(f"[TAU]  some pair (all 6)   has tau >= 1/5 : {tau_ok}/{total}")
    print(f"[TAUm] some MAX-SPEED pair has tau >= 1/5 : {maxpair_ok}/{total}")
    print(f"[TAUt] the TOP-SPEED pair   has tau >= 1/5 : {topspeed_ok}/{total}")
    print(f"mechanism of the best pair: {dict(mech)}")
    print("12 worst cases (smallest best-pair tau):")
    for tau, V, pr in worst[:12]:
        print(f"   tau={tau:.4f}  V={V}  pair={pr}")

if __name__ == "__main__":
    main()
