#!/usr/bin/env python3
"""
E2_generalize_n5.py — Reviewer Step 2: does the pair-collapse reduction generalize
beyond n = 4?

Conventions (match scripts/verify_tau_statement_k3.py exactly):
  V: n distinct positive speeds, sorted ascending, gcd(V) = 1.
  Target loneliness: 1/(n+1).
  For a pair (p, q): N = v_p + v_q; effective speeds
     eff = (v_p, v_i for i not in {p,q})   [residues mod N are implicit in ||.||_N]
  tau(p,q) = max over k in Z_N of min_{w in eff} ||w k||_N   (integer units).
  (TAU-n):  some pair has (n+1) * tau(p,q) >= N.

Per (n, speed bound B) battery:
  - pass rates: any pair / some pair containing the max speed / the top-speed
    (two largest) pair / the (min,max) pair;
  - which index-pairs certify (census over certifying pairs);
  - mechanism census of the best pair: j-gon (j | N, j <= n+1, no eff speed
    = 0 mod j) vs plain;
  - internal mirror stats: certifying pairs where two OTHER speeds sum to
    0 mod N (the "iteration available" signal — 4-effective collapses to
    3-effective on the same grid);
  - double-collapse certification: two disjoint pairs with EQUAL sum N
    (both mirrored on Z_N) — does the resulting 3-effective grid problem
    meet the target?
  - failure dumps + tightest sets.

All arithmetic exact (integers only).
"""
from math import gcd
from itertools import combinations
from collections import Counter

def tau_grid(N, eff):
    """max over k in Z_N of min_w ||w*k||_N ; returns (best, argmax_k)."""
    best, argk = 0, 0
    for k in range(1, N):
        m = N
        for w in eff:
            y = (w * k) % N
            d = y if y <= N - y else N - y
            if d < m:
                m = d
                if m == 0:
                    break
        if m == 0:
            continue
        if m > best:
            best, argk = m, k
    return best, argk

def jgon_label(N, eff, nmax):
    """smallest j in [2, n+1] with j | N and all eff speeds not 0 mod j, else None."""
    for j in range(2, nmax + 1):
        if N % j == 0 and all(w % j != 0 for w in eff):
            return j
    return None

def battery(n, B, dump_limit=12):
    T = n + 1
    total = 0
    ok_any = ok_max = ok_top = ok_minmax = 0
    mech = Counter()
    cert_pair_census = Counter()
    only_nonmax = 0            # sets certified but by NO pair containing v_max
    mirror_cert = 0            # sets whose SOME certifying pair has internal mirror
    mirror_cert_needed = 0     # sets where ONLY mirror pairs certify
    dbl_total = dbl_cert = 0   # sets with two disjoint equal-sum pairs; certified by it
    worst = []                 # (best_ratio, V)
    failures = []
    tight_exact = 0

    for V in combinations(range(1, B + 1), n):
        g = 0
        for x in V:
            g = gcd(g, x)
        if g != 1:
            continue
        total += 1
        pairs = list(combinations(range(n), 2))
        pairdata = []
        for (p, q) in pairs:
            N = V[p] + V[q]
            others = [V[i] for i in range(n) if i not in (p, q)]
            eff = (V[p],) + tuple(others)
            best, argk = tau_grid(N, eff)
            cert = (T * best >= N)
            mirror = any((others[a] + others[b]) % N == 0
                         for a in range(len(others)) for b in range(a + 1, len(others)))
            pairdata.append(dict(pq=(p, q), N=N, tau=best, cert=cert, mirror=mirror,
                                 eff=eff, argk=argk))
        certs = [d for d in pairdata if d['cert']]
        best_ratio = max((T * d['tau']) / d['N'] for d in pairdata)
        worst.append((best_ratio, V))

        if certs:
            ok_any += 1
            # mechanism of the best pair
            bd = max(pairdata, key=lambda d: (T * d['tau']) / d['N'])
            j = jgon_label(bd['N'], bd['eff'], T)
            mech['jgon-j%d' % j if j else 'plain'] += 1
            if T * bd['tau'] == bd['N']:
                tight_exact += 1
            for d in certs:
                cert_pair_census[d['pq']] += 1
            maxcert = [d for d in certs if n - 1 in d['pq']]
            if maxcert:
                ok_max += 1
            else:
                only_nonmax += 1
            if (n - 2, n - 1) in [d['pq'] for d in certs]:
                ok_top += 1
            if (0, n - 1) in [d['pq'] for d in certs]:
                ok_minmax += 1
            if any(d['mirror'] for d in certs):
                mirror_cert += 1
                if all(d['mirror'] for d in certs):
                    mirror_cert_needed += 1
        else:
            failures.append(V)

        # double collapse: two disjoint pairs with EQUAL pair sum (both mirrored
        # on the common grid Z_N; the 5-runner min becomes a 3-effective min)
        dbl_pairs = []
        for i1, (p, q) in enumerate(pairs):
            for (y, z) in pairs[i1 + 1:]:
                if len({p, q, y, z}) == 4 and V[p] + V[q] == V[y] + V[z]:
                    dbl_pairs.append((p, q, y, z))
        if dbl_pairs:
            dbl_total += 1
            for (p, q, y, z) in dbl_pairs:
                N = V[p] + V[q]
                reps = (V[p], V[y]) + tuple(V[i] for i in range(n)
                                             if i not in (p, q, y, z))
                t3, _ = tau_grid(N, reps)
                if T * t3 >= N:
                    dbl_cert += 1
                    break   # count the set once

    worst.sort()
    print(f"\n===== n = {n}, speeds <= {B}, target 1/{T} =====")
    print(f"sets (gcd 1): {total}")
    print(f"[any pair]        (TAU-n) certified : {ok_any}/{total}")
    print(f"[max-speed pair]  some pair w/ v_max: {ok_max}/{total}"
          + (f"   (ONLY non-max pairs certify: {only_nonmax})" if only_nonmax else ""))
    print(f"[top-speed pair]  two largest       : {ok_top}/{total}")
    print(f"[(min,max) pair]                    : {ok_minmax}/{total}")
    print(f"mechanism of best pair : {dict(mech)}")
    print(f"best-pair ratio exactly 1 (tight)   : {tight_exact}")
    print(f"certifying pair index census        : {dict(sorted(cert_pair_census.items()))}")
    print(f"some certifying pair has internal mirror : {mirror_cert}/{ok_any}"
          f"   (needed, i.e. ALL certifying pairs mirrored: {mirror_cert_needed})")
    print(f"two disjoint equal-sum pairs present: {dbl_total}/{total};"
          f" certified by the double collapse: {dbl_cert}/{dbl_total}")
    print(f"FAILURES (no pair certifies): {len(failures)}")
    for V in failures[:dump_limit]:
        print(f"   FAIL V={V}")
    print("tightest sets (smallest best-pair ratio):")
    for r, V in worst[:dump_limit]:
        print(f"   ratio={r:.4f}  V={V}")
    return dict(total=total, ok_any=ok_any, ok_max=ok_max, ok_top=ok_top,
                failures=failures)

if __name__ == "__main__":
    battery(3, 30)
    battery(4, 16)
    battery(5, 16)
    battery(6, 14)
