#!/usr/bin/env python3
"""T12e — Kill-level statistics for the (1K,3U) cell at p^2 (full re-run of
the T12c search with complete diagnostics). For every filter-passing residue
triple we record whether it dies
  (b) 'abstract'  -- some class has an EMPTY good set: no (v2, v3) covers
      that class even with free translations (per-class impossibility), or
  (c) 'coupling'  -- every class is individually coverable, but the joint
      intersection over classes is empty.
Also verify every found covering by direct bad_set evaluation (fresh check).
"""
import json
from itertools import combinations_with_replacement
from math import gcd

OUT = {}


def dist_mod(x, N):
    r = x % N
    return min(r, N - r)


def bad_set(w, N):
    return frozenset(k for k in range(N) if 6 * dist_mod(w * k, N) < N)


def search_cell_stats(p):
    N = p * p
    k = p // 6
    Bk = set(r for r in range(p) if min(r, p - r) <= k)
    nonkernel = [r for r in range(p) if r not in Bk]
    S = [frozenset(s for s in range(p) if 6 * dist_mod(rp + p * s, N) < N)
         for rp in range(p)]
    units = list(range(1, p))
    inv = {u: pow(u, -1, p) for u in units}
    balls = {u: frozenset((u * r) % p for r in Bk) for u in units}

    passing = []
    for triple in combinations_with_replacement(units, 3):
        if p % 6 == 1:
            cov = balls[triple[0]] | balls[triple[1]] | balls[triple[2]]
            if not (set(nonkernel) <= cov):
                continue
        elif p % 6 == 5:
            ok = True
            for r in nonkernel:
                if sum(1 for u in triple if (u * r) % p in Bk) > 1:
                    ok = False
                    break
            if not ok:
                continue
        passing.append(triple)

    n_abstract = n_coupling = n_alive = 0
    abstract_classes = {}
    alive_examples = []
    for triple in passing:
        u1, u2, u3 = triple
        good = None
        died = None
        for r in nonkernel:
            T1, T2, T3 = (u1 * r) // p, (u2 * r) // p, (u3 * r) // p
            rp1, rp2, rp3 = (u1 * r) % p, (u2 * r) % p, (u3 * r) % p
            A1 = frozenset((inv[u1] * (s - T1)) % p for s in S[rp1])
            B2 = frozenset((inv[u2] * s) % p for s in S[rp2])
            B3 = frozenset((inv[u3] * s) % p for s in S[rp3])
            invr = inv[r]
            step2 = (-inv[u2] * r) % p
            base2 = (-inv[u2] * T2) % p
            A1m = 0
            for q in A1:
                A1m |= 1 << q
            B3list = sorted(B3)
            full = (1 << p) - 1
            Gr = {}
            for v2 in range(p):
                t2 = (base2 + step2 * v2) % p
                B2m = 0
                for q in B2:
                    B2m |= 1 << ((q + t2) % p)
                D = A1m | B2m
                comp = full & ~D
                if comp == 0:
                    Gr[v2] = set(range(p))
                    continue
                if comp.bit_count() > len(B3):
                    continue
                cand = full
                x = comp
                while x and cand:
                    lb = (x & -x).bit_length() - 1
                    m3 = 0
                    for q in B3list:
                        m3 |= 1 << ((lb - q) % p)
                    cand &= m3
                    x &= x - 1
                if cand:
                    v3s = set()
                    t = cand
                    while t:
                        lb = (t & -t).bit_length() - 1
                        v3s.add(((-u3 * lb - T3) * invr) % p)
                        t &= t - 1
                    Gr[v2] = v3s
            if not Gr:
                m_r = sum(1 for u in triple if (u * r) % p in Bk)
                died = ("abstract", r, m_r)
                abstract_classes[m_r] = abstract_classes.get(m_r, 0) + 1
                break
            if good is None:
                good = {v2: set(v3s) for v2, v3s in Gr.items()}
            else:
                good = {v2: (good[v2] & Gr[v2]) for v2 in good if v2 in Gr}
                good = {v2: s for v2, s in good.items() if s}
                if not good:
                    died = ("coupling", r, None)
                    break
        if died is None:
            found = False
            for v2, v3s in (good or {}).items():
                for v3 in v3s:
                    u_1 = u1                      # v1 = 0
                    u_2 = u2 + v2 * p
                    u_3 = u3 + v3 * p
                    fam = [p, u_1, u_2, u_3]
                    if len(set(fam)) < 4:
                        continue
                    sets = [bad_set(w, N) for w in fam]
                    u = frozenset()
                    for sset in sets:
                        u |= sset
                    if len(u) == N:
                        alive_examples.append((triple, v2, v3, fam))
                        found = True
            if found:
                n_alive += 1
            else:
                n_coupling += 1     # all survivors were degenerate collisions
        elif died[0] == "abstract":
            n_abstract += 1
        else:
            n_coupling += 1
    return dict(p=p, n_triples=len(passing), n_abstract=n_abstract,
                n_coupling=n_coupling, n_alive=n_alive,
                abstract_classes=abstract_classes,
                alive_examples=alive_examples[:4],
                verified_coverings=[ex[3] for ex in alive_examples[:4]])


PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 47]
stats = []
for p in PRIMES:
    res = search_cell_stats(p)
    stats.append(res)
    print(f"p={p:>3} ({p % 6} mod 6): triples={res['n_triples']:>4}, "
          f"abstract-kills={res['n_abstract']:>4}, "
          f"coupling-kills={res['n_coupling']:>4}, "
          f"ALIVE={res['n_alive']}")
    if res["abstract_classes"]:
        print(f"        abstract kills by class-m-value: "
              f"{res['abstract_classes']}")
    if res["verified_coverings"]:
        print(f"        verified covering families (fresh bad_set check): "
              f"{res['verified_coverings']}")

with open("/home/z/my-project/scripts/out_T12e.json", "w") as f:
    json.dump(stats, f, indent=1, default=str)
print("saved scripts/out_T12e.json")
