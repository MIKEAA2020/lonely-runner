#!/usr/bin/env python3
"""T12c — Full structured (1K,3U) cell search at N = p^2.

Question: does there exist a covering
    B_{p a} U B_{u1} U B_{u2} U B_{u3} = Z_{p^2}
with kernel velocity p.a (gcd(a,p)=1) and three units u_i in Z_{p^2}^*?

Normalization (proved): multiply all velocities by a unit c = a^{-1} + t.p
with the right t; the kernel becomes exactly p and v_1 = 0.  So enumerate
(u'_1, v_1=0), (u'_2, v_2), (u'_3, v_3) with u'_i = u_i mod p, v_i the
p-digit.  Config space per p: (p-1) * [p(p-1)]^3 / (normalization) ->
(p-1)^3 * p^2 free parameters, searched via per-class good sets.

Per class r in Z_p \\ B_k (non-kernel classes), the units must cover the
whole fiber {q p + r}:
    Arc_i(r) = { q : (u'_i q + T_i(r) + v_i r) mod p in S(r'_i) },
    T_i(r) = floor(u'_i r / p),  r'_i = u'_i r mod p,
    S(r')  = { s in [0,p) : 6*dist(r' + p s, p^2) < p^2 }.
As v_i varies, Arc_i(r) = B_i + t_i is a full translate family
(t_i = -u'^-1_i (T_i + v_i r), bijective since r != 0 on non-kernel
classes).  Good set G_r = {(v_2,v_3) covering class r}; config alive iff
intersection over all non-kernel classes is nonempty.

Mod-p filters (necessary, proved from the capacity profiles):
  p = 6k+1: every non-kernel class r needs m(r) >= 1, i.e.
            Z_p \\ B_k  <=  u_1'^-1 B_k U u_2'^-1 B_k U u_3'^-1 B_k
            (the four-ball covering condition).
  p = 6k+5: every non-kernel class r needs m(r) <= 1 (avoidance).

Anchors: p=3 (N=9) MUST find the committed covering B_1 U B_2 U B_3 U B_4;
p=5,7,11,13 (N=25,49,121,169) MUST be empty (committed).
"""
import json
from itertools import combinations_with_replacement

OUT = {}


def dist_mod(x, N):
    r = x % N
    return min(r, N - r)


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = b"\x00" * len(sieve[i * i:: i])
    return [i for i in range(2, n + 1) if sieve[i]]


def search_cell(p, verbose=True):
    N = p * p
    k = p // 6
    Bk = set(r for r in range(p) if min(r, p - r) <= k)
    nonkernel = [r for r in range(p) if r not in Bk]
    # S(r') sets
    S = []
    for rp in range(p):
        S.append(frozenset(s for s in range(p)
                           if 6 * dist_mod(rp + p * s, N) < N))
    units = list(range(1, p))
    inv = {u: pow(u, -1, p) for u in units}
    balls = {u: frozenset((u * r) % p for r in Bk) for u in units}

    # ---- mod-p filter over residue triples (with repetition allowed) ----
    passing = []
    for triple in combinations_with_replacement(units, 3):
        if p % 6 == 1:
            cov = balls[triple[0]] | balls[triple[1]] | balls[triple[2]]
            if not (set(nonkernel) <= cov):
                continue
        else:  # p % 6 in {5, 3}
            ok = True
            for r in nonkernel:
                m = sum(1 for u in triple if (u * r) % p in Bk)
                if p % 6 == 5 and m > 1:
                    ok = False
                    break
                if p == 3 and m > 0:  # p=3: capacity = 3 = p forces m = 0
                    ok = False
                    break
            if not ok:
                continue
        passing.append(triple)
    if verbose:
        print(f"p={p:>3} (k={k}, {p % 6} mod 6): {len(nonkernel)} non-kernel "
              f"classes; {len(passing)} residue triples pass the mod-p filter")

    # ---- per-triple v-search via good sets ----
    coverings = []
    diag = []
    for triple in passing:
        u1, u2, u3 = triple
        good = None  # dict v2 -> set(v3)
        killed_at = None
        for r in nonkernel:
            # fixed shapes
            T1, T2, T3 = (u1 * r) // p, (u2 * r) // p, (u3 * r) // p
            rp1, rp2, rp3 = (u1 * r) % p, (u2 * r) % p, (u3 * r) % p
            A1 = frozenset((inv[u1] * (s - T1)) % p for s in S[rp1])
            B2 = frozenset((inv[u2] * s) % p for s in S[rp2])
            B3 = frozenset((inv[u3] * s) % p for s in S[rp3])
            # v -> t maps: t_i = -inv[u_i] * (T_i + v_i * r)
            # v_i = (-u_i * t_i - T_i) * inv[r]
            invr = inv[r % p] if r % p else None
            Gr = {}
            # iterate v2 directly (t2 affine in v2)
            step2 = (-inv[u2] * r) % p
            base2 = (-inv[u2] * T2) % p
            step3 = (-inv[u3] * r) % p
            base3 = (-inv[u3] * T3) % p
            A1m = 0
            for q in A1:
                A1m |= 1 << q
            B3list = sorted(B3)
            full = (1 << p) - 1
            for v2 in range(p):
                t2 = (base2 + step2 * v2) % p
                B2m = 0
                for q in B2:
                    B2m |= 1 << ((q + t2) % p)
                D = A1m | B2m
                comp = full & ~D
                if comp == 0:
                    # class covered by A1 and arc2 alone; all v3 allowed
                    Gr[v2] = set(range(p))
                    continue
                if comp.bit_count() > len(B3):
                    continue
                # valid t3: comp \u2286 B3 + t3
                cand = full
                x = comp
                ok = True
                while x:
                    lb = (x & -x).bit_length() - 1
                    m3 = 0
                    for q in B3list:
                        m3 |= 1 << ((lb - q) % p)
                    cand &= m3
                    x &= x - 1
                    if cand == 0:
                        ok = False
                        break
                if cand:
                    v3s = set()
                    t = cand
                    while t:
                        lb = (t & -t).bit_length() - 1
                        v3 = ((-u3 * lb - T3) * invr) % p
                        v3s.add(v3)
                        t &= t - 1
                    Gr[v2] = v3s
            if not Gr:
                killed_at = (r, "empty good set")
                break
            if good is None:
                good = {v2: set(v3s) for v2, v3s in Gr.items()}
            else:
                good = {v2: (good[v2] & Gr[v2]) for v2 in good
                        if v2 in Gr}
                good = {v2: s for v2, s in good.items() if s}
                if not good:
                    killed_at = (r, "intersection empty")
                    break
        if killed_at is None and good:
            for v2, v3s in good.items():
                for v3 in v3s:
                    coverings.append((triple, v2, v3))
        elif len(diag) < 6:
            diag.append((triple, killed_at))
    return dict(p=p, k=k, c6=p % 6, n_triples=len(passing),
                coverings=coverings, diagnostics=diag)


PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
results = []
for p in PRIMES:
    res = search_cell(p)
    results.append(res)
    ncov = len(res["coverings"])
    print(f"   -> coverings found: {ncov}"
          + (f"  {res['coverings'][:3]}" if ncov else ""))
    if not ncov and res["diagnostics"]:
        print(f"   -> first kills: {res['diagnostics'][:3]}")

print("\n=== Summary (1K,3U) cell at p^2 ===")
for res in results:
    p = res["p"]
    status = "ALIVE" if res["coverings"] else "dead"
    print(f"p={p:>3} N={p*p:>5}: {status:>5}  "
          f"(filter-passing triples: {res['n_triples']})")

with open("/home/z/my-project/scripts/out_T12c.json", "w") as f:
    json.dump(results, f, indent=1, default=str)
print("saved scripts/out_T12c.json")
