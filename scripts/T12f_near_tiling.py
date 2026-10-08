#!/usr/bin/env python3
"""T12f — The compact near-tiling inequality, cell-relevant sizes.

LEMMA (target, mirrors the user's AP near-tiling lemma):
For prime p = +-1 mod 6, k = floor(p/6), arcs
    I (difference 1), d2*J, d3*K   with d2, d3 not in {+-1}, d2 != +-d3,
sizes n_i in {2k, 2k+1, 2k+2} (the cell's arc sizes), positions arbitrary:
    I u d2J u d3K = Z_p  with overlap eps = sum(n_i) - p <= 2
has NO solutions once p >= 5  (empirically: for all p tested).

Note: unrestricted sizes would admit degenerate huge-arc families
(one arc of size p-C, two small APs); the cell never produces them
(arc sizes are pinned to ~p/3 by the fiber capacity profile).
"""
import json

OUT = {}


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = b"\x00" * len(sieve[i * i:: i])
    return [i for i in range(2, n + 1) if sieve[i]]


def shape_mask(d, n, p):
    m = 0
    for s in range(n):
        m |= 1 << ((d * s) % p)
    return m


def rot(m, t, p):
    """cyclic shift of a p-bit mask by t."""
    return ((m << t) | (m >> (p - t))) & ((1 << p) - 1)


def near_tiling_search(p, eps=2, P_MAX_UNUSED=None):
    k = p // 6
    full = (1 << p) - 1
    sizes = [n for n in (2 * k, 2 * k + 1, 2 * k + 2) if 1 <= n < p]
    # ordered size triples with sum in [p, p+eps]
    triples = []
    for n1 in sizes:
        for n2 in sizes:
            for n3 in sizes:
                s = n1 + n2 + n3
                if p <= s <= p + eps:
                    triples.append((n1, n2, n3))
    if not triples:
        return 0, None
    reps = [r for r in range(2, (p + 1) // 2)]   # d2, d3 reps, != 1
    sh = {}
    for r in reps:
        for n in sizes:
            sh[(r, n)] = shape_mask(r, n, p)
    sols = 0
    example = None
    for (n1, n2, n3) in triples:
        I = shape_mask(1, n1, p)
        for r2 in reps:
            B2 = sh[(r2, n2)]
            for t2 in range(p):
                B2t = rot(B2, t2, p)
                if (I & B2t).bit_count() > eps:
                    continue
                D = I | B2t
                comp = full & ~D
                cb = comp.bit_count()
                if cb == 0 or cb > n3:
                    continue
                for r3 in reps:
                    if r3 == r2:
                        continue
                    B3 = sh[(r3, n3)]
                    cand = full
                    x = comp
                    while x and cand:
                        lb = (x & -x).bit_length() - 1
                        m3 = 0
                        b = B3
                        while b:
                            qb = (b & -b).bit_length() - 1
                            m3 |= 1 << ((lb - qb) % p)
                            b &= b - 1
                        cand &= m3
                        x &= x - 1
                    if cand:
                        sols += 1
                        if example is None:
                            t3 = (cand & -cand).bit_length() - 1
                            example = (n1, n2, n3, r2, r3, t2, t3)
    return sols, example


rows = []
for p in primes_upto(155):
    if p % 6 not in (1, 5) or p < 5:
        continue
    sols, ex = near_tiling_search(p)
    rows.append((p, sols, ex))
    print(f"p={p:>4} ({p % 6} mod 6, k={p // 6}): near-tilings "
          f"(sizes ~p/3, eps<=2, +-distinct): {sols}"
          + (f"  e.g. {ex}" if ex else ""))

survivors = [(p, s, e) for p, s, e in rows if s > 0]
print("\nsurvivors:", survivors or "NONE — the inequality holds throughout")
with open("/home/z/my-project/scripts/out_T12f.json", "w") as f:
    json.dump(rows, f, indent=1, default=str)
print("saved scripts/out_T12f.json")
