#!/usr/bin/env python3
"""T12b — Four-ball landscape = extended rigid-zone classification, plus
validation of the count-reduction chain.

Prime-level question (equivalent to the pure-unit cell at p^2 via committed
Theorem 4, and the necessary condition for the (1K,3U) cell at p = 1 mod 6):
how many +/--distinct dilates of B_k = {x : |x|_p <= k}, k = floor(p/6), are
needed to cover Z_p?  Rigid zone (n=5, T=6) committed as {7,13,17,19,37},
verified for N <= 79; finiteness beyond is open (worklog line "provable
rigid-zone finiteness").

Reduction chain to validate (proved separately, p = 6k+1):
  4 +/--distinct dilates cover Z_p
    <=> the 4-fold middle-arc intersection A_0 & v1^-1 A_0 & v2^-1 A_0 &
        v3^-1 A_0 is empty
    => A_0 (size 4k) is covered by three dilates mu_j B_k
    => each |mu_j B cap A_0| >= 2k-2, so |B cap mu_j B| <= 3
    => rho(mu_j) = #{i in [1,k] : ||i mu_j||_p <= k} <= 2   (nonzero part).
With min_nu rho(nu) = floor(k/5) (T12a), any 4-covering at p = 6k+1 with
k >= 15 (p >= 97) is impossible -- conditional on the floor lemma.
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


def rho(nu, p, k):
    c = 0
    for i in range(1, k + 1):
        v = (i * nu) % p
        if v <= k or v >= p - k:
            c += 1
    return c


def analyze_prime(p):
    k = p // 6
    B = [r for r in range(p) if min(r, p - r) <= k]      # includes 0
    Bstar = [b for b in B if b != 0]
    full = (1 << p) - 1
    reps = list(range(1, (p - 1) // 2 + 1))              # +/- representatives
    mask = {}
    for lam in reps:
        m = 0
        for b in B:
            m |= 1 << ((lam * b) % p)
        mask[lam] = m
    # valid_lam[x] = set of reps lam with x in lam B  (x != 0)
    inv = {b: pow(b, -1, p) for b in Bstar}
    valid_lam = [None] * p
    for x in range(1, p):
        s = 0
        for b in Bstar:
            lam = (x * inv[b]) % p
            r = min(lam, p - lam)
            s |= 1 << r
        valid_lam[x] = s
    base = mask[1]
    others = [r for r in reps if r != 1]

    sol3 = None
    sol4 = None
    n3 = n4 = 0
    # 3-cover and 4-cover search: WLOG the family contains 1.
    for i2 in range(len(others)):
        c2 = base | mask[others[i2]]
        for i3 in range(i2 + 1, len(others)):
            c3 = c2 | mask[others[i3]]
            if c3 == full:
                if sol3 is None:
                    sol3 = (others[i2], others[i3])
                n3 += 1
                continue
            comp = full & ~c3
            if comp.bit_count() > len(B):
                continue
            # lambda_4 candidates: intersection of valid_lam[x]
            cand = full
            x = comp
            while x:
                lb = (x & -x).bit_length() - 1
                cand &= valid_lam[lb]
                x &= x - 1
                if cand == 0:
                    break
            if cand:
                # exclude 1, lam2, lam3 themselves (+ their +/-, same rep)
                cand &= ~((1 << 1) | (1 << others[i2]) | (1 << others[i3]))
                # also clear bit 0 (not a rep)
                cand &= ~1
                while cand:
                    lb = (cand & -cand).bit_length() - 1
                    if sol4 is None:
                        sol4 = (others[i2], others[i3], lb)
                    n4 += 1
                    cand &= cand - 1
    if sol3 is not None:
        mincover, family = 3, [1, sol3[0], sol3[1]]
    elif sol4 is not None:
        mincover, family = 4, [1, sol4[0], sol4[1], sol4[2]]
    else:
        mincover, family = 5, None
    return dict(p=p, k=k, c6=p % 6, min_dilates=mincover, family=family,
                n3=n3, n4_pairs=n4)


P_MAX = 700
rows = []
for p in primes_upto(P_MAX):
    if p < 7:
        continue
    rows.append(analyze_prime(p))

rigid = [r["p"] for r in rows if r["min_dilates"] <= 4]
print("=== Extended rigid zone (min # +/--distinct dilates of B_k covering Z_p) ===")
for r in rows:
    if r["min_dilates"] <= 4 or r["p"] < 60:
        print(f"p={r['p']:>4} k={r['k']:>3} ({r['c6']} mod 6): min dilates = "
              f"{r['min_dilates']}" +
              (f"  e.g. {r['family']}" if r["family"] else ""))
print("rigid zone (min <= 4) up to", P_MAX, ":", rigid)
print("new beyond committed {7,13,17,19,37}:",
      [p for p in rigid if p not in (7, 13, 17, 19, 37)])
OUT["rigidzone"] = rows

# ---- validate the reduction chain at the rigid primes, p = 1 mod 6 ----
print("\n=== Reduction-chain validation (p = 6k+1 rigid primes) ===")
chain = []
for r in rows:
    p = r["p"]
    if p % 6 != 1 or r["family"] is None:
        continue
    k = r["k"]
    fam = r["family"]
    # the three non-1 dilates must cover the middle arc A_0 = (k, p-k)
    A0 = [x for x in range(k + 1, p - k)]
    cov = set()
    for lam in fam[1:]:
        cov |= {((lam * b) % p) for b in range(-k, k + 1) if b }
    covers_A0 = set(A0) <= cov
    rhos = [rho(lam, p, k) for lam in fam[1:]]
    inA0 = [sum(1 for b in range(1, k + 1) if k < (lam * b) % p < p - k)
            for lam in fam[1:]]
    chain.append(dict(p=p, family=fam, covers_A0=covers_A0,
                      rho_of_units=rhos, units_in_A0=inA0))
    print(f"p={p}: family {fam}: units cover A_0: {covers_A0}, "
          f"rho(lambda_j) = {rhos} (need <= 2), "
          f"#elements of lambda_j B inside A_0 (nonzero) = {inA0} "
          f"(need >= 2k-3 = {2*k-3})")
OUT["chain_validation"] = chain

# ---- rho floor recap vs the kill threshold ----
print("\n=== Floor vs threshold ===")
for r in rows:
    if r["p"] % 6 == 1 and r["p"] >= 43:
        k = r["k"]
        floor = k // 5
        verdict = "KILLED by floor" if floor >= 3 else "floor insufficient"
        if r["min_dilates"] <= 4 and floor >= 3:
            verdict += "  << CONTRADICTION"
        if r["p"] < 120:
            print(f"p={r['p']:>4} k={k:>3}: floor(k/5)={floor}, "
                  f"min_dilates={r['min_dilates']} -> {verdict}")

with open("/home/z/my-project/scripts/out_T12b.json", "w") as f:
    json.dump(OUT, f, indent=1, default=str)
print("\nsaved scripts/out_T12b.json")
