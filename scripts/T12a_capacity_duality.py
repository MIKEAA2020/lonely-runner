#!/usr/bin/env python3
"""T12a — Capacity profiles, duality, and the return-count landscape.

Context (turn T-12, user's direct attack on the (1K,3U) cell at p^2):
The user's message claims, at p = 6k+5, "L(r') <= 2k+1 always", and concludes
the (1K,3U) cell is dead for p = 5 mod 6 by pure capacity.  The committed
Lemma 1 fiber cap is cap(p) = 2*floor(p/6)+1+[p=5 mod 6], which equals 2k+2
at p = 6k+5 -- contradicting that claim.  This script settles it by direct
enumeration, then verifies the m-duality structure and computes the
return-count landscape rho(nu) used by the four-ball analysis.

Locked convention (nua_proof_verify.py): bad_set(w,N) = {k : 6*min(wk mod N,
N - wk mod N) < N}, strict, symmetric.
"""
import json
from math import gcd
from random import Random

T = 6
OUT = {}


def dist_mod(x, N):
    r = x % N
    return min(r, N - r)


def bad_set(w, N):
    return frozenset(k for k in range(N) if T * dist_mod(w * k, N) < N)


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = b"\x00" * len(sieve[i * i:: i])
    return [i for i in range(2, n + 1) if sieve[i]]


# ---------------------------------------------------------------- Part A
# Fiber capacity profile at N = p^2.
# For a unit u with residue u' = u mod p and class r in Z_p (k = qp + r),
# the fiber of B_u over class r depends only on r' := u' * r mod p:
#   L(r') = |{ s in [0,p) : 6*dist(r' + p*s, p^2) < p^2 }|.
def fiber_L(rp, p):
    N = p * p
    return sum(1 for s in range(p) if T * dist_mod(rp + p * s, N) < N)


def ball(p):
    k = p // 6
    return frozenset(r for r in range(p) if min(r, p - r) <= k)


partA = []
for p in primes_upto(130):
    if p < 5:
        continue
    k = p // 6
    Bk = ball(p)
    L = [fiber_L(rp, p) for rp in range(p)]
    cap_committed = 2 * (p // 6) + 1 + (1 if p % 6 == 5 else 0)
    # my corrected profile
    if p % 6 == 1:
        off = 2 * k
    else:
        off = 2 * k + 2
    mine_ok = all(L[rp] == (2 * k + 1 if rp in Bk else off) for rp in range(p))
    # user's p=1 mod 6 table: 2k+1 on B_k, 2k on middle  (same as mine there)
    # user's p=5 mod 6 claim: L <= 2k+1 everywhere
    user5_ok = all(x <= 2 * k + 1 for x in L) if p % 6 == 5 else None
    # sanity: fiber sizes also match direct bad_set of a unit at p^2 in class r
    partA.append(dict(p=p, k=k, c=p % 6, Lmin=min(L), Lmax=max(L),
                      cap=cap_committed, cap_ok=(max(L) == cap_committed),
                      corrected_profile_ok=mine_ok,
                      user_p5_claim_ok=user5_ok,
                      multiset=sorted(set(L))))

# direct cross-check of the fiber formula on a real unit at p^2
check_fib = []
rng = Random(1)
for p in [5, 7, 11, 13, 17, 19, 29, 31]:
    N = p * p
    for _ in range(6):
        u = rng.randrange(1, N)
        if gcd(u, p) != 1:
            continue
        up = u % p
        B = bad_set(u, N)
        for r in list(range(p))[: p]:
            fiber = frozenset(q for q in range(p) if (q * p + r) in B)
            rp = (up * r) % p
            if len(fiber) != fiber_L(rp, p):
                raise AssertionError(f"fiber formula mismatch p={p} u={u} r={r}")
    check_fib.append(p)
OUT["A_fiber_formula_crosscheck_moduli"] = check_fib

# ---------------------------------------------------------------- Part B
# Duality.  capacity(r) = sum_i L(u'_i * r mod p).  Verify:
#   p = 6k+1: capacity - p = m(r) - 1,  m = #{i : u'_i r in B_k}
#   p = 6k+5: capacity - p = 1 - m(r)
# and the resulting necessary conditions on non-kernel classes.
partB = []
for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    k = p // 6
    Bk = ball(p)
    L = [fiber_L(rp, p) for rp in range(p)]
    units = [u for u in range(1, p)]
    rng = Random(p)
    trials = 4000 if p < 30 else 8000
    bad_cases = 0
    max_cap, min_cap = 0, 10 ** 9
    for _ in range(trials):
        us = rng.sample(units, 3)
        for r in range(p):
            if r in Bk:
                continue
            cap = sum(L[(u * r) % p] for u in us)
            m = sum(1 for u in us if (u * r) % p in Bk)
            if p % 6 == 1:
                if cap - p != m - 1:
                    bad_cases += 1
            else:
                if cap - p != 1 - m:
                    bad_cases += 1
            max_cap = max(max_cap, cap)
            min_cap = min(min_cap, cap)
    partB.append(dict(p=p, k=k, c=p % 6, mismatches=bad_cases,
                      cap_min=min_cap, cap_max=max_cap, p_is=p))
OUT["B_duality"] = partB

# ---------------------------------------------------------------- Part C
# Return-count landscape: rho(nu) = #{ i in [1,k] : ||i*nu||_p <= k },
# p = 6k+1.  The four-ball reduction (proved separately) forces all three
# witness multipliers to have rho <= 2, so min_nu rho(nu) is the key curve.
partC = []
for p in primes_upto(700):
    if p < 43 or p % 6 != 1:
        continue
    k = p // 6
    best, args = None, []
    for nu in range(1, p):
        c = 0
        for i in range(1, k + 1):
            v = (i * nu) % p
            if v <= k or v >= p - k:
                c += 1
        if best is None or c < best:
            best, args = c, [nu]
        elif c == best:
            args.append(nu)
    partC.append(dict(p=p, k=k, min_rho=best, argmin=args[:8],
                      n_argmin=len(args)))
OUT["C_return_count_landscape"] = partC

# ---------------------------------------------------------------- report
print("=== Part A: fiber capacity profiles at N = p^2 (direct enumeration) ===")
print(f"{'p':>4} {'k':>3} {'mod6':>4} {'Lmin':>4} {'Lmax':>4} {'cap':>4} "
      f"{'capOK':>5} {'corrOK':>6} {'userP5OK':>8}  profile")
for row in partA:
    print(f"{row['p']:>4} {row['k']:>3} {row['c']:>4} {row['Lmin']:>4} "
          f"{row['Lmax']:>4} {row['cap']:>4} {str(row['cap_ok']):>5} "
          f"{str(row['corrected_profile_ok']):>6} "
          f"{str(row['user_p5_claim_ok']):>8}  {row['multiset']}")

n5 = [r for r in partA if r["c"] == 5]
print("\np=5 mod 6 primes where the user's claim 'L <= 2k+1 always' holds:",
      [r["p"] for r in n5 if r["user_p5_claim_ok"]] or "NONE")
print("p=5 mod 6 primes where corrected profile (2k+1 on B_k, 2k+2 off) holds:",
      [r["p"] for r in n5 if r["corrected_profile_ok"]])
print("committed Lemma-1 cap 2*floor(p/6)+1+[p=5 mod 6] equals max L:",
      all(r["cap_ok"] for r in partA))

print("\n=== Part B: duality capacity checks (random unit triples) ===")
for row in partB:
    print(f"p={row['p']:>3} (={row['c']} mod 6): mismatches={row['mismatches']}, "
          f"cap range [{row['cap_min']},{row['cap_max']}] vs p={row['p_is']} "
          f"-> excess range [{row['cap_min']-row['p_is']},{row['cap_max']-row['p_is']}]")

print("\n=== Part C: min_nu rho(nu), p = 6k+1, p >= 43 ===")
for row in partC:
    print(f"p={row['p']:>4} k={row['k']:>3} min_rho={row['min_rho']:>3} "
          f"argmin={row['argmin'][:5]} (#{row['n_argmin']} minimizers)")
lows = [r for r in partC if r["min_rho"] <= 2]
print("primes with min_rho <= 2 (four-ball not killed by the count condition):",
      [(r["p"], r["min_rho"]) for r in lows] or "NONE for p >= 43")

with open("/home/z/my-project/scripts/out_T12a.json", "w") as f:
    json.dump(OUT, f, indent=1)
print("\nsaved scripts/out_T12a.json")
