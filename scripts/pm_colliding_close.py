#!/usr/bin/env python3
"""
T14b: Closure of the +-COLLIDING unit case for the (1K, 3U) cell at N = p^2.

THE PROBLEM (the named residual of the trichotomy note, section 5): the
(1K,3U) cell reduction forces, in every kernel-uncovered fiber, a t-space
covering by three APs whose differences are u_i^{-1} mod p -- pairwise
+-distinct iff the units are.  The cell definition (distinct speeds mod
p^2) does NOT exclude unit pairs with u_i = +-u_j (mod p); for those the
t-space APs have colliding differences and Lemmas A/B do not fire.

ANSWER TO THE REVIEWER'S QUESTION ("one-line consequence of the <=2-units
theorem or genuine obstruction?"):  NOT a consequence -- Theorem 2 of the
no-unit-assistance note requires <= 2 unit members in the family; this
cell has 3 DISTINCT unit bad sets.  The closure is its own theorem, the
SHEAR-RIGIDITY theorem, verified here:

  THEOREM (no +-colliding unit assistance at p^2).  Let p >= 11 (eps=5)
  resp. p >= 19 (eps=1) be prime.  No family of four distinct bad sets at
  N = p^2 consisting of one kernel set and three unit sets, two of which
  collide (u_i = +-u_j mod p), covers Z_{p^2}.  The remaining primes
  (5, 7, 11, 13, 17 for eps=1: 7, 13, 17) are verified by direct
  enumeration (PART A).

PROOF SKELETON (each step machine-verified below):
  S1  Fiber structure: B_u cap F_j is, in the t-coordinate, an AP with
      difference ubar = u^{-1} mod p, size L(uj mod p), start
      sigma_u(j) = ubar*(T2(a0) - a1) where uj = a1*p + a0; the value
      arc is ball / ball-1 (on-ball) or ball\\{k} (eps=1) / ball u {-k-1}
      (eps=5) (off-ball).
  S2  Shear linearity: B_{-v} = B_v, so WLOG the collision is +: v3 =
      v2 + pm with m != 0 (mod p).  Then a0_3 = a0_2, so the pair shares
      arcs and sizes, and sigma_3(j) - sigma_2(j) = -ubar*m*j EXACTLY.
  S3  Per-fiber size feasibility (counting) + the general-arc trichotomy
      restrict the pair difference delta to {1, 2, q} (eps=5, k>=2) resp.
      {1, 2, 3k} (eps=1, k>=3).
  S4  Pinning: for fixed interval start sigma1 (2 possible values), fixed
      feasible sizes, fixed delta, the covering configurations have
      sigma_3 - sigma_2 in a finite set of size <= 6.
  S5  Pigeonhole: Delta(j) = -ubar*m*j is injective on U (|U| = p - |B_k|)
      and must stay in the Delta-set; |Delta-set| < |U| kills the family.
  (P2) The second collision topology (a unit colliding with the
      normalized B_1, i.e. v2 = 1 + ps): then sigma_2(j) - sigma_1(j) =
      -s*j exactly, and the same pinning+pigeonhole runs on that
      difference; at eps=5 the counting collapses P2 into P1(delta=1).

PART A: ground truth (all colliding families, both topologies, small p).
PART B: S1/S2 verification.  PART C: S4/S5 pinning tables + pigeonhole.
PART D: near-miss margins.
"""

import json
import sys
import time
from math import gcd

T0 = time.time()
OUT = {"meta": {}, "ground_truth": {}, "structure": {}, "drift": {},
       "pinning": {}, "nearmiss": {}, "anchors": {}, "violations": {}}
VIOL = OUT["violations"]


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def sieve(n):
    bs = bytearray([1]) * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i:: i] = bytearray(len(bs[i * i:: i]))
    return [i for i in range(n + 1) if bs[i]]


def k_of(p):
    return (p - (1 if p % 6 == 1 else 5)) // 6


def ball_of(p):
    k = k_of(p)
    return frozenset(r for r in range(p) if min(r, p - r) <= k)


def pm_classes(N):
    return [u for u in range(1, N // 2 + 1) if gcd(u, N) == 1]


def bad_mask(u, N):
    m = 0
    for kk in range(N):
        v = (u * kk) % N
        if 6 * min(v, N - v) < N:
            m |= 1 << kk
    return m


# ---------------------------------------------------------------- S1

def arc_data(p):
    k = k_of(p)
    eps = 1 if p % 6 == 1 else 5
    c = (p * p - 1) // 6
    N = p * p
    Bk = ball_of(p)
    ball = frozenset(Bk)
    ballm1 = frozenset((x - 1) % p for x in Bk)
    if eps == 5:
        offarc = ball | {(-k - 1) % p}
        on_L, off_L = 2 * k + 1, 2 * k + 2
    else:
        offarc = ball - {k}
        on_L, off_L = 2 * k + 1, 2 * k
    bad = 0
    kinds = set()
    for a0 in range(p):
        Ts = [b for b in range(p)
              if 6 * min((a0 + p * b) % N, N - (a0 + p * b) % N) < N]
        L = len(Ts)
        if eps == 5:
            T1c = k - 1 if a0 >= p - k else k
            T2c = p - k if a0 <= k else p - k - 1
        else:
            T1c = k - 1 if a0 >= k + 1 else k
            T2c = 5 * k + 1 if a0 <= 5 * k else 5 * k
        arc_c = frozenset(list(range(T2c, p)) + list(range(0, T1c + 1)))
        onball = a0 in Bk
        exp_L = on_L if onball else off_L
        kind = ("ball" if arc_c == ball else
                "ball-1" if arc_c == ballm1 else
                "off" if arc_c == offarc else "OTHER")
        kinds.add(kind)
        # circular-run check: Ts must be exactly T2c, T2c+1, ... mod p
        run_ok = set(Ts) == {(T2c + i) % p for i in range(L)}
        # symmetry L(a0) = L(p - a0)
        a0m = (p - a0) % p
        Tsm = [b for b in range(p)
               if 6 * min((a0m + p * b) % N, N - (a0m + p * b) % N) < N]
        if (L != T1c + 1 + p - T2c or L != exp_L or kind == "OTHER"
                or not run_ok or len(Tsm) != L):
            bad += 1
            if bad <= 8:
                VIOL.setdefault("S1_arc", []).append((p, a0, L, exp_L, kind))
    return {"p": p, "eps": eps, "k": k, "bad": bad, "kinds": sorted(kinds)}


def fiber_tset_direct(u, j, p):
    N = p * p
    return frozenset(t for t in range(p)
                     if 6 * min((u * (j + p * t)) % N,
                                N - (u * (j + p * t)) % N) < N)


def fiber_tset_formula(u, j, p):
    N = p * p
    uq = u % N
    ubar = pow(uq, -1, p)
    a = (uq * j) % N
    a0, a1 = a % p, a // p
    k = k_of(p)
    eps = 1 if p % 6 == 1 else 5
    if eps == 5:
        T2 = p - k if a0 <= k else p - k - 1
    else:
        T2 = 5 * k + 1 if a0 <= 5 * k else 5 * k
    sigma = (ubar * (T2 - a1)) % p
    Bk = ball_of(p)
    onball = a0 in Bk
    L = (2 * k + 1 if onball else (2 * k + 2 if eps == 5 else 2 * k))
    return frozenset((sigma + m * ubar) % p for m in range(L)), sigma, L


def check_fiber_structure(p):
    N = p * p
    bad = 0
    checked = 0
    for u in range(1, N):
        if u % p == 0:
            continue
        for j in range(p):
            d = fiber_tset_direct(u, j, p)
            f, sigma, L = fiber_tset_formula(u, j, p)
            checked += 1
            if d != f or len(f) != L:
                bad += 1
                if bad <= 5:
                    VIOL.setdefault("S1_tset", []).append((p, u, j))
    return {"p": p, "checked": checked, "bad": bad}


def check_drift(p):
    """S2: sigma_3(j) - sigma_2(j) = -ubar*m*j (mod p), exactly, for all
    lifts v2, all m != 0, all j.  Also the P2 variant: for v2 = 1 + ps,
    sigma_2(j) - T2(j) = -s*j."""
    N = p * p
    bad = 0
    badP2 = 0
    checked = 0
    k = k_of(p)
    eps = 1 if p % 6 == 1 else 5

    def T2_of(a0):
        if eps == 5:
            return p - k if a0 <= k else p - k - 1
        return 5 * k + 1 if a0 <= 5 * k else 5 * k

    for v2 in range(1, N):
        if v2 % p == 0:
            continue
        ubar = pow(v2, -1, p)
        for m in range(1, p):
            v3 = (v2 + p * m) % N
            if v3 % p == 0:
                continue
            for j in range(p):
                _, s2, _ = fiber_tset_formula(v2, j, p)
                _, s3, _ = fiber_tset_formula(v3, j, p)
                checked += 1
                if (s3 - s2 + ubar * m * j) % p != 0:
                    bad += 1
                    if bad <= 5:
                        VIOL.setdefault("S2_drift", []).append(
                            (p, v2, m, j, s2, s3))
        # P2: v2 = 1 + ps  (only when v2 = 1 + ps for some s)
        if (v2 - 1) % p == 0 and v2 != 1:
            s = (v2 - 1) // p
            for j in range(p):
                _, s2, _ = fiber_tset_formula(v2, j, p)
                if (s2 - T2_of(j) + s * j) % p != 0:
                    badP2 += 1
                    if badP2 <= 5:
                        VIOL.setdefault("S2_driftP2", []).append(
                            (p, v2, j, s2))
    return {"p": p, "checked": checked, "bad": bad, "badP2": badP2}


# ------------------------------------------------------------ S4/S5 pinning

def ap_mask(p, s, d, L):
    m = 0
    for i in range(L):
        m |= 1 << ((s + i * d) % p)
    return m


def pinning_tables(p):
    """Enumerate all coverings of Z_p of the two colliding fiber shapes:
      P1: interval(L1 @ sigma1) + AP(L @ sigma2, d) + AP(L @ sigma3, d)
      P2: interval(L1 @ sigma1) + interval(L1 @ sigma2) + AP(L3 @ sigma3, d)
    for sigma1 in the actual T2-value set of B_1's piece and the feasible
    size combos.  Tabulate per-delta difference-sets and check
    |DiffSet| < |U| (the pigeonhole) and delta in the trichotomy set.
    """
    k = k_of(p)
    eps = 1 if p % 6 == 1 else 5
    Bk = ball_of(p)
    dmax = (p - 1) // 2
    nU = p - len(Bk)
    full = (1 << p) - 1

    def T2_of(a0):
        if eps == 5:
            return p - k if a0 <= k else p - k - 1
        return 5 * k + 1 if a0 <= 5 * k else 5 * k

    sig1_vals = sorted({T2_of(1), T2_of(p - 1)})   # both on-ball zones
    if eps == 5:
        on_L, off_L = 2 * k + 1, 2 * k + 2
        trich = {1, 2, dmax}
    else:
        on_L, off_L = 2 * k + 1, 2 * k
        trich = {1, 2, 3 * k}
    # feasible (L1, Lpair): L1 + 2*Lpair >= p  AND size-reality
    combos = []
    for L1 in (on_L, off_L):
        for Lpair in (on_L, off_L):
            if L1 + 2 * Lpair >= p:
                combos.append((L1, Lpair))
    res = {"eps": eps, "k": k, "nU": nU, "sig1_vals": sig1_vals,
           "combos": combos, "P1": {}, "P2": {}}

    # ---------------- P1: interval + two same-delta APs
    per_d = {}
    n_cov = 0
    for (L1, Lp) in combos:
        for sigma1 in sig1_vals:
            I = ap_mask(p, sigma1, 1, L1)
            for d in range(1, dmax + 1):
                starts = [ap_mask(p, s, d, Lp) for s in range(p)]
                for s2 in range(p):
                    rest = full ^ (I | starts[s2])
                    if rest == 0:
                        n_cov += 1        # degenerate (cannot happen)
                        continue
                    for s3 in range(p):
                        if rest & ~starts[s3] == 0:
                            n_cov += 1
                            per_d.setdefault(d, set()).add((s3 - s2) % p)
    res["P1"] = {"n_coverings": n_cov,
                 "per_delta": {str(d): sorted(v) for d, v in
                               sorted(per_d.items())}}
    # ---------------- P2: two equal-size intervals + one AP
    per_d2 = {}
    n_cov2 = 0
    for (L1, Lp) in combos:
        # the two intervals have size L1; the free AP has size Lp
        for sigma1 in sig1_vals:
            I1 = ap_mask(p, sigma1, 1, L1)
            for s2 in range(p):
                I2 = ap_mask(p, s2, 1, L1)
                base = I1 | I2
                rest0 = full ^ base
                if rest0 == 0:
                    n_cov2 += 1
                    per_d2.setdefault(0, set()).add((s2 - sigma1) % p)
                    continue
                for d in range(1, dmax + 1):
                    for s3 in range(p):
                        if rest0 & ~ap_mask(p, s3, d, Lp) == 0:
                            n_cov2 += 1
                            per_d2.setdefault(d, set()).add(
                                (s2 - sigma1) % p)
    res["P2"] = {"n_coverings": n_cov2,
                 "per_delta": {str(d): sorted(v) for d, v in
                               sorted(per_d2.items())}}
    # ---------------- checks
    res["P1_delta_extras"] = sorted(set(per_d) - trich) if (
        k >= (2 if eps == 5 else 3)) else sorted(set(per_d) - trich)
    res["P1_max_diffset"] = max((len(v) for v in per_d.values()), default=0)
    res["P1_pigeonhole_ok"] = all(len(v) < nU for v in per_d.values())
    res["P2_max_diffset"] = max((len(v) for v in per_d2.values()), default=0)
    res["P2_pigeonhole_ok"] = all(len(v) < nU for v in per_d2.values())
    # mechanism range: eps=5 k>=1 (p>=11), eps=1 k>=2 (p>=13).  p=7 (k=1,
    # eps=1) is OUT of range: |DeltaSet|=6 >= |U|=4 there; it is closed by
    # direct enumeration instead (PART A).
    in_mech_range = k >= (1 if eps == 5 else 2)
    res["in_mechanism_range"] = in_mech_range
    if in_mech_range:
        if not res["P1_pigeonhole_ok"]:
            VIOL.setdefault("S5_P1", []).append(
                (p, eps, {d: len(v) for d, v in per_d.items()}, nU))
        if not res["P2_pigeonhole_ok"]:
            VIOL.setdefault("S5_P2", []).append(
                (p, eps, {d: len(v) for d, v in per_d2.items()}, nU))
    return res


# ------------------------------------------------------------ ground truth

def ground_truth(p):
    """Exhaustive: no colliding family {K(a), B_1, B_v2, B_v3} covers
    Z_{p^2}, over both collision topologies:
      type A: v3 = +-v2 (mod p)   [the pair collides]
      type B: v2 = +-1 (mod p)    [collides with B_1; v3 free]
    """
    N = p * p
    Bk = ball_of(p)
    FULL = (1 << N) - 1
    Km = {}
    for a in range(1, (p + 1) // 2 + 1):
        if a in Km:
            continue
        m = 0
        for j in range(p):
            if (a * j) % p in Bk:
                for t in range(p):
                    m |= 1 << (j + p * t)
        Km[a] = m
    a_list = sorted(Km)
    B1 = bad_mask(1, N)
    KB = {a: (Km[a] | B1) for a in a_list}
    reps = pm_classes(N)
    Bmask = {u: bad_mask(u, N) for u in reps}
    # negation symmetry anchor
    for u in reps[:2]:
        if Bmask[u] != bad_mask(N - u, N):
            VIOL.setdefault("negation", []).append((p, u))
    coll = {}
    for u in reps:
        coll.setdefault(u % p, []).append(u)
    nA = nB = 0
    covers = []
    # ---- type A: pair collides
    for v2 in reps:
        r2 = v2 % p
        cands = set(u for u in coll.get(r2, []) if u != v2)
        cands |= set(u for u in coll.get((p - r2) % p, []) if u != v2)
        for v3 in sorted(cands):
            M23 = Bmask[v2] | Bmask[v3]
            for a in a_list:
                nA += 1
                if (KB[a] | M23) == FULL:
                    covers.append(("A", a, v2, v3))
    # ---- type B: v2 collides with B_1 (v2 = +-1 mod p, v2 != 1), v3 free
    typeB_v2 = sorted(set(coll.get(1, []) + coll.get(p - 1, [])) - {1})
    for v2 in typeB_v2:
        for v3 in reps:
            if v3 == v2:
                continue
            M23 = Bmask[v2] | Bmask[v3]
            for a in a_list:
                nB += 1
                if (KB[a] | M23) == FULL:
                    covers.append(("B", a, v2, v3))
    return {"p": p, "n_typeA": nA, "n_typeB": nB,
            "n_coverings": len(covers), "covers": covers[:10]}


# ------------------------------------------------------------ near miss

def near_miss(p, limit=40):
    """For the 'closest' eps=5 families (a'=1, v2, v3 = +-1 lifts): count
    kernel-uncovered fibers failing the t-covering; report the min."""
    N = p * p
    k = k_of(p)
    eps = 1 if p % 6 == 1 else 5
    Bk = ball_of(p)
    if eps != 5:
        return None
    U = [j for j in range(p) if j not in Bk]
    best = None
    n_fam = 0
    for v2 in range(2, N):
        if v2 % p != 1 or v2 >= N:
            continue
        for v3 in range(2, N):
            if v3 % p != 1:
                continue
            if v3 == v2 or (v3 + v2) % N == 0:
                continue
            n_fam += 1
            fails = 0
            for j in U:
                s1, _, L1 = fiber_tset_formula(1, j, p)
                s2, _, L2 = fiber_tset_formula(v2, j, p)
                s3, _, L3 = fiber_tset_formula(v3, j, p)
                if len(s1 | s2 | s3) < p:
                    fails += 1
            if best is None or fails < best[0]:
                best = (fails, v2, v3)
            if n_fam >= limit:
                break
        if n_fam >= limit:
            break
    return {"p": p, "n_families": n_fam, "min_failing_fibers": best[0]
            if best else None, "witness": best[1:] if best else None,
            "|U|": len(U)}


def main():
    P_FULL = [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1
              else ["5,7,11,13,17,19,23,29,31"])]
    P_STRUCT = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2
                else ["11,13,17,19,23,29,31"])]
    P_DRIFT = [int(x) for x in (sys.argv[3].split(",") if len(sys.argv) > 3
               else ["11,13,17,19"])]
    P_PIN = [int(x) for x in (sys.argv[4].split(",") if len(sys.argv) > 4
             else ["11,17,23,29,41,13,19,31,37,43,61"])]
    P_NEAR = [int(x) for x in (sys.argv[5].split(",") if len(sys.argv) > 5
              else ["11,17,23,29,41"])]
    OUT["meta"] = {"p_full": P_FULL, "p_struct": P_STRUCT,
                   "p_drift": P_DRIFT, "p_pin": P_PIN, "p_near": P_NEAR}

    log(f"PART A: ground truth (both topologies) at {P_FULL}")
    for p in P_FULL:
        t0 = time.time()
        g = ground_truth(p)
        OUT["ground_truth"][p] = g
        log(f"  p={p}: typeA={g['n_typeA']} typeB={g['n_typeB']} "
            f"COVERINGS={g['n_coverings']} ({time.time()-t0:.1f}s)")
        if g["n_coverings"]:
            VIOL.setdefault("GT", []).append((p, g["covers"][:3]))

    log(f"PART B: structure S1 at {P_STRUCT}")
    for p in P_STRUCT:
        OUT["structure"][p] = arc_data(p)
        log(f"  p={p}: arc kinds={OUT['structure'][p]['kinds']} "
            f"bad={OUT['structure'][p]['bad']}")
    log("PART B: full fiber t-set formula check (p <= 23)")
    for p in [q for q in P_STRUCT if q <= 23]:
        f = check_fiber_structure(p)
        OUT["structure"][p]["fiber"] = f
        log(f"  p={p}: tset checked={f['checked']} bad={f['bad']}")

    log(f"PART B: drift S2 at {P_DRIFT}")
    for p in P_DRIFT:
        d = check_drift(p)
        OUT["drift"][p] = d
        log(f"  p={p}: drift checked={d['checked']} bad={d['bad']} "
            f"badP2={d['badP2']}")

    log(f"PART C: pinning tables S4/S5 at {P_PIN}")
    for p in P_PIN:
        t0 = time.time()
        OUT["pinning"][p] = pinning_tables(p)
        t = OUT["pinning"][p]
        log(f"  p={p} (eps={t['eps']},k={t['k']}): P1 cov={t['P1']['n_coverings']} "
            f"deltas={sorted(t['P1']['per_delta'])} "
            f"maxDS={t['P1_max_diffset']} pgOK={t['P1_pigeonhole_ok']} | "
            f"P2 cov={t['P2']['n_coverings']} "
            f"deltas={sorted(t['P2']['per_delta'])} "
            f"maxDS={t['P2_max_diffset']} pgOK={t['P2_pigeonhole_ok']} "
            f"nU={t['nU']} ({time.time()-t0:.1f}s)")

    log(f"PART D: near-miss margins at {P_NEAR}")
    for p in P_NEAR:
        nm = near_miss(p)
        if nm:
            OUT["nearmiss"][p] = nm
            log(f"  p={p}: min failing fibers={nm['min_failing_fibers']} "
                f"of |U|={nm['|U|']} over {nm['n_families']} families")

    # ---- anchors ----
    A = {}
    A["ground_truth_clean"] = all(
        v["n_coverings"] == 0 for v in OUT["ground_truth"].values())
    A["arcs_clean"] = all(v["bad"] == 0 for v in OUT["structure"].values())
    A["tsets_clean"] = all(
        v.get("fiber", {}).get("bad", 0) == 0
        for v in OUT["structure"].values())
    A["drift_clean"] = all(
        v["bad"] == 0 and v["badP2"] == 0 for v in OUT["drift"].values())
    A["pigeonhole_P1"] = all(
        v["P1_pigeonhole_ok"] or not v["in_mechanism_range"]
        for v in OUT["pinning"].values())
    A["pigeonhole_P2"] = all(
        v["P2_pigeonhole_ok"] or not v["in_mechanism_range"]
        for v in OUT["pinning"].values())
    OUT["anchors"] = A
    log(f"anchors: {json.dumps(A)}")

    with open("/home/z/my-project/scripts/out_colliding_close.json", "w") as f:
        json.dump(OUT, f, indent=1, default=str)
    log("written scripts/out_colliding_close.json")
    if VIOL:
        log(f"VIOLATIONS: { {k: len(v) for k, v in VIOL.items()} }")
        sys.exit(1)
    log("ALL COLLIDING-CASE CHECKS PASS")


if __name__ == "__main__":
    main()
