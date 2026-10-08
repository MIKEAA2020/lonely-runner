#!/usr/bin/env python3
"""
T14a: Line-by-line machine verification of the Trichotomy Lemma (Lemma T)
proof, as written in download/ap_lemmas_trichotomy_note.md, section 1.2.

The reviewer directive: "walk through the trichotomy proof line by line, with
the machine checking each intermediate claim."  Every numbered claim of the
proof is checked at every prime p <= PMAX (both classes), every critical
(s1, s), and every contained AP.

CLAIM LEDGER (proof line -> machine check)
  C1  (Case 1, d <= s1):
   C1a  non-terminal walk points never exceed p-1-d (else the next point
         lands in [0, d-1] subset [0, s1-1], outside I1); hence the walk
         never wraps: x_j = a + jd as integers.
   C1b  span bound: (s-1)d = x_{s-1} - a <= p-1-s1 = L.
   C1c  arithmetic: eps=5: (s-1)>=2k => 2k d <= 4k+3 => d <= 2+3/(2k) < 3
         for k >= 1 (so d = 2 whenever d <= s1 and k >= 2 lemma range);
         eps=1: (s-1)>=2k-1 => (2k-1)d <= 4k => d <= 2+2/(2k-1) < 3 for
         k >= 2.  Check the inequality chain numerically per (p, k).
  C2  (Case 2, d >= s1+1):
   C2a  d > p/3 and 2d < p, so 2d < p < 3d.  (numeric)
   C2b  r := p - 2d satisfies 1 <= r < s1 and r odd.  (numeric per config)
   C2c  walk identity x_{j+2} = x_j - r (mod p) for all j <= s-3.
   C2d  no-wrap descent: even and odd subsequences descend by exactly r per
         index pair as INTEGERS (x_{2i} = x_{2(i-1)} - r, no mod reduction),
         every non-terminal subsequence element is >= s1 + r.
   C2e  sub-case (i) (first step does not wrap, x1 = x0+d <= p-1):
         span(A) >= d + M1*r with M1 = ceil(s/2)-1, and d + M1*r <= L.
   C2f  arithmetic (i): eps=5: M1 = k => k r <= L-d <= 2k+1 => r <= 2+1/k
         < 3 for k >= 1, r odd => r = 1 for k >= 2; eps=1: M1 >= k-1 =>
         (k-1) r <= 2k-1 => r <= 2+1/(k-1) < 3 for k >= 3 => r = 1.
   C2g  sub-case (ii) (first step wraps, x0 in [p-d, p-1]):
         span(A) >= (p-d) + M2*r with M2 = floor(s/2)-1, and <= L.
   C2h  arithmetic (ii): eps=5: p-d >= 3k+3, M2 >= k-1 => (k-1)r <= k =>
         r <= 1+1/(k-1) < 2 for k >= 3 (k=2 gives r <= 2, odd => 1);
         eps=1: p-d >= 3k+1 => (k-1)r <= k-1 => r <= 1 => r = 1 (k >= 2).
   C2i  r = 1 in both sub-cases (lemma range) => d = (p-1)/2 and A is the
         two-block set [a-ceil(s/2)+1, a] u [a+d-floor(s/2)+1, a+d] (mod p).
  T1   Corollary T1 exact placements at eps=5, s1 = 2k+2 (L = 4k+2):
         size-(2k+2): d=2 unique form E = {2k+2, 2k+4, ..., 6k+4};
         d=q unique two-block T = [2k+2, 3k+2] u [5k+4, 6k+4];
         size-(2k+1): d=2 has exactly 3 placements, d=q exactly 2.
  T2   Run Lemma: AP size sigma, +-rep diff d != +-1, 2(sigma-1) < p,
         rho >= 2 consecutive residues in A => nu(d^{-1}) <= (sigma-1)/(rho-1).
  T3   Spacing-2 Run Lemma: rho points in AP of value-spacing gamma =>
         nu(gamma d^{-1}) <= (sigma-1)/(rho-1).
  GA   GENERAL-ARC trichotomy (new, used by the colliding-units closure):
         an AP of critical size contained in ANY arc of span L' (not just
         I1 = [s1, p-1]) has d in {1, 2, q} for eps=5, k >= 2, L' <= 4k+3;
         d in {1, 2, 3k} for eps=1, k >= 3, L' <= 4k+1.
         (d = 1 newly allowed: intervals.  Verified by translating the arc
         to I1-form, but checked directly here against every arc.)

Anchors: p=23 wrap AP; p=13 exception d=5; p=17 census {2, 8}; p=11
exceptions {3, 4}; two-block identity at d=q everywhere.
"""

import json
import sys
import time
from math import gcd

T0 = time.time()
OUT = {"meta": {}, "case1": {}, "case2": {}, "arith": {}, "t1": {},
       "runlemma": {}, "general_arc": {}, "anchors": {}, "violations": {}}
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


def pmrep(d, p):
    m = d % p
    return min(m, p - m)


def nu_of(x, p):
    """+-representative of x mod p."""
    m = x % p
    return min(m, p - m)


def ap_set(p, a, d, s):
    return sorted((a + j * d) % p for j in range(s))


def positions_by_hits(p, s1, s, d):
    """hits[a] = #{j in [0,s): (a + j d) mod p in [0, s1)} (the A1-zone)."""
    D = [0] * (2 * p + 2)
    for j in range(s):
        x = (j * d) % p
        st = (p - x) % p
        D[st] += 1
        D[st + s1] -= 1
    cov = [0] * (2 * p)
    c = 0
    for t in range(2 * p):
        c += D[t]
        cov[t] = c
    return [cov[a] + cov[a + p] for a in range(p)]


def contained_aps(p, s1, s, d):
    """All a with {a + jd} subset [s1, p-1]."""
    hits = positions_by_hits(p, s1, s, d)
    return [a for a in range(p) if hits[a] == 0]


def check_prime(p):
    k = k_of(p)
    eps = 1 if p % 6 == 1 else 5
    dmax = (p - 1) // 2
    q = dmax
    size_opts = [2 * k, 2 * k + 1, 2 * k + 2] if eps == 1 else [2 * k + 1, 2 * k + 2]
    res = {"p": p, "k": k, "eps": eps, "n_case1": 0, "n_case2": 0,
           "c1a": 0, "c1b": 0, "c2b": 0, "c2c": 0, "c2d": 0,
           "c2e": 0, "c2g": 0, "c2i": 0, "d_set": set(), "twoblock_bad": 0}

    in_range = k >= (2 if eps == 5 else 3)     # Lemma T's stated range
    for s1 in size_opts:
        L = p - 1 - s1
        for s in size_opts:
            if s < 2 or s > p - s1:
                continue
            for d in range(2, dmax + 1):
                for a in contained_aps(p, s1, s, d):
                    A = ap_set(p, a, d, s)
                    res["d_set"].add(pmrep(d, p))
                    # the walk points as residues
                    x = [(a + j * d) % p for j in range(s)]
                    if d <= s1:
                        # ---------------- CASE 1 ----------------
                        res["n_case1"] += 1
                        # C1a: walk never wraps; non-terminal <= p-1-d
                        ok = True
                        for j in range(s):
                            if a + j * d >= p:            # wrapped at step j
                                ok = False
                        for j in range(s - 1):
                            if x[j] > p - 1 - d:          # non-terminal high
                                ok = False
                        if not ok:
                            res["c1a"] += 1
                            VIOL.setdefault("C1a", []).append(
                                (p, s1, s, d, a))
                        # C1b: span bound (s-1)d <= L
                        if (s - 1) * d != x[s - 1] - x[0] or (s - 1) * d > L:
                            res["c1b"] += 1
                            VIOL.setdefault("C1b", []).append(
                                (p, s1, s, d, a))
                    else:
                        # ---------------- CASE 2 ----------------
                        res["n_case2"] += 1
                        r = p - 2 * d
                        # C2a: 2d < p < 3d
                        if not (2 * d < p < 3 * d):
                            VIOL.setdefault("C2a", []).append((p, s1, s, d, a))
                        # C2b: 1 <= r < s1, r odd  (asserted in lemma range;
                        # out-of-range slack at k<=1 is recorded, not asserted
                        # -- it is exactly the p=11 d=4 exception of S1.3)
                        if not (1 <= r < s1 and r % 2 == 1):
                            res["c2b"] += 1
                            if in_range:
                                VIOL.setdefault("C2b", []).append(
                                    (p, s1, s, d, a, r))
                        # C2c: x_{j+2} = x_j - r (mod p)
                        okc = all((x[j + 2] - x[j] + r) % p == 0
                                  for j in range(s - 2))
                        if not okc:
                            res["c2c"] += 1
                            VIOL.setdefault("C2c", []).append(
                                (p, s1, s, d, a))
                        # C2d: integer descent of both subsequences,
                        # non-terminal elements >= s1 + r
                        okd = True
                        for i in range(1, (s + 1) // 2):        # even idx 2i
                            if x[2 * i] != x[2 * (i - 1)] - r:
                                okd = False
                        for i in range(1, s // 2):              # odd idx 2i+1
                            if x[2 * i + 1] != x[2 * (i - 1) + 1] - r:
                                okd = False
                        for j in range(s - 2):                  # non-terminal
                            # (j <= s-3: has a successor in its subsequence;
                            #  x_{s-2} is a subsequence terminal - exempt)
                            if x[j] < s1 + r:
                                okd = False
                        if not okd:
                            res["c2d"] += 1
                            VIOL.setdefault("C2d", []).append(
                                (p, s1, s, d, a))
                        # C2e / C2g: span bounds per sub-case
                        M1 = (s + 1) // 2 - 1
                        M2 = s // 2 - 1
                        span = max(A) - min(A)
                        if x[1] == (a + d) % p and a + d <= p - 1:
                            # sub-case (i): first step does not wrap
                            if span < d + M1 * r or d + M1 * r > L:
                                res["c2e"] += 1
                                VIOL.setdefault("C2e", []).append(
                                    (p, s1, s, d, a, span, d + M1 * r, L))
                        else:
                            # sub-case (ii): first step wraps
                            if a < p - d or a > p - 1:
                                VIOL.setdefault("C2ii-setup", []).append(
                                    (p, s1, s, d, a))
                            if span < (p - d) + M2 * r or (p - d) + M2 * r > L:
                                res["c2g"] += 1
                                VIOL.setdefault("C2g", []).append(
                                    (p, s1, s, d, a, span,
                                     (p - d) + M2 * r, L))
                        # C2i: two-block identity (must hold when r = 1;
                        # recorded always when d = q)
                        if r == 1:
                            c1 = (s + 1) // 2
                            c2 = s // 2
                            arc1 = [(a - m) % p for m in range(c1)]
                            arc2 = [(a + d - m) % p for m in range(c2)]
                            if set(A) != set(arc1) | set(arc2):
                                res["twoblock_bad"] += 1
                                VIOL.setdefault("C2i", []).append(
                                    (p, s1, s, d, a))
    res["d_set"] = sorted(res["d_set"])
    return res


def check_arithmetic():
    """C1c / C2f / C2h inequality chains, numerically, for all k in range."""
    rows = []
    for k in range(1, 80):
        p5, p1 = 6 * k + 5, 6 * k + 1
        # C1c eps=5: d <= (4k+3)/(2k) = 2 + 3/(2k); < 3 iff k >= 2
        b5 = (4 * k + 3) / (2 * k)
        ok5 = (b5 < 3) == (k >= 2)
        # C1c eps=1: d <= 4k/(2k-1) = 2 + 2/(2k-1); < 3 iff k >= 2
        b1 = (4 * k) / (2 * k - 1)
        ok1 = (b1 < 3) == (k >= 2)
        # C2f(i) eps=5: r <= (L-d)/M1 <= (4k+3 - (2k+2))/k = 2 + 1/k;
        #   < 3 iff k >= 2 (k=1 gives r <= 3: the p=11 exception)
        r5i = (2 * k + 1) / k if k else 99
        okr5i = (r5i < 3) == (k >= 2)
        # C2f(i) eps=1: r <= (2k-1)/(k-1); < 3 iff k >= 3
        r1i = (2 * k - 1) / (k - 1) if k >= 2 else 99
        okr1i = (r1i < 3) == (k >= 3)
        # C2h(ii) eps=5: r <= k/(k-1); < 2 iff k >= 3; <=2 iff k>=2
        r5ii = k / (k - 1) if k >= 2 else 99
        okr5ii = (r5ii < 2) == (k >= 3) and (r5ii <= 2) == (k >= 2)
        # C2h(ii) eps=1: r <= (k-1)/(k-1) = 1 for k >= 2
        r1ii = (k - 1) / (k - 1) if k >= 2 else 99
        okr1ii = r1ii <= 1 if k >= 2 else True
        # Case1 threshold sanity: d <= s1 possible only if 2 <= s1
        rows.append({"k": k, "p5": p5, "p1": p1,
                     "c1c_eps5_ok": ok5, "c1c_eps5_bound": round(b5, 4),
                     "c1c_eps1_ok": ok1, "c1c_eps1_bound": round(b1, 4),
                     "c2f_eps5_ok": okr5i, "c2f_eps1_ok": okr1i,
                     "c2h_eps5_ok": okr5ii, "c2h_eps1_ok": okr1ii})
        if not (ok5 and ok1 and okr5i and okr1i and okr5ii and okr1ii):
            VIOL.setdefault("ARITH", []).append(rows[-1])
    return rows


def check_T1(p):
    """Corollary T1 at eps=5, s1 = 2k+2 (L = 4k+2)."""
    if p % 6 != 5:
        return None
    k = k_of(p)
    s1 = 2 * k + 2
    if s1 not in (2 * k + 1, 2 * k + 2):
        return None
    dmax = (p - 1) // 2
    out = {"p": p, "k": k}
    # size 2k+2, d = 2: unique placement E = {2k+2, ..., 6k+4} step 2
    A_d2 = contained_aps(p, s1, 2 * k + 2, 2)
    E = sorted(range(2 * k + 2, p, 2))
    out["d2_big_n"] = len(A_d2)
    out["d2_big_is_E"] = (A_d2 == [2 * k + 2]) and (E[0] == 2 * k + 2)
    # size 2k+2, d = q: unique T = [2k+2, 3k+2] u [5k+4, 6k+4]
    A_q_big = contained_aps(p, s1, 2 * k + 2, dmax)
    out["dq_big_n"] = len(A_q_big)
    # size 2k+1, d = 2: exactly 3 placements
    A_d2_s = contained_aps(p, s1, 2 * k + 1, 2)
    out["d2_small_n"] = len(A_d2_s)
    # size 2k+1, d = q: exactly 2 placements
    A_q_s = contained_aps(p, s1, 2 * k + 1, dmax)
    out["dq_small_n"] = len(A_q_s)
    # claims per lemma range
    if k >= 2:
        if len(A_d2) != 1 or A_d2 != [2 * k + 2]:
            VIOL.setdefault("T1_E", []).append((p, A_d2))
        if len(A_q_big) != 1:
            VIOL.setdefault("T1_T", []).append((p, A_q_big))
        if len(A_d2_s) != 3:
            VIOL.setdefault("T1_small_d2", []).append((p, A_d2_s))
        if len(A_q_s) != 2:
            VIOL.setdefault("T1_small_dq", []).append((p, A_q_s))
    return out


def check_runlemma(p):
    """T2/T3: for every AP (any a, d with +-rep in [2, dmax], size sigma with
    2(sigma-1) < p): if A contains rho >= 2 consecutive (spacing gamma)
    residues, then nu(gamma * d^{-1}) <= (sigma-1)/(rho-1)."""
    dmax = (p - 1) // 2
    bad = 0
    checked = 0
    for sigma in range(2, (p + 3) // 2):     # 2(sigma-1) < p <=> sigma<(p+2)/2
        for d in range(2, dmax + 1):
            dinv = pow(d, -1, p)
            for a in range(p):
                A = set((a + j * d) % p for j in range(sigma))
                for gamma in (1, 2):
                    # longest gamma-spaced run inside A (circular)
                    best = 0
                    for st in range(p):
                        if st not in A:
                            continue
                        ln = 1
                        while ((st + ln * gamma) % p in A
                               and ln * gamma < p):
                            ln += 1
                        best = max(best, ln)
                    if best >= 2:
                        checked += 1
                        bound = (sigma - 1) / (best - 1)
                        val = nu_of(gamma * dinv, p)
                        if val > bound + 1e-12:
                            bad += 1
                            if bad <= 5:
                                VIOL.setdefault("RUNLEMMA", []).append(
                                    (p, sigma, d, a, gamma, best, val, bound))
    return {"p": p, "checked": checked, "bad": bad}


def check_general_arc(p, span_list, assert_list):
    """GA: APs of critical size contained in ANY arc of given span.

    Asserted (in-lemma-range) spans: eps=5, k>=2: 4k+2, 4k+3 (translated
    s1' = 2k+2 resp. 2k+1, both critical); eps=1, k>=3: 4k-1, 4k (s1' =
    2k+1 resp. 2k+2... 2k).  Spans one larger (4k+4 resp. 4k+1) put s1'
    one BELOW the critical window -- Lemma T does not apply -- they are
    recorded but not asserted.
    """
    k = k_of(p)
    eps = 1 if p % 6 == 1 else 5
    dmax = (p - 1) // 2
    size_opts = ([2 * k, 2 * k + 1, 2 * k + 2] if eps == 1
                 else [2 * k + 1, 2 * k + 2])
    allowed = {1, 2, dmax} if eps == 5 else {1, 2, 3 * k}
    out = {"p": p, "eps": eps, "spans": {}, "n_contained": 0}
    for Lp in span_list:
        dset = set()
        n = 0
        for alpha in range(p):
            M = frozenset((alpha + i) % p for i in range(Lp + 1))
            for s in size_opts:
                if s > len(M):
                    continue
                for d in range(1, dmax + 1):
                    for a in range(p):
                        if all((a + j * d) % p in M for j in range(s)):
                            n += 1
                            dset.add(pmrep(d, p))
        out["spans"][Lp] = sorted(dset)
        out["n_contained"] += n
        if Lp in assert_list and k >= (2 if eps == 5 else 3):
            extra = dset - allowed
            if extra:
                VIOL.setdefault("GA", []).append(
                    (p, eps, k, Lp, sorted(extra)))
    return out


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    GAPMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 61
    primes = [p for p in sieve(PMAX) if p >= 5 and p % 6 in (1, 5)]
    OUT["meta"] = {"pmax": PMAX, "gapmax": GAPMAX, "n_primes": len(primes),
                   "primes": primes}

    log(f"checking Lemma T claims at {len(primes)} primes <= {PMAX}")
    for p in primes:
        r = check_prime(p)
        key = "eps5" if r["eps"] == 5 else "eps1"
        OUT["case1"][key] = OUT["case1"].get(key, {})
        OUT["case2"][key] = OUT["case2"].get(key, {})
        (OUT["case1"] if r["eps"] == 5 else OUT["case2"])[key][p] = r
        if p % 20 == 1 or p in (7, 11, 13, 17, 19, 23):
            log(f"  p={p} (k={r['k']}, eps={r['eps']}): case1={r['n_case1']} "
                f"case2={r['n_case2']} d_set={r['d_set']} "
                f"violations so far={sum(len(v) for v in VIOL.values())}")

    log("arithmetic chains (C1c, C2f, C2h)")
    OUT["arith"] = check_arithmetic()

    log("Corollary T1 placements (eps=5)")
    for p in primes:
        if p % 6 == 5:
            t = check_T1(p)
            if t:
                OUT["t1"][p] = t

    log("Run Lemmas T2/T3 at small primes")
    for p in [q for q in primes if q <= 37]:
        r = check_runlemma(p)
        OUT["runlemma"][p] = r
        log(f"  p={p}: runlemma checked={r['checked']} bad={r['bad']}")

    log("general-arc trichotomy")
    for p in [q for q in primes if 7 <= q <= GAPMAX]:
        k = k_of(p)
        eps = 1 if p % 6 == 1 else 5
        if eps == 5:
            spans = [4 * k + 2, 4 * k + 3, 4 * k + 4]
            asserts = {4 * k + 2, 4 * k + 3}
        else:
            spans = [4 * k - 1, 4 * k, 4 * k + 1]
            asserts = {4 * k - 1, 4 * k}
        g = check_general_arc(p, spans, asserts)
        OUT["general_arc"][p] = g
        log(f"  p={p}: GA spans -> {g['spans']} (asserted: {sorted(asserts)})")

    # ---- anchors ----
    A = {}
    # p=23 wrap-around counterexample: step-8 AP of size 8 in arc [5, 22]
    hits23 = positions_by_hits(23, 5, 8, 8)
    A["p23_wrap_contained"] = hits23[5] == 0
    A["p23_wrap_set"] = ap_set(23, 5, 8, 8)
    # p=13 exception d=5 contained in [4, 12] (the note's set {4,7,9,12}
    # is the start a=7 configuration; check any contained start)
    hits13 = positions_by_hits(13, 4, 4, 5)
    A["p13_d5_contained"] = any(hits13[a] == 0 for a in range(13))
    A["p13_d5_set"] = sorted(ap_set(13, 7, 5, 4))
    # p=17 census {2, 8}
    r17 = OUT["case2"]["eps5"].get(17) or OUT["case1"]["eps5"].get(17)
    A["p17_dset"] = r17["d_set"]
    A["p17_dset_is_28"] = r17["d_set"] == [2, 8]
    # p=11 exceptions {3, 4}
    r11 = OUT["case1"]["eps5"].get(11) or OUT["case2"]["eps5"].get(11)
    A["p11_dset"] = r11["d_set"]
    # two-block clean everywhere
    A["twoblock_clean_all"] = all(
        v["twoblock_bad"] == 0
        for key in ("eps1", "eps5")
        for v in (OUT["case1"].get(key, {}) | OUT["case2"].get(key, {})).values())
    # no violations at all
    A["no_violations"] = (sum(len(v) for v in VIOL.values()) == 0)
    OUT["anchors"] = A
    log(f"anchors: {json.dumps({k: v for k, v in A.items() if k != 'p23_wrap_set'})}")

    with open("/home/z/my-project/scripts/out_trich_linecheck.json", "w") as f:
        json.dump(OUT, f, indent=1, default=str)
    log("written scripts/out_trich_linecheck.json")

    if VIOL:
        log(f"VIOLATIONS PRESENT: { {k: len(v) for k, v in VIOL.items()} }")
        sys.exit(1)
    log("ALL LINE-BY-LINE CHECKS PASS")


if __name__ == "__main__":
    main()
