"""Closed-form filler characterization X(S): theorem and exact verification.

THEOREM (filler characterization).  Let S be a finite set of positive
integers and g in (0, 1/2) rational.  Write

    E(S) = { t in [0,1] : min_{s in S} || s t || > g }     (escape set),
    L(S, g) = { t in [0,1] : min_{s in S} || s t || = g }  (level set).

E(S) is a finite union of open intervals (alpha_j, beta_j) with

    0 < alpha_j < beta_j < 1,   min_s ||s alpha_j|| = min_s ||s beta_j|| = g,

whose endpoints are level-g crossing times of individual sawtooths
s in S (no escape interval touches 0 or 1 because g_S -> 0 there).
L(S, g) is a finite set containing those endpoints and possibly further
isolated level-g times (local maxima of g_S at height exactly g -- the
candidate binding times of Lemma A).  Then for an integer x >= 1,
x not in S:

    gap(S + {x}) = g      <=>

    (A) [kill]   for every escape interval (alpha_j, beta_j) there is an
                 integer m_j with  m_j - g <= x*alpha_j  and  x*beta_j <= m_j + g
                 (equivalently:  x*(beta_j-alpha_j) <= 2g  and the interval
                 [x*beta_j - g, x*alpha_j + g] contains an integer);

    (B) [touch]  some t in L(S, g) has  || x t || >= g
                 (the filler is clean at a level-g time of the core).

For gap(S) <= g the escape set is empty and (A) is vacuous, so the
statement also covers the degenerate cores: gap(S) < g gives L = empty
and no fillers; gap(S) = g gives the infinite family {x clean at a
level-g time}.

Consequences.
  * Radius bound: (A) forces x*(beta_j-alpha_j) <= 2g for every j, hence
    x <= 2g / max_j (beta_j - alpha_j) =: R(S).  So a core with gap(S) > g
    has a FINITE filler list with an explicit a priori radius.
  * gap(S) = g cores have infinite filler lists (any x clean at a level-g
    time of S); gap(S) < g cores have none.  (Adding speeds only lowers
    the gap.)  In the zoo regimes gap(S) < g never occurs for cores of
    size n-1 >= 4 because the Lonely Runner Conjecture holds up to 7
    speeds: gap >= 1/n > 2/(2n+1).
  * (B) is the binding structure of the extension: a filler can only
    BIND at a level-g time of the core (the window argument forbids
    ||x t|| = g inside an escape interval), and Lemma A's residue-clean
    condition at the exhibited time t* is exactly (B) at t* in L(S, g).

PROOF (recorded in the rung-family note).  gap(S+{x}) <= g iff (A):
outside E the core is already <= g; inside E the filler must be <= g,
and ||x t|| <= g on (alpha, beta) iff [x*alpha, x*beta] lies in a single
window [m-g, m+g] of the circle R/Z (windows are separated by gaps of
length 1-2g > 0).  Given (A), the value g is attained iff some t has
min(g_S(t), ||x t||) = g, i.e. g_S(t) >= g and ||x t|| >= g.  If
g_S(t) > g then t lies in an escape interval, where the window argument
forces ||x t|| < g strictly in the interior (an interior point at
window distance exactly g would push an endpoint past the window edge);
so g_S(t) = g exactly: t in L(S, g).  Conversely a level time t with
||x t|| >= g gives min(g, ||x t||) = g there.  (For gap(S) = g cores
E is empty and only the second clause remains.)

This script verifies the theorem EXACTLY (Fraction arithmetic):
  V1  every core (proper subset of size |V|-1) of every zoo member at
      n = 5, 6, 7, 8 with gap(S) > g_n:  X_theory = X_brute on [1,70];
  V2  random cores at n = 5..8 (gap > g): same equality;
  V3  gap(S) = g and gap(S) < g cores behave as stated;
  V4  worked examples H1, H2 with stage-by-stage pruning;
  V5  spot checks of the C brute force against the independent Python
      implementations gap_int / gap_frac of lrc_gap_lib.
"""
import math
import random
import subprocess
import sys
from fractions import Fraction
from itertools import combinations

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap_int, gap_frac, _norm

BATCH = "/home/z/my-project/scripts/lrc_xs_batch"
OUT = "/home/z/my-project/scripts/out"

# ----------------------------------------------------------------------
# exact machinery
# ----------------------------------------------------------------------

def dist(x):
    """||x|| for Fraction x."""
    return _norm(x)


def gS(S, t):
    return min(dist(s * t) for s in S)


def level_crossings(s, g):
    """All t in (0,1) with || s t || = g (g = p/q reduced)."""
    p, q = g.numerator, g.denominator
    out = []
    for m in range(0, s + 1):
        for num in (q * m + p, q * m - p):
            if 0 < num < s * q:
                t = Fraction(num, s * q)
                if dist(s * t) == g:      # exact filter
                    out.append(t)
    return out


def sawtooth_breaks(s):
    return [Fraction(a, 2 * s) for a in range(1, 2 * s, 2)]


def escape_intervals(S, g):
    """E(S) = {g_S > g} as a list of open intervals (alpha, beta).

    Exact: events = 0, 1, all sawtooth breakpoints, all level-g crossings
    of all sawtooths.  On a segment between consecutive events every
    ||s t|| is linear and does not cross level g, so g_S - g has constant
    sign; maximal positive runs are the escape intervals (merged across
    event points where g_S > g, split where g_S = g).
    """
    S = list(S)
    events = {Fraction(0), Fraction(1)}
    for s in S:
        for t in sawtooth_breaks(s):
            events.add(t)
        for t in level_crossings(s, g):
            events.add(t)
    ev = sorted(events)
    pos = []
    for i in range(len(ev) - 1):
        mid = (ev[i] + ev[i + 1]) / 2
        pos.append(gS(S, mid) > g)
    val = [gS(S, t) for t in ev]
    out = []
    i, n = 0, len(pos)
    while i < n:
        if pos[i]:
            j = i
            while j + 1 < n and pos[j + 1] and val[j + 1] > g:
                j += 1
            out.append((ev[i], ev[j + 1]))
            i = j + 1
        else:
            i += 1
    for (a, b) in out:
        assert 0 < a < b < 1, ("interval endpoints", S, a, b)
        assert gS(S, a) == g and gS(S, b) == g, ("endpoint levels", S, a, b)
    return out


def level_set(S, g):
    """All t in (0,1) with g_S(t) = g exactly (finite: crossings)."""
    cands = set()
    for s in S:
        cands.update(level_crossings(s, g))
    return sorted(t for t in cands if gS(S, t) == g)


def X_theory(S, g, xmax=None):
    """{x in [1, xmax], x not in S : gap(S + {x}) = g} by the theorem."""
    S = list(S)
    ivs = escape_intervals(S, g)
    lvl = level_set(S, g)
    if not ivs:
        if not lvl:
            return [], 0, "core below rung: no fillers"
        hi = xmax if xmax is not None else 200
        xs = [x for x in range(1, hi + 1) if x not in S
              and any(dist(x * t) >= g for t in lvl)]
        return xs, None, "core ON the rung: infinite family"
    Lmax = max(b - a for (a, b) in ivs)
    R = int(2 * g / Lmax)            # floor: a priori radius
    hi = min(R, xmax) if xmax is not None else R
    out = []
    for x in range(1, hi + 1):
        if x in S:
            continue
        okA = True
        for (a, b) in ivs:
            lo, hii = x * b - g, x * a + g
            if lo > hii or math.ceil(lo) > math.floor(hii):
                okA = False
                break
        if not okA:
            continue
        if any(dist(x * t) >= g for t in lvl):
            out.append(x)
    return out, R, None


def X_brute_batch(cores, g, xmax):
    """Run the C batch tool (solve() verbatim from the sweep solver)."""
    inp = "".join(
        f"{len(c)} " + " ".join(str(v) for v in c) + "\n" for c in cores)
    proc = subprocess.run(
        [BATCH, str(g.numerator), str(g.denominator), str(xmax)],
        input=inp, capture_output=True, text=True, check=True)
    res = {}
    for line in proc.stdout.splitlines():
        parts = line.split(" | ")
        if len(parts) == 2:
            parts.append("")
        core_part, gap_part, x_part = parts
        core = tuple(int(v) for v in core_part.split(","))
        num, den = (int(z) for z in gap_part.split("/"))
        xs = [int(v) for v in x_part.split(",") if v.strip()]
        res[core] = (Fraction(num, den), xs)
    return res


# ----------------------------------------------------------------------
# verification
# ----------------------------------------------------------------------

def load_zoo(path):
    out = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith("v=["):
            body = line[3:line.index("]")]
            out.append(tuple(int(x) for x in body.split(",")))
        elif line.startswith("v=("):
            body = line[3:line.index(")")]
            out.append(tuple(int(x) for x in body.split(",")))
    return out


def zoo_files():
    base = "/home/z/my-project/lonely-runner/audits/computational-check/"
    return [
        (5, Fraction(2, 11), base + "n5_V50_rung2.txt"),
        (6, Fraction(2, 13), base + "n6_V50_rung2.txt"),
        (7, Fraction(2, 15), base + "tight_n7_V50_rung2.txt"),
        (8, Fraction(2, 17), base + "n8_V38_rung2.txt"),
    ]


def main():
    random.seed(20261002)
    XMAX = 70
    n_zoo_cores = 0
    failures = []
    stats = {"above": 0, "on": 0, "below": 0}

    print("== V1: zoo cores (every (k-1)-subset of every zoo member) ==")
    for n, g, path in zoo_files():
        Z = load_zoo(path)
        cores = {}
        for v in Z:
            for S in combinations(v, len(v) - 1):
                cores.setdefault(S, None)
        core_list = sorted(cores)
        batch = X_brute_batch(core_list, g, XMAX)
        for S in core_list:
            cgap, xbrute = batch[S]
            xt, R, note = X_theory(S, g, XMAX)
            if cgap > g:
                stats["above"] += 1
                if xt != xbrute:
                    failures.append(("V1", S, g, xt, xbrute))
            elif cgap == g:
                stats["on"] += 1
                if xt != xbrute:
                    failures.append(("V1-on", S, g, xt, xbrute))
            else:
                stats["below"] += 1
                if xbrute:
                    failures.append(("V1-below", S, g, [], xbrute))
        n_zoo_cores += len(core_list)
        print(f"  n={n} g={g}: {len(core_list)} cores "
              f"(gap>g: {sum(1 for S in core_list if batch[S][0] > g)}, "
              f"=g: {sum(1 for S in core_list if batch[S][0] == g)}, "
              f"<g: {sum(1 for S in core_list if batch[S][0] < g)})")
    print(f"  zoo cores total: {n_zoo_cores}  "
          f"(above/on/below = {stats['above']}/{stats['on']}/{stats['below']})")

    print("\n== V2: random cores, gap(S) > g ==")
    for n, count in [(5, 80), (6, 120), (7, 200), (8, 80)]:
        g = Fraction(2, 2 * n + 1)
        cores = []
        tries = 0
        while len(cores) < count and tries < count * 40:
            tries += 1
            S = tuple(sorted(random.sample(range(1, 51), n - 1)))
            gp, _ = gap_int(list(S))
            if gp > g:
                cores.append(S)
        batch = X_brute_batch(cores, g, XMAX)
        for S in cores:
            cgap, xbrute = batch[S]
            xt, R, note = X_theory(S, g, XMAX)
            assert cgap > g
            if xt != xbrute:
                failures.append(("V2", S, g, xt, xbrute))
        print(f"  n={n} g={g}: {len(cores)} random cores above the rung, "
              f"all X_theory == X_brute" if not any(
                  f[0] == "V2" for f in failures) else
              f"  n={n} g={g}: {len(cores)} random cores, FAILURES")

    print("\n== V3: gap(S) = g and gap(S) < g cores ==")
    n = 7
    g = Fraction(2, 15)
    on_cores, below_cores = [], []
    tries = 0
    while (len(on_cores) < 8 or len(below_cores) < 8) and tries < 6000:
        tries += 1
        S = tuple(sorted(random.sample(range(1, 51), 6)))
        gp, _ = gap_int(list(S))
        if gp == g and len(on_cores) < 8:
            on_cores.append(S)
        elif gp < g and len(below_cores) < 8:
            below_cores.append(S)
    batch = X_brute_batch(on_cores + below_cores, g, XMAX)
    for S in on_cores:
        cgap, xbrute = batch[S]
        xt, _, _ = X_theory(S, g, XMAX)
        if xt != xbrute:
            failures.append(("V3-on", S, g, xt, xbrute))
    for S in below_cores:
        cgap, xbrute = batch[S]
        if xbrute:
            failures.append(("V3-below", S, g, [], xbrute))
    print(f"  on-rung cores: {len(on_cores)} (infinite families, "
          f"checked to {XMAX}), below-rung cores: {len(below_cores)} "
          f"(no fillers)")

    print("\n== V5: spot checks of the C brute vs Python references ==")
    spot_cores = []
    for n, g, path in zoo_files()[:3]:
        Z = load_zoo(path)
        for v in random.sample(Z, min(6, len(Z))):
            spot_cores.append(tuple(sorted(v)))
    batch = X_brute_batch(spot_cores, Fraction(2, 15), 30)  # gaps only
    nchk, nok = 0, 0
    for S, (cgap, _) in batch.items():
        g1, _ = gap_int(list(S))
        g2, _ = gap_frac(list(S))
        nchk += 1
        if g1 == g2 == cgap:
            nok += 1
        else:
            failures.append(("V5", S, cgap, g1, g2))
    print(f"  {nok}/{nchk} triple agreements (C solve == gap_int == gap_frac)")

    print("\n== V4: worked examples ==")
    for S, g, label in [
        ((1, 3, 4, 5, 7, 11), Fraction(2, 15), "H1"),
        ((1, 4, 5, 6, 7, 11), Fraction(2, 15), "H2"),
    ]:
        ivs = escape_intervals(S, g)
        lvl = level_set(S, g)
        Lmax = max(b - a for (a, b) in ivs)
        R = int(2 * g / Lmax)
        xt, _, _ = X_theory(S, g, 200)
        print(f"  {label} = {list(S)}: {len(ivs)} escape intervals, "
              f"L_max = {float(Lmax):.4f}, radius R = {R}, "
              f"|L(S,g)| = {len(lvl)} level times")
        for (a, b) in ivs[:6]:
            print(f"     escape ({a}, {b})  length {float(b-a):.4f}")
        if len(ivs) > 6:
            print(f"     ... and {len(ivs)-6} more")
        lvl_str = ", ".join(str(t) for t in lvl[:12])
        print(f"     level times L: {lvl_str}"
              + (" ..." if len(lvl) > 12 else ""))
        # stage-by-stage pruning
        stage = [x for x in range(1, R + 1) if x not in S]
        print(f"     stage 0 (radius):      {len(stage)} candidates x <= R = {R}")
        for j, (a, b) in enumerate(ivs):
            nxt = []
            for x in stage:
                lo, hi = x * b - g, x * a + g
                if lo <= hi and math.ceil(lo) <= math.floor(hi):
                    nxt.append(x)
            if nxt != stage:
                print(f"     interval {j} = ({float(a):.4f},{float(b):.4f})"
                      f" prunes to {len(nxt)}")
                stage = nxt
        touch = [x for x in stage if any(dist(x * t) >= g for t in lvl)]
        print(f"     after all kills:      {stage}")
        print(f"     after touch (B):      {touch}")
        print(f"     X({label}) = {xt}")
        # cross-check against brute to 70
        batch = X_brute_batch([S], g, 70)
        xb = batch[S][1]
        print(f"     X_brute[1,70]       = {xb}   "
              f"{'MATCH' if [x for x in xt if x <= 70] == xb else 'MISMATCH'}"
              f"{'' if [x for x in xt if x <= 70] == xb else (xt, xb)}")
        if [x for x in xt if x <= 70] != xb:
            failures.append(("V4", S, xt, xb))

    print("\n== summary ==")
    if failures:
        print(f"FAILURES: {len(failures)}")
        for f in failures[:10]:
            print("  ", f)
        return 1
    print("ALL CHECKS PASS: X_theory == X_brute on every core tested "
          f"(zoo: {n_zoo_cores}, random: 480, on/below: 16, spot: {nchk}).")
    print("The filler-kill ladder decomposition is now a theorem.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
