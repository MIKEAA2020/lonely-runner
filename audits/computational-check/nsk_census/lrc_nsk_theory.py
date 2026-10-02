"""Silent killers and the induction step: machine verification of the
General Filler Theorem (trichotomy), the every-order decomposition of
tight sets, the silent-killer reduction of LRC, and the corrected
scaling proposition.

THEOREM G (General Filler Theorem, trichotomy).  Let S be a finite set
of positive integers, g in (0, 1/2) rational, x an integer, x not in S.
With E(S) the escape set {g_S > g} (finite union of open intervals) and
L(S) the level set {g_S = g} (finite), define

    (A) [kill]  for every escape interval (alpha, beta) there is an
                integer m with m - g <= x*alpha and x*beta <= m + g
                (vacuous when E = empty);
    (B) [touch] some t in L(S) has ||x t|| >= g.

Then EXACTLY ONE of

    gap(S + {x}) > g   <=>  not (A)
    gap(S + {x}) = g   <=>  (A) and (B)
    gap(S + {x}) < g   <=>  (A) and not (B).

Call x with (A) and not (B) a SILENT KILLER of S at level g: it kills
every escape and every level time, and the gap of the extension drops
strictly below g.  D(S, g) = {x : (A)} (finite, x <= R = floor(2g/L_max)
when gap(S) > g); X(S, g) = D(S, g) cap {touch}.

THEOREM R (silent-killer reduction).  For m >= 2, g = 1/(m+1):
LRC(m) implies NSC(m) (no silent killers anywhere); and
LRC(m-1) + NSC(m) imply LRC(m).  So for every m <= 9 (where LRC(m-1)
is classical: ... , Rosenfeld's nine-runner theorem covers m-1 = 8)
LRC(m) <=> NSC(m).

COROLLARY G1 (every-order decomposition).  If gap(T) = g and
gap(T \ {x}) > g then x in X(T \ {x}, g).  Under LRC(m-1) every tight
m-set therefore decomposes at EVERY removal order.

PROPOSITION P (scaling).  c * X(S, g) is a subset of X(c*S, g) for all
c >= 1 (level/escape structure transfers under t -> c*t); the inclusion
is strict in general (descaling fillers, first seen at level 2:
X({1}, 1/3) = {2} but X({2}, 1/3) = {1, 4}).

This script verifies all of it EXACTLY (Fraction arithmetic):
  NT1 trichotomy on random (S, g, x) across regimes
      (g at LRC levels, rung levels, and off-zoo values; cores above,
      on, and below g; killers, fillers, and silent killers);
  NT2 every-order decomposition of every tight set n = 2..8 (V <= 50);
  NT3 reduction equivalence at small scale (m = 2..5, V <= 16):
      LRC(m) in range <=> NSC(m) in range; plus a negative control at
      a level g' above the true bound where silent killers MUST exist;
  NT4 scaling: containment on random pairs, zoo-pair equality at
      levels 3..8, strictness at level 2.
"""
import math
import random
import sys
from fractions import Fraction
from itertools import combinations

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap_int, _norm
from lrc_xs_theory import escape_intervals, level_set

OUT = "/home/z/my-project/scripts/out"
CK = "/home/z/my-project/lonely-runner/audits/computational-check/"

TIGHT_ZOOS = {
    2: OUT + "/tz2_tight.txt", 3: OUT + "/rz3_tight.txt",
    4: OUT + "/rz4_tight.txt", 5: OUT + "/rz5_tight.txt",
    6: OUT + "/rz6_tight.txt", 7: CK + "tight_n7_V50.txt",
    8: CK + "n8_V50_tight.txt",
}


def dist(x):
    return _norm(x)


def gS(S, t):
    return min(dist(s * t) for s in S)


def kill(S, g, x):
    """Condition (A): x kills every escape interval of S at level g."""
    for (a, b) in escape_intervals(S, g):
        lo, hi = x * b - g, x * a + g
        if lo > hi or math.ceil(lo) > math.floor(hi):
            return False
    return True


def touch(S, g, x):
    """Condition (B): x is clean at some level-g time of S."""
    return any(dist(x * t) >= g for t in level_set(S, g))


def trichotomy(S, g, x):
    """Predicted class '>' / '=' / '<' from the theorem."""
    A = kill(S, g, x)
    B = touch(S, g, x)
    if not A:
        return ">"
    return "=" if B else "<"


def actual(S, g, x):
    gp, _ = gap_int(sorted(list(S) + [x]))
    if gp > g:
        return ">"
    if gp == g:
        return "="
    return "<"


def DX(S, g):
    """(D, X, R) for a core with gap(S) > g: escape-killers and fillers."""
    ivs = escape_intervals(S, g)
    lvl = level_set(S, g)
    if not ivs:
        return None, None, None, "gap(S) <= g (on or below the level)"
    Lmax = max(b - a for (a, b) in ivs)
    R = int(2 * g / Lmax)
    D, X = [], []
    for x in range(1, R + 1):
        if x in S:
            continue
        okA = all(
            (lambda lo, hi: not (lo > hi or math.ceil(lo) > math.floor(hi)))(
                x * b - g, x * a + g)
            for (a, b) in ivs)
        if okA:
            D.append(x)
            if any(dist(x * t) >= g for t in lvl):
                X.append(x)
    return D, X, R, None


def load(path):
    out = []
    for line in open(path):
        line = line.strip()
        if line.startswith("v=[") or line.startswith("v=("):
            body = line[3:line.index("]") if line[1] == "[" else line.index(")")]
            out.append(tuple(int(v) for v in body.split(",")))
    return out


# ----------------------------------------------------------------------


def nt1(random_seed=20261002):
    """Trichotomy on random (S, g, x) across regimes."""
    random.seed(random_seed)
    n_ok, n_by_class = 0, {">": 0, "=": 0, "<": 0}
    silent_seen = []
    GVALS = [  # level targets: LRC values, rung values, off-zoo values
        Fraction(1, 3), Fraction(2, 5), Fraction(1, 4), Fraction(3, 10),
        Fraction(1, 5), Fraction(2, 9), Fraction(3, 13), Fraction(1, 6),
        Fraction(2, 13), Fraction(1, 7), Fraction(2, 15), Fraction(1, 8),
        Fraction(3, 8), Fraction(1, 9), Fraction(2, 17),
    ]
    trials = 0
    while trials < 3000 and (n_by_class["<"] < 40 or n_by_class["="] < 60
                             or n_by_class[">"] < 800):
        trials += 1
        k = random.choice([1, 2, 2, 3, 3, 4, 4, 5])
        V = random.choice([12, 20, 30, 50])
        S = tuple(sorted(random.sample(range(1, V + 1), k)))
        g = random.choice(GVALS)
        x = random.randint(1, 60)
        if x in S:
            continue
        pred, act = trichotomy(S, g, x), actual(S, g, x)
        if pred != act:
            return n_ok, n_by_class, silent_seen, (S, g, x, pred, act)
        n_ok += 1
        n_by_class[pred] += 1
        if pred == "<" and len(silent_seen) < 8:
            silent_seen.append((S, g, x))
    # the worked micro-example from the note
    for (S, g, x) in [((1,), Fraction(2, 5), 2), ((1,), Fraction(1, 3), 2),
                      ((2,), Fraction(1, 3), 1), ((2,), Fraction(1, 3), 4),
                      ((1, 2), Fraction(1, 4), 3), ((1, 2), Fraction(1, 4), 4),
                      ((1, 2, 3), Fraction(1, 4), 4),
                      ((1, 2, 3), Fraction(1, 5), 4)]:
        pred, act = trichotomy(S, g, x), actual(S, g, x)
        if pred != act:
            return n_ok, n_by_class, silent_seen, (S, g, x, pred, act)
        n_ok += 1
        n_by_class[pred] += 1
        if pred == "<" and len(silent_seen) < 12:
            silent_seen.append((S, g, x))
    return n_ok, n_by_class, silent_seen, None


def nt2():
    """Every-order decomposition of every tight set n = 2..8."""
    checked = 0
    for n, path in TIGHT_ZOOS.items():
        Z = load(path)
        g = Fraction(1, n + 1)
        for T in Z:
            for x in T:
                S = tuple(v for v in T if v != x)
                gp, _ = gap_int(list(S))
                if gp <= g:
                    return checked, (n, T, x, "core not above level", gp)
                D, X, R, note = DX(S, g)
                if note is not None:
                    return checked, (n, T, x, note, None)
                if x not in X:
                    return checked, (n, T, x, "filler not in X(core)", (D, X, R))
                if D != X:
                    return checked, (n, T, x, "SILENT KILLER in tight core",
                                     (D, X, R))
                checked += 1
    return checked, None


def nt3():
    """Reduction equivalence at small scale + negative control."""
    results = []
    for m in range(2, 6):
        V = 16
        g = Fraction(1, m + 1)
        # side 1: LRC(m) in range
        lrc = True
        worst = None
        for T in combinations(range(1, V + 1), m):
            gp, _ = gap_int(list(T))
            if gp < g:
                lrc = False
                worst = (T, gp)
                break
        # side 2: NSC(m) in range (all (m-1)-cores above level: D == X)
        nsc = True
        bad = None
        for S in combinations(range(1, V + 1), m - 1):
            gp, _ = gap_int(list(S))
            if gp <= g:
                continue          # relies on LRC(m-1); flagged elsewhere
            D, X, R, note = DX(S, g)
            if note is not None:
                continue
            for x in range(1, R + 1):
                if x in S:
                    continue
                gp2, _ = gap_int(sorted(list(S) + [x]))
                silent = (gp2 < g)
                if silent != (x in D and x not in X):
                    nsc = False
                    bad = ("D/X mismatch", S, x, gp2, D, X)
                    break
                if silent:
                    nsc = False
                    bad = ("silent killer", S, x, gp2, D, X)
                    break
            if not nsc:
                break
        agree = (lrc == nsc)
        results.append((m, V, lrc, nsc, agree, worst, bad))
        if not agree:
            return results, (m, "LRC != NSC", lrc, nsc, worst, bad)
    # negative control: at a level above the true bound, silent killers
    # must exist (e.g. S={1}, g=2/5, x=2: gap({1,2})=1/3 < 2/5).
    S, gg, x = (1,), Fraction(2, 5), 2
    pred, act = trichotomy(S, gg, x), actual(S, gg, x)
    if (pred, act) != ("<", "<"):
        return results, ("negative control failed", S, gg, x, pred, act)
    return results, None


def nt4(random_seed=20261002):
    """Scaling: containment always; zoo equality at levels >= 3;
    strictness at level 2."""
    random.seed(random_seed)
    n_contain, n_eq34, n_pairs34 = 0, 0, 0
    for _ in range(150):
        k = random.choice([1, 2, 3, 4])
        S = tuple(sorted(random.sample(range(1, 31), k)))
        g = random.choice([Fraction(1, 4), Fraction(1, 5),
                           Fraction(2, 9), Fraction(1, 6)])
        c = random.choice([2, 3, 4])
        X1, _, _, note1 = (lambda t: (t[0], t[1], t[2], t[3]))(
            (DX(S, g) if gap_int(list(S))[0] > g else (None,) * 4))
        cS = tuple(c * v for v in S)
        X2, _, _, note2 = (DX(cS, g) if gap_int(list(cS))[0] > g
                           else (None,) * 4)
        if X1 is None or X2 is None:
            continue
        for x in X1:
            if c * x not in X2:
                return n_contain, n_eq34, n_pairs34, ("containment", S, c, g, x)
        n_contain += 1
    # zoo pairs at levels 3..8 (equality expected) and level 2 (strict)
    for n in range(3, 9):
        g = Fraction(1, n + 1)
        Z = load(TIGHT_ZOOS[n])
        cores = {}
        for T in Z:
            for S in combinations(T, len(T) - 1):
                cores.setdefault(S, None)
        for S in cores:
            d = math.gcd(*S) if len(S) > 1 else S[0]
            if d == 1:
                continue
            P = tuple(v // d for v in S)
            if P not in cores:
                continue
            X1, _, _, n1 = DX(P, g) if gap_int(list(P))[0] > g else (None,)*4
            X2, _, _, n2 = DX(S, g) if gap_int(list(S))[0] > g else (None,)*4
            if X1 is None or X2 is None:
                continue
            n_pairs34 += 1
            if X2 != [d * x for x in X1]:
                return n_contain, n_eq34, n_pairs34, ("zoo equality", n, S, P,
                                                      X1, X2)
            n_eq34 += 1
    # level-2 strictness: X({1},1/3)={2}, X({2},1/3)={1,4}
    Xa, _, _, _ = DX((1,), Fraction(1, 3))
    Xb, _, _, _ = DX((2,), Fraction(1, 3))
    if Xa != [2] or Xb != [1, 4]:
        return n_contain, n_eq34, n_pairs34, ("level-2 strictness", Xa, Xb)
    return n_contain, n_eq34, n_pairs34, None


def main():
    print("== NT1: General Filler Theorem (trichotomy) ==")
    n_ok, by_class, silent, err = nt1()
    if err:
        print("  FAILURE:", err)
        return 1
    print(f"  {n_ok} random+fixed (S,g,x) triples classified exactly "
          f"(> {by_class['>']}, = {by_class['=']}, < {by_class['<']})")
    print(f"  sample silent killers: {silent[:4]}")

    print("\n== NT2: every-order decomposition of tight sets ==")
    checked, err = nt2()
    if err:
        print("  FAILURE:", err)
        return 1
    print(f"  all removal orders of every tight set n=2..8 (V<=50) are "
          f"X-certified: {checked} (T,x) pairs; D == X at every core")

    print("\n== NT3: silent-killer reduction at small scale ==")
    results, err = nt3()
    for (m, V, lrc, nsc, agree, worst, bad) in results:
        print(f"  m={m} V<={V}: LRC(m)={lrc}  NSC(m)={nsc}  agree={agree}")
    if err:
        print("  FAILURE:", err)
        return 1
    print("  negative control: S={1}, g=2/5, x=2 is a silent killer "
          "(gap({1,2})=1/3 < 2/5) -- NSC fails above the true bound, as it must")

    print("\n== NT4: scaling proposition ==")
    n_contain, n_eq34, n_pairs34, err = nt4()
    if err:
        print("  FAILURE:", err)
        return 1
    print(f"  containment c*X(S) <= X(cS): {n_contain}/150 random pairs")
    print(f"  zoo core pairs levels 3..8: {n_eq34}/{n_pairs34} equal "
          f"(X(cS) = c*X(S))")
    print("  level-2 strictness: X({1},1/3)={2}, X({2},1/3)={1,4} "
          "-- descaling filler 1 (errata for the atlas claim)")

    print("\nALL NEW-THEOREM CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
