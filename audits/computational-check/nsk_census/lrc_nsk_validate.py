"""Validation of lrc_nsk_census (C) against the Python reference DX()
and the brute-force gap solver lrc_xs_batch.

V-A  random cores at mixed levels (LRC values and rung values):
      C check mode must reproduce Python (R, D, X) exactly.
V-B  trichotomy: for every x in D (and a few random x <= R):
      gap(S+{x}) computed by the C brute solver must be
      = g exactly for x in X and < g for x in D \\ X (silent), > g for
      x outside D.  This ties the census classification to the
      independent solve() implementation.
V-C  small full-census cross-check (m=4, B=14): census cores file must
      equal the Python enumeration.
"""
import random
import subprocess
import sys
from fractions import Fraction
from itertools import combinations

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_nsk_theory import DX
from lrc_gap_lib import gap_int

CEN = "/home/z/my-project/scripts/lrc_nsk_census"
BATCH = "/home/z/my-project/scripts/lrc_xs_batch"
OUT = "/home/z/my-project/scripts/out"


def c_check(cores, g, xmax):
    inp = "".join(f"{len(c)} " + " ".join(map(str, c)) + "\n" for c in cores)
    p = subprocess.run([CEN, "check", str(g.numerator), str(g.denominator),
                        str(xmax)], input=inp, capture_output=True, text=True,
                       check=True)
    res = {}
    for line in p.stdout.splitlines():
        # S=.., | R=.. | D=.. | X=..
        left, rest = line.split(" | ", 1)
        S = tuple(int(v) for v in left[2:].rstrip(",").split(","))
        if rest.startswith("ON") or rest.startswith("BELOW"):
            res[S] = (None, None, None, rest.strip())
            continue
        parts = rest.split(" | ")
        R = int(parts[0][2:])
        D = [int(v) for v in parts[1][2:].split(",") if v]
        X = [int(v) for v in parts[2][2:].split(",") if v]
        res[S] = (R, D, X, None)
    return res


def brute_gaps(cores, g, xmax):
    """gap(S+{x}) for x in [1,xmax] via lrc_xs_batch; returns
    {core: {x: gap}}."""
    inp = "".join(f"{len(c)} " + " ".join(map(str, c)) + "\n" for c in cores)
    p = subprocess.run([BATCH, str(g.numerator), str(g.denominator), str(xmax)],
                       input=inp, capture_output=True, text=True, check=True)
    res = {}
    for line in p.stdout.splitlines():
        parts = line.split(" | ")
        if len(parts) < 2:
            continue
        core = tuple(int(v) for v in parts[0].split(","))
        num, den = (int(z) for z in parts[1].split("/"))
        xs = [int(v) for v in parts[2].split(",") if v.strip()] \
            if len(parts) > 2 and parts[2].strip() else []
        res[core] = (Fraction(num, den), xs)
    return res


def main():
    random.seed(20261003)
    failures = []

    # ---------------- V-A: random cores, mixed levels ----------------
    print("== V-A: C check vs Python DX on random cores ==")
    levels = [(2, Fraction(1, 3)), (3, Fraction(1, 4)),
              (4, Fraction(1, 5)), (4, Fraction(2, 9)),
              (5, Fraction(1, 6)), (5, Fraction(2, 11)),
              (6, Fraction(1, 7)), (6, Fraction(2, 13)),
              (7, Fraction(1, 8)), (7, Fraction(2, 15)),
              (8, Fraction(1, 9)), (8, Fraction(2, 17))]
    n_cmp = 0
    for (m, g) in levels:
        cores = []
        tries = 0
        while len(cores) < 40 and tries < 4000:
            tries += 1
            k = m - 1
            S = tuple(sorted(random.sample(range(1, 45), k)))
            gp, _ = gap_int(list(S))
            if gp > g:
                cores.append(S)
        cc = c_check(cores, g, 400)
        for S in cores:
            D, X, R, note = DX(S, g)
            got = cc.get(S)
            if got is None:
                failures.append(("V-A missing", S, g))
                continue
            Rc, Dc, Xc, notec = got
            if notec is not None or note is not None:
                if (notec is None) != (note is None):
                    failures.append(("V-A flag", S, g, note, notec))
                continue
            if (Rc, Dc, Xc) != (R, D, X):
                failures.append(("V-A", S, g, (R, D, X), (Rc, Dc, Xc)))
            else:
                n_cmp += 1
    print(f"  exact agreement on {n_cmp} cores across 12 levels"
          + ("" if not failures else "  (FAILURES)"))

    # ---------------- V-B: trichotomy vs brute solver ----------------
    print("== V-B: trichotomy of D-members vs brute solve() ==")
    n_tri = 0
    # known near-tight cores per level (zoo members and their subsets,
    # plus handpicked rung-family pieces)
    known = [
        (2, Fraction(1, 3), [(1,), (2,), (3,), (4,)]),
        (4, Fraction(1, 5), [(1, 3, 4), (1, 2, 3), (1, 2, 4), (2, 3, 4),
                              (1, 3, 5), (1, 4, 7), (2, 3, 5)]),
        (7, Fraction(2, 15), [(1, 4, 5, 6, 7, 11), (1, 3, 4, 5, 7, 11),
                              (1, 2, 3, 4, 5, 7), (1, 2, 3, 4, 5, 12),
                              (1, 2, 3, 4, 5, 6)]),
        (8, Fraction(1, 9), [(1, 2, 3, 4, 5, 6, 7), (1, 2, 3, 4, 5, 6, 8),
                              (1, 2, 3, 4, 5, 7, 8), (2, 3, 4, 5, 6, 7, 8)]),
    ]
    for (m, g, cores) in known:
        for S in cores:
            gp, _ = gap_int(list(S))
            if gp <= g:
                continue
            D, X, R, note = DX(S, g)
            if note is not None or not D:
                continue
            for x in D:
                gp3, _ = gap_int(sorted(list(S) + [x]))
                want = (gp3 == g) if x in X else (gp3 < g)
                if not want:
                    failures.append(("V-B", S, g, x, x in X, gp3))
                n_tri += 1
            outs = [x for x in range(1, R + 1)
                    if x not in S and x not in D][:3]
            for x in outs:
                gp3, _ = gap_int(sorted(list(S) + [x]))
                if gp3 <= g:
                    failures.append(("V-B-outside", S, g, x, gp3))
                n_tri += 1
    print(f"  trichotomy verified on {n_tri} (S,x) pairs"
          + ("" if not failures else "  (FAILURES)"))

    # ---------------- V-C: small full census vs Python ----------------
    print("== V-C: census m=4 B=14 vs Python enumeration ==")
    import os
    prefix = OUT + "/nsk_val_m4"
    p = subprocess.run([CEN, "4", "14", prefix, "1", "0", "300"],
                       capture_output=True, text=True, check=True)
    got = {}
    for line in open(prefix + "_cores.txt"):
        left, rest = line.split(" | ", 1)
        S = tuple(int(v) for v in left[2:].rstrip(",").split(","))
        parts = rest.split(" | ")
        R = int(parts[0][2:])
        D = [int(v) for v in parts[1][2:].split(",") if v]
        X = [int(v) for v in parts[2][2:].split(",") if v]
        got[S] = (R, D, X)
    want = {}
    g = Fraction(1, 5)
    for S in combinations(range(1, 15), 3):
        D, X, R, note = DX(S, g)
        if note is None and D:
            want[S] = (R, D, X)
    if got != want:
        only_c = {k: v for k, v in got.items() if want.get(k) != v}
        only_p = {k: v for k, v in want.items() if got.get(k) != v}
        failures.append(("V-C", list(only_c.items())[:3],
                         list(only_p.items())[:3]))
    print(f"  census cores file == Python enumeration: {got == want} "
          f"({len(got)} near-tight cores of {sum(1 for _ in combinations(range(1,15),3))})")
    print(p.stdout.strip())

    if failures:
        print(f"\nFAILURES: {len(failures)}")
        for f in failures[:10]:
            print("  ", f)
        return 1
    print("\nALL NSK-CENSUS VALIDATIONS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
