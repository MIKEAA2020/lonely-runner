"""Analysis of the n=8, V=50 hunt (536,878,650 vectors, v4 filtered solver).

  H1  accounting: sum(spectrum) + sum(filtered) == C(50,8)
  H2  sample validation: every sampled record re-verified against the
      independent Python reference (gap, argmax, feasibility)
  H3  tight census: closed under scalars of the consecutive family?
  H4  rung-2 census (21 vectors): primitive classes, binding pairs,
      generalized Lemma A (grids d = 17e), new members vs V=38
  H5  rung-k zoos k = 3..6: family members {1..7, 8k} present, Lemma A
      on grids d = (8k+1)e
  H6  inter-rung window (1/9, 2/17): EMPTY (the hunt's target question)
  H7  ladder: every new rung-2 member = core + filler with filler in
      X(core) (the closed-form filler characterization, applied at n=8)
  H8  spectrum bottom table
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap_int
from lrc_xs_theory import X_theory, escape_intervals, level_set

OUT = "/home/z/my-project/scripts/out"
N = 8
G2 = Fraction(2, 17)          # rung-2 value at n = 8
TOT = 536878650               # C(50,8)


def parse_vec(line):
    body = line[line.index("(") + 1:line.index(")")]
    return tuple(int(x) for x in body.split(","))


def load(prefix, kind):
    out = []
    for i in range(8):
        try:
            for line in open(f"{OUT}/{prefix}_s{i}_{kind}.txt"):
                line = line.strip()
                if line.startswith("v=("):
                    out.append(parse_vec(line))
        except FileNotFoundError:
            pass
    return out


def binders(line):
    b = line[line.index("bind=(") + 6:line.index(")", line.index("bind=("))]
    return tuple(int(x) for x in b.split(","))


def load_bind(prefix, kind):
    out = []
    for i in range(8):
        try:
            for line in open(f"{OUT}/{prefix}_s{i}_{kind}.txt"):
                line = line.strip()
                if line.startswith("v=("):
                    out.append((parse_vec(line), binders(line),
                                 line.split(" t=")[1].split()[0]))
        except FileNotFoundError:
            pass
    return out


def main():
    ok = True

    # ---- H1 accounting ----
    spec_tot = filtered = fullsolved = 0
    spectrum = {}
    for i in range(8):
        for line in open(f"{OUT}/n8_V50_s{i}_spectrum.txt"):
            p = line.split()
            if p[0] == "filtered":
                filtered += int(p[1])
                continue
            num, den, c = int(p[0]), int(p[1]), int(p[2])
            g = gcd(num, den)
            spectrum[(num // g, den // g)] = \
                spectrum.get((num // g, den // g), 0) + c
            spec_tot += c
    print("H1 accounting:")
    print(f"  spectrum total {spec_tot} + filtered {filtered} = "
          f"{spec_tot + filtered}  (C(50,8) = {TOT}) "
          f"{'OK' if spec_tot + filtered == TOT else 'FAIL'}")
    ok &= (spec_tot + filtered == TOT)

    # ---- H2 sample validation ----
    nsamp = nok = 0
    for i in range(8):
        for line in open(f"{OUT}/n8_V50_s{i}_sample.txt"):
            line = line.strip()
            if not line.startswith("v=("):
                continue
            v = parse_vec(line)
            t = Fraction(*map(int,
                              line.split(" t=")[1].split()[0].split("/")))
            gapv = Fraction(*map(int,
                                 line.split(" gap=")[1].split()[0].split("/")))
            feas = int(line.split("feasible=")[1])
            g, ts = gap_int(list(v))
            nsamp += 1
            if g == gapv and (t in ts) and \
               (feas == (g * (N + 1) >= 1)):
                nok += 1
            else:
                print(f"  SAMPLE MISMATCH: {line}")
    print(f"H2 sample validation: {nok}/{nsamp} records match the "
          f"independent Python reference exactly")
    ok &= (nok == nsamp)

    # ---- H3 tight census ----
    tight = load("n8_V50", "tight")
    tset = set(tight)
    scal_consec = {tuple(c * j for j in range(1, 9)) for c in range(1, 7)}
    print(f"H3 tight census: {len(tight)} vectors at gap 1/9")
    for v in sorted(tset):
        is_scal = v in scal_consec
        print(f"   {v}  "
              f"{'consecutive family x %d' % v[0] if is_scal else 'NEW!'}")
    ok &= (tset == {s for s in scal_consec if max(s) <= 50})

    # ---- H4 rung-2 census ----
    r2 = load_bind("n8_V50", "rung2")
    r2v = [x[0] for x in r2]
    print(f"\nH4 rung-2 census: {len(r2v)} vectors at gap 2/17")
    prim = {}
    for v in sorted(set(r2v)):
        c = gcd(*v)
        p = tuple(x // c for x in v)
        prim.setdefault(p, []).append(c)
    for p in sorted(prim):
        print(f"   {p}  x {sorted(prim[p])}")
    # generalized Lemma A check on grids d = 17e + new members
    new = []
    lemmaA_ok = True
    for v, b, tstr in r2:
        t = Fraction(*map(int, tstr.split("/")))
        d, a = t.denominator, t.numerator
        e = d // 17
        if d % 17 or e == 0:
            lemmaA_ok = False
            continue
        res = {u: (a * u) % d for u in v}
        if min(min(r, d - r) for r in res.values()) != 2 * e:
            lemmaA_ok = False
        if not any(res[u] == 2 * e for u in b):
            lemmaA_ok = False
        if not any(res[u] == d - 2 * e for u in b):
            lemmaA_ok = False
        if not any((b[i] + b[j]) % d == 0
                   for i in range(len(b)) for j in range(len(b))
                   if i != j):
            lemmaA_ok = False
        if max(v) > 38:
            new.append((v, b, tstr))
    print(f"   Lemma A (grid d = 17e, residues >= 2e, binders at +-2e, "
          f"binder pair sum 0 mod d): "
          f"{'holds for all %d vectors' % len(r2) if lemmaA_ok else 'FAIL'}")
    ok &= lemmaA_ok
    print(f"   new members (v_8 > 38): {len(new)}")
    for v, b, tstr in new:
        print(f"     {v} bind={b} t*={tstr}")

    # ---- H5 rung-k zoos ----
    print(f"\nH5 rung-k zoos at n = 8 (value k/(8k+1)):")
    for k in range(3, 7):
        rk = load_bind("n8_V50", f"rung{k}")
        target = Fraction(k, 8 * k + 1)
        fam = tuple(list(range(1, 8)) + [8 * k])
        # generalized Lemma A on grids d = (8k+1)e
        la_ok = True
        for v, b, tstr in rk:
            t = Fraction(*map(int, tstr.split("/")))
            d, a = t.denominator, t.numerator
            D0 = 8 * k + 1
            if d % D0:
                la_ok = False
                continue
            e = d // D0
            res = {u: (a * u) % d for u in v}
            if min(min(r, d - r) for r in res.values()) != k * e:
                la_ok = False
            if not any(res[u] == k * e for u in b):
                la_ok = False
            if not any(res[u] == d - k * e for u in b):
                la_ok = False
            if not any((b[i] + b[j]) % d == 0
                       for i in range(len(b)) for j in range(len(b))
                       if i != j):
                la_ok = False
        prim = {}
        for v, b, tstr in rk:
            c = gcd(*v)
            prim.setdefault(tuple(x // c for x in v), []).append(c)
        print(f"   k={k}: {len(rk)} vectors, value {target}, "
              f"{len(prim)} primitive classes, family member {fam} "
              f"{'present' if fam in [x[0] for x in rk] else 'MISSING!'}, "
              f"Lemma A(mod {D0}e): {'OK' if la_ok else 'FAIL'}")
        for p in sorted(prim):
            print(f"      {p} x {sorted(prim[p])}")
        ok &= la_ok and (fam in [x[0] for x in rk])

    # ---- H6 inter-rung ----
    inter = load("n8_V50", "interrung")
    print(f"\nH6 inter-rung window (1/9, 2/17): "
          f"{'EMPTY at V=50' if not inter else len(inter)} "
          f"(vs 4 vectors at n=7, V=50; 2 at n=6)")
    ok &= (not inter)

    # ---- H7 ladder: new rung-2 members via X(S) ----
    print(f"\nH7 ladder via the filler characterization (n=8, g=2/17):")
    zoo8 = set(tuple(v) for v in r2v)
    for v, b, tstr in sorted(new):
        explained = []
        for S in combinations(v, 7):
            xt, R, note = X_theory(S, G2, 60)
            filler = tuple(x for x in v if x not in S)
            if len(filler) == 1 and filler[0] in xt:
                explained.append((S, filler[0], R))
        if explained:
            S, x, R = explained[0]
            print(f"   {v}: core {S} + filler {x} in X(core) "
                  f"(R={R}) -- {len(explained)} supporting cores")
        else:
            print(f"   {v}: NOT YET explained by a 7-core "
                  f"(sporadic / deeper structure)")

    # ---- H8 spectrum bottom ----
    print(f"\nH8 spectrum bottom at V=50 (exact for gap <= 2/15):")
    bot = sorted(((Fraction(num, den), c) for (num, den), c
                  in spectrum.items() if Fraction(num, den) <= Fraction(2, 15)),
                 key=lambda z: (z[0], -z[1]))
    for val, c in bot:
        mark = ""
        for k in range(1, 7):
            if val == Fraction(k, 8 * k + 1):
                mark = f"   <- rung-{k} value"
        print(f"   {val}  {c}{mark}")

    print(f"\nSUMMARY: {'ALL CHECKS PASS' if ok else 'SOME CHECKS FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
