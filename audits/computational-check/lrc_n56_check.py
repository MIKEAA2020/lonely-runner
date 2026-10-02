"""Necessary-condition verification for the rung-2 zoos at n = 5, 6
(values 2/11, 2/13), mirroring the n = 7 / n = 8 analysis."""
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap

def check(path, n):
    MOD = 2 * n + 1
    TARGET = Fraction(2, MOD)
    vecs = []
    for line in open(path):
        line = line.strip()
        if line:
            body = line.split("(")[1].split(")")[0] if "(" in line \
                else line.split("[")[1].split("]")[0]
            vecs.append(tuple(int(x) for x in body.split(",")))
    ok = True
    grids, pairset = set(), set()
    for v in vecs:
        g, ts = gap(list(v))
        assert g == TARGET, (v, g)
        for t in ts:
            a, d = t.numerator, t.denominator
            e, r = divmod(d, MOD)
            if r != 0:
                ok = False
                print("GRID FAIL", v, t)
                continue
            dist2 = 2 * e
            res = {u: (u * a) % d for u in v}
            bind = [u for u in v if min(res[u], d - res[u]) == dist2]
            bad = [u for u in v if min(res[u], d - res[u]) < dist2]
            pairs = [(u, w) for u in bind for w in bind
                     if w > u and (u + w) % d == 0]
            if bad or not pairs:
                ok = False
                print("NEC FAIL", v, t, bad, pairs)
            grids.add(e)
            for p in pairs:
                pairset.add(tuple(x % MOD for x in p))
    print(f"n={n}: {len(vecs)} vectors, necessary condition "
          f"{'HOLDS' if ok else 'FAILS'}, grid levels e={sorted(grids)}, "
          f"binding pair classes mod {MOD}: {sorted(pairset)}")
    # primitive classes
    prim = {}
    for v in vecs:
        d0 = gcd(*v)
        prim.setdefault(tuple(x // d0 for x in v), []).append(d0)
    print(f"  primitive classes: {len(prim)}")
    for base, mults in sorted(prim.items()):
        print(f"    {list(base)} x {sorted(mults)}")
    return ok

check("/home/z/my-project/scripts/out/n5_V50_rung2.txt", 5)
print()
check("/home/z/my-project/scripts/out/n6_V50_rung2.txt", 6)
