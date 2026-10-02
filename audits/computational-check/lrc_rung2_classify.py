"""Rung-2 classification at n = 7 (item (a)): the 29 vectors with gap 2/15.

General necessary condition (Lemma 6.1 mechanism + value 2/15):
  gap(V) = 2/15 attained at t*  ==>
    t* = a/d reduced with 15 | d (write d = 15e, e >= 1),
    t* is a sum-crossing: some binding pair v,w in V has (v+w) t* in Z,
    v + w = 0 mod d  (in particular v + w = 0 mod 15),
    v a = +2e mod d and w a = -2e mod d (binders sit at +-2e = +-2d/15),
    every u in V has dist(u a mod d, d) >= 2e  (i.e. ua not in 0,+-1..+-2e-1).

For every census vector we compute, exactly:
  - gap re-verification (must be 2/15),
  - ALL argmax times t* (they come in mirror pairs t <-> 1-t),
  - the grid level e = d/15 of each argmax,
  - binding speeds and binding pairs at each argmax,
  - residue conditions above (checked, not assumed),
  - primitivity class (gcd) and scalar multiple structure,
  - Schur triples (u+v=w inside V) -- the sum-rich census,
  - killer cores: inclusion-minimal subsets S with gap(S) <= 2/15.
"""
import sys, json
from fractions import Fraction
from itertools import combinations
from math import gcd

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap, dist_at

TARGET = Fraction(2, 15)
MOD = 15

def schur_triples(v):
    S = set(v)
    return [(a, b, a + b) for a in v for b in v if b >= a and a + b in S]

def analyze_argmax(v, t):
    """Full residue picture of one argmax t (Fraction, reduced a/d)."""
    a, d = t.numerator, t.denominator
    e, r = divmod(d, MOD)
    assert r == 0, (v, t, "15 must divide the argmax denominator")
    dist2 = 2 * e                              # binder distance in units of 1/d
    res = {u: (u * a) % d for u in v}
    bind = sorted(u for u in v if min(res[u], d - res[u]) == dist2)
    bad = [u for u in v if min(res[u], d - res[u]) < dist2]
    pairs = [(u, w) for u in bind for w in bind
             if w > u and (u + w) % d == 0]
    return dict(t=str(t), a=a, d=d, e=e,
                binders=bind, pairs=pairs,
                residues={u: res[u] for u in sorted(res)},
                bad=bad)

def minimal_cores(v, thr=TARGET):
    """Inclusion-minimal subsets with gap <= thr, smallest sizes first."""
    hits = []
    for r in range(1, len(v) + 1):
        for sub in combinations(v, r):
            g, _ = gap(list(sub))
            if g <= thr:
                hits.append(sub)
    mini = [h for h in hits if not any(o != h and set(o) < set(h) for o in hits)]
    return mini

def main():
    src = ("/home/z/my-project/lonely-runner/audits/computational-check/"
           "tight_n7_V50_rung2.txt")
    vecs = []
    for line in open(src):
        line = line.strip()
        if not line:
            continue
        vecs.append([int(x) for x in
                     line.split("[")[1].split("]")[0].split(",")])
    assert len(vecs) == 29, len(vecs)

    rows = []
    nec_ok = True
    for v in vecs:
        g, ts = gap(v)
        assert g == TARGET, (v, g)
        d0 = gcd(*v)
        amaxes = [analyze_argmax(v, t) for t in ts]
        for am in amaxes:
            if am["bad"] or not am["pairs"]:
                nec_ok = False
        cores = minimal_cores(v)
        rows.append(dict(v=v, gcd=d0, amaxes=amaxes,
                         schur=[list(map(int, s)) for s in schur_triples(v)],
                         cores=[list(c) for c in cores],
                         core_gaps=[str(gap(list(c))[0]) for c in cores]))

    # ---- aggregate views -------------------------------------------------
    print("== necessary condition (15 | d, v+w = 0 mod d, +-2e, residues) ==")
    print("all 29 satisfy:", nec_ok)

    print("\n== grid levels e = d/15 of the argmaxes ==")
    lev = {}
    for r in rows:
        for am in r["amaxes"]:
            lev.setdefault(am["e"], set()).add(tuple(am["pairs"]))
    for e in sorted(lev):
        print(f"e={e} (d={15*e}): binding pair residue classes "
              f"{sorted(lev[e])}")

    print("\n== binding pairs observed (mod 15 classes in brackets) ==")
    pair_classes = {}
    for r in rows:
        key = tuple(sorted({tuple(p) for am in r["amaxes"] for p in am["pairs"]}))
        pair_classes.setdefault(key, []).append((r["v"], r["gcd"]))
    for k, vs in sorted(pair_classes.items()):
        print(f"pairs {k}: {len(vs)} vectors")

    print("\n== primitivity classes (base x multipliers) ==")
    prim = {}
    for r in rows:
        base = tuple(x // r["gcd"] for x in r["v"])
        prim.setdefault(base, []).append(r["gcd"])
    for base, mults in sorted(prim.items()):
        print(f"{list(base)}  x {sorted(mults)}")

    print("\n== killer cores (minimal subsets with gap <= 2/15) ==")
    core_set = {}
    for r in rows:
        for c, cg in zip(r["cores"], r["core_gaps"]):
            core_set.setdefault(tuple(c), set()).add(cg)
    for c, gaps_ in sorted(core_set.items(), key=lambda kv: (len(kv[0]), kv[0])):
        print(f"core {list(c)} gap={sorted(gaps_)} in "
              f"{sum(1 for r in rows if set(c) <= set(r['v']))} zoo vectors")

    print("\n== sum-richness: Schur triples ==")
    for r in rows:
        print(f"v={r['v']} gcd={r['gcd']} schur={r['schur']}")

    print("\n== per-argmax detail ==")
    for r in rows:
        for am in r["amaxes"]:
            print(f"v={r['v']} t*={am['t']} e={am['e']} binders={am['binders']}"
                  f" pairs={am['pairs']}")

    out = "/home/z/my-project/scripts/out/rung2_n7_classification.json"
    with open(out, "w") as f:
        json.dump(rows, f, indent=1)
    print(f"\nwrote {out}")
    return 0 if nec_ok else 1

if __name__ == "__main__":
    sys.exit(main())
