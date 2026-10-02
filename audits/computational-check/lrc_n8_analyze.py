"""n = 8 zoo analysis and validation (item (b): does the tight-set zoo stay
closed at n = 8, V = 38, value 2/17?).

Layers:
  V1  merge shard outputs; reconcile counts (sum of shard counts = C(V,8);
      sum of tight/rung2/interrung = totals).
  V2  Python reference re-verification of EVERY census vector (gap exact,
      argmax on the d-grid, 17 | d, binding pair v+w = 0 mod d, binders at
      +-2d/17, residue cleanliness) -- the generalized necessary condition.
  V3  primitivity classes, scalar multiplicities, binding-pair classes,
      hub-core decomposition against the n = 6 / n = 7 ladders.
  V4  prediction check: zoo members must extend a lower-rung core
      (subset structure); report which zoo_n structures lift.
  V5  spectrum bottom (merged histogram) + V30 -> V38 reconciliation.
"""
import sys, json, glob
from fractions import Fraction
from itertools import combinations
from math import gcd

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap, dist_at

N, MOD = 8, 17
VLO, VHI = 30, 38

def parse_vecs(path):
    out = []
    try:
        for line in open(path):
            line = line.strip()
            if line:
                body = line.split("(")[1].split(")")[0] if "(" in line \
                    else line.split("[")[1].split("]")[0]
                out.append(tuple(int(x) for x in body.split(",")))
    except FileNotFoundError:
        pass
    return out

def analyze_argmax(v, t):
    a, d = t.numerator, t.denominator
    e, r = divmod(d, MOD)
    if r != 0:
        return None
    dist2 = 2 * e
    res = {u: (u * a) % d for u in v}
    bind = sorted(u for u in v if min(res[u], d - res[u]) == dist2)
    bad = [u for u in v if min(res[u], d - res[u]) < dist2]
    pairs = [(u, w) for u in bind for w in bind if w > u and (u + w) % d == 0]
    return dict(t=str(t), a=a, d=d, e=e, binders=bind, pairs=pairs, bad=bad)

def main():
    out = "/home/z/my-project/scripts/out"
    # ---- V1: merge + reconcile -------------------------------------------
    tot = dict(total=0, feasible=0, tight=0, rung2=0, interrung=0,
               vf=0, secs=0.0)
    for s in sorted(glob.glob(f"{out}/n8_V38_s*_summary.txt")):
        for line in open(s):
            kv = dict(p.split("=") for p in line.split())
            tot["total"] += int(kv["count"])
            tot["feasible"] += int(kv["feasible"])
            tot["tight"] += int(kv["tight"])
            tot["rung2"] += int(kv["rung2"])
            tot["interrung"] += int(kv["interrung"])
            tot["vf"] += int(kv["verify_failures"])
            tot["secs"] += float(kv["seconds"])
    print("== V1 shard reconciliation ==")
    print(tot)
    assert tot["total"] == 48903492, tot["total"]
    assert tot["vf"] == 0

    zoo = parse_vecs(f"{out}/n8_V38_rung2_merged.txt")
    tight = parse_vecs(f"{out}/n8_V38_tight_merged.txt")
    inter = parse_vecs(f"{out}/n8_V38_inter_merged.txt")
    assert len(zoo) == tot["rung2"], (len(zoo), tot["rung2"])
    assert len(tight) == tot["tight"]
    assert len(inter) == tot["interrung"]

    # ---- V2: reference re-verification -----------------------------------
    print("\n== V2 reference re-verification of census vectors ==")
    nec_ok = True
    rows = []
    for v in zoo:
        g, ts = gap(list(v))
        assert g == Fraction(2, 17), (v, g)
        amaxes = [analyze_argmax(v, t) for t in ts]
        amaxes = [am for am in amaxes if am]
        for am in amaxes:
            if am["bad"] or not am["pairs"]:
                nec_ok = False
                print("NEC FAIL", v, am)
        d0 = gcd(*v)
        rows.append(dict(v=list(v), gcd=d0, amaxes=amaxes))
    for v in tight:
        g, ts = gap(list(v))
        assert g == Fraction(1, 9), (v, g)
    for v in inter:
        g, _ = gap(list(v))
        assert Fraction(1, 9) < g < Fraction(2, 17), (v, g)
    print(f"all {len(zoo)} rung2 + {len(tight)} tight + {len(inter)} "
          f"interrung re-verified; necessary condition: {nec_ok}")

    # ---- V3: structure ----------------------------------------------------
    print("\n== V3 primitivity classes ==")
    prim = {}
    for r in rows:
        base = tuple(x // r["gcd"] for x in r["v"])
        prim.setdefault(base, []).append(r["gcd"])
    for base, mults in sorted(prim.items()):
        print(f"{list(base)}  x {sorted(mults)}")

    print("\n== V3 binding-pair classes (at e=1 grids d=17) ==")
    pcls = {}
    for r in rows:
        key = tuple(sorted({tuple(p) for am in r["amaxes"] if am["e"] == 1
                            for p in am["pairs"]}))
        pcls.setdefault(key, []).append(tuple(x // r["gcd"] for x in r["v"]))
    for k, vs in sorted(pcls.items()):
        print(f"e=1 pairs {k}: {sorted(set(vs))}")

    print("\n== V3 grid levels ==")
    lev = sorted({am["e"] for r in rows for am in r["amaxes"]})
    print("e values:", lev, "(d = 17e)")

    # ---- V4: ladder / hub structure --------------------------------------
    print("\n== V4 hub-core decomposition (which lower zoos lift) ==")
    zoo6 = set(parse_vecs("/home/z/my-project/scripts/out/n6_V50_rung2.txt"))
    zoo7 = set(parse_vecs("/home/z/my-project/lonely-runner/audits/"
                          "computational-check/tight_n7_V50_rung2.txt"))
    zoo5 = set(parse_vecs("/home/z/my-project/scripts/out/n5_V50_rung2.txt"))
    for name, zn in (("zoo5", zoo5), ("zoo6", zoo6), ("zoo7", zoo7)):
        hits = 0
        for r in rows:
            v = set(r["v"])
            for core in zn:
                c = set(core)
                if c <= v:
                    hits += 1
                    break
        print(f"zoo members extending a {name} core: {hits}/{len(rows)}")
    # 7-subset cores of the n=8 zoo (zoo members minus one element)
    print("\n== V4 7-subset cores of the n=8 zoo with multiple extensions ==")
    core_ext = {}
    for v in set(zoo):
        for S in combinations(v, 7):
            core_ext.setdefault(S, 0)
    ext_of = {}
    for v in zoo:
        for S in combinations(v, 7):
            ext_of.setdefault(S, set()).add(v)
    for S, exts in sorted(ext_of.items(), key=lambda kv: -len(kv[1])):
        if len(exts) > 1 or True:
            gS, tsS = gap(list(S))
            inz6 = S in zoo6
            inz7 = S in zoo7
            tag = ("zoo6" if inz6 else ("zoo7" if inz7 else ""))
            print(f"core {list(S)} gap={gS} {tag} exts="
                  f"{sorted(x[0] if x[0] not in S else x[-1] for x in exts)}")

    # ---- V5: spectrum bottom + reconciliation ----------------------------
    print("\n== V5 spectrum bottom (V=38) vs V=30 ==")
    hist = {}
    for s in glob.glob(f"{out}/n8_V38_s*_spectrum.txt"):
        for line in open(s):
            num, den, cnt = line.split()
            key = Fraction(int(num), int(den))
            hist[key] = hist.get(key, 0) + int(cnt)
    hist30 = {}
    try:
        for line in open(f"{out}/n8_V30_spectrum.txt"):
            pass  # run mode has no spectrum; use shard mode only
    except FileNotFoundError:
        pass
    bottom = sorted(hist.items(), key=lambda kv: kv[0])[:15]
    for val, cnt in bottom:
        print(f"  {val}: {cnt}")
    # V=30 comparison via zoo files
    zoo30 = parse_vecs(f"{out}/n8_V30_rung2.txt")
    zoo30_v38 = [v for v in zoo if max(v) <= 30]
    print(f"rung2: V=30 census {len(zoo30)}, V<=30 members at V=38 "
          f"{len(zoo30_v38)} (match: {set(zoo30) == set(zoo30_v38)})")
    tight30 = parse_vecs(f"{out}/n8_V30_tight.txt")
    tight30_v38 = [v for v in tight if max(v) <= 30]
    print(f"tight: V=30 census {len(tight30)}, V<=30 members at V=38 "
          f"{len(tight30_v38)} (match: {set(tight30) == set(tight30_v38)})")

    json.dump(dict(rows=rows, tot=tot,
                   prim={str(k): v for k, v in prim.items()}),
              open(f"{out}/rung2_n8_classification.json", "w"), indent=1)
    print("\nwrote rung2_n8_classification.json")
    return 0 if nec_ok else 1

if __name__ == "__main__":
    sys.exit(main())
