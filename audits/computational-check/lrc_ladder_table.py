"""Ladder table at fixed V = 38: rung-1 (tight), rung-2, inter-rung counts
and primitive rung-2 classes for n = 5..8, plus cross-checks against the
committed censuses."""
import sys
from math import gcd

sys.path.insert(0, "/home/z/my-project/scripts")

def load(path):
    out = []
    try:
        for line in open(path):
            line = line.strip()
            if line:
                body = line.split("(")[1].split(")")[0] if "(" in line \
                    else line.split("[")[1].split("]")[0]
                out.append(tuple(int(x) for x in body.split(",")))
    except FileNotFoundError:
        print("missing", path)
    return out

def report(n, V, zoo, tight, inter):
    zv = [v for v in zoo if max(v) <= V]
    tv = [v for v in tight if max(v) <= V]
    iv = [v for v in inter if max(v) <= V]
    prim = {}
    for v in zv:
        d0 = gcd(*v)
        prim.setdefault(tuple(x // d0 for x in v), []).append(d0)
    print(f"n={n} V<={V}: tight={len(tv)} rung2={len(zv)} "
          f"(primitive classes: {len(prim)}) interrung={len(iv)}")
    return zv, tv, iv, prim

out = "/home/z/my-project/scripts/out"
repo = "/home/z/my-project/lonely-runner/audits/computational-check"

zoo5 = load(f"{out}/n5_V50_rung2.txt")
t5 = load(f"{out}/n5_V50_tight.txt")
i5 = load(f"{out}/n5_V50_interrung.txt")

zoo6 = load(f"{out}/n6_V50_rung2.txt")
t6 = load(f"{out}/n6_V50_tight.txt") if False else load(f"{out}/n6_V50_s0_tight.txt") + \
     load(f"{out}/n6_V50_s1_tight.txt") + load(f"{out}/n6_V50_s2_tight.txt") + \
     load(f"{out}/n6_V50_s3_tight.txt")
i6 = load(f"{out}/n6_V50_interrung.txt")

zoo7 = load(f"{repo}/tight_n7_V50_rung2.txt")
t7 = load(f"{repo}/tight_n7_V50.txt")
i7 = load(f"{repo}/inter_rung_n7_V50.txt")

zoo8 = load(f"{out}/n8_V38_rung2_merged.txt")
t8 = load(f"{out}/n8_V38_tight_merged.txt")
i8 = load(f"{out}/n8_V38_inter_merged.txt")

print("== ladder at V <= 38 ==")
report(5, 38, zoo5, t5, i5)
report(6, 38, zoo6, t6, i6)
report(7, 38, zoo7, t7, i7)
report(8, 38, zoo8, t8, i8)

print("\n== ladder at V <= 50 (n = 5..7), V <= 38 (n = 8) ==")
for n, zoo in ((5, zoo5), (6, zoo6), (7, zoo7), (8, zoo8)):
    prim = {}
    for v in zoo:
        d0 = gcd(*v)
        prim.setdefault(tuple(x // d0 for x in v), []).append(d0)
    print(f"n={n}: rung2 vectors={len(zoo)}, primitive classes={len(prim)}")

print("\n== cross-checks ==")
# v3 n=7 V=38 rung2 vs committed V=50 census filtered to 38
z7v3 = load(f"{out}/n7_V38v3_rung2.txt")
z7_38 = [v for v in zoo7 if max(v) <= 38]
print(f"n=7 V=38 rung2: v3={len(z7v3)}, committed-filtered={len(z7_38)}, "
      f"match={set(z7v3) == set(z7_38)}")
t7v3 = load(f"{out}/n7_V38v3_tight.txt")
t7_38 = [v for v in t7 if max(v) <= 38]
print(f"n=7 V=38 tight: v3={len(t7v3)}, committed-filtered={len(t7_38)}, "
      f"match={set(t7v3) == set(t7_38)}")
i7v3 = load(f"{out}/n7_V38v3_interrung.txt")
i7_38 = [v for v in i7 if max(v) <= 38]
print(f"n=7 V=38 interrung: v3={len(i7v3)}, committed-filtered={len(i7_38)}, "
      f"match={set(i7v3) == set(i7_38)}")

# V=30 reconciliation for n=8
z8_30 = load(f"{out}/n8_V30_rung2.txt")
z8_30b = [v for v in zoo8 if max(v) <= 30]
print(f"n=8 V=30 rung2: V30-run={len(z8_30)}, V38-filtered={len(z8_30b)}, "
      f"match={set(z8_30) == set(z8_30b)}")

# direct lifts: zoo_n member + one speed in zoo_{n+1}
print("\n== direct lifts (zoo member + one speed = next zoo) ==")
for nm, zn, znext in (("5->6", zoo5, zoo6), ("6->7", zoo6, zoo7),
                      ("7->8", zoo7, zoo8)):
    lifts = []
    zset = set(znext)
    for v in zn:
        for x in range(1, 51):
            w = tuple(sorted(v + (x,))) if x not in v else None
            if w and w in zset:
                lifts.append((v, x, w))
    print(f"{nm}: {len(lifts)} direct lifts")
    for v, x, w in lifts:
        print(f"   {list(v)} + {x} = {list(w)}")
