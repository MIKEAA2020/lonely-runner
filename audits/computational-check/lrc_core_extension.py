"""Core-extension analysis of the n=7 rung-2 zoo.

Zoo members are inclusion-minimal with gap <= 2/15 (verified), so the
classification runs through 6-subsets ("cores") and their admissible
extensions: X(S) = {x not in S : gap(S + x) == 2/15}.

For each 6-subset S of each zoo member we compute:
  - gap(S) and its argmaxes ("escape times" of the core),
  - X(S) restricted to [1,50],
and then decompose the zoo as  {S + x : x in X(S)} over supporting cores.

Mechanistic reading (Lemma 6.2 transplant): at an escape time t' of S the
core has min > 2/15; a filler x is admissible only if ||x t'|| <= 2/15 at
EVERY escape time -- the filler must kill all of the core's escapes.
"""
import sys
from fractions import Fraction
from itertools import combinations

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap, dist_at

TARGET = Fraction(2, 15)
VMAX = 50

def zoo():
    src = ("/home/z/my-project/lonely-runner/audits/computational-check/"
           "tight_n7_V50_rung2.txt")
    out = []
    for line in open(src):
        line = line.strip()
        if line:
            out.append(tuple(int(x) for x in
                             line.split("[")[1].split("]")[0].split(",")))
    return out

def main():
    Z = zoo()
    zset = set(Z)

    # gap of every 6-subset of every zoo member; find zoo-supporting cores
    core_info = {}
    for v in Z:
        for S in combinations(v, 6):
            if S in core_info:
                continue
            g, ts = gap(list(S))
            ext = []
            for x in range(1, VMAX + 1):
                if x in S:
                    continue
                gx, _ = gap(list(S) + [x])
                if gx == TARGET:
                    ext.append(x)
            core_info[S] = (g, ts, ext)

    print("== 6-subsets of zoo members: gap, argmaxes, extensions X(S) ==")
    for S in sorted(core_info, key=lambda s: (core_info[s][0], s)):
        g, ts, ext = core_info[S]
        mark = "SUPPORTS" if ext else "         "
        print(f"{mark} S={list(S)} gap={g} argmax={[str(t) for t in ts]} "
              f"X={ext}")

    # supporting cores generate the zoo?
    print("\n== zoo decomposition over supporting cores ==")
    gen = {}
    for S, (g, ts, ext) in core_info.items():
        if ext:
            for x in ext:
                W = tuple(sorted(S + (x,)))
                if W in zset:
                    gen.setdefault(S, []).append(x)
    for S in sorted(gen):
        print(f"core {list(S)} (gap={core_info[S][0]}) -> "
              f"{len(gen[S])} zoo members via x in {sorted(gen[S])}")

    covered = {tuple(sorted(S + (x,))) for S in gen for x in gen[S]}
    print(f"\nzoo members generated: {len(covered)} / {len(zset)}")
    missing = zset - covered
    if missing:
        print("NOT generated (need deeper cores):")
        for m in sorted(missing):
            print("  ", list(m))

    # escape-time picture for the two main supporting cores
    print("\n== escape times of key cores (times where min_S > 2/15) ==")
    for S in sorted(gen, key=lambda s: -len(gen[s]))[:6]:
        g, ts, _ = core_info[S]
        print(f"core {list(S)} gap={g} argmaxes={[str(t) for t in ts]}")
        for x in gen[S][:8]:
            # where does the filler sit at the core's argmaxes?
            pos = [str(dist_at(x, t)) for t in ts]
            print(f"   filler x={x:2d}: ||x t'|| at core argmaxes = {pos}")

    return 0

if __name__ == "__main__":
    sys.exit(main())
