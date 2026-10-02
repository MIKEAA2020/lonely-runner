"""Validation layers for the n = 8, V = 38 sharded sweep.

  S1  spectrum reconciliation: the merged reduced histogram must sum to
      C(38,8) = 48,903,492.
  S2  sample cross-check: every 99733rd vector recorded by each shard
      (full record: v, k, t, gap, feasible) is re-computed with the
      independent Python reference and must match exactly (gap value,
      argmax value, feasibility; k is a function of v and the argmax).
  S3  census totals: tight/rung2/interrung lists have the sizes the
      summaries claim, and every listed vector re-verifies exactly.
"""
import sys, glob
from fractions import Fraction

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap

def parse_line(line):
    d = {}
    for part in line.split():
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d

def main():
    out = "/home/z/my-project/scripts/out"
    # S1: spectrum sums to C(38,8)
    tot = 0
    for s in glob.glob(f"{out}/n8_V38_s*_spectrum.txt"):
        for line in open(s):
            num, den, cnt = line.split()
            tot += int(cnt)
    print(f"S1 spectrum total = {tot} (expect 48903492): "
          f"{'OK' if tot == 48903492 else 'FAIL'}")

    # S2: sample cross-check against the Python reference
    n_checked = n_ok = 0
    for s in sorted(glob.glob(f"{out}/n8_V38_s*_sample.txt")):
        for line in open(s):
            d = parse_line(line.strip())
            if "v" not in d:
                continue
            v = [int(x) for x in d["v"].strip("()").split(",")]
            krec = [int(x) for x in d["k"].strip("()").split(",")]
            ta, tb = map(int, d["t"].split("/"))
            num, den = map(int, d["gap"].split("/"))
            feas = int(d["feasible"])
            g, ts = gap(v)
            n_checked += 1
            ok = (g == Fraction(num, den)
                  and Fraction(ta, tb) in ts
                  and feas == (1 if g >= Fraction(1, 9) else 0)
                  and krec == [ta * x // tb for x in v])
            n_ok += ok
            if not ok:
                print("SAMPLE MISMATCH", d, g, ts)
    print(f"S2 sample cross-check: {n_ok}/{n_checked} OK "
          f"({100.0 * n_ok / max(n_checked, 1):.1f}%)")

    # S3: census sizes + re-verification
    for name, expect, target in (("tight", 4, Fraction(1, 9)),
                                 ("rung2", 14, Fraction(2, 17)),
                                 ("inter", 0, None)):
        vecs = []
        for s in glob.glob(f"{out}/n8_V38_s*_{name}.txt"):
            for line in open(s):
                line = line.strip()
                if line:
                    vecs.append([int(x) for x in
                                 line.split("(")[1].split(")")[0].split(",")])
        ok = len(vecs) == expect
        if target is not None:
            for v in vecs:
                g, _ = gap(v)
                if g != target:
                    ok = False
                    print("CENSUS VALUE MISMATCH", v, g)
        print(f"S3 census {name}: {len(vecs)} (expect {expect}): "
              f"{'OK' if ok else 'FAIL'}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
