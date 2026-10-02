"""Analysis of the no-silent-kill (NSK) census: merge shards, verify
every X-filler against the independent gap references, classify the
productive cores, reconcile with the Task-3 productivity census, and
write the report + JSON.

Census claim being certified (Theorem R, silent-killer reduction):
for every level m and every (m-1)-core S in the universe,
    D(S, 1/(m+1)) = X(S, 1/(m+1)),
i.e. every escape-killer touches.  Consequence: no m-set T with
T \\ {max(T)} inside the universe has gap < 1/(m+1).
"""
import json
import subprocess
import sys
from fractions import Fraction
from math import gcd
from functools import reduce

sys.path.insert(0, "/home/z/my-project/scripts")
from lrc_gap_lib import gap_int

OUT = "/home/z/my-project/scripts/out"
BATCH = "/home/z/my-project/scripts/lrc_xs_batch"

TIERS = [
    (2, 500, ["nsk_m2_B500"]),
    (3, 300, ["nsk_m3_B300"]),
    (4, 200, ["nsk_m4_B200"]),
    (5, 150, ["nsk_m5_B150"]),
    (6, 120, ["nsk_m6_s0", "nsk_m6_s1"]),
    (7, 80, ["nsk_m7_s0", "nsk_m7_s1", "nsk_m7_s2", "nsk_m7_s3"]),
    (8, 60, ["nsk_m8_s0", "nsk_m8_s1", "nsk_m8_s2", "nsk_m8_s3",
             "nsk_m8_s4", "nsk_m8_s5"]),
    (9, 50, ["nsk_m9_s0", "nsk_m9_s1", "nsk_m9_s2", "nsk_m9_s3",
             "nsk_m9_s4", "nsk_m9_s5", "nsk_m9_s6", "nsk_m9_s7"]),
    (10, 40, ["nsk_m10_s0", "nsk_m10_s1", "nsk_m10_s2", "nsk_m10_s3",
              "nsk_m10_s4", "nsk_m10_s5"]),
]

# known productive seeds from the Task-3 productivity census (zoo cores)
KNOWN_SEEDS = {
    4: [((1, 3, 4), {2, 7}), ((3, 4, 7), {1})],
    5: [((1, 3, 4, 5), {2, 9}), ((1, 4, 5, 9), {3})],
    7: [((1, 2, 3, 4, 5, 12), {7}), ((1, 4, 5, 6, 7, 11), {13})],
}


def load_cores(prefix):
    out = {}
    for line in open(f"{OUT}/{prefix}_cores.txt"):
        left, rest = line.split(" | ", 1)
        S = tuple(int(v) for v in left[2:].rstrip(",").split(","))
        parts = rest.split(" | ")
        R = int(parts[0][2:])
        D = [int(v) for v in parts[1][2:].split(",") if v]
        X = [int(v) for v in parts[2][2:].split(",")
             if v and v != "-1"]
        out[S] = (R, D, X)
    return out


def load_summary(prefix):
    d = {}
    for line in open(f"{OUT}/{prefix}_summary.txt"):
        for kv in line.split():
            if "=" in kv:
                k, v = kv.split("=", 1)
                try:
                    d[k] = int(v)
                except ValueError:
                    d[k] = v
    return d


def brute_gaps(sets):
    inp = "".join(f"{len(t)} " + " ".join(map(str, t)) + "\n" for t in sets)
    p = subprocess.run([BATCH, "1", "2", "0"], input=inp, capture_output=True,
                       text=True, check=True)
    res = {}
    for line in p.stdout.splitlines():
        parts = line.split(" | ")
        if len(parts) < 2:
            continue
        T = tuple(int(v) for v in parts[0].split(","))
        num, den = (int(z) for z in parts[1].split("/"))
        res[T] = Fraction(num, den)
    return res


def fam_minus_one(m):
    """c * ([m] minus one element) family, as a classifier."""
    def classify(S, m):
        d = reduce(gcd, S)
        P = tuple(v // d for v in S)
        base = set(range(1, m + 1))
        for j in range(1, m + 1):
            if P == tuple(sorted(base - {j})):
                return ("ladder", d, j)
        return (None, d, None)
    return classify


def main():
    report = {"tiers": {}, "productive": {}, "reconciliation": {},
              "verification": {}}
    failures = []
    total_cores = 0
    total_silent = 0

    print("=" * 72)
    print("NSK CENSUS ANALYSIS  (D(S, 1/(m+1)) = X(S, 1/(m+1)) ?)")
    print("=" * 72)
    for (m, B, prefixes) in TIERS:
        g = Fraction(1, m + 1)
        tallies = []
        cores = {}
        for p in prefixes:
            s = load_summary(p)
            tallies.append(s)
            for S, v in load_cores(p).items():
                if S in cores and cores[S] != v:
                    failures.append(("shard conflict", m, S, cores[S], v))
                cores[S] = v
            total_silent += s.get("SILENT", 0)
            flagtxt = open(f"{OUT}/{p}_flags.txt").read().strip()
            if flagtxt:
                failures.append(("flags nonempty", m, p, flagtxt[:200]))
        tot = sum(t.get("cores", 0) for t in tallies)
        total_cores += tot
        nD = len(cores)
        sumD = sum(len(v[1]) for v in cores.values())
        sumX = sum(len(v[2]) for v in cores.values())
        maxR = max((v[0] for v in cores.values()), default=0)
        exp = 1
        for i in range(m - 1):
            exp *= (B - i)
        for i in range(m - 1):
            exp //= (i + 1)
        ok_count = (tot == exp)
        print(f"\nm={m}  B={B}  g=1/{m+1}")
        print(f"  cores enumerated: {tot:,} (expected C({B},{m-1}) = "
              f"{exp:,})  {'OK' if ok_count else 'MISMATCH'}")
        print(f"  skip_certified={sum(t.get('skip_certified',0) for t in tallies):,}"
              f"  full={sum(t.get('full',0) for t in tallies):,}")
        print(f"  SILENT={sum(t.get('SILENT',0) for t in tallies)}   "
              f"productive cores={nD}   sumD={sumD}  sumX={sumX}   "
              f"max_R={maxR}")
        if sumD != sumX:
            failures.append(("D != X totals", m, sumD, sumX))
        report["tiers"][str(m)] = {
            "B": B, "g": f"1/{m+1}", "cores": tot, "expected": exp,
            "productive": nD, "sumD": sumD, "sumX": sumX, "max_R": maxR,
            "silent": sum(t.get("SILENT", 0) for t in tallies),
        }
        report["productive"][str(m)] = {
            ",".join(map(str, S)): v for S, v in sorted(cores.items())
        }

        # classification: ladder family vs known seeds vs other
        classify = fam_minus_one(m)
        fam, seeds_ok, other = 0, 0, []
        for S, (R, D, X) in cores.items():
            kind, d, j = classify(S, m)
            if kind == "ladder":
                fam += 1
                continue
            if m in KNOWN_SEEDS:
                d = reduce(gcd, S)
                P = tuple(v // d for v in S)
                if any(P == ks for (ks, _) in KNOWN_SEEDS[m]):
                    seeds_ok += 1
                    continue
            other.append((S, (R, D, X)))
        print(f"  productive classification: ladder(c*[m]\\j)={fam}, "
              f"known mod-seeds={seeds_ok}, OTHER={len(other)}")
        if other:
            for (S, v) in other[:6]:
                print(f"    OTHER: {S} -> R={v[0]} D={v[1]} X={v[2]}")
        report["tiers"][str(m)]["ladder"] = fam
        report["tiers"][str(m)]["seeds"] = seeds_ok
        report["tiers"][str(m)]["other"] = len(other)

        # reconciliation with Task-3
        if m in KNOWN_SEEDS:
            for (S, Xwant) in KNOWN_SEEDS[m]:
                got = cores.get(S)
                if got is None or set(got[2]) != Xwant:
                    failures.append(("seed reconcile", m, S, Xwant, got))
                else:
                    print(f"  reconcile {S} -> X={sorted(Xwant)} OK")

    # ---------------- verification of all X-fillers ----------------
    print("\n" + "=" * 72)
    print("VERIFICATION: every X-filler gives gap exactly 1/(m+1)")
    print("=" * 72)
    n_ver = 0
    sample_neg = []
    for (m, B, prefixes) in TIERS:
        g = Fraction(1, m + 1)
        cores = {}
        for p in prefixes:
            cores.update(load_cores(p))
        sets_ = []
        for S, (R, D, X) in cores.items():
            for x in X:
                sets_.append(tuple(sorted(S + (x,))))
            negs = [x for x in range(1, R + 1)
                    if x not in S and x not in D][:2]
            for x in negs:
                sample_neg.append((m, g, tuple(sorted(S + (x,)))))
        gaps = brute_gaps(sets_)
        for S, (R, D, X) in cores.items():
            for x in X:
                T = tuple(sorted(S + (x,)))
                if gaps.get(T) != g:
                    failures.append(("X-filler gap", m, S, x, gaps.get(T)))
                n_ver += 1
    print(f"  {n_ver} X-fillers verified exact (C brute lrc_xs_batch)")
    # dual spot check with gap_int on a sample across tiers
    import random
    random.seed(20261004)
    spot = random.sample(sample_neg, min(120, len(sample_neg)))
    n_spot = 0
    for (m, g, T) in spot:
        gp, _ = gap_int(list(T))
        if gp <= g:
            failures.append(("negative sample", m, T, gp))
        n_spot += 1
    print(f"  {n_spot} negative samples (x not in D) all have gap > g "
          "(gap_int dual check)")
    report["verification"] = {"x_fillers": n_ver, "negative_samples": n_spot}

    print("\n" + "=" * 72)
    print(f"TOTALS: {total_cores:,} cores enumerated; "
          f"SILENT KILLERS: {total_silent}")
    print("=" * 72)
    report["totals"] = {"cores": total_cores, "silent": total_silent}

    if failures:
        print(f"\nFAILURES: {len(failures)}")
        for f in failures[:12]:
            print("  ", f)
        return 1, report
    print("\nALL NSK CENSUS CLAIMS VERIFIED")
    with open(f"{OUT}/nsk_census_report.json", "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)
    print(f"wrote {OUT}/nsk_census_report.json")
    return 0, report


if __name__ == "__main__":
    rc, rep = main()
    sys.exit(rc)
