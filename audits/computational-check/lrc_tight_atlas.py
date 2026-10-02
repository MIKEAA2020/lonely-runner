"""Tight-set atlas: canonical X-chains, seed taxonomy, productivity census,
X-predictions beyond the census horizon, and scaling covariance.

Levels: an m-speed set sits at level m; its LRC-tight value is 1/(m+1).
Zoos at level m (V <= 50): tight (gap 1/(m+1)), rung-2 (2/(2m+1)),
rung-3 (3/(3m+1)), interrung (strictly between tight and rung-2).

A. CANONICAL CHAINS. For every tight n-set T: repeatedly remove the filler
   whose core has minimal gap; record the ascending gap sequence down to a
   single speed. Every step is X-certified (gap(U) = target at each level).

B. SEED TAXONOMY. Every census core S (level n-1) is classified by zoo
   membership at level n-1: tight / rung2 / rung3 / interrung / composite.

C. PRODUCTIVITY CENSUS. For every zoo member S at level n-1 (n = 4..8):
   X(S, 1/(n+1)) -- which members spawn tight children, and which children.

D. PREDICTIONS. For cores with X-elements beyond V = 50: verify
   gap(S + x) = 1/(n+1) exactly (falsifiable out-of-range check).

E. SCALING COVARIANCE. X(c*S, g) == c * X(S, g) on every scaled core pair.
"""
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations

sys.path.insert(0, "/home/z/my-project/scripts")
sys.path.insert(0, "/home/z/my-project/lonely-runner/audits/computational-check")
from lrc_gap_lib import gap_int
from lrc_xs_theory import X_theory, X_brute_batch

OUT = "/home/z/my-project/scripts/out"
CK = "/home/z/my-project/lonely-runner/audits/computational-check/"

ZOOS = {  # level -> {label: path}
    3: {"tight": OUT + "/rz3_tight.txt", "rung2": OUT + "/rz3_rung2.txt",
        "rung3": OUT + "/rz3_rung3.txt"},
    4: {"tight": OUT + "/rz4_tight.txt", "rung2": OUT + "/rz4_rung2.txt",
        "rung3": OUT + "/rz4_rung3.txt"},
    5: {"tight": OUT + "/rz5_tight.txt", "rung2": OUT + "/rz5_rung2.txt",
        "rung3": OUT + "/rz5_rung3.txt"},
    6: {"tight": OUT + "/rz6_tight.txt", "rung2": OUT + "/rz6_rung2.txt",
        "rung3": OUT + "/rz6_rung3.txt", "interrung": OUT + "/rz6_interrung.txt"},
    7: {"tight": CK + "tight_n7_V50.txt",
        "rung2": CK + "tight_n7_V50_rung2.txt"},
}
TIGHT_ZOOS = {
    2: OUT + "/tz2_tight.txt", 3: OUT + "/rz3_tight.txt",
    4: OUT + "/rz4_tight.txt", 5: OUT + "/rz5_tight.txt",
    6: OUT + "/rz6_tight.txt", 7: CK + "tight_n7_V50.txt",
    8: CK + "n8_V50_tight.txt",
}


def load(path):
    out = []
    for line in open(path):
        line = line.strip()
        if line.startswith("v=(") or line.startswith("v=["):
            close = ")" if line.startswith("v=(") else "]"
            body = line[3:line.index(close)]
            out.append(tuple(sorted(int(x) for x in body.split(","))))
    return sorted(set(out))


def batch_gaps(sets, target, xmax=0):
    if not sets:
        return {}
    inp = "".join(f"{len(c)} " + " ".join(map(str, c)) + "\n" for c in sets)
    proc = subprocess.run if False else None
    import subprocess as sp
    proc = sp.run(["/home/z/my-project/scripts/lrc_xs_batch",
                   str(target.numerator), str(target.denominator), str(xmax)],
                  input=inp, capture_output=True, text=True, check=True)
    res = {}
    for line in proc.stdout.splitlines():
        parts = line.split(" | ")
        if len(parts) >= 2:
            core = tuple(int(v) for v in parts[0].split(","))
            num, den = (int(z) for z in parts[1].split("/"))
            xs = [int(v) for v in parts[2].split(",")] if len(parts) > 2 and parts[2].strip() else []
            res[core] = (Fraction(num, den), xs)
    return res


def main():
    rep = {"chains": {}, "taxonomy": {}, "productivity": {},
           "predictions": {}, "scaling": {}, "chain_spectrum": {}}
    all_chain_gaps = Counter()
    failures = []

    # ---------------- A: canonical chains -------------------------------
    print("== A. canonical min-gap X-chains for every tight set ==")
    for n in range(2, 9):
        Z = load(TIGHT_ZOOS[n])
        chains = {}
        gap_cache = {}

        def gaps_of(cores):
            need = [c for c in cores if c not in gap_cache]
            if need:
                b = batch_gaps(need, Fraction(1, n + 1), 0)
                for c, (gp, _) in b.items():
                    gap_cache[c] = gp
            return {c: gap_cache[c] for c in cores}

        for T in Z:
            chain, cur, tgt = [], T, None
            gp_T, _ = gap_int(list(T))
            tgt = gp_T
            while len(cur) > 1:
                cores = [tuple(sorted(v for v in cur if v != x)) for x in cur]
                # dedupe (equal speeds impossible; cores distinct)
                gm = gaps_of(cores)
                S = min(cores, key=lambda c: (gm[c], c))   # canonical: min gap
                chain.append({"set": list(cur), "gap": str(tgt),
                              "removed": [x for x in cur
                                          if tuple(sorted(v for v in cur if v != x)) == S][0],
                              "core_gap": str(gm[S])})
                cur, tgt = S, gm[S]
            chain.append({"set": list(cur), "gap": "1/2", "removed": None,
                          "core_gap": None})
            chains[T] = chain
            for step in chain[:-1]:
                all_chain_gaps[step["core_gap"]] += 1
        lens = Counter(len(c) for c in chains.values())
        seeds = Counter(c[-1]["set"][0] for c in chains.values())
        rep["chains"][str(n)] = {
            "zoo": len(Z),
            "chain_length_dist": {str(k): v for k, v in lens.items()},
            "seeds": {str(k): v for k, v in seeds.items()},
            "example": chains[Z[0]] if Z else None,
        }
        print(f"  n={n}: {len(Z)} tight sets; chain lengths {dict(lens)}; "
              f"seed speeds {dict(seeds)}")
        # spot-print the mod-family chains at n=4,5,7
        if n in (4, 5, 7):
            for T in Z:
                if T[0] == 1 and len(set(T)) == len(T):
                    seq = [f"{s['set']}@{s['gap']}" for s in chains[T]]
                    print(f"    chain: {' -> '.join(seq)}")

    rep["chain_spectrum"] = {str(k): v for k, v in
                             sorted(all_chain_gaps.items(),
                                    key=lambda kv: Fraction(kv[0]))}
    print(f"  chain gap spectrum (union over all levels): "
          f"{ {str(k): v for k, v in sorted(all_chain_gaps.items(), key=lambda kv: Fraction(kv[0]))} }")

    # ---------------- B: seed taxonomy ----------------------------------
    print("\n== B. seed taxonomy: census cores by zoo membership ==")
    decomp = json.load(open(OUT + "/tight_decomp.json"))
    for n in range(3, 9):
        cores = [tuple(int(v) for v in k.split(","))
                 for k in decomp["levels"][str(n)]["core_details"]]
        zoos_m = ZOOS.get(n - 1, {})
        members = {}
        for label, path in zoos_m.items():
            for S in load(path):
                members.setdefault(S, []).append(label)
        tax = Counter()
        tax_detail = defaultdict(list)
        for S in cores:
            lab = "+".join(sorted(members.get(S, []))) if S in members else "composite"
            tax[lab] += 1
            tax_detail[lab].append(S)
        rep["taxonomy"][str(n)] = dict(tax)
        print(f"  n={n}: census cores by level-{n-1} zoo membership: {dict(tax)}")
        for lab in sorted(tax_detail):
            if lab not in ("tight", "rung2", "rung3", "tight+rung2"):
                ex = tax_detail[lab][:4]
                print(f"     {lab}: e.g. {ex}")

    # ---------------- C: productivity census ----------------------------
    print("\n== C. productivity census: zoo members spawning tight children ==")
    for n in range(4, 9):
        g = Fraction(1, n + 1)
        zoos_m = ZOOS.get(n - 1, {})
        prod = {}
        for label, path in sorted(zoos_m.items()):
            hits = []
            for S in load(path):
                xt, R, note = X_theory(S, g, None)
                if xt:
                    hits.append((S, xt, R))
            if hits:
                prod[label] = hits
        rep["productivity"][str(n)] = {
            label: [{"seed": list(S), "X": xt, "R": R} for S, xt, R in hits]
            for label, hits in prod.items()}
        for label, hits in prod.items():
            exts = []
            for S, xt, R in hits:
                for x in xt:
                    T = tuple(sorted(S + (x,)))
                    gp, _ = gap_int(list(T))
                    ok = (gp == g)
                    if not ok:
                        failures.append(("C-ext", n, T, str(gp)))
                    exts.append((list(S), x, str(gp), ok))
            print(f"  n={n} from level-{n-1} {label} seeds: {len(hits)} "
                  f"productive, extensions: "
                  f"{[(e[0], e[1]) for e in exts if e[3]][:8]}"
                  f"{' ...' if len(exts) > 8 else ''}")
        if not prod:
            print(f"  n={n}: no zoo seeds with nonempty X at 1/{n+1}")

    # ---------------- D: predictions beyond V=50 ------------------------
    print("\n== D. X-predictions beyond the census horizon (V > 50) ==")
    n_pred = n_ok = 0
    preds = []
    for n in range(3, 9):
        g = Fraction(1, n + 1)
        for k, info in decomp["levels"][str(n)]["core_details"].items():
            S = tuple(int(v) for v in k.split(","))
            for x in info["X"]:
                if x > 50:
                    T = tuple(sorted(S + (x,)))
                    gp, _ = gap_int(list(T))
                    n_pred += 1
                    if gp == g:
                        n_ok += 1
                    else:
                        failures.append(("D", n, T, str(gp)))
                    preds.append((n, list(S), x, str(gp)))
    rep["predictions"] = {"count": n_pred, "verified_tight": n_ok,
                          "sample": preds[:10]}
    print(f"  predicted tight sets with max speed > 50: {n_pred}; "
          f"verified gap = target: {n_ok}; failures: {n_pred - n_ok}")
    for p in preds[:6]:
        print(f"    n={p[0]}: {p[1]} + {p[2]}  -> gap {p[3]}")

    # ---------------- E: scaling covariance ------------------------------
    print("\n== E. scaling covariance X(c*S) == c*X(S) ==")
    n_pairs = n_ok = 0
    bad = []
    for n in range(3, 9):
        g = Fraction(1, n + 1)
        det = decomp["levels"][str(n)]["core_details"]
        # primitive cores = those with gcd 1; scaled partners = c*primitive
        prims = {}
        for k, info in det.items():
            S = tuple(int(v) for v in k.split(","))
            from math import gcd
            from functools import reduce
            d = reduce(gcd, S)
            if d == 1:
                prims[S] = info["X"]
        for k, info in det.items():
            S = tuple(int(v) for v in k.split(","))
            from math import gcd
            from functools import reduce
            d = reduce(gcd, S)
            if d == 1:
                continue
            P = tuple(v // d for v in S)
            if P in prims:
                n_pairs += 1
                if info["X"] == [d * x for x in prims[P]]:
                    n_ok += 1
                else:
                    bad.append((n, S, info["X"], [d * x for x in prims[P]]))
    rep["scaling"] = {"pairs": n_pairs, "covariant": n_ok, "violations": bad[:5]}
    print(f"  scaled core pairs: {n_pairs}; X(c*S) == c*X(S): {n_ok}; "
          f"violations: {n_pairs - n_ok}")
    for b in bad[:5]:
        print("    VIOLATION:", b)

    # ---------------- verdict --------------------------------------------
    print("\n" + "=" * 70)
    if failures:
        print(f"FAILURES: {len(failures)}")
        for f in failures[:10]:
            print("  ", f)
    else:
        print("ALL ATLAS CHECKS PASS")
    rep["failures"] = [list(map(str, f)) for f in failures[:20]]
    with open(OUT + "/tight_atlas.json", "w") as f:
        json.dump(rep, f, indent=1, default=str)
    print("wrote", OUT + "/tight_atlas.json")


if __name__ == "__main__":
    import subprocess
    sys.exit(main())
