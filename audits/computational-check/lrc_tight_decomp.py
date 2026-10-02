"""Tight-set decomposition experiment: point the X(S) lemma at the tight sets.

For every tight n-set T (gap(T) = 1/(n+1), max speed <= V_n) and every
x in T, let S = T \ {x} (an (n-1)-core).  By the filler theorem,

    gap(S + {x}) = 1/(n+1)   <=>   x in X(S, 1/(n+1))   [kill + touch].

This script:
  T1  validates the theorem on every (T, x) pair of every tight zoo
      n = 2..8 (theory membership == brute-force membership);
  T2  classifies every core: tight@n-1 (gap = 1/n), rung-k@n-1
      (gap = k/(k(n-1)+1)), or other (exact value recorded);
  T3  naive-ladder test: which tight sets have a tight core;
  T4  X-completeness at every core: zoo extensions over S ==
      X(S, 1/(n+1)) cut to [1, V_n]  (the zoo is X-closed);
  T5  core spectrum and productive tight sets (which Z_{n-1} members
      appear as tight cores);
  T6  ladder chains: T -> tight core -> ... -> seed, for every zoo member.

Output: report to stdout + JSON to out/tight_decomp.json.
"""
import json
import subprocess
import sys
from collections import Counter, defaultdict
from fractions import Fraction

sys.path.insert(0, "/home/z/my-project/scripts")
sys.path.insert(0, "/home/z/my-project/lonely-runner/audits/computational-check")
from lrc_gap_lib import gap_int
from lrc_xs_theory import X_theory, escape_intervals, level_set, X_brute_batch

OUT = "/home/z/my-project/scripts/out"
BATCH = "/home/z/my-project/scripts/lrc_xs_batch"

ZOOS = {
    2: (50, OUT + "/tz2_tight.txt"),
    3: (50, OUT + "/tz3_tight.txt"),
    4: (50, OUT + "/tz4_tight.txt"),
    5: (50, OUT + "/tz5_tight.txt"),
    6: (50, OUT + "/tz6_tight.txt"),
    7: (50, "/home/z/my-project/lonely-runner/audits/computational-check/"
            "tight_n7_V50.txt"),
    8: (50, "/home/z/my-project/lonely-runner/audits/computational-check/"
            "n8_V50_tight.txt"),
}
XMAX = 60          # brute-force filler horizon (covers V_n = 50)


def load_zoo(path):
    out = []
    for line in open(path):
        line = line.strip()
        if line.startswith("v=("):
            body = line[3:line.index(")")]
            out.append(tuple(int(x) for x in body.split(",")))
    return sorted(set(out))


def rung_k(gap, m):
    """If gap = k/(k*m+1) for integer k >= 2 with m speeds, return k."""
    if gap is None:
        return None
    for k in range(2, 40):
        if gap == Fraction(k, k * m + 1):
            return k
    return None


def main():
    report = {"levels": {}}
    all_fail = []

    for n, (Vn, path) in sorted(ZOOS.items()):
        g = Fraction(1, n + 1)
        Z = load_zoo(path)
        assert all(len(T) == n for T in Z)
        # sanity: every member really is tight
        for T in Z:
            gp, _ = gap_int(list(T))
            assert gp == g, (T, gp, g)

        # ---- collect cores -------------------------------------------------
        pairs = []          # (T, x, S)
        cores = {}
        for T in Z:
            for x in T:
                S = tuple(sorted(v for v in T if v != x))
                pairs.append((T, x, S))
                cores.setdefault(S, [])

        # ---- brute X lists (C batch) --------------------------------------
        core_list = sorted(cores)
        batch = X_brute_batch(core_list, g, XMAX)

        # ---- T1 + T4 ------------------------------------------------------
        level = {
            "n": n, "V": Vn, "gap_target": f"1/{n+1}", "zoo_size": len(Z),
            "pairs": len(pairs), "cores": len(core_list),
        }
        core_info = {}
        for S in core_list:
            cgap, xbrute = batch[S]
            assert cgap > g, ("core below target?!", S, cgap)
            xt, R, note = X_theory(S, g, None)
            # T1: theory == brute on [1, min(R, XMAX)] and brute empty above R
            lo = min(R if R is not None else XMAX, XMAX)
            xt_lo = [x for x in xt if x <= lo]
            xb_lo = [x for x in xbrute if x <= lo]
            ok1 = xt_lo == xb_lo and not [x for x in xbrute if x > lo]
            # classify core
            cls = None
            if cgap == Fraction(1, n):
                cls = "tight"
            else:
                k = rung_k(cgap, n - 1)
                if k:
                    cls = f"rung{k}"
                else:
                    cls = f"other:{cgap}"
            core_info[S] = {
                "gap": f"{cgap}", "class": cls, "R": R,
                "X_theory": xt, "X_brute": xbrute, "t1_ok": ok1,
            }
            if not ok1:
                all_fail.append(("T1", n, S, xt_lo, xb_lo,
                                 [x for x in xbrute if x > lo]))

        # T4: X-completeness per core against the zoo
        zoo_by_core = defaultdict(set)
        for (T, x, S) in pairs:
            zoo_by_core[S].add(x)
        t4_fail = []
        for S in core_list:
            xs = core_info[S]["X_theory"]
            pred = set(x for x in xs if x <= Vn)
            actual = zoo_by_core[S]
            if pred != actual:
                t4_fail.append((S, sorted(pred - actual),
                                sorted(actual - pred)))
                all_fail.append(("T4", n, S, sorted(pred - actual),
                                 sorted(actual - pred)))
        level["t4_xclosure_failures"] = len(t4_fail)

        # ---- T2/T3/T5 -----------------------------------------------------
        cls_counter = Counter(core_info[S]["class"] for S in core_list)
        level["core_classes"] = dict(cls_counter)

        naive = [T for T in Z
                 if any(core_info[tuple(sorted(v for v in T if v != x))]["class"]
                        == "tight" for x in T)]
        level["zoo_with_tight_core"] = len(naive)
        level["zoo_without_tight_core"] = len(Z) - len(naive)

        # productive tight sets: members of Z_{n-1} appearing as tight cores
        productive = sorted(S for S in core_list
                            if core_info[S]["class"] == "tight")
        level["tight_cores_appearing"] = len(productive)

        # core gap spectrum (per decomposition pair)
        spectrum = Counter()
        for (T, x, S) in pairs:
            spectrum[core_info[S]["gap"]] += 1
        level["core_gap_spectrum"] = {str(k): v for k, v in
                                      sorted(spectrum.items(),
                                             key=lambda kv: Fraction(kv[0]))}

        # ---- T6: ladder chains via tight cores ---------------------------
        prev_zoo = set()
        if n - 1 in report["levels"]:
            prev_zoo = set(load_zoo(ZOOS[n - 1][1]))
        chains = {}
        for T in Z:
            chain = [T]
            cur = T
            while len(cur) > 2:
                nxt = None
                for x in cur:
                    S = tuple(sorted(v for v in cur if v != x))
                    if core_info.get(S, {}).get("class") == "tight" or (
                            n > 2 and S in prev_zoo and len(S) == len(cur) - 1
                            and _is_tight(S, len(S))):
                        nxt = S
                        break
                if nxt is None:
                    break
                chain.append(nxt)
                cur = nxt
            chains[T] = chain
        seeds = Counter(c[-1] for c in chains.values())
        level["chain_lengths"] = dict(Counter(len(c) for c in chains.values()))
        level["seeds"] = {str(k): v for k, v in seeds.items()}
        # chains that do NOT reach a seed of size 2 are partial
        level["chains_reaching_2"] = sum(1 for c in chains.values()
                                         if len(c[-1]) == 2)

        # per-core detail for the report
        level["core_details"] = {
            ",".join(map(str, S)): {
                "gap": core_info[S]["gap"], "class": core_info[S]["class"],
                "R": core_info[S]["R"],
                "X": core_info[S]["X_theory"],
                "zoo_ext": sorted(zoo_by_core[S]),
            } for S in core_list
        }
        report["levels"][str(n)] = level

        # ---- console summary ----------------------------------------------
        print(f"== n={n}  target g=1/{n+1}  V<= {Vn}  zoo={len(Z)} "
              f"pairs={len(pairs)} distinct cores={len(core_list)}")
        print(f"   core classes: {dict(cls_counter)}")
        print(f"   naive ladder: {len(naive)}/{len(Z)} tight sets have a "
              f"tight core; {len(Z)-len(naive)} do not")
        print(f"   T1 theory-vs-brute failures: "
              f"{sum(1 for S in core_list if not core_info[S]['t1_ok'])}"
              f"   T4 X-closure failures: {len(t4_fail)}")
        print(f"   core gap spectrum: "
              f"{ {str(k): v for k, v in sorted(spectrum.items(), key=lambda kv: Fraction(kv[0]))} }")
        print(f"   chains reach n=2 seed: {level['chains_reaching_2']}/{len(Z)}"
              f"  lengths {level['chain_lengths']}")
        if n >= 4 and t4_fail:
            for f in t4_fail[:5]:
                print("   T4 FAIL:", f)
        print()

    # ---- global verdict ----------------------------------------------------
    print("=" * 70)
    if all_fail:
        print(f"FAILURES: {len(all_fail)}")
        for f in all_fail[:15]:
            print("  ", f)
        report["verdict"] = "FAILURES PRESENT"
    else:
        report["verdict"] = ("ALL PASS: X(S,1/(n+1)) certified on every "
                             "decomposition of every tight set n=2..8; "
                             "zoo is X-closed at every core.")
    print(report["verdict"])

    with open(OUT + "/tight_decomp.json", "w") as f:
        json.dump(report, f, indent=1, default=str)
    print("wrote", OUT + "/tight_decomp.json")


def _is_tight(S, m):
    gp, _ = gap_int(list(S))
    return gp == Fraction(1, m + 1)


if __name__ == "__main__":
    sys.exit(main())
