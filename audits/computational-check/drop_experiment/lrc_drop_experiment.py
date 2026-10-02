"""THE DROP EXPERIMENT (the advisor's question, mapped to the X(S) framework).

Question.  Let S be a core with gap g = gap(S), and let x be a filler that is
killed by S (trivially: at g = gap(S) every x kills, since S has no escapes
above its own gap) but NOT clean at any level-g time (i.e. x does not touch
the argmax of S: ||x t*|| < g for every argmax time t* of S).
What is gap(S u {x})?  Is the drop bounded?  Is there a formula?
Does it depend only on S, or on x too?

Theorem D (drop formula), verified here machine-exactly:
  gap(S u {x}) <= gap(S), with equality iff touch.
  In a drop case, every t_new in Argmax(S u {x}) is exactly one of
    (C) CROSSING : ||x t_new|| = q(t_new) = gap(S u {x}),  q = min_{v in S} ||v t||
    (Q) Q-CORNER : q(t_new) = gap(S u {x}) < ||x t_new||, and t_new is then
                   automatically a local max of q  [ proved: p > q on a nbhd,
                   h = min(p,q) = q there, global max => local max of q ].
    (P) P-CORNER : ||x t_new|| = 1/2 = gap(S u {x})  -- impossible in a drop
                   case (would force no drop).
  Hence  gap(S u {x}) = max( Q*, C* )  with
    Q* = max{ q(t) : t a local-max candidate of q, ||x t|| >= q(t) }   (S's OWN ladder)
    C* = best crossing value                                             (mixed arithmetic).

Records, per (S, x):
  touch / drop; g_S, g_new, drop; classification of every new argmax time;
  Q* (local-max restricted) and Q_any (all candidates); LRC margin
  g_new - 1/(k+2); closed-form flags of g_new (tight 1/m, rung a/(an+1));
  binding elements and binding-pair sums mod the new-argmax denominator.

Stages:  complete (k=1..4 exhaustive small universes)  |  zoo (X-relevant
cores from the NSK census files, k=2..7)  |  analyze.
Idempotent: (S,x) pairs already in the JSONL are skipped.
Run: python3 lrc_drop_experiment.py --stage complete --shard 0/4 ... etc.
"""
import argparse
import json
import os
import sys
from fractions import Fraction
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc_gap_lib import gap_int, gap_frac, _norm, _candidates_int  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RECS = os.path.join(OUT, "drop_records.jsonl")
NSK = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "download", "lonely-runner-nsk-census")


def fstr(fr):
    return f"{fr.numerator}/{fr.denominator}" if fr.denominator != 1 else str(fr.numerator)


# ---------------------------------------------------------------- core cache
def _full_grid(S):
    """ALL breakpoints of q = min_{v in S} ||v t||: profile breakpoints
    j/(2v) (peaks AND valleys) plus pairwise profile crossings a/(v_i+-v_j).
    Between consecutive grid points every profile is linear and the binding
    set is constant, so q is linear there; hence every local max of q is a
    grid point.  (The gap solver's candidate list omits valleys -- fine for
    the global max, NOT for local-max detection.)"""
    pts = {Fraction(0, 1), Fraction(1, 1)}
    for v in S:
        for j in range(0, 2 * v + 1):
            pts.add(Fraction(j, 2 * v))
    n = len(S)
    for i in range(n):
        for j in range(i + 1, n):
            s_ = S[i] + S[j]
            for a in range(0, s_ + 1):
                pts.add(Fraction(a, s_))
            d_ = S[j] - S[i]
            if d_ > 0:
                for a in range(0, d_ + 1):
                    pts.add(Fraction(a, d_))
    return sorted(pts)


class CoreInfo:
    """Precomputed per-core data: gap, argmaxes, full breakpoint grid with q
    values and local-max flags (q linear between consecutive grid points)."""

    def __init__(self, S):
        self.S = tuple(S)
        self.k = len(S)
        g, ts = gap_int(list(S))
        self.g = g
        self.argmax = ts
        cands = _full_grid(self.S)
        qs = [min(_norm(v * t) for v in S) for t in cands]
        lmax = []
        for i in range(len(cands)):
            lo = qs[i - 1] if i > 0 else Fraction(0)
            hi = qs[i + 1] if i + 1 < len(cands) else Fraction(0)
            # q linear between consecutive grid points => grid local max is a
            # true local max (non-constant pieces => adjacent values differ)
            lmax.append(qs[i] >= lo and qs[i] >= hi and qs[i] > 0)
        self.cands = cands
        self.qs = qs
        self.lmax = lmax
        self.candset = set(cands)


_CACHE = {}


def core_info(S):
    key = tuple(sorted(S))
    if key not in _CACHE:
        _CACHE[key] = CoreInfo(key)
    return _CACHE[key]


# ---------------------------------------------------------------- one pair
def classify(v2, ci, x):
    """Full record for the pair (S, x).  v2 = sorted(S + [x])."""
    S, k, g_S = ci.S, ci.k, ci.g
    g_new, ts_new = gap_int(v2)
    # touch: some argmax time of S cleared by x
    touch = any(_norm(x * t) >= g_S for t in ci.argmax)
    rec = {
        "S": list(S), "x": x, "k": k, "aug": k + 1,
        "gS": fstr(g_S), "g_new": fstr(g_new),
        "drop": fstr(g_S - g_new), "touch": touch,
    }
    # Q* and Q_any (S's own ladder cleared by x)
    qstar = Fraction(0)
    qany = Fraction(0)
    for t, qv, lm in zip(ci.cands, ci.qs, ci.lmax):
        if _norm(x * t) >= qv:
            if lm and qv > qstar:
                qstar = qv
            if qv > qany:
                qany = qv
    rec["Qstar"] = fstr(qstar)
    rec["Qany"] = fstr(qany)
    # LRC margin for the augmented set (size k+1): g_new - 1/(k+2)
    floor = Fraction(1, k + 2)
    margin = g_new - floor
    rec["floor"] = fstr(floor)
    rec["margin"] = fstr(margin)
    rec["margin_flag"] = ("floor" if margin == 0 else
                          "above" if margin > 0 else "SILENT")
    # closed-form flags of g_new
    a, b = g_new.numerator, g_new.denominator
    rec["g_new_form"] = ("tight" if a == 1 and b >= 2 else
                         "rung" if a >= 2 and (b - 1) % a == 0 and (b - 1) // a >= 2
                         else "other")
    if rec["g_new_form"] == "rung":
        rec["rung_n"] = (b - 1) // a
    # new-argmax classification
    cls_list = []
    for t in ts_new:
        p = _norm(x * t)
        q = min(_norm(v * t) for v in S)
        binders = [v for v in v2 if _norm(v * t) == g_new]
        entry = {"t": fstr(t), "p": fstr(p), "q": fstr(q),
                 "binders": binders,
                 "in_cand_S": t in ci.candset}
        if p == q:
            entry["cls"] = "crossing"
            # Lemma-A-like residue structure at the new argmax
            den = t.denominator
            bp = []
            for i in range(len(binders)):
                for j in range(i + 1, len(binders)):
                    s_, d_ = binders[i] + binders[j], abs(binders[j] - binders[i])
                    bp.append((s_ % den == 0) or (d_ % den == 0 and d_ > 0))
            entry["bpair_mod_den"] = bp
        elif q == g_new:
            entry["cls"] = "qcorner"
            entry["is_cand_S"] = t in ci.candset
        elif p == g_new:
            entry["cls"] = "pcorner"
        else:
            entry["cls"] = "???"
        cls_list.append(entry)
    rec["amax_new"] = cls_list
    rec["cls"] = (cls_list[0]["cls"] if len({e["cls"] for e in cls_list}) == 1
                  else "mixed")
    return rec, g_new


def check(rec, ci, g_new):
    """Machine-verify Theorem D on this record.  Returns list of violations."""
    bad = []
    g_S = ci.g
    # T1: touch <=> drop 0
    if rec["touch"] != (g_new == g_S):
        bad.append("T1")
    if not rec["touch"]:
        for e in rec["amax_new"]:
            if e["cls"] == "qcorner":
                # q-binding => t must be a candidate of S (local max of q)
                if not e["in_cand_S"]:
                    bad.append("Qcand")
                # Q* must equal g_new when a q-corner attains the max
                if Fraction(*map(int, rec["Qstar"].split("/"))) != g_new:
                    bad.append("Qstar")
            elif e["cls"] == "pcorner":
                if g_new != Fraction(1, 2):
                    bad.append("Pcorner")
            elif e["cls"] == "???":
                bad.append("cls")
        if Fraction(*map(int, rec["Qany"].split("/"))) > g_new:
            bad.append("Qany>g")
    return bad


# ---------------------------------------------------------------- universes
def load_done():
    done = set()
    if os.path.exists(RECS):
        with open(RECS) as f:
            for line in f:
                try:
                    r = json.loads(line)
                    done.add((tuple(r["S"]), r["x"]))
                except Exception:
                    pass
    return done


def pairs_complete():
    out = []
    for v in range(1, 25):                                   # k = 1
        out.extend(([v], x) for x in range(1, 41) if x != v)
    from itertools import combinations
    for S in combinations(range(1, 15), 2):                  # k = 2
        out.extend((list(S), x) for x in range(1, 41) if x not in S)
    for S in combinations(range(1, 12), 3):                  # k = 3
        out.extend((list(S), x) for x in range(1, 37) if x not in S)
    for S in combinations(range(1, 11), 4):                  # k = 4
        out.extend((list(S), x) for x in range(1, 31) if x not in S)
    return out


def pairs_zoo():
    """X-relevant cores (productive cores of the NSK census), primitive,
    elements <= 60 (<= 80 for k=7); x in [1..60] (extensions for k=7)."""
    out, seen = [], set()
    import re
    for m in range(3, 9):  # files cores_m3..cores_m8  (core size k = m-1)
        fn = os.path.join(NSK, f"cores_m{m}_B{{}}.txt")
        import glob
        g = glob.glob(os.path.join(NSK, f"cores_m{m}_B*.txt"))
        if not g:
            continue
        for line in open(g[0]):
            mm = re.match(r"S=([0-9,]+),", line.strip())
            if not mm:
                continue
            S = tuple(int(t) for t in mm.group(1).split(","))
            k = len(S)
            if k <= 1 or k >= 8:
                continue
            cap = 60 if k <= 6 else 80
            if max(S) > cap:
                continue
            if gcd(*S) != 1:
                continue
            if S in seen:
                continue
            seen.add(S)
            if k <= 6:
                xs = [x for x in range(1, 61) if x not in S]
            else:
                xs = [x for x in range(max(S) + 1, 81) if x not in S]
            out.extend((list(S), x) for x in xs)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["complete", "zoo", "count"])
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--consensus", type=int, default=97,
                    help="consensus-check every Nth pair with gap_frac")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    if args.stage == "count":
        pc, pz = pairs_complete(), pairs_zoo()
        print(f"complete: {len(pc)} pairs | zoo: {len(pz)} pairs")
        return

    pairs = (pairs_complete() if args.stage == "complete" else pairs_zoo())
    si, sn = map(int, args.shard.split("/"))
    pairs = [p for idx, p in enumerate(pairs) if idx % sn == si]
    done = load_done()
    pairs = [p for p in pairs if (tuple(sorted(p[0])), p[1]) not in done]
    print(f"stage={args.stage} shard={args.shard}: {len(pairs)} pairs to do "
          f"({len(done)} already done)", flush=True)

    violations, n_done, n_touch, n_drop = [], 0, 0, 0
    tCHECK = 0
    with open(RECS, "a") as f:
        for S, x in pairs:
            ci = core_info(S)
            v2 = sorted(S + [x])
            rec, g_new = classify(v2, ci, x)
            bad = check(rec, ci, g_new)
            if bad:
                violations.append((S, x, bad))
            n_done += 1
            if rec["touch"]:
                n_touch += 1
            else:
                n_drop += 1
            # consensus with the independent Fraction implementation
            tCHECK += 1
            if tCHECK % args.consensus == 0 or (not rec["touch"] and
                                                 rec["margin_flag"] == "floor"):
                g2, t2 = gap_frac(v2)
                t1 = [Fraction(*map(int, e["t"].split("/"))) for e in rec["amax_new"]]
                if g2 != g_new or t2 != t1:
                    violations.append((S, x, ["CONSENSUS"]))
            f.write(json.dumps(rec) + "\n")
            if n_done % 500 == 0:
                f.flush()
                print(f"  {n_done}/{len(pairs)} touch={n_touch} drop={n_drop}",
                      flush=True)
    print(f"DONE stage={args.stage} shard={args.shard}: {n_done} pairs, "
          f"touch={n_touch}, drop={n_drop}")
    if violations:
        print(f"VIOLATIONS ({len(violations)}):")
        for v in violations[:20]:
            print("  ", v)
        sys.exit(1)
    print("ALL THEOREM-D CHECKS PASS")


if __name__ == "__main__":
    main()
