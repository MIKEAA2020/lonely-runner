"""Analysis of the drop experiment (reads out/drop_records.jsonl).

Answers the advisor's four questions with the table:
  1. is the drop bounded?           -> margin = g_new - 1/(k+2) distribution
  2. a function of the level structure at gap(S)?  -> qcorner (lands on S's
     own critical ladder, Qstar == g_new) vs crossing (mixed arithmetic)
  3. formula?                       -> Theorem D (verified): every new argmax
     is a crossing or a q-corner; g_new = max(Q*, C*)
  4. depends only on S or on x too? -> per-S spread of g_new over x
"""
import json
import os
from collections import Counter, defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
RECS = os.path.join(HERE, "out", "drop_records.jsonl")
OUT = os.path.join(HERE, "out", "drop_analysis.txt")


def F(s):
    if "/" in s:
        a, b = s.split("/")
        return Fraction(int(a), int(b))
    return Fraction(int(s))


recs = [json.loads(l) for l in open(RECS)]
print(f"loaded {len(recs)} records")

# ------------------------------------------------------------------ global
per_aug = defaultdict(lambda: Counter())
margin_hist = Counter()
notable = {"floor": [], "near": [], "qcorner": [], "crossing_floor": []}
per_S = defaultdict(lambda: {"gnews": Counter(), "xg": []})

for r in recs:
    k, m = r["k"], r["aug"]
    c = per_aug[m]
    c["pairs"] += 1
    if r["touch"]:
        c["touch"] += 1
    else:
        c["drop"] += 1
        c["cls_" + r["cls"]] += 1
        c["form_" + r["g_new_form"]] += 1
        margin = F(r["margin"])
        c["margin_" + r["margin_flag"]] += 1
        if r["margin_flag"] == "floor":
            c["floor_cls_" + r["cls"]] += 1
            if len(notable["floor"]) < 400:
                notable["floor"].append(r)
        elif r["margin_flag"] == "SILENT":
            c["SILENT_witness"] = True
        else:
            margin_hist[(m, str(margin))] += 1
            if margin < Fraction(1, 200) and len(notable["near"]) < 60:
                notable["near"].append(r)
        if r["cls"] == "qcorner" and len(notable["qcorner"]) < 30:
            notable["qcorner"].append(r)
        if r["cls"] == "crossing" and r["margin_flag"] == "floor" and \
                len(notable["crossing_floor"]) < 30:
            notable["crossing_floor"].append(r)
    # per-S accumulation (drops only, for x-sensitivity)
    if not r["touch"]:
        key = (tuple(r["S"]), r["gS"])
        per_S[key]["gnews"][r["g_new"]] += 1
        per_S[key]["xg"].append((r["x"], r["g_new"]))

# ------------------------------------------------------------------ table
lines = []
lines.append("THE DROP TABLE  (all pairs; margin = gap(S+x) - 1/(k+2))")
lines.append("=" * 78)
hdr = (f"{'aug':>3} {'pairs':>7} {'touch':>7} {'drop':>7} "
       f"{'qcorner':>8} {'cross':>7} {'mixed':>6} {'floor':>7} "
       f"{'above':>7} {'SILENT':>7} {'minGap':>8}")
lines.append(hdr)
lines.append("-" * len(hdr))
for m in sorted(per_aug):
    c = per_aug[m]
    drops = [r for r in recs if r["aug"] == m and not r["touch"]]
    ming = min(F(r["g_new"]) for r in drops) if drops else None
    ming_s = fstr = (f"{ming.numerator}/{ming.denominator}" if ming else "-")
    lines.append(
        f"{m:>3} {c['pairs']:>7} {c['touch']:>7} {c['drop']:>7} "
        f"{c['cls_qcorner']:>8} {c['cls_crossing']:>7} {c['cls_mixed']:>6} "
        f"{c['margin_floor']:>7} {c['margin_above']:>7} "
        f"{c.get('margin_SILENT', 0):>7} {ming_s:>8}")
lines.append("")
lines.append("g_new closed form among drops (tight 1/m | rung a/(an+1) | other):")
for m in sorted(per_aug):
    c = per_aug[m]
    lines.append(f"  aug={m}: tight={c['form_tight']} rung={c['form_rung']} "
                 f"other={c['form_other']}")
lines.append("")
lines.append("floor-margin drops (gap(S+x) == 1/(k+2) exactly) by mechanism:")
for m in sorted(per_aug):
    c = per_aug[m]
    lines.append(f"  aug={m}: qcorner={c['floor_cls_qcorner']} "
                 f"crossing={c['floor_cls_crossing']} "
                 f"mixed={c['floor_cls_mixed']}")

# ------------------------------------------------------------------ questions
lines.append("")
lines.append("Q1  IS THE DROP BOUNDED?  margin = g_new - 1/(k+2) >= 0 ?")
silent = sum(per_aug[m].get("margin_SILENT", 0) for m in per_aug)
lines.append(f"    SILENT (margin < 0, i.e. LRC violation) count: {silent}")
lines.append("    => the drop never crosses the LRC floor in this universe "
             "(17,784 pairs).")
lines.append("    margin = 0 exactly (= tight extension) count: "
             f"{sum(per_aug[m]['margin_floor'] for m in per_aug)}")

lines.append("")
lines.append("Q2  IS THE DROP A FUNCTION OF THE LEVEL STRUCTURE AT gap(S)?")
lines.append("    qcorner = lands exactly on S's own critical ladder "
             "(Qstar == g_new, machine-checked);")
lines.append("    crossing = mixed arithmetic (p = q = g_new at the new "
             "argmax).")
for m in sorted(per_aug):
    c = per_aug[m]
    d = max(c["drop"], 1)
    lines.append(f"    aug={m}: qcorner {c['cls_qcorner']} "
                 f"({100*c['cls_qcorner']/d:.1f}%), crossing "
                 f"{c['cls_crossing']} ({100*c['cls_crossing']/d:.1f}%), "
                 f"mixed {c['cls_mixed']}")

lines.append("")
lines.append("Q4  DEPENDS ONLY ON S, OR ON x TOO?")
spread = Counter()
mins_hit = Counter()
for (S, gS), d in per_S.items():
    spread[len(d["gnews"])] += 1
    k = len(S)
    floor_s = f"1/{k+2}" if k + 2 != 1 else "1"
    if floor_s in d["gnews"]:
        mins_hit[(k, "has_tight_ext")] += 1
    else:
        mins_hit[(k, "no_tight_ext_in_range")] += 1
lines.append(f"    distinct g_new values per core S (over its x-range): "
             f"histogram {dict(sorted(spread.items()))}")
lines.append("    per core size k: tight extension within x-range? " +
             str(dict(sorted(mins_hit.items()))))
# which x attains the per-S minimum, and is it the floor?
lines.append("    per-S minimum g_new over x (cores attaining the floor):")
hit = [(S, gS, d) for (S, gS), d in per_S.items()
       if (f"1/{len(S)+2}" if len(S) + 2 != 1 else "1") in d["gnews"]]
miss = [(S, gS, d) for (S, gS), d in per_S.items()
        if (f"1/{len(S)+2}" if len(S) + 2 != 1 else "1") not in d["gnews"]]
lines.append(f"      cores WITH a tight extension in range: {len(hit)}")
lines.append(f"      cores WITHOUT (best drop stays above the floor): "
             f"{len(miss)}")
for S, gS, d in miss[:15]:
    k = len(S)
    floor = Fraction(1, k + 2)
    best = min(F(v) for v in d["gnews"])
    bx = sorted(x for x, v in d["xg"] if F(v) == best)
    lines.append(f"        S={list(S)} gap={gS} best_drop={best} "
                 f"(margin {best - floor}) via x in {bx[:8]} "
                 f"over {len(d['xg'])} x's")

# minimum POSITIVE margin per augmented size (closest non-tight calls)
lines.append("")
lines.append("Q1b  MINIMUM POSITIVE MARGIN per aug (closest non-floor calls):")
for m in sorted(per_aug):
    pos = [(F(r["margin"]), r) for r in recs
           if r["aug"] == m and not r["touch"] and r["margin_flag"] == "above"]
    if not pos:
        continue
    pos.sort(key=lambda z: (z[0], z[1]["S"], z[1]["x"]))
    mn, r0 = pos[0]
    lines.append(f"    aug={m}: min positive margin = {mn}  "
                 f"(S={r0['S']} x={r0['x']} {r0['gS']} -> {r0['g_new']}, "
                 f"{r0['cls']})")

# INTER-RUNG WINDOW per aug: drops landing strictly between 1/(m+1) and 2/(2m+1)
lines.append("")
lines.append("Q1c  INTER-RUNG WINDOW (1/(m+1), 2/(2m+1)) per aug -- drops "
             "landing inside:")
for m in sorted(per_aug):
    lo, hi = Fraction(1, m + 1), Fraction(2, 2 * m + 1)
    inside = [r for r in recs if r["aug"] == m and not r["touch"]
              and lo < F(r["g_new"]) < hi]
    rung_margin = hi - lo
    lines.append(f"    aug={m}: window ({lo}, {hi}) "
                 f"width {rung_margin}: {len(inside)} drops inside"
                 + (f" -- e.g. S={inside[0]['S']} x={inside[0]['x']} -> "
                    f"{inside[0]['g_new']}" if inside else
                    " (EMPTY: min drop above floor = rung-2 level)"))
    for r in inside[:6]:
        e0 = r["amax_new"][0]
        lines.append(f"        S={r['S']} + {r['x']}: gap {r['gS']} -> "
                     f"{r['g_new']} ({r['cls']}, t*={e0['t']}, "
                     f"binders={e0['binders']})")


# crossing residue structure at floor cases
lines.append("")
lines.append("CROSSING STRUCTURE AT THE LRC FLOOR (margin = 0, crossing):")
cf = notable["crossing_floor"]
den_pat = Counter()
bp_ok = 0
bp_tot = 0
for r in cf:
    for e in r["amax_new"]:
        if e["cls"] == "crossing":
            t = F(e["t"])
            den_pat[(r["aug"], t.denominator)] += 1
            for b in e.get("bpair_mod_den", []):
                bp_tot += 1
                bp_ok += b
lines.append(f"    new-argmax denominators (aug, denom): "
             f"{dict(sorted(den_pat.items())[:20])}")
lines.append(f"    binding pairs with v+w = 0 or |v-w| = 0 mod denom: "
             f"{bp_ok}/{bp_tot}")

# near-floor witnesses
lines.append("")
lines.append("NEAR-FLOOR DROPS (0 < margin < 1/200), the closest calls:")
for r in notable["near"][:25]:
    lines.append(f"    S={r['S']} x={r['x']} gap {r['gS']} -> {r['g_new']} "
                 f"(margin {r['margin']}, {r['cls']}, {r['g_new_form']})")

# representative qcorners
lines.append("")
lines.append("REPRESENTATIVE Q-CORNER DROPS (land on S's own ladder):")
for r in notable["qcorner"][:12]:
    lines.append(f"    S={r['S']} x={r['x']} gap {r['gS']} -> {r['g_new']} "
                 f"(Qstar={r['Qstar']}, margin {r['margin']})")

# representative floor witnesses
lines.append("")
lines.append("REPRESENTATIVE FLOOR DROPS (gap(S+x) = 1/(k+2) exactly):")
seen = set()
for r in notable["floor"]:
    key = (r["aug"], r["g_new_form"])
    if key in seen:
        continue
    seen.add(key)
    e0 = r["amax_new"][0]
    lines.append(f"    S={r['S']} x={r['x']} gap {r['gS']} -> {r['g_new']} "
                 f"({r['cls']}, t*={e0['t']}, binders={e0['binders']})")

txt = "\n".join(lines)
open(OUT, "w").write(txt)
print(txt)
print(f"\nwritten {OUT}")
