"""C-HALF PROBE: dissect every crossing in the drop experiment.

Advisor directive (this session): attack the C-half via the residue law.
  C-half conjecture: for S with gap(S) >= 1/(k+1), x killed-not-clean, and
  t* = a/b the new argmax at a crossing,  max(Q*, C*) >= 1/(k+2).
Equivalently, with m := C* * b (integer: C* = ||v a/b|| = dist_b(va)/b):
  b <= m(k+2),   i.e.  deficit d := m(k+2) - b >= 0.

Laws probed (see report):
  DEFICIT   d >= 0; joint (m, d) table; d=0 => m=1 (tight); rung d=m-1;
            d=1 inhabitants (m=2 rung-2, m=3 inter-rung).
  GRID      all argmaxes of the core share ONE denominator b_S, and
            gap(S) = m_S / b_S (m_S integer).
  OVERSHOOT no-touch (drop) => max over core argmaxes j/b_S of
            dist_{b_S}(x j) <= m_S - 1   (x is near a multiple of b_S).
  M-LAW     m vs l*m_S where l = (x - eps)/b_S, eps = signed residue.
  BINDER    the crossing's core binder v: in the core's own argmax binder
            set?  which side (rise/fall)?  sum-law vs diff-law, j.
  GAPS      residues {u a mod b} u in V, plus 0: all in [m, b-m]; binder
            residues exactly m and b-m; interior gaps average (b-2m)/k;
            C-half <=> average interior gap <= m.  Max gap G vs m.
  COVERING  exact count of i in Z_b with min_u dist_b(i u) <= m (the
            covering bound says b <= (k+1)(2m+1) -- factor 2 off target).
  C1/C2     crossings at core candidate times (value = core ladder value,
            Q* carries too) vs genuinely new times (in_cand_S false) --
            the C-half's real content is C2.

Run: python3 lrc_chalf_probe.py [--sample-every N]
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lrc_gap_lib import gap_int, _norm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RECS = os.path.join(OUT, "drop_records.jsonl")


def frac(s):
    if "/" in s:
        n, d = s.split("/")
        return Fraction(int(n), int(d))
    return Fraction(int(s))


def dist_b(r, b):
    r %= b
    return r if r < b - r else b - r


# ---------------------------------------------------------------- core cache
_CACHE = {}


def core_info(S):
    key = tuple(sorted(S))
    if key in _CACHE:
        return _CACHE[key]
    g, ts = gap_int(list(key))
    dens = sorted({t.denominator for t in ts})
    unified = len(dens) == 1
    b_S = dens[-1]
    m_S = None
    if unified:
        val = g * b_S
        assert val.denominator == 1, ("m_S not integer", key, g, b_S)
        m_S = val.numerator
    amax_binders = set()
    for t in ts:
        for v in key:
            if _norm(v * t) == g:
                amax_binders.add(v)
    # argmax indices j (numerators) for unified cores
    js = sorted({t.numerator for t in ts}) if unified else None
    ci = dict(S=key, g=g, ts=ts, dens=dens, unified=unified, b_S=b_S,
              m_S=m_S, amax_binders=amax_binders, js=js)
    _CACHE[key] = ci
    return ci


def signed_residue(x, b):
    e = x % b
    return e if 2 * e <= b else e - b


# ---------------------------------------------------------------- one entry
def dissect(rec, ci, entry, rows, errs, cover_every):
    S, x, k = ci["S"], rec["x"], rec["k"]
    V = tuple(sorted(S + (x,)))
    t = frac(entry["t"])
    a, b = t.numerator, t.denominator
    g_new = frac(rec["g_new"])
    m_val = g_new * b
    if m_val.denominator != 1:
        errs.append(("m_not_int", S, x, t))
        return
    m = m_val.numerator
    d = m * (k + 2) - b
    binders = entry["binders"]
    if x not in binders:
        errs.append(("x_not_binder", S, x, t))
        return
    # residues of ALL elements of V at t*
    res = {u: (u * a) % b for u in V}
    ok_range = all(m <= r <= b - m for r in res.values())
    if not ok_range:
        errs.append(("residue_range", S, x, t))
    r_x = res[x]
    side_x = "rise" if r_x == m and b != 2 * m else (
        "fall" if r_x == b - m else "?%d" % r_x)
    # partners and laws
    partners = []
    for v in binders:
        if v == x:
            continue
        if (x + v) % b == 0:
            law, jj = "sum", (x + v) // b
        elif x != v and abs(x - v) % b == 0:
            law, jj = "diff", abs(x - v) // b
        else:
            law, jj = "none", 0
            errs.append(("pair_law", S, x, t, v))
        partners.append((v, law, jj, v in ci["amax_binders"]))
    # gap structure of {0} u residues
    pts = sorted(set([0] + [r for r in res.values()]))
    gaps = [pts[i + 1] - pts[i] for i in range(len(pts) - 1)] + [b - pts[-1] + pts[0]]
    G = max(gaps)
    n_inter = k  # interior gaps between the k+1 residues (r_min..r_max)
    sum_inter = (b - 2 * m)  # = sum of gaps strictly between m and b-m
    avg_inter = Fraction(b - 2 * m, k) if k else None
    # exact covering count (i in Z_b with min_u dist_b(i u) <= m)
    cov = None
    if cover_every is not None and b <= 300:
        cov = 0
        for i in range(b):
            c = min(dist_b(i * u, b) for u in V)
            if c <= m:
                cov += 1
    # overshoot vs core grid
    ov = None
    if ci["unified"]:
        b_S, m_S = ci["b_S"], ci["m_S"]
        maxd = max(dist_b(x * j, b_S) for j in ci["js"])
        ov = dict(b_S=b_S, m_S=m_S, maxdist=maxd,
                  eps=signed_residue(x, b_S),
                  l=(x - signed_residue(x, b_S)) // b_S)
    rows.append(dict(
        S=list(S), x=x, k=k, aug=rec["aug"], gS=rec["gS"],
        g_new=rec["g_new"], form=rec.get("g_new_form"),
        t=entry["t"], a=a, b=b, m=m, d=int(d),
        binders=binders, partners=[[v, law, jj, inb]
                                   for v, law, jj, inb in partners],
        r_x=r_x, side_x=side_x, ok_range=ok_range,
        in_cand_S=entry.get("in_cand_S"),
        Qstar=rec["Qstar"], touch=rec["touch"],
        G=int(G), ngaps=len(gaps), gaps=gaps if d <= 2 else None,
        avg_inter_num=b - 2 * m, n_inter_gaps=n_inter,
        gaps_gt_m=sum(1 for g in gaps if g > m) if d <= 40 else None,
        cov=cov, gcds=[gcd(u, b) for u in V],
        overshoot=ov,
        m_minus_l_mS=(m - ov["l"] * ov["m_S"]) if ov else None,
    ))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-every", type=int, default=None,
                    help="compute covering count only every Nth entry")
    args = ap.parse_args()
    cover_every = args.sample_every

    rows, errs = [], []
    n_rec = n_drop = n_touch = 0
    n_cross_ent = n_q_ent = 0
    c1 = c2 = 0
    with open(RECS) as f:
        for line in f:
            rec = json.loads(line)
            n_rec += 1
            ci = core_info(rec["S"])
            if rec["touch"]:
                n_touch += 1
                continue
            n_drop += 1
            for entry in rec["amax_new"]:
                if entry["cls"] == "crossing":
                    n_cross_ent += 1
                    if entry.get("in_cand_S"):
                        c1 += 1
                    else:
                        c2 += 1
                    dissect(rec, ci, entry, rows, errs, cover_every)
                elif entry["cls"] == "qcorner":
                    n_q_ent += 1

    # ------------------------------------------------------------ aggregate
    R = []
    R.append("C-HALF PROBE  (crossing dissection of the drop experiment)")
    R.append("=" * 72)
    R.append(f"records={n_rec} drops={n_drop} touch={n_touch} "
             f"crossing_entries={n_cross_ent} (C1={c1} core-cand, C2={c2} new) "
             f"qcorner_entries={n_q_ent}")
    R.append(f"errors: {len(errs)}  {errs[:10]}")
    R.append("")

    # DEFICIT law
    neg = [r for r in rows if r["d"] < 0]
    R.append("[DEFICIT] d = m(k+2) - b  (C-half: d >= 0)")
    R.append(f"  violations (d<0): {len(neg)}")
    by_aug = defaultdict(list)
    for r in rows:
        by_aug[r["aug"]].append(r["d"])
    for aug in sorted(by_aug):
        ds = by_aug[aug]
        R.append(f"  aug={aug}: n={len(ds)} min_d={min(ds)} max_d={max(ds)} "
                 f"d=0:{sum(1 for x in ds if x == 0)} "
                 f"d=1:{sum(1 for x in ds if x == 1)}")
    d0 = [r for r in rows if r["d"] == 0]
    R.append(f"  d=0 cases: {len(d0)}; of these m=1: "
             f"{sum(1 for r in d0 if r['m'] == 1)}, m>1: "
             f"{Counter(r['m'] for r in d0 if r['m'] > 1)}")
    d1 = [r for r in rows if r["d"] == 1]
    R.append(f"  d=1 cases: {len(d1)}; (m,aug) histogram: "
             f"{sorted(Counter((r['m'], r['aug']) for r in d1).items())[:20]}")
    # rung-form check
    rung = [r for r in rows if r["form"] == "rung"]
    rung_ok = sum(1 for r in rung if r["d"] == r["m"] - 1)
    R.append(f"  rung-form entries: {len(rung)}, of which d == m-1: {rung_ok}")
    md = sorted(set((r["m"], r["d"]) for r in rows))
    R.append(f"  observed (m,d) pairs: {md[:40]}{' ...' if len(md) > 40 else ''}")
    R.append("")

    # GRID law (over ALL cores seen, incl. touch records)
    uni = sum(1 for c in _CACHE.values() if c["unified"])
    R.append(f"[GRID] cores={len(_CACHE)} unified-argmax-denominator={uni} "
             f"exceptions={[c['S'] for c in _CACHE.values() if not c['unified']][:12]}")
    R.append("")

    # OVERSHOOT law
    ov_rows = [r for r in rows if r["overshoot"]]
    bad_ov = [r for r in ov_rows if r["overshoot"]["maxdist"] > r["overshoot"]["m_S"] - 1]
    R.append("[OVERSHOOT] drop => max_j dist_{b_S}(x j) <= m_S - 1")
    R.append(f"  unified-core crossing rows: {len(ov_rows)}, violations: {len(bad_ov)}")
    if bad_ov:
        for r in bad_ov[:6]:
            R.append(f"    VIOLATION S={r['S']} x={r['x']} "
                     f"b_S={r['overshoot']['b_S']} m_S={r['overshoot']['m_S']} "
                     f"maxdist={r['overshoot']['maxdist']}")
    eps_hist = Counter(abs(r["overshoot"]["eps"]) for r in ov_rows)
    R.append(f"  |eps| histogram (signed x mod b_S): {sorted(eps_hist.items())[:15]}")
    mS_hist = Counter(r["overshoot"]["m_S"] for r in ov_rows)
    R.append(f"  m_S histogram: {sorted(mS_hist.items())}")
    bad_mS = [r for r in ov_rows if abs(r["overshoot"]["eps"]) > r["overshoot"]["m_S"] - 1]
    R.append(f"  |eps| > m_S-1 violations: {len(bad_mS)}")
    if bad_mS:
        for r in bad_mS[:6]:
            R.append(f"    VIOLATION S={r['S']} x={r['x']} eps={r['overshoot']['eps']} "
                     f"m_S={r['overshoot']['m_S']} b_S={r['overshoot']['b_S']} "
                     f"maxdist={r['overshoot']['maxdist']}")
    R.append("")

    # M-LAW
    mm = Counter(r["m_minus_l_mS"] for r in ov_rows)
    R.append("[M-LAW] m - l*m_S  (l = (x-eps)/b_S)")
    R.append(f"  histogram: {sorted(mm.items())[:25]}")
    R.append("")

    # BINDER law
    allp = [p for r in rows for p in r["partners"]]
    R.append("[BINDER] crossing partner pairs")
    R.append(f"  total partner pairs: {len(allp)}")
    R.append(f"  law: {Counter(p[1] for p in allp)}")
    R.append(f"  j histogram: {sorted(Counter(p[2] for p in allp).items())[:15]}")
    R.append(f"  partner in core argmax binders: "
             f"{sum(1 for p in allp if p[3])}/{len(allp)}")
    R.append(f"  side of x: {Counter(r['side_x'] for r in rows)}")
    R.append("")

    # GAPS law
    R.append("[GAPS] residue structure at the crossing")
    big = [(r["G"] - r["m"], r) for r in rows if r["gaps_gt_m"] is not None]
    big.sort(key=lambda z: -z[0])
    R.append(f"  max(G - m) over rows: {big[0][0] if big else '-'} "
             f"(S={big[0][1]['S']} x={big[0][1]['x']} m={big[0][1]['m']} "
             f"b={big[0][1]['b']})" if big else "")
    gm = Counter(r["gaps_gt_m"] for r in rows if r["gaps_gt_m"] is not None)
    R.append(f"  #interior gaps > m, histogram: {sorted(gm.items())}")
    # avg interior gap vs m  <=> d>=0 (identity check)
    R.append(f"  (b-2m)/k > m count (= d<0): {sum(1 for r in rows if r['avg_inter_num'] > r['m'] * r['k'])}")
    R.append("")

    # COVERING
    covs = [r for r in rows if r["cov"] is not None]
    if covs:
        R.append("[COVERING] exact |{i in Z_b: min_u dist_b(iu) <= m}| (target b)")
        ratios = [(r["cov"], r["b"], r["k"], r["m"]) for r in covs]
        full = sum(1 for c, bb, _, _ in ratios if c >= bb)
        R.append(f"  computed for {len(covs)} rows; cover-all: {full}")
        short = [(bb - c, r) for r, (c, bb, _, _) in zip(covs, ratios)]
        R.append(f"  max shortfall b - cov: {max(bb - c for c, bb, _, _ in ratios)}")
    R.append("")

    # gcd structure
    gmulti = sum(1 for r in rows if any(g > 1 for g in r["gcds"]))
    R.append(f"[GCD] rows with some gcd(u,b)>1: {gmulti}/{len(rows)}")

    # C1/C2
    R.append("")
    R.append("[C1/C2] crossing entries at core-candidate times vs new times")
    for lab, sel in (("C1 (core cand)", True), ("C2 (new)", False)):
        rs = [r for r in rows if r["in_cand_S"] == sel]
        if not rs:
            continue
        R.append(f"  {lab}: n={len(rs)} min_d={min(r['d'] for r in rs)} "
                 f"d=0:{sum(1 for r in rs if r['d'] == 0)} "
                 f"d=1:{sum(1 for r in rs if r['d'] == 1)}")
        qcar = sum(1 for r in rs if frac(r["Qstar"]) == frac(r["g_new"]))
        R.append(f"    Qstar == g_new (ladder value): {qcar}/{len(rs)}")

    report = "\n".join(R)
    print(report)
    with open(os.path.join(OUT, "chalf_probe_report.txt"), "w") as f:
        f.write(report + "\n")
    with open(os.path.join(OUT, "chalf_rows.json"), "w") as f:
        json.dump(rows, f)
    print(f"\nwrote {len(rows)} rows -> out/chalf_rows.json")


if __name__ == "__main__":
    main()
