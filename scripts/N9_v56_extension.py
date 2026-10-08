#!/usr/bin/env python3
"""
N9_v56_extension.py — THE V=56 DATA POINT: the reviewer's decisive
experiment ("Running one more V-step (V=56) would distinguish 'stable
coverage with growing classification burden' from 'fundamental decay'").
Extends BOTH decay series (same battery AND augmented battery) one step
past the committed V=48 point.  HONEST READING ENFORCED (reviewer):
V=48's augmented 80.21% was a MARGINAL pass of the 80% bar by 0.21
points, and four data points cannot distinguish deceleration from
roughly-linear residual growth with noise — so this run prints the
residual series and its increments explicitly, and the decision logic
reports the reading the numbers support, not a confidence story.

BATTERY: byte-identical to the committed N6 (same battery: P4 at
{7,13,17}), N7 (augmented battery: P4 at {7,13,17,19,37}) and N8
(V=48, both batteries) runs — the ONLY change is range, not battery:
  - tables / capability / signature-capability maps for ALL N <= 111
    (pair sums reach 2*56-1 = 111);
  - capability scan of the genuinely new moduli 96..111 (odd composites
    99, 105, 111; primes 97, 101, 103, 107, 109 are inside the committed
    (47,150] scan and are re-derived as controls, NOT new scan);
  - odd-composite structure table + independent set-based cross-check in
    the new range;
  - flag=>cert soundness sweeps extended to 96..111 (random) with
    exhaustive checks at the new-range endpoints 99 and 111.

HARD CONTROLS (asserts, all vs committed logs):
  - n=4 positive control (v<=16): 1651/1661/1744, OPEN 1, open set
    (3,5,8,13)  [committed n5/n6 record]
  - same battery: V=16 4179/4193/4311 OPEN 0; V=24 35778/36362/41656
    OPEN 0 with V-tier 4878+416; V=32 156461/161596/196751 OPEN 0;
    V=40 480454/500280/641166 OPEN 0; V=48 1200777/1253258/1665356
    OPEN 0  [scripts/n6_run.log, scripts/n8_run.log]
  - augmented battery: V=24 38611/38924/41656; V=32 173337/176376/196751;
    V=40 526607/540415/641166; V=48 1295478/1335857/1665356, OPEN 0,
    self-check unflagged allinv @19/37 == 0 everywhere
    [scripts/n7_run.log, scripts/n8_run.log]
  - capability: n=4 {7,11,13}; n=5 N<=47 {7,13,17,19,37} with counts
    20/68/16/64/32; capable in 48..95 == []  [committed N6/N8]
  - sigcap: 48..79 == committed N6 dict; 80..95 == committed N8 dict
  - V=48 unflagged decompositions: same 1041235 = 358745 P4-target +
    532700 @N<=47 + 137484 @48..79 + 12306 @80..95; augmented 682490 =
    532700 + 137484 + 12306  [scripts/n8_run.log]
  - V=56 corpus size: 3712576 primitive 5-sets (Mobius count; the same
    formula reproduces the committed V=48 total 1665356 exactly)

Deterministic.  Runtime ~24 min total, split into SIX FOREGROUND STAGES
(this environment reaps background processes between tool calls, so each
stage must run in the foreground under the 10-minute budget; per-corpus
results persist in n9_state.json):
  stage 1: prelude + soundness sweeps + controls (n=4, V=16..40 same)
  stage 2: V=48 same control + augmented controls V=24..40
  stage 3: SAME battery PRIMARY V=56
  stage 4: V=48 augmented control
  stage 5: AUGMENTED battery PRIMARY V=56
  stage 6: trend summary + the honest reading (requires all tags)
Run: python3 -u N9_v56_extension.py <1..6|all>   (append per-stage output
to n9_run.log)
"""
from math import gcd
from itertools import combinations, combinations_with_replacement
from collections import Counter
import json
import random
import sys
import time

RIGID_PRIMES = (7, 13, 17, 19, 37)              # committed n=5 rigid zone
SAME_P4 = (7, 13, 17)                           # N6 battery
AUG_P4 = (7, 13, 17, 19, 37)                    # N7 battery
N4_P4 = (7, 11, 13)                             # n=4 battery

# ------------------------------------------------------------------ tables

def build_tables(N, T):
    full = (1 << N) - 1
    masks = [0] * N
    sizes = [0] * N
    orders = [0] * N
    invs = [False] * N
    for w in range(N):
        m = 0
        for k in range(N):
            r = (w * k) % N
            d = r if r <= N - r else N - r
            if T * d < N:
                m |= 1 << k
        masks[w] = m
        sizes[w] = bin(m).count('1')
        orders[w] = N // gcd(w, N)
        invs[w] = (w != 0) and gcd(w, N) == 1
    return dict(masks=masks, sizes=sizes, orders=orders, invs=invs, full=full)

def check_size_formula(N, T, tab):
    """Lemma S5: |B_w| = (N/ord(w)) * (2*floor((ord(w)-1)/T) + 1)."""
    bad = []
    for w in range(N):
        j = tab['orders'][w]
        if (N // j) * (2 * ((j - 1) // T) + 1) != tab['sizes'][w]:
            bad.append((N, w, j))
    return bad

# ------------------------------------------------------- covering scanners

def capable_units(N, tab, k):
    """All scale-normalized unit k-multisets (1, b2, .., bk) covering Z_N."""
    masks, full = tab['masks'], tab['full']
    units = [w for w in range(1, N) if tab['invs'][w]]
    covers = []
    m1 = masks[1]
    idx = list(range(len(units)))
    if k == 3:
        for i in idx:
            ma = m1 | masks[units[i]]
            for jj in idx[i:]:
                if ma | masks[units[jj]] == full:
                    covers.append((1, units[i], units[jj]))
    elif k == 4:
        for i in idx:
            ma = m1 | masks[units[i]]
            for jj in idx[i:]:
                mb = ma | masks[units[jj]]
                for kk in idx[jj:]:
                    if mb | masks[units[kk]] == full:
                        covers.append((1, units[i], units[jj], units[kk]))
    return covers

def signature_capability(N, tab):
    """order-signature -> example covering 4-multiset (nonzero residues)."""
    masks, full, orders = tab['masks'], tab['full'], tab['orders']
    cap = {}
    res = list(range(1, N))
    idx = list(range(len(res)))
    for i in idx:
        mi = masks[res[i]]
        oi = orders[res[i]]
        for jj in idx[i:]:
            mij = mi | masks[res[jj]]
            for kk in idx[jj:]:
                mijk = mij | masks[res[kk]]
                for ll in idx[kk:]:
                    if mijk | masks[res[ll]] == full:
                        sig = tuple(sorted((oi, orders[res[jj]],
                                            orders[res[kk]],
                                            orders[res[ll]])))
                        if sig not in cap:
                            cap[sig] = (res[i], res[jj], res[kk], res[ll])
    return cap

# --------------------------------------------- independent set-based check

def capable_setbased(N, T=6):
    """Independent (frozenset, no bitmask) unit 4-tuple capability scan."""
    units = [w for w in range(1, N) if gcd(w, N) == 1]

    def bad_set(w):
        return frozenset(k for k in range(N) if T * cls(w * k, N) < N)

    B = {w: bad_set(w) for w in units}
    B1 = B[1]
    full = frozenset(range(N))
    covers = []
    n = len(units)
    for i in range(n):
        sb = B1 | B[units[i]]
        for jj in range(i, n):
            sc = sb | B[units[jj]]
            for kk in range(jj, n):
                if sc | B[units[kk]] == full:
                    covers.append((1, units[i], units[jj], units[kk]))
    return covers

# ---------------------------------------------- odd-composite structure

def odd_composite_structure(N, T=6):
    """Why is (or isn't) an odd composite unit-covering-capable?  Non-unit
    classes can only be reached through j <= h with gcd(j, N) > 1."""
    h = (N - 1) // T
    nonunit = {cls(x, N) for x in range(1, N) if gcd(x, N) > 1}
    unitc = {cls(x, N) for x in range(1, N) if gcd(x, N) == 1}
    js_u = [j for j in range(1, h + 1) if gcd(j, N) == 1]
    js_n = [j for j in range(1, h + 1) if gcd(j, N) > 1]
    reach = set()
    for j in js_n:
        for x in range(1, N):
            if gcd(x, N) == 1:
                reach.add(cls(j * x, N))
    return dict(h=h, n_nonunit_cls=len(nonunit), n_unit_cls=len(unitc),
                unit_slots_per_set=len(js_u), nonunit_slots_per_set=len(js_n),
                nonunit_reachable=len(reach & nonunit))

# ------------------------------------------------------------- class tools

def cls(x, N):
    r = x % N
    return r if r <= N - r else N - r

def cls_inv(w, N):
    return cls(pow(w, -1, N), N)

def flag_classification(N, rs, T, p4set):
    """P4: proved covering classifications. Returns flag iff NON-covering is
    proved via the class reduction (Lemma C5). Caller guarantees units."""
    if N not in p4set:
        return None
    h = (N - 1) // T
    m = (N - 1) // 2
    covered = set()
    for r in rs:
        X = cls_inv(r, N)
        for j in range(1, h + 1):
            covered.add(cls(X * j, N))
    return 'P4:cls%d' % N if len(covered) < m else None

def flag_D3(N, rs, n, T, invs):
    """P5 Lemma D3: +/--coincidence + invertible + counting => no covering."""
    if not all(invs[r] for r in rs):
        return None
    pm = len({(r if r <= N - r else N - r) for r in rs})
    if pm > n - 2:
        return None
    h = (N - 1) // T
    if pm * (2 * h + 1) - (pm - 1) < N:
        return 'P5:D3'
    return None

def battery_flags(N, rs, n, T, tabs, capable, sigcap, p4set):
    """All battery flags for a pair with residues rs (tier logic identical
    to committed N6/N7; only the P4 modulus set is parameterized)."""
    tab = tabs[N]
    invs = tab['invs']
    fl = []
    for j in range(2, T + 1):
        if N % j == 0 and all(r % j for r in rs):
            fl.append('P1:jgon%d' % j)
            break
    allinv = all(invs[r] for r in rs)
    if allinv:
        f = flag_classification(N, rs, T, p4set)
        if f:
            fl.append(f)
    f = flag_D3(N, rs, n, T, invs)
    if f:
        fl.append(f)
    if sum(tab['sizes'][r] for r in rs) - (n - 2) < N:
        fl.append('P2:count')
    if allinv and not capable[N]:
        fl.append('V1:nocovmod')
    if sigcap is not None and all(r != 0 for r in rs):
        sig = tuple(sorted(tab['orders'][r] for r in rs))
        if sig not in sigcap.get(N, {}):
            fl.append('V2:sigcap')
    return fl

# ---------------------------------------------------------- soundness parts

def direct_grid_max(N, eff):
    """max over k in 1..N-1 of min_w ||w*k||_N (independent method)."""
    best = 0
    for k in range(1, N):
        m = N
        for w in eff:
            r = (w * k) % N
            d = r if r <= N - r else N - r
            if d < m:
                m = d
        if m > best:
            best = m
    return best

def soundness_sweep(tabs, capable, n, T, Nrange, p4set, exhaustive=True,
                    n_random=0, seed=1, sigcap=None):
    """Every flag implies non-covering, over residue k-multisets."""
    k = n - 1
    viol = []
    checked = 0
    rng = random.Random(seed)
    for N in Nrange:
        tab = tabs[N]
        masks, full = tab['masks'], tab['full']
        if exhaustive:
            for rs in combinations_with_replacement(range(N), k):
                mm = 0
                for r in rs:
                    mm |= masks[r]
                fl = battery_flags(N, rs, n, T, tabs, capable, sigcap, p4set)
                if fl and mm == full:
                    viol.append((N, rs, fl))
                checked += 1
        if n_random:
            for _ in range(n_random):
                rs = tuple(sorted(rng.randrange(N) for _ in range(k)))
                mm = 0
                for r in rs:
                    mm |= masks[r]
                fl = battery_flags(N, rs, n, T, tabs, capable, sigcap, p4set)
                if fl and mm == full:
                    viol.append((N, rs, fl))
                checked += 1
    return viol, checked

def classification_exactness(tabs, T, n, p4set):
    """At each P4 modulus: predicate(non-cover) == mask(non-cover) for ALL
    unit k-multisets."""
    k = n - 1
    mism = []
    for N in p4set:
        tab = tabs[N]
        masks, full = tab['masks'], tab['full']
        units = [w for w in range(1, N) if tab['invs'][w]]
        for rs in combinations_with_replacement(units, k):
            mm = 0
            for r in rs:
                mm |= masks[r]
            noncover_mask = (mm != full)
            noncover_pred = (flag_classification(N, rs, T, p4set) is not None)
            if noncover_mask != noncover_pred:
                mism.append((N, rs))
    return mism

# ------------------------------------------------------------ corpus run

def corpus_run(n, VMAX, T, tabs, capable, sigcap, p4set, open_cap=30):
    k = n - 1
    total = 0
    closed = Counter()          # strict / P5only / Vonly / OPEN
    jgon_plain = Counter()
    pairstats = Counter()
    cert_hist = Counter()
    ground_fail = []
    open_sets = []
    open_pair_stats = Counter()
    all_unflagged = Counter()
    p1_witness_fail = 0
    m4_mismatch = 0
    open_plain = 0
    danger_pairs = Counter()
    cls1937_sets = 0
    missed1937 = 0     # self-check: unflagged allinv pairs at 19/37 (must
                       # be 0 under the augmented battery)
    v1_closed = 0
    v2_only_closed = 0
    rng = random.Random(7)
    for V in combinations(range(1, VMAX + 1), n):
        g = 0
        for x in V:
            g = gcd(g, x)
        if g != 1:
            continue
        total += 1
        certs = 0
        got = {'s': False, 'f': False, 'v': False, 'j': False,
               'v1': False, 'v2': False, 'c1937': False}
        odetail = []
        set_unflagged = []
        for (p, q) in combinations(range(n), 2):
            N = V[p] + V[q]
            tab = tabs[N]
            others = [V[i] for i in range(n) if i != p and i != q]
            eff = (V[p],) + tuple(others)
            rs = tuple(w % N for w in eff)
            mm = tab['masks'][rs[0]]
            for r in rs[1:]:
                mm |= tab['masks'][r]
            cert = (mm != tab['full'])
            if not cert:
                continue
            certs += 1
            fl = battery_flags(N, rs, n, T, tabs, capable, sigcap, p4set)
            allinv = all(tab['invs'][r] for r in rs)
            if N in RIGID_PRIMES:
                danger_pairs[(N, allinv)] += 1
            for f in fl:
                pairstats[f] += 1
            for f in fl:
                if f.startswith('P1:jgon'):
                    j = int(f[7:])
                    if (mm >> (N // j)) & 1:
                        p1_witness_fail += 1
            if any(f.startswith(('P1', 'P2', 'P4')) for f in fl):
                got['s'] = True
            if any(not f.startswith('V') for f in fl):
                got['f'] = True
            if any(f.startswith('V') for f in fl):
                got['v'] = True
            if any(f == 'V1:nocovmod' for f in fl):
                got['v1'] = True
            if any(f == 'V2:sigcap' for f in fl):
                got['v2'] = True
            if any(f.startswith('P1') for f in fl):
                got['j'] = True
            if not fl:
                sig = tuple(sorted(tab['orders'][r] for r in rs))
                key = (N, allinv, capable[N] if allinv else None, sig)
                set_unflagged.append(key)
                all_unflagged[key] += 1
                if allinv and N in (19, 37):
                    got['c1937'] = True
                    missed1937 += 1
                if len(odetail) < 4:
                    odetail.append((V[p], V[q], N, rs, sig))
            if rng.random() < 0.002:
                best = direct_grid_max(N, eff)
                if (T * best >= N) != cert:
                    m4_mismatch += 1
        if certs == 0:
            ground_fail.append(V)
        cert_hist[certs] += 1
        jgon_plain['jgon' if got['j'] else 'plain'] += 1
        if got['c1937']:
            cls1937_sets += 1
        if got['s']:
            closed['strict'] += 1
        elif got['f']:
            closed['P5only'] += 1
        elif got['v']:
            closed['Vonly'] += 1
            if got['v1']:
                v1_closed += 1
            else:
                v2_only_closed += 1
        else:
            closed['OPEN'] += 1
            for key in set_unflagged:
                open_pair_stats[key] += 1
            if not got['j']:
                open_plain += 1
            if len(open_sets) < open_cap:
                open_sets.append((V, odetail))
    return dict(total=total, closed=closed, jgon_plain=jgon_plain,
                pairstats=pairstats, cert_hist=cert_hist,
                ground_fail=ground_fail, open_sets=open_sets,
                open_pair_stats=open_pair_stats, all_unflagged=all_unflagged,
                open_plain=open_plain, danger_pairs=danger_pairs,
                cls1937_sets=cls1937_sets, missed1937=missed1937,
                v1_closed=v1_closed, v2_only_closed=v2_only_closed,
                p1_witness_fail=p1_witness_fail, m4_mismatch=m4_mismatch)

def report(res, n, VMAX, tag, battery):
    t = res['total']
    c = res['closed']
    strict = c['strict']
    full = c['strict'] + c['P5only']
    anyv = full + c['Vonly']
    print("\n---- corpus %s (n=%d, v<=%d, %d sets) [%s battery] ----"
          % (tag, n, VMAX, t, battery))
    print("  closed by PROVED ports only (P1+P2+P4):   %6d  (%.2f%%)"
          % (strict, 100.0 * strict / t))
    print("  closed by PROVED incl. new P5 (D3):        %6d  (%.2f%%)"
          % (full, 100.0 * full / t))
    print("  closed by PROVED+VERIFIED (+V1/V2):        %6d  (%.2f%%)"
          % (anyv, 100.0 * anyv / t))
    print("  OPEN:                                      %6d  (%.2f%%)"
          % (c['OPEN'], 100.0 * c['OPEN'] / t))
    print("  V-tier split: with-V1-pair %d, V2-only %d"
          % (res['v1_closed'], res['v2_only_closed']))
    print("  j-gon/plain (any-pair) split: %s" % dict(res['jgon_plain']))
    print("  certifying-pairs-per-set histogram: %s"
          % dict(sorted(res['cert_hist'].items())))
    print("  per-flag pair counts: %s"
          % dict(sorted(res['pairstats'].items())))
    print("  rigid-prime certifying pairs (N,allinv)->count: %s"
          % dict(sorted(res['danger_pairs'].items(), key=str)))
    print("  cls1937 candidate sets: %d"
          % res['cls1937_sets'])
    print("  SELF-CHECK unflagged allinv pairs at 19/37: %d"
          % res['missed1937'])
    print("  ground-truth failures (sets w/o cert pair): %d"
          % len(res['ground_fail']))
    print("  soundness: P1 witness fails %d, m4 mismatches %d"
          % (res['p1_witness_fail'], res['m4_mismatch']))
    return strict, full, anyv, c['OPEN']

# ------------------------------------------------------- staged execution
# Background processes are reaped between tool calls in this environment,
# so the run is split into six foreground stages, each well under the
# 10-minute budget, persisting per-corpus results to n9_state.json.

STATE_FILE = '/home/z/my-project/scripts/n9_state.json'

def load_state():
    try:
        with open(STATE_FILE) as f:
            return json.load(f)
    except Exception:
        return {}

def unflagged_buckets(res, augmented):
    """Kernel-world decomposition of the unflagged certifying pairs."""
    au = res['all_unflagged']
    tot = sum(au.values())
    by = Counter()
    for (N, allinv, cap, sig), c in au.items():
        if allinv and N in (19, 37) and not augmented:
            by['P4-target: allinv @ 19/37'] += c
        elif allinv:
            by['allinv @ capable %d' % N] += c
        elif N >= 96:
            by['kernel @ 96..111'] += c
        elif N >= 80:
            by['kernel @ 80..95'] += c
        elif N >= 48:
            by['kernel @ 48..79'] += c
        else:
            by['kernel @ N<=47'] += c
    return tot, by

def unflagged_decomposition(res, label, augmented):
    tot, by = unflagged_buckets(res, augmented)
    print("  [%s] unflagged certifying pairs corpus-wide: %d" % (label, tot))
    for k, v in sorted(by.items()):
        print("      %-34s %7d" % (k, v))
    return tot, by

def save_tag(tag, res, augmented):
    tot, by = unflagged_buckets(res, augmented)
    st = load_state()
    st[tag] = dict(
        total=res['total'], strict=res['closed']['strict'],
        p5only=res['closed']['P5only'], vonly=res['closed']['Vonly'],
        opensets=res['closed']['OPEN'], v1=res['v1_closed'],
        v2=res['v2_only_closed'], missed1937=res['missed1937'],
        cls1937=res['cls1937_sets'], open_plain=res['open_plain'],
        ground_fail=len(res['ground_fail']),
        p1_fail=res['p1_witness_fail'], m4=res['m4_mismatch'],
        unflag_total=tot, unflag_buckets=dict(by),
        jgon=res['jgon_plain'].get('jgon', 0),
        plain=res['jgon_plain'].get('plain', 0))
    with open(STATE_FILE, 'w') as f:
        json.dump(st, f, indent=1, sort_keys=True)
    print("  [state] tag '%s' saved (%d sets)" % (tag, res['total']))

def run_save(n, VMAX, T, tabs, cap, sigcap, p4set, tag, label, battery,
             augmented):
    res = corpus_run(n, VMAX, T, tabs, cap, sigcap, p4set)
    report(res, n, VMAX, label, battery)
    save_tag(tag, res, augmented)
    return res

# ---------------------------------------------------------------- prelude

def prelude(with_soundness):
    print("[tables] building n=4 (N<=31) and n=5 (N<=111) tables...")
    tabs4 = {N: build_tables(N, 5) for N in range(3, 32)}
    tabs5 = {N: build_tables(N, 6) for N in range(3, 112)}
    bad4 = [b for N in range(3, 32) for b in check_size_formula(N, 5, tabs4[N])]
    bad5 = [b for N in range(3, 112) for b in check_size_formula(N, 6, tabs5[N])]
    print("[check] Lemma S5 exact-size violations: n4 %d, n5(N<=111) %d"
          % (len(bad4), len(bad5)))
    assert not bad4 and not bad5

    t1 = time.time()
    cov4 = {N: capable_units(N, tabs4[N], 3) for N in range(3, 32)}
    cap4 = {N: bool(v) for N, v in cov4.items()}
    cov5 = {N: capable_units(N, tabs5[N], 4) for N in range(3, 112)}
    cap5 = {N: bool(v) for N, v in cov5.items()}
    cap5_list = sorted(N for N in cap5 if cap5[N])
    ctrl_n4 = sorted(N for N in range(3, 32) if cap4[N])
    print("[n=4] unit-triple capable moduli (T=5): %s" % ctrl_n4)
    assert ctrl_n4 == [7, 11, 13], "n=4 capability control FAILED"
    ctrl_n5_47 = sorted(N for N in range(3, 48) if cap5[N])
    print("[n=5] unit 4-tuple capable moduli, N<=47 (control): %s" % ctrl_n5_47)
    assert ctrl_n5_47 == [7, 13, 17, 19, 37], "n=5 capability control FAILED"
    cap_48_95 = [N for N in cap5_list if 48 <= N <= 95]
    print("[n=5] capable moduli in 48..95 (control, committed []): %s"
          % cap_48_95)
    assert cap_48_95 == [], "conflict with committed N6/N8 (48..95 safe)"
    for N in (83, 89, 97, 101, 103, 107, 109):
        assert not cap5[N], \
            "N=%d capable — conflicts with committed (47,150] prime scan" % N
    print("[n=5] capable primes 83, 89, 97, 101, 103, 107, 109 (inside "
          "committed (47,150] scan): all safe")
    new_cap = [N for N in cap5_list if N >= 96]
    print("[n=5] *** capable moduli in 96..111 (NEW RANGE): %s ***" % new_cap)
    print("[n=5] covering-configuration counts per capable modulus: %s"
          % {N: len(cov5[N]) for N in cap5_list})
    assert {N: len(cov5[N]) for N in cap5_list if N <= 47} == \
        {7: 20, 13: 68, 17: 16, 19: 64, 37: 32}, "committed counts FAILED"

    new_odd_composites = [99, 105, 111]
    print("\n[odd composites 96..111] structure vs capability (T=6):")
    print("      N    h   nonunitCls reach  unitSlots/4set  capable  #cov")
    for N in new_odd_composites:
        st = odd_composite_structure(N)
        print("     %3d  %2d  %5d  %5d  %2d+%2d  %8s  %4d"
              % (N, st['h'], st['n_nonunit_cls'], st['nonunit_reachable'],
                 st['unit_slots_per_set'], st['nonunit_slots_per_set'],
                 cap5[N], len(cov5[N])))

    print("\n[cross-check] independent set-based capability scan (new range):")
    for N in new_odd_composites + [97, 101, 103, 107, 109]:
        cov_sb = capable_setbased(N)
        ok = (len(cov_sb) == len(cov5[N]))
        print("      N=%3d: set-based %d vs mask-based %d coverings  %s"
              % (N, len(cov_sb), len(cov5[N]), "AGREE" if ok else "MISMATCH"))
        assert ok, "set-based/mask-based capability mismatch at N=%d" % N
    print("      (elapsed %.1fs)" % (time.time() - t1))

    t1 = time.time()
    print("\n[sigcap] building order-signature capability maps, N=3..111...")
    sigcap = {}
    for N in range(3, 112):
        sigcap[N] = signature_capability(N, tabs5[N])
    nsig_new = {N: len(sigcap[N]) for N in range(96, 112) if sigcap[N]}
    print("[n=5] capable signatures, N in 96..111 (NEW): %s" % nsig_new)
    ctrl_sig_48_79 = {N: len(sigcap[N]) for N in range(48, 80) if sigcap[N]}
    assert ctrl_sig_48_79 == {48: 13, 49: 2, 50: 2, 51: 2, 52: 9, 54: 4,
                              56: 22, 57: 2, 58: 1, 60: 5, 62: 1, 63: 9,
                              64: 10, 65: 3, 66: 2, 68: 3, 69: 1, 70: 13,
                              72: 19, 74: 2, 75: 2, 76: 4, 77: 3, 78: 11}, \
        "sigcap control (48..79) vs committed N6 FAILED"
    ctrl_sig_80_95 = {N: len(sigcap[N]) for N in range(80, 96) if sigcap[N]}
    assert ctrl_sig_80_95 == {80: 14, 81: 3, 82: 1, 84: 22, 85: 1, 86: 1,
                              87: 1, 88: 10, 90: 6, 91: 6, 92: 2, 93: 1,
                              94: 1, 95: 1}, \
        "sigcap control (80..95) vs committed N8 FAILED"
    print("[n=5] sigcap controls 48..79 (N6) and 80..95 (N8): EXACT")
    print("      (elapsed %.1fs)" % (time.time() - t1))

    if with_soundness:
        t1 = time.time()
        mism5 = classification_exactness(tabs5, 6, 5, AUG_P4)
        print("\n[check] P4 classification exactness (augmented set) "
              "mismatches: %d" % len(mism5))
        assert not mism5
        v5a, c5a = soundness_sweep(tabs5, cap5, 5, 6, range(96, 112), AUG_P4,
                                    exhaustive=False, n_random=4000,
                                    sigcap=sigcap)
        v5b, c5b = soundness_sweep(tabs5, cap5, 5, 6, [99], AUG_P4,
                                    exhaustive=True, sigcap=sigcap)
        v5c, c5c = soundness_sweep(tabs5, cap5, 5, 6, [111], AUG_P4,
                                    exhaustive=True, sigcap=sigcap)
        print("[check] flag=>cert sweep (NEW range): 96..111 random %d/%d "
              "viol; N=99 exhaustive %d/%d viol; N=111 exhaustive %d/%d viol"
              % (len(v5a), c5a, len(v5b), c5b, len(v5c), c5c))
        assert not v5a and not v5b and not v5c
        print("      (elapsed %.1fs)" % (time.time() - t1))
    return tabs4, tabs5, cap4, cap5, sigcap

# ----------------------------------------------------------------- stages

def stage1():
    t0 = time.time()
    print("=" * 72)
    print("N9 STAGE 1/6 — prelude + soundness + controls (n=4, V=16..40)")
    print("=" * 72)
    tabs4, tabs5, cap4, cap5, sigcap = prelude(True)

    print("\n===== PART A: n=4 positive control (v<=16) =====")
    res4 = corpus_run(4, 16, 5, tabs4, cap4, None, N4_P4)
    s4, f4, v4c, o4 = report(res4, 4, 16, "control", "n=4")
    assert (s4, f4, v4c, o4, res4['total']) == (1651, 1661, 1744, 1, 1745), \
        "n=4 CONTROL FAILED: got %s" % ((s4, f4, v4c, o4, res4['total']),)
    assert res4['open_sets'][0][0] == (3, 5, 8, 13), \
        "n=4 control open set mismatch"
    save_tag('n4', res4, False)
    print("  [CONTROL PASS] n=4 reproduces committed: 1651/1661/1744 of "
          "1745, open set = (3,5,8,13)")

    print("\n===== PART B: SAME battery controls V=16..40 =====")
    res16 = run_save(5, 16, 6, tabs5, cap5, sigcap, SAME_P4, "V16s",
                     "control-V16", "same", False)
    assert (res16['closed']['strict'],
            res16['closed']['strict'] + res16['closed']['P5only'],
            res16['closed']['OPEN'], res16['total']) == \
        (4179, 4193, 0, 4311), "V=16 same-battery CONTROL FAILED"
    print("  [CONTROL PASS] V=16 same: 4179/4193/4311, OPEN 0")

    res24 = run_save(5, 24, 6, tabs5, cap5, sigcap, SAME_P4, "V24s",
                     "control-V24", "same", False)
    f24 = res24['closed']['strict'] + res24['closed']['P5only']
    assert (res24['closed']['strict'], f24, res24['closed']['OPEN'],
            res24['total']) == (35778, 36362, 0, 41656), \
        "V=24 same-battery CONTROL FAILED"
    assert f24 + res24['v1_closed'] == 41240, "V=24 +V1 tier FAILED"
    assert f24 + res24['v1_closed'] + res24['v2_only_closed'] == 41656
    unflagged_decomposition(res24, "V=24 same", augmented=False)
    print("  [CONTROL PASS] V=24 same: 35778/36362/41240/41656, OPEN 0 "
          "(V-tier: %d V1-closed + %d V2-only)"
          % (res24['v1_closed'], res24['v2_only_closed']))

    res32 = run_save(5, 32, 6, tabs5, cap5, sigcap, SAME_P4, "V32s",
                     "control-V32", "same", False)
    f32 = res32['closed']['strict'] + res32['closed']['P5only']
    assert (res32['closed']['strict'], f32, res32['closed']['OPEN'],
            res32['total']) == (156461, 161596, 0, 196751), \
        "V=32 same-battery CONTROL FAILED"
    print("  [CONTROL PASS] V=32 same: 156461/161596/196751, OPEN 0")

    res40s = run_save(5, 40, 6, tabs5, cap5, sigcap, SAME_P4, "V40s",
                      "control-V40", "same", False)
    f40s = res40s['closed']['strict'] + res40s['closed']['P5only']
    assert (res40s['closed']['strict'], f40s, res40s['closed']['OPEN'],
            res40s['total']) == (480454, 500280, 0, 641166), \
        "V=40 same-battery CONTROL FAILED"
    tot40s, by40s = unflagged_decomposition(res40s, "V=40 same",
                                            augmented=False)
    assert tot40s == 534098, "V=40 same unflagged total CONTROL FAILED"
    print("  [CONTROL PASS] V=40 same: 480454/500280/641166, OPEN 0")
    print("\n  stage 1 elapsed: %.1fs" % (time.time() - t0))

def stage2():
    t0 = time.time()
    print("=" * 72)
    print("N9 STAGE 2/6 — same battery V=48 control + augmented V=24..40")
    print("=" * 72)
    tabs4, tabs5, cap4, cap5, sigcap = prelude(False)

    print("\n===== PART C: SAME battery CONTROL V=48 (committed N8) =====")
    res48s = run_save(5, 48, 6, tabs5, cap5, sigcap, SAME_P4, "V48s",
                      "control-V48", "same", False)
    f48s = res48s['closed']['strict'] + res48s['closed']['P5only']
    assert (res48s['closed']['strict'], f48s, res48s['closed']['OPEN'],
            res48s['total']) == (1200777, 1253258, 0, 1665356), \
        "V=48 same-battery CONTROL FAILED"
    tot48s, by48s = unflagged_decomposition(res48s, "V=48 same",
                                            augmented=False)
    assert by48s['P4-target: allinv @ 19/37'] == 358745
    assert by48s['kernel @ N<=47'] == 532700
    assert by48s['kernel @ 48..79'] == 137484
    assert by48s['kernel @ 80..95'] == 12306
    assert tot48s == 1041235, "V=48 same unflagged total CONTROL FAILED"
    print("  [CONTROL PASS] V=48 same: 1200777/1253258/1665356, OPEN 0; "
          "unflagged decomposition exact (358745/532700/137484/12306)")

    print("\n===== PART D: AUGMENTED battery controls V=24..40 =====")
    res24a = run_save(5, 24, 6, tabs5, cap5, sigcap, AUG_P4, "V24a",
                      "control-V24", "AUGMENTED", True)
    f24a = res24a['closed']['strict'] + res24a['closed']['P5only']
    assert (res24a['closed']['strict'], f24a, res24a['closed']['OPEN'],
            res24a['total']) == (38611, 38924, 0, 41656), \
        "V=24 augmented CONTROL FAILED"
    assert res24a['missed1937'] == 0
    print("  [CONTROL PASS] V=24 augmented: 38611/38924/41656, OPEN 0")

    res32a = run_save(5, 32, 6, tabs5, cap5, sigcap, AUG_P4, "V32a",
                      "control-V32", "AUGMENTED", True)
    f32a = res32a['closed']['strict'] + res32a['closed']['P5only']
    assert (res32a['closed']['strict'], f32a, res32a['closed']['OPEN'],
            res32a['total']) == (173337, 176376, 0, 196751), \
        "V=32 augmented CONTROL FAILED"
    assert res32a['missed1937'] == 0
    print("  [CONTROL PASS] V=32 augmented: 173337/176376/196751, OPEN 0")

    res40a = run_save(5, 40, 6, tabs5, cap5, sigcap, AUG_P4, "V40a",
                      "control-V40", "AUGMENTED", True)
    f40a = res40a['closed']['strict'] + res40a['closed']['P5only']
    assert (res40a['closed']['strict'], f40a, res40a['closed']['OPEN'],
            res40a['total']) == (526607, 540415, 0, 641166), \
        "V=40 augmented CONTROL FAILED"
    assert res40a['missed1937'] == 0
    tot40a, by40a = unflagged_decomposition(res40a, "V=40 augmented",
                                            augmented=True)
    assert tot40a == 338126, "V=40 augmented unflagged total CONTROL FAILED"
    print("  [CONTROL PASS] V=40 augmented: 526607/540415/641166, OPEN 0")
    print("\n  stage 2 elapsed: %.1fs" % (time.time() - t0))

def stage3():
    t0 = time.time()
    print("=" * 72)
    print("N9 STAGE 3/6 — SAME battery PRIMARY V=56")
    print("=" * 72)
    tabs4, tabs5, cap4, cap5, sigcap = prelude(False)
    print("\n===== PART C2: SAME battery PRIMARY V=56 =====")
    res56s = run_save(5, 56, 6, tabs5, cap5, sigcap, SAME_P4, "V56s",
                      "primary-V56", "same", False)
    assert res56s['total'] == 3712576, \
        "V=56 corpus size != Mobius count 3712576"
    unflagged_decomposition(res56s, "V=56 same", augmented=False)
    print("\n  stage 3 elapsed: %.1fs" % (time.time() - t0))

def stage4():
    t0 = time.time()
    print("=" * 72)
    print("N9 STAGE 4/6 — AUGMENTED battery CONTROL V=48 (committed N8)")
    print("=" * 72)
    tabs4, tabs5, cap4, cap5, sigcap = prelude(False)
    print("\n===== PART E: AUGMENTED battery CONTROL V=48 =====")
    res48a = run_save(5, 48, 6, tabs5, cap5, sigcap, AUG_P4, "V48a",
                      "control-V48", "AUGMENTED", True)
    f48a = res48a['closed']['strict'] + res48a['closed']['P5only']
    assert (res48a['closed']['strict'], f48a, res48a['closed']['OPEN'],
            res48a['total']) == (1295478, 1335857, 0, 1665356), \
        "V=48 augmented CONTROL FAILED"
    assert res48a['missed1937'] == 0
    tot48a, by48a = unflagged_decomposition(res48a, "V=48 augmented",
                                            augmented=True)
    assert by48a['kernel @ N<=47'] == 532700
    assert by48a['kernel @ 48..79'] == 137484
    assert by48a['kernel @ 80..95'] == 12306
    assert tot48a == 682490, "V=48 augmented unflagged total CONTROL FAILED"
    print("  [CONTROL PASS] V=48 augmented: 1295478/1335857/1665356, OPEN 0; "
          "unflagged decomposition exact (532700/137484/12306)")
    print("\n  stage 4 elapsed: %.1fs" % (time.time() - t0))

def stage5():
    t0 = time.time()
    print("=" * 72)
    print("N9 STAGE 5/6 — AUGMENTED battery PRIMARY V=56")
    print("=" * 72)
    tabs4, tabs5, cap4, cap5, sigcap = prelude(False)
    print("\n===== PART E2: AUGMENTED battery PRIMARY V=56 =====")
    res56a = run_save(5, 56, 6, tabs5, cap5, sigcap, AUG_P4, "V56a",
                      "primary-V56", "AUGMENTED", True)
    assert res56a['total'] == 3712576, \
        "V=56 corpus size != Mobius count 3712576"
    assert res56a['missed1937'] == 0
    unflagged_decomposition(res56a, "V=56 augmented", augmented=True)
    print("\n  stage 5 elapsed: %.1fs" % (time.time() - t0))

def stage6():
    t0 = time.time()
    st = load_state()
    need = ['n4', 'V16s', 'V24s', 'V32s', 'V40s', 'V48s', 'V56s',
            'V24a', 'V32a', 'V40a', 'V48a', 'V56a']
    missing = [tg for tg in need if tg not in st]
    assert not missing, "missing tags in state: %s" % missing
    for tg in need:
        r = st[tg]
        assert r['ground_fail'] == 0 and r['p1_fail'] == 0 and r['m4'] == 0, \
            "soundness failure in tag %s" % tg
    print("=" * 72)
    print("N9 STAGE 6/6 — TREND SUMMARY AND THE HONEST READING")
    print("=" * 72)
    n4 = st['n4']
    print("\n  n=4 control: strict %d / +D3 %d / +V %d / OPEN %d of %d "
          "(open set (3,5,8,13))"
          % (n4['strict'], n4['strict'] + n4['p5only'],
             n4['strict'] + n4['p5only'] + n4['vonly'], n4['opensets'],
             n4['total']))
    print("\n===== TREND: SAME battery (P4 at {7,13,17}) =====")
    print("  V    sets    ports      +D3      +V(any)   OPEN")
    for tg, VV in (('V16s', 16), ('V24s', 24), ('V32s', 32), ('V40s', 40),
                   ('V48s', 48), ('V56s', 56)):
        r = st[tg]
        s = r['strict']
        fu = s + r['p5only']
        av = fu + r['vonly']
        print("  %2d  %6d  %6d %6.2f%%  %6.2f%%  %6.2f%%  %4d"
              % (VV, r['total'], s, 100.0 * s / r['total'],
                 100.0 * fu / r['total'], 100.0 * av / r['total'],
                 r['opensets']))
    print("\n===== TREND: AUGMENTED battery (P4 at {7,13,17,19,37}) =====")
    print("  V    sets    ports      +D3      +V(any)   OPEN")
    aug = []
    for tg, VV in (('V24a', 24), ('V32a', 32), ('V40a', 40), ('V48a', 48),
                   ('V56a', 56)):
        r = st[tg]
        s = r['strict']
        fu = s + r['p5only']
        av = fu + r['vonly']
        aug.append((VV, r['total'], s, fu))
        print("  %2d  %6d  %6d %6.2f%%  %6.2f%%  %6.2f%%  %4d"
              % (VV, r['total'], s, 100.0 * s / r['total'],
                 100.0 * fu / r['total'], 100.0 * av / r['total'],
                 r['opensets']))
    print("\n  V=56 unflagged decomposition (same battery): %s"
          % st['V56s']['unflag_buckets'])
    print("  V=56 unflagged decomposition (augmented):    %s"
          % st['V56a']['unflag_buckets'])

    aug_d3 = [100.0 * fu / t for (VV, t, s, fu) in aug]
    resid = [100.0 - x for x in aug_d3]
    incs = [resid[i + 1] - resid[i] for i in range(len(resid) - 1)]
    aug56_d3 = aug_d3[-1]
    aug56_ports = 100.0 * aug[-1][2] / aug[-1][1]
    print("\n  V=56 augmented: ports %.2f%%, +D3 %.2f%%"
          % (aug56_ports, aug56_d3))
    print("  augmented +D3 coverage series: %s"
          % " -> ".join("%.2f" % x for x in aug_d3))
    print("  residual series: %s"
          % " -> ".join("%.2f" % x for x in resid))
    print("  residual increments: %s"
          % ", ".join("%+.2f" % x for x in incs))
    print("  (V=48 margin over the 80%% bar was +0.21 pts; mean residual "
          "increment is %.2f pts per 8 V)" % (sum(incs) / len(incs)))
    if aug56_d3 >= 80.0:
        print("\n  READING: augmented +D3 at V=56 is STILL ABOVE 80%% "
              "(margin %+.2f pts)." % (aug56_d3 - 80.0))
        print("  The trend does not establish that coverage stabilizes; "
              "five data")
        print("  points still cannot distinguish deceleration from "
              "roughly-linear")
        print("  residual growth with noise. The residual is entirely "
              "kernel-world.")
    elif aug56_d3 >= 60.0:
        print("\n  READING: augmented +D3 at V=56 is BELOW the 80%% "
              "continuation")
        print("  bar (margin %+.2f pts). Within the tested range the "
              "reformulation" % (aug56_d3 - 80.0))
        print("  remains productive; the honest reading is fundamental "
              "decay /")
        print("  unresolved asymptotics. The residual is entirely "
              "kernel-world;")
        print("  closing it requires the kernel-world structure theorem.")
    else:
        print("\n  READING: augmented +D3 at V=56 is BELOW the 60%% pivot "
              "level —")
        print("  pivot territory per the original decision rule.")
    print("\n  stage 6 elapsed: %.1fs" % (time.time() - t0))

STAGES = {'1': stage1, '2': stage2, '3': stage3, '4': stage4,
          '5': stage5, '6': stage6}

if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else ''
    if which == 'all':
        for k in sorted(STAGES):
            STAGES[k]()
    elif which in STAGES:
        STAGES[which]()
    else:
        raise SystemExit(
            "usage: N9_v56_extension.py <1..6|all>\n"
            "  staged execution: this environment reaps background "
            "processes\n  between tool calls, so each stage runs in the "
            "foreground under the\n  10-minute budget; per-corpus results "
            "persist in n9_state.json;\n  stage 6 prints the trend summary "
            "and the honest reading once all\n  tags are present.")
