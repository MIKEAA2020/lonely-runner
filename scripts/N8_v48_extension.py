#!/usr/bin/env python3
"""
N8_v48_extension.py — the optional V=48 data point (reviewer: "if compute
allows"), completing BOTH decay series (same battery AND augmented battery)
before the paper is written.  Decision rule (reviewer, verbatim intent):
augmented coverage at V=48 still above 80%  -> decay slower than
extrapolated, paper framing CONFIDENT;  below 78% -> on the extrapolated
track, framing CAUTIOUS.  Either way the paper is writable.

BATTERY: byte-identical to the committed N6 (same battery: P4 at
{7,13,17}) and N7 (augmented battery: P4 at {7,13,17,19,37}) runs — the
ONLY change is range, not battery:
  - tables / capability / signature-capability maps for ALL N <= 95
    (pair sums reach 2*48-1 = 95);
  - capability scan of the genuinely new moduli 80..95 (odd composites
    81, 85, 87, 91, 93, 95; primes 83, 89 are inside the committed
    (47,150] scan and are re-derived as controls, NOT new scan);
  - odd-composite structure table + independent set-based cross-check in
    the new range;
  - flag=>cert soundness sweeps extended to 80..95 (random) with
    exhaustive checks at the new-range endpoints 81 and 95.

HARD CONTROLS (asserts, all vs committed logs):
  - n=4 positive control (v<=16): 1651/1661/1744, OPEN 1, open set
    (3,5,8,13)  [committed n5/n6 record]
  - same battery: V=16 4179/4193/4311 OPEN 0; V=24 35778/36362/41656
    OPEN 0 with V-tier 4878+416; V=32 156461/161596/196751 OPEN 0;
    V=40 480454/500280/641166 OPEN 0  [scripts/n6_run.log]
  - augmented battery: V=24 38611/38924/41656; V=32 173337/176376/196751;
    V=40 526607/540415/641166, OPEN 0, self-check unflagged allinv @19/37
    == 0 everywhere  [scripts/n7_run.log]
  - capability: n=4 {7,11,13}; n=5 N<=47 {7,13,17,19,37} with counts
    20/68/16/64/32; capable in 48..79 == []  [committed N6]

Deterministic.  Runtime ~7 min.  Run: python3 -u N8_v48_extension.py
"""
from math import gcd
from itertools import combinations, combinations_with_replacement
from collections import Counter
import random
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

def unflagged_decomposition(res, label, augmented):
    au = res['all_unflagged']
    tot = sum(au.values())
    by = Counter()
    for (N, allinv, cap, sig), c in au.items():
        if allinv and N in (19, 37) and not augmented:
            by['P4-target: allinv @ 19/37'] += c
        elif allinv:
            by['allinv @ capable %d' % N] += c
        elif N >= 80:
            by['kernel @ new composite N>=80'] += c
        elif N >= 48:
            by['kernel @ new composite N>=48'] += c
        else:
            by['kernel @ N<=47'] += c
    print("  [%s] unflagged certifying pairs corpus-wide: %d" % (label, tot))
    for k, v in sorted(by.items()):
        print("      %-34s %7d" % (k, v))
    return tot

# ------------------------------------------------------------------- main

def main():
    t0 = time.time()
    print("=" * 72)
    print("N8 — THE V=48 DATA POINT (both batteries; the reviewer's optional")
    print("     extrapolation check before the paper)")
    print("=" * 72)

    # ---------------- tables (N <= 95: pair sums reach 2*48-1 = 95)
    print("[tables] building n=4 (N<=31) and n=5 (N<=95) tables...")
    tabs4 = {N: build_tables(N, 5) for N in range(3, 32)}
    tabs5 = {N: build_tables(N, 6) for N in range(3, 96)}
    bad4 = [b for N in range(3, 32) for b in check_size_formula(N, 5, tabs4[N])]
    bad5 = [b for N in range(3, 96) for b in check_size_formula(N, 6, tabs5[N])]
    print("[check] Lemma S5 exact-size violations: n4 %d, n5(N<=95) %d"
          % (len(bad4), len(bad5)))
    assert not bad4 and not bad5

    # ---------------- capability maps (ALL N <= 95)
    t1 = time.time()
    cov4 = {N: capable_units(N, tabs4[N], 3) for N in range(3, 32)}
    cap4 = {N: bool(v) for N, v in cov4.items()}
    cov5 = {N: capable_units(N, tabs5[N], 4) for N in range(3, 96)}
    cap5 = {N: bool(v) for N, v in cov5.items()}
    cap5_list = sorted(N for N in cap5 if cap5[N])
    ctrl_n4 = sorted(N for N in range(3, 32) if cap4[N])
    print("[n=4] unit-triple capable moduli (T=5): %s" % ctrl_n4)
    assert ctrl_n4 == [7, 11, 13], "n=4 capability control FAILED"
    ctrl_n5_47 = sorted(N for N in range(3, 48) if cap5[N])
    print("[n=5] unit 4-tuple capable moduli, N<=47 (control): %s" % ctrl_n5_47)
    assert ctrl_n5_47 == [7, 13, 17, 19, 37], "n=5 capability control FAILED"
    cap_48_79 = [N for N in cap5_list if 48 <= N <= 79]
    print("[n=5] capable moduli in 48..79 (control, committed []): %s"
          % cap_48_79)
    assert cap_48_79 == [], "conflict with committed N6 (48..79 safe)"
    for N in (83, 89):
        assert not cap5[N], \
            "N=%d capable — conflicts with committed (47,150] prime scan" % N
    print("[n=5] capable primes 83, 89 (inside committed (47,150] scan): "
          "safe, safe")
    new_cap = [N for N in cap5_list if N >= 80]
    print("[n=5] *** capable moduli in 80..95 (NEW RANGE): %s ***" % new_cap)
    print("[n=5] covering-configuration counts per capable modulus: %s"
          % {N: len(cov5[N]) for N in cap5_list})
    assert {N: len(cov5[N]) for N in cap5_list if N <= 47} == \
        {7: 20, 13: 68, 17: 16, 19: 64, 37: 32}, "committed counts FAILED"

    # ---- odd-composite structural table for the new range
    new_odd_composites = [81, 85, 87, 91, 93, 95]
    print("\n[odd composites 81..95] structure vs capability (T=6):")
    print("      N    h   nonunitCls reach  unitSlots/4set  capable  #cov")
    for N in new_odd_composites:
        st = odd_composite_structure(N)
        print("     %3d  %2d  %5d  %5d  %2d+%2d  %8s  %4d"
              % (N, st['h'], st['n_nonunit_cls'], st['nonunit_reachable'],
                 st['unit_slots_per_set'], st['nonunit_slots_per_set'],
                 cap5[N], len(cov5[N])))

    # ---- independent set-based cross-check at the new-range danger moduli
    print("\n[cross-check] independent set-based capability scan (new range):")
    for N in new_odd_composites + [83, 89]:
        cov_sb = capable_setbased(N)
        ok = (len(cov_sb) == len(cov5[N]))
        print("      N=%3d: set-based %d vs mask-based %d coverings  %s"
              % (N, len(cov_sb), len(cov5[N]), "AGREE" if ok else "MISMATCH"))
        assert ok, "set-based/mask-based capability mismatch at N=%d" % N
    print("      (elapsed %.1fs)" % (time.time() - t1))

    # ---------------- signature capability map (all N <= 95)
    t1 = time.time()
    print("\n[sigcap] building order-signature capability maps, N=3..95...")
    sigcap = {}
    for N in range(3, 96):
        sigcap[N] = signature_capability(N, tabs5[N])
    nsig_new = {N: len(sigcap[N]) for N in range(80, 96) if sigcap[N]}
    print("[n=5] capable signatures, N in 80..95 (NEW): %s" % nsig_new)
    ctrl_sig_48_79 = {N: len(sigcap[N]) for N in range(48, 80) if sigcap[N]}
    assert ctrl_sig_48_79 == {48: 13, 49: 2, 50: 2, 51: 2, 52: 9, 54: 4,
                              56: 22, 57: 2, 58: 1, 60: 5, 62: 1, 63: 9,
                              64: 10, 65: 3, 66: 2, 68: 3, 69: 1, 70: 13,
                              72: 19, 74: 2, 75: 2, 76: 4, 77: 3, 78: 11}, \
        "sigcap control (48..79) vs committed N6 FAILED"
    print("[n=5] sigcap control 48..79 vs committed N6: EXACT")
    print("      (elapsed %.1fs)" % (time.time() - t1))

    # ---------------- soundness sweeps
    t1 = time.time()
    mism5 = classification_exactness(tabs5, 6, 5, AUG_P4)
    print("\n[check] P4 classification exactness (augmented set) mismatches: %d"
          % len(mism5))
    assert not mism5
    # committed ranges are committed; new range: random + endpoint exhaustive
    v5a, c5a = soundness_sweep(tabs5, cap5, 5, 6, range(80, 96), AUG_P4,
                                exhaustive=False, n_random=4000,
                                sigcap=sigcap)
    v5b, c5b = soundness_sweep(tabs5, cap5, 5, 6, [81], AUG_P4,
                                exhaustive=True, sigcap=sigcap)
    v5c, c5c = soundness_sweep(tabs5, cap5, 5, 6, [95], AUG_P4,
                                exhaustive=True, sigcap=sigcap)
    print("[check] flag=>cert sweep (NEW range): 80..95 random %d/%d viol; "
          "N=81 exhaustive %d/%d viol; N=95 exhaustive %d/%d viol"
          % (len(v5a), c5a, len(v5b), c5b, len(v5c), c5c))
    assert not v5a and not v5b and not v5c
    print("      (elapsed %.1fs)" % (time.time() - t1))

    # ---------------- Part A: n=4 positive control
    print("\n===== PART A: n=4 positive control (v<=16) =====")
    res4 = corpus_run(4, 16, 5, tabs4, cap4, None, N4_P4)
    s4, f4, v4c, o4 = report(res4, 4, 16, "control", "n=4")
    assert (s4, f4, v4c, o4, res4['total']) == (1651, 1661, 1744, 1, 1745), \
        "n=4 CONTROL FAILED: got %s" % ((s4, f4, v4c, o4, res4['total']),)
    assert res4['open_sets'][0][0] == (3, 5, 8, 13), \
        "n=4 control open set mismatch"
    print("  [CONTROL PASS] n=4 reproduces committed: 1651/1661/1744 of "
          "1745, open set = (3,5,8,13)")

    # ---------------- Part B: SAME battery (P4 at {7,13,17}) — controls + V=48
    print("\n===== PART B: SAME battery (P4 at {7,13,17}) =====")
    res16 = corpus_run(5, 16, 6, tabs5, cap5, sigcap, SAME_P4)
    s16, f16, v16, o16 = report(res16, 5, 16, "control-V16", "same")
    assert (s16, f16, o16, res16['total']) == (4179, 4193, 0, 4311), \
        "V=16 same-battery CONTROL FAILED"
    print("  [CONTROL PASS] V=16 same: 4179/4193/4311, OPEN 0")

    res24 = corpus_run(5, 24, 6, tabs5, cap5, sigcap, SAME_P4)
    s24, f24, v24, o24 = report(res24, 5, 24, "control-V24", "same")
    assert (s24, f24, o24, res24['total']) == (35778, 36362, 0, 41656), \
        "V=24 same-battery CONTROL FAILED"
    assert f24 + res24['v1_closed'] == 41240, "V=24 +V1 tier FAILED"
    assert f24 + res24['v1_closed'] + res24['v2_only_closed'] == 41656
    print("  [CONTROL PASS] V=24 same: 35778/36362/41240/41656, OPEN 0 "
          "(V-tier: %d V1-closed + %d V2-only)"
          % (res24['v1_closed'], res24['v2_only_closed']))
    unflagged_decomposition(res24, "V=24 same", augmented=False)

    res32 = corpus_run(5, 32, 6, tabs5, cap5, sigcap, SAME_P4)
    s32, f32, v32, o32 = report(res32, 5, 32, "control-V32", "same")
    assert (s32, f32, o32, res32['total']) == (156461, 161596, 0, 196751), \
        "V=32 same-battery CONTROL FAILED"
    print("  [CONTROL PASS] V=32 same: 156461/161596/196751, OPEN 0")

    res40s = corpus_run(5, 40, 6, tabs5, cap5, sigcap, SAME_P4)
    s40s, f40s, v40s, o40s = report(res40s, 5, 40, "control-V40", "same")
    assert (s40s, f40s, o40s, res40s['total']) == \
        (480454, 500280, 0, 641166), "V=40 same-battery CONTROL FAILED"
    print("  [CONTROL PASS] V=40 same: 480454/500280/641166, OPEN 0")
    unflagged_decomposition(res40s, "V=40 same", augmented=False)

    print("\n===== PART C: SAME battery PRIMARY V=48 =====")
    res48s = corpus_run(5, 48, 6, tabs5, cap5, sigcap, SAME_P4)
    s48s, f48s, v48s, o48s = report(res48s, 5, 48, "primary-V48", "same")
    unflagged_decomposition(res48s, "V=48 same", augmented=False)

    # ---------------- Part D: AUGMENTED battery (P4 at {7,13,17,19,37})
    print("\n===== PART D: AUGMENTED battery (P4 at {7,13,17,19,37}) =====")
    res24a = corpus_run(5, 24, 6, tabs5, cap5, sigcap, AUG_P4)
    s24a, f24a, v24a, o24a = report(res24a, 5, 24, "control-V24", "AUGMENTED")
    assert (s24a, f24a, o24a, res24a['total']) == (38611, 38924, 0, 41656), \
        "V=24 augmented CONTROL FAILED"
    assert res24a['missed1937'] == 0
    print("  [CONTROL PASS] V=24 augmented: 38611/38924/41656, OPEN 0")

    res32a = corpus_run(5, 32, 6, tabs5, cap5, sigcap, AUG_P4)
    s32a, f32a, v32a, o32a = report(res32a, 5, 32, "control-V32", "AUGMENTED")
    assert (s32a, f32a, o32a, res32a['total']) == \
        (173337, 176376, 0, 196751), "V=32 augmented CONTROL FAILED"
    assert res32a['missed1937'] == 0
    print("  [CONTROL PASS] V=32 augmented: 173337/176376/196751, OPEN 0")

    res40a = corpus_run(5, 40, 6, tabs5, cap5, sigcap, AUG_P4)
    s40a, f40a, v40a, o40a = report(res40a, 5, 40, "control-V40", "AUGMENTED")
    assert (s40a, f40a, o40a, res40a['total']) == \
        (526607, 540415, 0, 641166), "V=40 augmented CONTROL FAILED"
    assert res40a['missed1937'] == 0
    print("  [CONTROL PASS] V=40 augmented: 526607/540415/641166, OPEN 0")
    unflagged_decomposition(res40a, "V=40 augmented", augmented=True)

    print("\n===== PART E: AUGMENTED battery PRIMARY V=48 =====")
    res48a = corpus_run(5, 48, 6, tabs5, cap5, sigcap, AUG_P4)
    s48a, f48a, v48a, o48a = report(res48a, 5, 48, "primary-V48", "AUGMENTED")
    unflagged_decomposition(res48a, "V=48 augmented", augmented=True)
    assert res48a['missed1937'] == 0

    # ---------------- Part F: trends + the reviewer's decision rule
    print("\n" + "=" * 72)
    print("===== TREND: SAME battery (P4 at {7,13,17}) =====")
    print("=" * 72)
    print("  V    sets    ports      +D3      +V(any)   OPEN")
    for VV, r in ((16, res16), (24, res24), (32, res32), (40, res40s),
                  (48, res48s)):
        t = r['total']
        c = r['closed']
        st = c['strict']
        fu = st + c['P5only']
        av = fu + c['Vonly']
        print("  %2d  %6d  %6d %6.2f%%  %6.2f%%  %6.2f%%  %4d"
              % (VV, t, st, 100.0 * st / t, 100.0 * fu / t,
                 100.0 * av / t, c['OPEN']))

    print("\n" + "=" * 72)
    print("===== TREND: AUGMENTED battery (P4 at {7,13,17,19,37}) =====")
    print("=" * 72)
    print("  V    sets    ports      +D3      +V(any)   OPEN")
    aug_rows = []
    for VV, r in ((24, res24a), (32, res32a), (40, res40a), (48, res48a)):
        t = r['total']
        c = r['closed']
        st = c['strict']
        fu = st + c['P5only']
        av = fu + c['Vonly']
        aug_rows.append((VV, t, st, fu))
        print("  %2d  %6d  %6d %6.2f%%  %6.2f%%  %6.2f%%  %4d"
              % (VV, t, st, 100.0 * st / t, 100.0 * fu / t,
                 100.0 * av / t, c['OPEN']))

    # reviewer's decision rule on the V=48 augmented number
    aug48_d3 = 100.0 * aug_rows[-1][3] / aug_rows[-1][1]
    aug48_ports = 100.0 * aug_rows[-1][2] / aug_rows[-1][1]
    print("\n  V=48 augmented: ports %.2f%%, +D3 %.2f%%"
          % (aug48_ports, aug48_d3))
    print("  decay per 8 V-units (augmented +D3): "
          "%.2f -> %.2f -> %.2f pts"
          % (100.0 * aug_rows[0][3] / aug_rows[0][1]
             - 100.0 * aug_rows[1][3] / aug_rows[1][1],
             100.0 * aug_rows[1][3] / aug_rows[1][1]
             - 100.0 * aug_rows[2][3] / aug_rows[2][1],
             100.0 * aug_rows[2][3] / aug_rows[2][1] - aug48_d3))
    if aug48_d3 > 80.0:
        print("\n  DECISION (reviewer's rule): augmented +D3 at V=48 is "
              "ABOVE 80% ->\n  decay slower than extrapolated; paper "
              "framing CONFIDENT.")
    elif aug48_d3 < 78.0:
        print("\n  DECISION (reviewer's rule): augmented +D3 at V=48 is "
              "BELOW 78% ->\n  decay on the extrapolated track; paper "
              "framing CAUTIOUS.")
    else:
        print("\n  DECISION (reviewer's rule): augmented +D3 at V=48 is in "
              "the 78-80%\n  band -> borderline; report both readings, "
              "frame cautiously-optimistic.")
    print("\n  total elapsed: %.1fs" % (time.time() - t0))

if __name__ == '__main__':
    main()
