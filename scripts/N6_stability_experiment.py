#!/usr/bin/env python3
"""
N6_stability_experiment.py — the V=32/40 stability extension (the reviewer's
decisive experiment, amended pivot condition).

Run the SAME n=5 pipeline — battery IDENTICAL to the committed
scripts/N5_leverage_experiment.py (P1 j-gon / P2 counting / P4 at
{7,13,17} / P5 = Lemma D3 / V1 safe modulus / V2 safe signature) — at
V = 32 and V = 40, with hard controls at V = 16/24 and the n = 4 control.

Amended decision rule: continue the program if the rigid zone is provably
finite or the residual rate stays bounded as V grows.  Proved-coverage
stabilizes above 80% -> the reformulation has content at n=5, paper
writable; keeps dropping toward 60% -> base-rung-only, pivot to R2.

New machinery required by the range (NOT a battery change):
  - tables / capability / signature-capability maps for ALL N <= 79
    (pair sums reach 2*40-1 = 79);
  - the capability completion for 48..79 is genuinely open input: the
    committed "odd composites safe" argument (unit bad sets cannot reach
    non-unit classes when every j <= h has gcd(j,N)=1) fails from N = 49
    on (h = 8 >= 7), so odd composites 49..77 must be scanned, not assumed.
    Primes 53..79 are re-verified against the committed (47,150] scan.
  - independent set-based cross-check of capability at the danger moduli.
No prime scan beyond 150 (reviewer directive: fix the framing, not the
range).  No attack on the no-covering conjecture.

Diagnostics added (additive only; tier logic byte-identical to committed):
  - rigid-prime pair incidence per corpus;
  - cls1937 candidate sets (sets whose unflagged certifying pairs include
    an all-invertible pair at N in {19,37} — exactly the sets the K19/K37
    classifications would rescue; predicts the Phase-2 gain without
    changing this experiment's battery);
  - V1-only vs V2-only closure breakdown.
"""
from math import gcd
from itertools import combinations, combinations_with_replacement
from collections import Counter
import random
import time

CLS_MODULI = {5: (7, 11, 13), 6: (7, 13, 17)}   # SAME battery as committed
RIGID_PRIMES = (7, 13, 17, 19, 37)              # committed n=5 rigid zone

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

# ------------------------------------------------------------- class tools

def cls(x, N):
    r = x % N
    return r if r <= N - r else N - r

def cls_inv(w, N):
    return cls(pow(w, -1, N), N)

def flag_classification(N, rs, T):
    """P4: proved covering classifications. Returns flag iff NON-covering is
    proved via the class reduction (Lemma C port). Caller guarantees units."""
    if N not in CLS_MODULI[T]:
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

def battery_flags(N, rs, n, T, tabs, capable, sigcap=None):
    """All battery flags for a pair with residues rs. Soundness: every flag
    must imply certification (checked by the sweeps)."""
    tab = tabs[N]
    invs = tab['invs']
    fl = []
    for j in range(2, T + 1):
        if N % j == 0 and all(r % j for r in rs):
            fl.append('P1:jgon%d' % j)
            break
    allinv = all(invs[r] for r in rs)
    if allinv:
        f = flag_classification(N, rs, T)
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

def classification_exactness(tabs, T, n):
    """At each P4 modulus: predicate(non-cover) == mask(non-cover) for ALL
    unit k-multisets."""
    k = n - 1
    mism = []
    for N in CLS_MODULI[T]:
        tab = tabs[N]
        masks, full = tab['masks'], tab['full']
        units = [w for w in range(1, N) if tab['invs'][w]]
        for rs in combinations_with_replacement(units, k):
            mm = 0
            for r in rs:
                mm |= masks[r]
            noncover_mask = (mm != full)
            noncover_pred = (flag_classification(N, rs, T) is not None)
            if noncover_mask != noncover_pred:
                mism.append((N, rs))
    return mism

def soundness_sweep(tabs, capable, n, T, Nrange, exhaustive=True,
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
                fl = battery_flags(N, rs, n, T, tabs, capable, sigcap)
                if fl and mm == full:
                    viol.append((N, rs, fl))
                checked += 1
        if n_random:
            for _ in range(n_random):
                rs = tuple(sorted(rng.randrange(N) for _ in range(k)))
                mm = 0
                for r in rs:
                    mm |= masks[r]
                fl = battery_flags(N, rs, n, T, tabs, capable, sigcap)
                if fl and mm == full:
                    viol.append((N, rs, fl))
                checked += 1
    return viol, checked

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

# ------------------------------------------------------------ corpus run

def corpus_run(n, VMAX, T, tabs, capable, sigcap=None, open_cap=30):
    k = n - 1
    total = 0
    closed = Counter()          # strict / P5only / Vonly / OPEN
    jgon_plain = Counter()
    pairstats = Counter()
    cert_hist = Counter()
    ground_fail = []
    open_sets = []
    open_pair_stats = Counter()  # (N, allinv, capable, ordersig) -> count
    all_unflagged = Counter()    # same, over ALL sets (corpus-wide)
    p1_witness_fail = 0
    m4_mismatch = 0
    open_plain = 0
    # --- additive diagnostics (tier logic identical to committed) ---
    danger_pairs = Counter()     # (N, allinv) for certifying pairs at rigid primes
    cls1937_sets = 0             # sets w/ >=1 unflagged allinv cert pair at 19/37
    v1_closed = 0                # V-tier sets with a V1-flagged pair
    v2_only_closed = 0           # V-tier sets with no V1-flagged pair
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
            fl = battery_flags(N, rs, n, T, tabs, capable, sigcap)
            allinv = all(tab['invs'][r] for r in rs)
            if N in RIGID_PRIMES:
                danger_pairs[(N, allinv)] += 1
            for f in fl:
                pairstats[f] += 1
            for f in fl:
                if f.startswith('P1:jgon'):
                    j = int(f[7:])
                    if (mm >> (N // j)) & 1:
                        p1_witness_fail += 1   # constructive witness failed
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
                cls1937_sets=cls1937_sets, v1_closed=v1_closed,
                v2_only_closed=v2_only_closed,
                p1_witness_fail=p1_witness_fail, m4_mismatch=m4_mismatch)

def report(res, n, VMAX, tag):
    t = res['total']
    c = res['closed']
    strict = c['strict']
    full = c['strict'] + c['P5only']
    anyv = full + c['Vonly']
    print("\n---- corpus %s (n=%d, v<=%d, %d sets) ----" % (tag, n, VMAX, t))
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
    print("  cls1937 candidate sets (unflagged allinv pair at 19/37): %d"
          % res['cls1937_sets'])
    print("  ground-truth failures (sets w/o cert pair): %d"
          % len(res['ground_fail']))
    print("  soundness: P1 witness fails %d, m4 mismatches %d"
          % (res['p1_witness_fail'], res['m4_mismatch']))
    return strict, full, anyv, c['OPEN']

def unflagged_decomposition(res, label):
    """Decompose the corpus-wide unflagged certifying pairs by modulus class."""
    au = res['all_unflagged']
    tot = sum(au.values())
    by = Counter()
    for (N, allinv, cap, sig), c in au.items():
        if N in (19, 37) and allinv:
            by['P4-target: allinv @ 19/37'] += c
        elif allinv:
            by['allinv @ capable %d' % N] += c
        elif N >= 48:
            by['kernel @ new composite N>=48'] += c
        else:
            by['kernel @ N<=47'] += c
    print("  [%s] unflagged certifying pairs corpus-wide: %d" % (label, tot))
    for k, v in by.most_common():
        print("      %-34s %7d" % (k, v))
    return tot

# ------------------------------------------------------------------- main

def main():
    t0 = time.time()
    print("=" * 72)
    print("N6 STABILITY EXPERIMENT — V=32/40 extension of the n=5 pipeline")
    print("  battery IDENTICAL to committed N5 (P4 at {7,13,17})")
    print("=" * 72)

    # ---------------- tables
    print("[tables] building n=4 (N<=31) and n=5 (N<=79) tables...")
    tabs4 = {N: build_tables(N, 5) for N in range(3, 32)}
    tabs5 = {N: build_tables(N, 6) for N in range(3, 80)}
    bad4 = [b for N in range(3, 32) for b in check_size_formula(N, 5, tabs4[N])]
    bad5 = [b for N in range(3, 80) for b in check_size_formula(N, 6, tabs5[N])]
    print("[check] Lemma S5 exact-size violations: n4 %d, n5(N<=79) %d"
          % (len(bad4), len(bad5)))

    # ---------------- capability maps (ALL N <= 79)
    t1 = time.time()
    cov4 = {N: capable_units(N, tabs4[N], 3) for N in range(3, 32)}
    cap4 = {N: bool(v) for N, v in cov4.items()}
    cov5 = {N: capable_units(N, tabs5[N], 4) for N in range(3, 80)}
    cap5 = {N: bool(v) for N, v in cov5.items()}
    cap5_list = sorted(N for N in cap5 if cap5[N])
    # controls against the committed record
    ctrl_n4 = sorted(N for N in range(3, 32) if cap4[N])
    ctrl_n5_47 = sorted(N for N in range(3, 48) if cap5[N])
    print("[n=4] unit-triple capable moduli (T=5): %s" % ctrl_n4)
    assert ctrl_n4 == [7, 11, 13], "n=4 capability control FAILED"
    print("[n=5] unit 4-tuple capable moduli, N<=47 (control): %s" % ctrl_n5_47)
    assert ctrl_n5_47 == [7, 13, 17, 19, 37], "n=5 capability control FAILED"
    primes_53_79 = [N for N in (53, 59, 61, 67, 71, 73, 79) if cap5[N]]
    print("[n=5] capable primes in (47,79] (must be [] — committed scan): %s"
          % primes_53_79)
    assert primes_53_79 == [], "conflict with committed (47,150] prime scan"
    new_cap = [N for N in cap5_list if N > 47]
    print("[n=5] *** capable moduli in 48..79 (NEW RANGE, composites open): %s"
          " ***" % new_cap)
    print("[n=5] covering-configuration counts per capable modulus: %s"
          % {N: len(cov5[N]) for N in cap5_list})

    # ---- odd-composite structural table (explanation, not just verification)
    print("\n[odd composites 49..77] structure vs capability (T=6):")
    print("      N    h   nonunitCls reach  unitSlots/4set  capable  #cov")
    for N in [49, 51, 55, 57, 63, 65, 69, 75, 77]:
        st = odd_composite_structure(N)
        print("     %3d  %2d  %5d  %5d  %2d+%2d  %8s  %4d"
              % (N, st['h'], st['n_nonunit_cls'], st['nonunit_reachable'],
                 st['unit_slots_per_set'], st['nonunit_slots_per_set'],
                 cap5[N], len(cov5[N])))

    # ---- independent set-based cross-check at the danger moduli
    print("\n[cross-check] independent set-based capability scan:")
    for N in [19, 37, 49, 55, 63, 65, 69, 75, 77, 51, 57]:
        cov_sb = capable_setbased(N)
        ok = (len(cov_sb) == len(cov5[N]))
        print("      N=%3d: set-based %d vs mask-based %d coverings  %s"
              % (N, len(cov_sb), len(cov5[N]), "AGREE" if ok else "MISMATCH"))
        assert ok, "set-based/mask-based capability mismatch at N=%d" % N
    print("      (elapsed %.1fs)" % (time.time() - t1))

    # ---------------- signature capability map (all N <= 79)
    t1 = time.time()
    print("\n[sigcap] building order-signature capability maps, N=3..79...")
    sigcap = {}
    for N in range(3, 80):
        sigcap[N] = signature_capability(N, tabs5[N])
    nsig_new = {N: len(sigcap[N]) for N in range(48, 80) if sigcap[N]}
    nsig_old = {N: len(sigcap[N]) for N in range(5, 48) if sigcap[N]}
    print("[n=5] capable signatures, N in 5..47 (control, %d moduli): %s"
          % (len(nsig_old), nsig_old))
    print("[n=5] capable signatures, N in 48..79 (NEW): %s" % nsig_new)
    print("      (elapsed %.1fs)" % (time.time() - t1))

    # ---------------- soundness sweeps
    t1 = time.time()
    mism4 = classification_exactness(tabs4, 5, 4)
    mism5 = classification_exactness(tabs5, 6, 5)
    print("\n[check] P4 classification exactness mismatches: n4 %d, n5 %d"
          % (len(mism4), len(mism5)))
    v4, c4 = soundness_sweep(tabs4, cap4, 4, 5, range(3, 32))
    v5a, c5a = soundness_sweep(tabs5, cap5, 5, 6, range(3, 32), sigcap=sigcap)
    v5b, c5b = soundness_sweep(tabs5, cap5, 5, 6, range(32, 80),
                               exhaustive=False, n_random=4000,
                               sigcap=sigcap)
    v5c, c5c = soundness_sweep(tabs5, cap5, 5, 6, [49], exhaustive=True,
                               sigcap=sigcap)
    print("[check] flag=>cert sweep: n4 %d/%d viol; n5(3..31) %d/%d viol; "
          "n5(32..79 rand) %d/%d viol; n5(49 exhaustive) %d/%d viol"
          % (len(v4), c4, len(v5a), c5a, len(v5b), c5b, len(v5c), c5c))
    print("      (elapsed %.1fs)" % (time.time() - t1))

    # ---------------- Part A: n=4 positive control
    print("\n===== PART A: n=4 positive control (v<=16) =====")
    res4 = corpus_run(4, 16, 5, tabs4, cap4)
    s4, f4, v4c, o4 = report(res4, 4, 16, "control")
    assert (s4, res4['total'], o4) == (1651, 1745, 1), \
        "n=4 CONTROL FAILED: got %s" % ((s4, res4['total'], o4),)
    assert res4['open_sets'][0][0] == (3, 5, 8, 13), \
        "n=4 control open set mismatch"
    print("  [CONTROL PASS] n=4 reproduces committed: 1651/1745, "
          "open set = (3,5,8,13)")

    # ---------------- Part B: n=5 controls, then primaries
    print("\n===== PART B: n=5 control V=16 =====")
    res16 = corpus_run(5, 16, 6, tabs5, cap5, sigcap=sigcap)
    s16, f16, v16, o16 = report(res16, 5, 16, "control-V16")
    assert res16['total'] == 4311, "V=16 control total mismatch"
    print("  [control] committed matched-range: 4179 (96.94%%), 4193 (97.26%%)"
          " -> got %d (%.2f%%), %d (%.2f%%)"
          % (s16, 100.0 * s16 / 4311, f16, 100.0 * f16 / 4311))

    print("\n===== PART C: n=5 control V=24 (must reproduce committed exactly) =====")
    res24 = corpus_run(5, 24, 6, tabs5, cap5, sigcap=sigcap)
    s24, f24, v24, o24 = report(res24, 5, 24, "control-V24")
    # committed tiers: ports 35778; +D3 36362; +V1 41240 (= 36362 + V-tier
    # sets having a V1-flagged pair); +V2 41656; OPEN 0; total 41656.
    assert (s24, f24, o24, res24['total']) == (35778, 36362, 0, 41656), \
        "V=24 CONTROL FAILED: got %s" % ((s24, f24, o24, res24['total']),)
    assert f24 + res24['v1_closed'] == 41240, \
        "V=24 +V1 tier FAILED: %d + %d != 41240" % (f24, res24['v1_closed'])
    assert f24 + res24['v1_closed'] + res24['v2_only_closed'] == 41656
    print("  [CONTROL PASS] V=24 reproduces committed exactly: 35778/36362/"
          "41240/41656, OPEN 0 (V-tier: %d V1-closed + %d V2-only)"
          % (res24['v1_closed'], res24['v2_only_closed']))
    unflagged_decomposition(res24, "V=24")

    print("\n===== PART D: n=5 PRIMARY V=32 =====")
    res32 = corpus_run(5, 32, 6, tabs5, cap5, sigcap=sigcap)
    s32, f32, v32, o32 = report(res32, 5, 32, "primary-V32")
    unflagged_decomposition(res32, "V=32")

    print("\n===== PART E: n=5 PRIMARY V=40 =====")
    res40 = corpus_run(5, 40, 6, tabs5, cap5, sigcap=sigcap)
    s40, f40, v40, o40 = report(res40, 5, 40, "primary-V40")
    unflagged_decomposition(res40, "V=40")

    # ---------------- Part F: trend + verdict
    print("\n" + "=" * 72)
    print("===== TREND (n=5, SAME battery, primitive sets) =====")
    print("=" * 72)
    print("  V    sets    ports      +D3      +V(any)   OPEN")
    rows = []
    for VV, r in ((16, res16), (24, res24), (32, res32), (40, res40)):
        t = r['total']
        c = r['closed']
        st = c['strict']
        fu = st + c['P5only']
        av = fu + c['Vonly']
        rows.append((VV, t, st, fu, av, c['OPEN']))
        print("  %2d  %6d  %6d %6.2f%%  %6.2f%%  %6.2f%%  %4d"
              % (VV, t, st, 100.0 * st / t, 100.0 * fu / t,
                 100.0 * av / t, c['OPEN']))
    print("\n  matched-range note: V=16 row is the OVERLAP REGION only "
          "(4,311 sets);")
    print("  the V=24/32/40 rows carry the honest range-growth statement.")
    p40 = 100.0 * rows[3][2] / rows[3][1]
    f40p = 100.0 * rows[3][3] / rows[3][1]
    drop24_40 = 100.0 * rows[1][2] / rows[1][1] - p40
    print("\n  ports: V24 %.2f%% -> V32 %.2f%% -> V40 %.2f%%  "
          "(total drop 24->40: %.2f pts)"
          % (100.0 * rows[1][2] / rows[1][1],
             100.0 * rows[2][2] / rows[2][1], p40, drop24_40))
    print("  +D3 : V24 %.2f%% -> V32 %.2f%% -> V40 %.2f%%"
          % (100.0 * rows[1][3] / rows[1][1],
             100.0 * rows[2][3] / rows[2][1], f40p))
    print("  cls1937-rescuable sets (predicted K19/K37 gain): "
          "V24 %d, V32 %d, V40 %d"
          % (res24['cls1937_sets'], res32['cls1937_sets'],
             res40['cls1937_sets']))
    if p40 >= 80.0:
        print("\n  VERDICT (same battery): proved ports STABILIZE above 80%% "
              "at V=40 ->\n  per the amended condition the reformulation has "
              "content at n=5; paper writable.")
    elif p40 >= 60.0:
        print("\n  VERDICT (same battery): ports between 60 and 80 at V=40 — "
              "marginal;\n  augment with the K19/K37 classifications before "
              "deciding.")
    else:
        print("\n  VERDICT (same battery): ports below 60%% at V=40 — "
              "degrading toward the pivot threshold.")

    print("\n  total elapsed: %.1fs" % (time.time() - t0))

if __name__ == '__main__':
    main()
