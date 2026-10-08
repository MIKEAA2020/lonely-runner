#!/usr/bin/env python3
"""
N5_leverage_experiment.py — THE decisive experiment (reviewer rec. 1):
at n = 5, enumerate all primitive speed sets (v <= 24, the E2b base, 41656
sets), run the cyclic-lemma battery ported from n = 4, and measure the
fraction of sets closed by PROVED lemmas alone.

Decision rule (reviewer rec. 4): >= 80% closure -> continue the reformulation
program; < 60% -> pivot to R2 (character/kernel track) with the reformulation
as a bridge, not a program.

Framework (n speeds, T = n + 1, k = n - 1 effective residues per pair):
  pair (p, q), N = v_p + v_q; v_q = -v_p (mod N) shares the bad set of v_p.
  eff residues rs = (v_p, others) mod N.  Pair certifies (tau >= 1/T)
  iff some k has min_w ||w k||_N >= N/T
  iff the k bad sets  B_w = {k in Z_N : T*||w k||_N < N}  do NOT cover Z_N.
  N <= 47 < 64, so covering checks are bitmask ORs.

Battery (status per item; PROVED items have proofs in the note):
  P1  j-gon (Lemma J port, PROVED): j | N, 2 <= j <= T, all rs != 0 mod j
      => k = N/j certifies.  (n=4: j in 2..5; n=5: j in 2..6.)
  P2  counting (B0/B5 port + Lemma S5 exact sizes, PROVED):
      |B_w| = (N / ord(w)) * (2*floor((ord(w)-1)/T) + 1) exactly;
      sum of the k sizes - (k - 1) < N => no covering.
  P4  proved covering classifications at small moduli (PROVED):
      n=4: N in {7, 11, 13} (ported from L2, re-verified in-script);
      n=5: N in {7, 13, 17} (new this turn: class / domino / matching proofs).
      Predicate (Lemma C port): h = floor((N-1)/T); the classes of B_w are
      {cls(j*w^{-1}) : 1 <= j <= h}; covering iff they cover C_{(N-1)/2}.
  P5  Lemma D3 (NEW at n=5, PROVED; sound replacement of the dead B2 port):
      if the k residues span pm <= k-1 antipodal classes (a +/- coincidence),
      all invertible, and pm*(2h+1) - (pm-1) < N, then no covering.
  V1  per-modulus VERIFIED: all residues invertible and N not covering-capable
      (exhaustive scale-normalized unit enumeration at that N) => certifies.

Committed comparisons:
  - n=4 positive control (v <= 16): reproduce the committed coverage numbers
    (proved-only 1651/1745 = 94.6% per download/leverage_experiment.md §3).
  - n=5 ground truth: every set has >= 1 certifying pair (E2b: 41656/41656).

Also measured: the n=5 rigid zone (covering-capable moduli), the
order-signature capability map, 3-set coverings (B2-port falsification),
the plain-regime share, and full diagnostics of the sets left open.
"""
from math import gcd
from itertools import combinations, combinations_with_replacement
from collections import Counter
import random
import time

CLS_MODULI = {5: (7, 11, 13), 6: (7, 13, 17)}

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

def triple_cover_moduli(tabs, Nmax):
    """N -> example triple of DISTINCT residues whose bad sets cover Z_N
    (falsifier scan for the n=4 B2 coincidence port at n=5)."""
    out = {}
    for N in range(3, Nmax + 1):
        tab = tabs[N]
        masks, full = tab['masks'], tab['full']
        res = list(range(1, N))
        found = None
        for i in range(len(res)):
            ma = masks[res[i]]
            for jj in range(i + 1, len(res)):
                mb = ma | masks[res[jj]]
                for kk in range(jj + 1, len(res)):
                    if mb | masks[res[kk]] == full:
                        found = (res[i], res[jj], res[kk])
                        break
                if found:
                    break
            if found:
                break
        if found:
            out[N] = found
    return out

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

def lift_law_check(tabs5, trials=3000, seed=3):
    """(ma,mb,mc,md) covers mN' iff (a,b,c,d) covers N' — sample check."""
    rng = random.Random(seed)
    bad = 0
    done = 0
    for _ in range(trials):
        Np = rng.randrange(3, 16)
        m = rng.choice([2, 3, 4])
        N = m * Np
        if N > 47:
            continue
        done += 1
        tabN, tabNp = tabs5[N], tabs5[Np]
        ws = [rng.randrange(1, Np) for _ in range(4)]
        mm = 0
        for w in ws:
            mm |= tabN['masks'][(m * w) % N]
        base = 0
        for w in ws:
            base |= tabNp['masks'][w % Np]
        if (mm == tabN['full']) != (base == tabNp['full']):
            bad += 1
    return bad, done

def capable_prime_scan(Nmax):
    """Targeted extension: unit 4-tuple covering capability at primes in
    (47, Nmax] — for the rigid-zone open problem (n=4 analogue was the
    [4,60] scan behind the no-covering conjecture)."""
    out = {}
    for N in range(53, Nmax + 1):
        if any(N % p == 0 for p in range(2, int(N ** 0.5) + 1)):
            continue
        tab = build_tables(N, 6)
        cov = capable_units(N, tab, 4)
        if cov:
            out[N] = (len(cov), cov[0])
    return out

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
    rng = random.Random(7)
    for V in combinations(range(1, VMAX + 1), n):
        g = 0
        for x in V:
            g = gcd(g, x)
        if g != 1:
            continue
        total += 1
        certs = 0
        got = {'s': False, 'f': False, 'v': False, 'j': False}
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
            if any(f.startswith('P1') for f in fl):
                got['j'] = True
            if not fl:
                sig = tuple(sorted(tab['orders'][r] for r in rs))
                allinv = all(tab['invs'][r] for r in rs)
                key = (N, allinv, capable[N] if allinv else None, sig)
                set_unflagged.append(key)
                all_unflagged[key] += 1
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
        if got['s']:
            closed['strict'] += 1
        elif got['f']:
            closed['P5only'] += 1
        elif got['v']:
            closed['Vonly'] += 1
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
                open_plain=open_plain,
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
    print("  closed by PROVED+VERIFIED (+V1):           %6d  (%.2f%%)"
          % (anyv, 100.0 * anyv / t))
    print("  OPEN:                                      %6d  (%.2f%%)"
          % (c['OPEN'], 100.0 * c['OPEN'] / t))
    print("  j-gon/plain (any-pair) split: %s" % dict(res['jgon_plain']))
    print("  certifying-pairs-per-set histogram: %s"
          % dict(sorted(res['cert_hist'].items())))
    print("  per-flag pair counts: %s"
          % dict(sorted(res['pairstats'].items())))
    print("  ground-truth failures (sets w/o cert pair): %d"
          % len(res['ground_fail']))
    print("  soundness: P1 witness fails %d, m4 mismatches %d"
          % (res['p1_witness_fail'], res['m4_mismatch']))
    return strict, full, anyv, c['OPEN']

# ------------------------------------------------------------------- main

def main():
    t0 = time.time()
    print("=" * 72)
    print("N5 LEVERAGE EXPERIMENT — the decisive n=5 measurement")
    print("=" * 72)

    # ---------------- tables
    tabs4 = {N: build_tables(N, 5) for N in range(3, 32)}
    tabs5 = {N: build_tables(N, 6) for N in range(3, 48)}
    bad4 = [b for N in range(3, 32) for b in check_size_formula(N, 5, tabs4[N])]
    bad5 = [b for N in range(3, 48) for b in check_size_formula(N, 6, tabs5[N])]
    print("[check] Lemma S5 exact-size formula violations: n4 %d, n5 %d"
          % (len(bad4), len(bad5)))

    # ---------------- capability maps
    cov4 = {N: capable_units(N, tabs4[N], 3) for N in range(3, 32)}
    cap4 = {N: bool(v) for N, v in cov4.items()}
    print("[n=4] unit-triple covering-capable moduli (T=5): %s"
          % sorted(N for N in cap4 if cap4[N]))
    cov5 = {N: capable_units(N, tabs5[N], 4) for N in range(3, 48)}
    cap5 = {N: bool(v) for N, v in cov5.items()}
    cap5_list = sorted(N for N in cap5 if cap5[N])
    print("[n=5] unit 4-tuple covering-capable moduli (T=6): %s" % cap5_list)
    print("[n=5] covering configuration counts per capable modulus: %s"
          % {N: len(cov5[N]) for N in cap5_list})
    for N in cap5_list[:6]:
        print("        e.g. N=%d: %s" % (N, cov5[N][:3]))

    # hand-verified example from the note: eff (1,4,6,7) covers Z_17
    m17 = tabs5[17]['masks']
    u = m17[1] | m17[4] | m17[6] | m17[7]
    print("[check] Z_17 example (1,4,6,7) covers: %s"
          % (u == tabs5[17]['full']))

    # ---------------- 3-set coverings (B2 port falsification)
    tc5 = triple_cover_moduli(tabs5, 47)
    print("[n=5] 3-distinct-set covering moduli (B2-port world): %s"
          % {N: e for N, e in sorted(tc5.items())})

    # ---------------- signature capability map (n=5)
    t1 = time.time()
    sigcap = {}
    for N in range(3, 48):
        sigcap[N] = signature_capability(N, tabs5[N])
    nsig = {N: len(sigcap[N]) for N in range(5, 48)}
    print("[n=5] order-signature capability (nonzero residues), per N: %s" % nsig)
    anysig = sorted(N for N in range(5, 48) if sigcap[N])
    print("[n=5] N with ANY capable signature (unit + kernel world): %s"
          % anysig)
    print("[n=5] capable signatures at N in 7..24:")
    for N in range(7, 25):
        if sigcap[N]:
            print("        N=%d: %s" % (N, sorted(sigcap[N].keys())))
    print("        (elapsed %.1fs)" % (time.time() - t1))

    # ---------------- soundness sweeps
    mism4 = classification_exactness(tabs4, 5, 4)
    mism5 = classification_exactness(tabs5, 6, 5)
    print("[check] P4 classification exactness mismatches: n4 %d, n5 %d"
          % (len(mism4), len(mism5)))
    v4, c4 = soundness_sweep(tabs4, cap4, 4, 5, range(3, 32))
    v5a, c5a = soundness_sweep(tabs5, cap5, 5, 6, range(3, 32), sigcap=sigcap)
    v5b, c5b = soundness_sweep(tabs5, cap5, 5, 6, range(32, 48),
                               exhaustive=False, n_random=4000,
                               sigcap=sigcap)
    print("[check] flag=>cert sweep: n4 %d/%d viol, n5(3..31) %d/%d viol, "
          "n5(32..47 rand) %d/%d viol"
          % (len(v4), c4, len(v5a), c5a, len(v5b), c5b))
    lb, ld = lift_law_check(tabs5)
    print("[check] lift law sample: %d violations / %d trials" % (lb, ld))

    # ---------------- Part A: n=4 positive control
    print("\n===== PART A: n=4 positive control (v<=16) =====")
    res4 = corpus_run(4, 16, 5, tabs4, cap4)
    s4, f4, v4c, o4 = report(res4, 4, 16, "control")
    print("  committed comparison: proved-only 1651/1745 (94.6%%), "
          "full battery leaves 1 set open")
    if res4['open_sets']:
        print("  control OPEN set: V=%s (committed open set was (3,5,8,13))"
              % (res4['open_sets'][0][0],))

    # ---------------- Part B: n=5 corpora
    print("\n===== PART B: n=5 primary corpus (v<=24) =====")
    res5 = corpus_run(5, 24, 6, tabs5, cap5, sigcap=sigcap)
    s5, f5, v5c, o5 = report(res5, 5, 24, "primary")
    print("\n===== PART B2: n=5 comparability slice (v<=16) =====")
    res5b = corpus_run(5, 16, 6, tabs5, cap5, sigcap=sigcap)
    s5b, f5b, v5b2, o5b = report(res5b, 5, 16, "slice")

    # ---------------- Part C: diagnostics of the open tier
    print("\n===== PART C: open-tier diagnostics (primary corpus) =====")
    ops = res5['open_pair_stats']
    tot_op = sum(ops.values())
    print("  unflagged certifying pairs in OPEN sets: %d" % tot_op)
    print("  unflagged certifying pairs corpus-wide: %d"
          % sum(res5['all_unflagged'].values()))
    print("  OPEN sets with NO j-gon pair anywhere (plain): %d / %d"
          % (res5['open_plain'], res5['closed']['OPEN']))
    nHist = Counter()
    invHist = Counter()
    sigHist = Counter()
    for (N, allinv, cap, sig), c in ops.items():
        nHist[N] += c
        invHist[(allinv, cap)] += c
        sigHist[sig] += c
    print("  N histogram (OPEN sets only): %s" % dict(sorted(nHist.items())))
    print("  (all-invertible?, modulus capable?) -> count: %s"
          % dict(sorted(invHist.items(), key=str)))
    print("  top order signatures (OPEN sets only): %s"
          % sigHist.most_common(14))
    print("  top (N, sig) blocks (OPEN sets only): %s"
          % Counter({(kk[0], kk[3]): v for kk, v in ops.items()})
          .most_common(14))
    print("  first OPEN sets (with unflagged certifying pairs):")
    for (V, det) in res5['open_sets'][:12]:
        print("    V=%s" % (V,))
        for d in det:
            print("        pair=(%d,%d) N=%d rs=%s orders=%s" % d)

    # ---------------- Part D: decision metrics
    print("\n===== PART D: DECISION METRICS =====")
    t5 = res5['total']
    print("  n=4 control (proved-only, apples-to-apples): "
          "%d/%d = %.2f%% (committed: 94.6%%)" % (s4, res4['total'],
          100.0 * s4 / res4['total']))
    print("  n=5 PRIMARY (v<=24, %d sets):" % t5)
    print("     proved ports only (P1+P2+P4):   %6d = %5.2f%%"
          % (s5, 100.0 * s5 / t5))
    print("     proved + new lemma P5 (D3):      %6d = %5.2f%%"
          % (f5, 100.0 * f5 / t5))
    print("     proved + verified (V1):          %6d = %5.2f%%"
          % (v5c, 100.0 * v5c / t5))
    print("  reviewer thresholds: >= 80%% continue program, < 60%% pivot to R2")

    # ---------------- targeted extension: capable primes beyond the corpus
    t2 = time.time()
    ext = capable_prime_scan(150)
    print("[n=5 extension] unit-capable primes in (47,150]: %s"
          % {N: c for N, (c, e) in sorted(ext.items())})
    for N, (c, e) in sorted(ext.items())[:8]:
        print("        N=%d: %d coverings, e.g. %s" % (N, c, e))
    print("  (extension elapsed %.1fs; total %.1fs)"
          % (time.time() - t2, time.time() - t0))
    print("  total elapsed: %.1fs" % (time.time() - t0))

if __name__ == '__main__':
    main()
