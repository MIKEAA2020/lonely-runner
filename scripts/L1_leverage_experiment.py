#!/usr/bin/env python3
"""
L1_leverage_experiment.py — Reviewer directive 3: measure whether the cyclic-group
language admits NEW tools that close the gridonly tier (the 70 n=4 sets whose
(TAU) certificates are fully grid-specific, i.e. not at any classical breakpoint).

Cyclic language: pair (u,w), N = v_u + v_w, effective triple eff = (v_p, v_y, v_z)
residues on Z_N.  Pair certifies iff exists k with 5*||w k||_N >= N for all three,
i.e. iff the three BAD SETS  B_w = {k : 5*||wk||_N < N}  do NOT cover Z_N.
Failure = covering.  So (TAU)-4 is exactly a covering-avoidance problem on Z_N.

Parts:
  A. rebuild corpus (v<=16, gcd 1) + E3 tiers (classical/semi/gridonly/fail)
  B. gridonly sets: certifying pairs, bad-set structure (additive orders, gcds),
     argmax residues — the geometry the new tools must engage
  C. per-modulus covering classification: for N in 4..60, over ALL invertible
     triples (scale-normalized a=1), which N admit a covering?  Also 2-runner
     grid law check (B2 lemma) and mixed kernel-order coverage.
  D. battery coverage: which certifying pairs are closed by PROVED lemmas
     (B1 j-gon, B2 antipodal-coincidence [verified], B3 no-covering moduli),
     and what the leftover looks like.
"""
from math import gcd
from itertools import combinations
from collections import Counter, defaultdict

VMAX = 16

def distN(x, N):
    r = x % N
    return r if r <= N - r else N - r

def grid_max(N, eff):
    """max over k in 1..N-1 of min_w distN(w*k, N); returns (best, argmax list)."""
    best, args = -1, []
    for k in range(1, N):
        m = N
        for w in eff:
            d = distN(w * k, N)
            if d < m:
                m = d
        if m > best:
            best, args = m, [k]
        elif m == best:
            args.append(k)
    return best, args

def is_breakpoint(k, N, eff):
    for w in eff:
        if (2 * w * k) % N == 0:
            return True
    for (x, y) in combinations(eff, 2):
        if ((x + y) * k) % N == 0 or (abs(x - y) * k) % N == 0:
            return True
    return False

def unrestricted_opt3(eff):
    """Exact 3-runner unrestricted optimum (breakpoint enumeration), E3 method."""
    a, b, c = eff
    dens = sorted({d for d in (2*a, 2*b, 2*c, a+b, a+c, b+c,
                               abs(a-b), abs(a-c), abs(b-c)) if d > 0})
    bn, bd, args = 0, 1, []
    for d in dens:
        for m in range(d):
            val = min(distN(w * m, d) for w in (a, b, c))
            if val * bd > bn * d:
                bn, bd, args = val, d, [(d, m)]
            elif val * bd == bn * d:
                args.append((d, m))
    return bn, bd, args

# ---------------------------------------------------------------- Part A+B
def analyze_corpus():
    sets, gridonly = [], []
    tier = Counter()
    grid_pairs = []          # (V, pair, N, eff, m3, args, cert) for gridonly sets
    for V in combinations(range(1, VMAX + 1), 4):
        if gcd(*V) != 1:
            continue
        sets.append(V)
        per_pair, certs = [], []
        for (p, q) in combinations(range(4), 2):
            N = V[p] + V[q]
            others = [V[i] for i in range(4) if i not in (p, q)]
            eff = (V[p], others[0], others[1])
            m3, args = grid_max(N, eff)
            cert = 5 * m3 >= N
            if cert:
                certs.append((p, q))
            mu_n, mu_d, mu_args = unrestricted_opt3(eff)
            classical = any((N * m) % d == 0 for (d, m) in mu_args)
            stuck = any(w % N == 0 for w in eff)
            collide = len({w % N for w in eff}) < 3
            bp_max = max((min(distN(w * k, N) for w in eff)
                          for k in range(1, N) if is_breakpoint(k, N, eff)),
                         default=0)
            if stuck or collide:
                cls = 'degen'
            elif classical:
                cls = 'classical'
            elif cert and 5 * bp_max >= N:
                cls = 'semi'
            elif cert:
                cls = 'gridonly'
            else:
                cls = 'fail'
            per_pair.append(dict(pq=(V[p], V[q]), N=N, eff=eff, m3=m3,
                                 args=args, cert=cert, cls=cls,
                                 mu=(mu_n, mu_d)))
        if any(t['cls'] == 'classical' for t in per_pair):
            tier['classical'] += 1
        elif any(t['cls'] == 'semi' for t in per_pair if t['cert']):
            tier['semi'] += 1
        elif certs:
            tier['gridonly'] += 1
            gridonly.append(V)
            grid_pairs.append((V, per_pair))
        else:
            tier['fail'] += 1
    return sets, tier, gridonly, grid_pairs

# ---------------------------------------------------------------- Part C
def covering_scan(N, max_pairs=100000):
    """All invertible (b,c), b<=c, does A u bA u cA cover Z_N?  (a=1 fixed,
    scale-invariance.)  Returns list of covering (b,c)."""
    full = set(range(N))
    A = {k for k in range(N) if 5 * distN(k, N) < N}
    inv = [w for w in range(1, N) if gcd(w, N) == 1]
    covers = []
    for i, b in enumerate(inv):
        Bb = {k for k in range(N) if 5 * distN(b * k, N) < N}
        for c in inv[i:]:
            Bc = {k for k in range(N) if 5 * distN(c * k, N) < N}
            if (A | Bb | Bc) == full:
                covers.append((b, c))
    return covers

def two_runner_grid_law(N):
    """Check: for all a<b in Z_N^*, max_k min(||ak||,||bk||) >= N/5? (B2 lemma).
    Returns (violations, min_ratio_num, min_ratio_den)."""
    worst = (10**9, 1)
    viol = []
    for a in range(1, N):
        for b in range(a + 1, N):
            best = max(min(distN(a * k, N), distN(b * k, N))
                       for k in range(1, N))
            # ratio best/N vs 1/5 : 5*best >= N ?
            if 5 * best < N:
                viol.append((a, b, best))
            if best * worst[1] < worst[0] * N:
                worst = (best, N)
    return viol, worst

# ---------------------------------------------------------------- Part D
def battery(V, t):
    """Which proved lemmas close certifying pair t?  Returns flags."""
    N, eff, args, m3 = t['N'], t['eff'], t['args'], t['m3']
    out = []
    # B1 j-gon: j | N, j in {2,3,4,5}, all eff % j != 0  =>  m3 >= N/j >= N/5
    for j in (2, 3, 4, 5):
        if N % j == 0 and all(w % j != 0 for w in eff):
            out.append(f'B1:jgon{j}')
    # B2 antipodal coincidence: two eff residues coincide or are antipodal
    for (x, y) in combinations(eff, 2):
        if x % N == y % N or (x + y) % N == 0:
            out.append('B2:coincide')
            break
    # B3: covering impossible for this modulus class (filled by caller)
    return out

def main():
    print("== Part A: corpus and tiers (v <= %d) ==" % VMAX)
    sets, tier, gridonly, grid_pairs = analyze_corpus()
    print(f"sets: {len(sets)}   tier split: {dict(tier)}")
    print(f"gridonly sets: {len(gridonly)}")

    print("\n== Part B: gridonly certifying pairs - structure ==")
    order_counter = Counter()
    N_counter = Counter()
    cert_argmax_nonbp = 0
    cert_total = 0
    examples = []
    for (V, per_pair) in grid_pairs:
        for t in per_pair:
            if not t['cert']:
                continue
            cert_total += 1
            N, eff, args = t['N'], t['eff'], t['args']
            N_counter[N] += 1
            orders = tuple(sorted(N // gcd(w % N or N, N) if w % N else 1
                                  for w in eff))
            order_counter[orders] += 1
            if not any(is_breakpoint(k, N, eff) for k in args):
                cert_argmax_nonbp += 1
                if len(examples) < 14:
                    resid = [(w * args[0]) % N for w in eff]
                    examples.append((V, t['pq'], N, eff, m3 if (m3 := t['m3']) else 0,
                                     args[:3], [distN(r, N) for r in resid]))
    print(f"certifying pairs inside gridonly sets: {cert_total}")
    print(f"N distribution: {dict(sorted(N_counter.items()))}")
    print(f"additive-order signatures: {dict(order_counter.most_common(12))}")
    print(f"certifying argmaxes that are NOT breakpoints: {cert_argmax_nonbp}")
    print("\nexample gridonly certificates (V, pair, N, eff, m3, argmax k, dists):")
    for e in examples:
        print("  ", e)

    print("\n== Part C: per-modulus covering classification (invertible triples) ==")
    capable = {}
    for N in range(4, 61):
        cov = covering_scan(N)
        if cov:
            capable[N] = cov
    print(f"covering-capable moduli in [4,60]: {sorted(capable.keys())}")
    for N in sorted(capable):
        print(f"  N={N}: {len(capable[N])} covering (b,c) pairs (a=1); "
              f"e.g. {capable[N][:6]}")

    print("\n== Part C2: 2-runner grid law (B2 support), N in [3,40] ==")
    allviol = []
    worst_global = (10**9, 1)
    for N in range(3, 41):
        viol, worst = two_runner_grid_law(N)
        allviol += [(N,) + v for v in viol]
        if worst[0] * worst_global[1] < worst_global[0] * worst[1]:
            worst_global = worst
    print(f"violations of max_k min(||ak||,||bk||) >= N/5 : {len(allviol)}"
          + (f"  {allviol[:6]}" if allviol else ""))
    print(f"worst 2-runner grid ratio observed: {worst_global[0]}/{worst_global[1]}"
          f" = {worst_global[0]/worst_global[1]:.4f}")

    print("\n== Part D: battery coverage of gridonly certifying pairs ==")
    covering_capable = set(capable.keys())
    batt = Counter()
    leftover = []
    for (V, per_pair) in grid_pairs:
        closed_set = False
        for t in per_pair:
            if not t['cert']:
                continue
            N, eff = t['N'], t['eff']
            flags = battery(V, t)
            all_inv = all(gcd(w % N, N) == 1 for w in eff if w % N != 0) and \
                      all(w % N != 0 for w in eff)
            if all_inv and N not in covering_capable:
                flags.append('B3:nocover-modulus')
            if flags:
                closed_set = True
                batt[flags[0]] += 1
            else:
                leftover.append((V, t['pq'], N, eff, t['m3'], t['args'][:2]))
        batt['SET-CLOSED' if closed_set else 'SET-OPEN'] += 1
    print(f"battery tally over certifying pairs: {dict(batt)}")
    print(f"leftover (no lemma applies) certifying pairs: {len(leftover)}")
    for lv in leftover[:15]:
        print("   LEFTOVER:", lv)
    return gridonly, capable

if __name__ == "__main__":
    main()
