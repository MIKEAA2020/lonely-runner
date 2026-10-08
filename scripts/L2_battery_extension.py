#!/usr/bin/env python3
"""
L2_battery_extension.py — extended lemma battery for the gridonly tier, using
ONLY proved or per-modulus-verified cyclic statements:

  B1  j-gon (proved): j | N, j in {2,3,4,5}, all eff % j != 0  => certifies
      at k = N/j with m3 >= N/j >= N/5.
  B2  antipodal coincidence (2-runner grid law: VERIFIED N<=80, proof pending):
      two eff residues equal or antipodal mod N => certifies.
  B3  no-covering modulus (per-modulus verified [4,60]): all eff invertible
      and N not in {7,11,13} => covering impossible => certifies.
  B3' covering classification at N in {7,11,13} (PROVED, code-verified):
      covering occurs iff the antipodal-class condition fails:
        N=7 : {class(a),class(b),class(c)} = all of C_3
        N=11: classes-as-C_5-subsets in {{0,1,3},{0,2,3},{0,2,4}} (a=1)
        N=13: {classes} = parity class {X, Xg^2, Xg^4} of C_6
      => pair certifies iff condition FAILS.
  B5  subgroup counting (proved): an eff residue of additive order 5 makes
      B_w = 5Z_N (size N/5); union bound kills covering => certifies.
      Order 4: killed when N/4 + 2(2h+1) - 2 < N (checked per N).
  B6  order-6+ kernel residues: B_w is a coset-union; NOT closed (residual).

Also: verify the Z_7/11/13 classification claims independently (counts and
membership), verify B5 arithmetic per N, extend the 2-runner law to N<=80,
and list the sets left open by the full battery with full structure for the
residual-gap analysis.
"""
from math import gcd
from itertools import combinations
from collections import Counter

VMAX = 16

def distN(x, N):
    r = x % N
    return r if r <= N - r else N - r

def grid_max(N, eff):
    best, args = -1, []
    for k in range(1, N):
        m = min(distN(w * k, N) for w in eff)
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

def h_of(N):
    return max(0, (N - 1) // 5) if N % 5 else N // 5 - 1

# ---------------- classification verifications ----------------
def bad_set(w, N):
    return {k for k in range(N) if 5 * distN(w * k, N) < N}

def verify_classifications():
    print("== verify proved covering classifications ==")
    for N in (7, 11, 13):
        full = set(range(N))
        inv = [w for w in range(1, N) if gcd(w, N) == 1]
        # classes: x ~ -x ; class rep = min(x, N-x); class group via multiplication
        def cls(x):
            r = x % N
            return min(r, N - r)
        covers = []
        for b in inv:
            for c in inv:
                if c < b:
                    continue
                if bad_set(1, N) | bad_set(b, N) | bad_set(c, N) == full:
                    covers.append((b, c))
        # independent class-based prediction
        m = (N - 1) // 2
        # class group: identify 2's class order and build class multiplication
        # classes labeled 1..m (rep r); class product: rep of (r*s mod N)
        def cprod(x, y):
            return cls(x * y)
        # enumerate class triples directly
        pred = []
        classes = list(range(1, m + 1))
        for B in classes:
            for C in classes:
                if C < B:
                    continue
                # dilates a=1, b=?, c=? have classes X_w = class(w)^-1...
                # prediction via domino/tiling logic implemented brute-force:
                # bad classes of dilate w = {class(w^-1 * j) : 1<=j<=h}
                pass
        print(f"  N={N}: scan covering pairs = {covers}")
        if N == 7:
            # prediction: all three classes distinct
            pred = []
            for b in inv:
                for c in inv:
                    if c < b:
                        continue
                    ca, cb, cc = cls(1), cls(b), cls(c)
                    if len({ca, cb, cc}) == 3:
                        pred.append((b, c))
            print(f"       predicted (3 distinct classes) = {pred}  "
                  f"match={pred == covers}")
        if N == 11:
            # C_5 domino logic: classes of dilate w: {X, X.g} in C_5, g=class(2).
            # covering iff union of three dominos covers all 5 classes.
            g = cls(2)
            def cdiv(x, y):  # x / y in class group of order 5
                # brute force: find z with cprod(z, y) == x
                for z in range(1, m + 1):
                    if cprod(z, y) == x:
                        return z
            def domino(w):
                X = cdiv(cls(1), cls(w))  # class(w)^-1
                return {X, cprod(X, g)}
            pred = []
            for b in inv:
                for c in inv:
                    if c < b:
                        continue
                    u = domino(1) | domino(b) | domino(c)
                    if len(u) == 5:
                        pred.append((b, c))
            print(f"       predicted (C_5 domino cover)  = {pred}  "
                  f"match={pred == covers}")
        if N == 13:
            # C_6 tiling: three dominos {X, X.g} tile C_6 iff starts form a
            # coset of the order-3 subgroup <g^2>: {X, Xg^2, Xg^4}.
            g = cls(2)
            def cdiv(x, y):
                for z in range(1, m + 1):
                    if cprod(z, y) == x:
                        return z
            def gpow(x, e):
                r = x
                for _ in range(e):
                    r = cprod(r, g)
                return r
            pred = []
            for b in inv:
                for c in inv:
                    if c < b:
                        continue
                    Xa = cdiv(cls(1), cls(1))
                    Xb = cdiv(cls(1), cls(b))
                    Xc = cdiv(cls(1), cls(c))
                    tri = {Xa, Xb, Xc}
                    ok = False
                    for X in range(1, m + 1):
                        coset = {X, gpow(X, 2), gpow(X, 4)}
                        if tri == coset:
                            ok = True
                    if ok:
                        pred.append((b, c))
            print(f"       predicted (C_6 parity coset)  = {pred}  "
                  f"match={pred == covers}")

def verify_B5():
    print("\n== verify B5 subgroup-counting arithmetic ==")
    ok = True
    for N in range(5, 61):
        h = h_of(N)
        for j in (2, 3, 4, 5):
            if N % j:
                continue
            # order-j residue exists iff N/j >= 1; sizes: N/j + 2*(2h+1) - 2
            tot = N // j + 2 * (2 * h + 1) - 2
            if j == 5 and tot >= N:
                print(f"  B5 FAILS N={N} j=5 tot={tot}")
                ok = False
            if j == 4 and tot < N:
                pass
    print("  order-5 counting always kills covering:", ok)
    kills4 = [N for N in range(8, 61) if N % 4 == 0 and
              N // 4 + 2 * (2 * h_of(N) + 1) - 2 < N]
    print(f"  order-4 counting kills for N in {kills4}")

def two_runner_law(N):
    worst = None
    viol = []
    for a in range(1, N):
        for b in range(a + 1, N):
            best = max(min(distN(a * k, N), distN(b * k, N))
                       for k in range(1, N))
            if 5 * best < N:
                viol.append((a, b, best))
            if worst is None or best * worst[1] < worst[0] * N:
                worst = (best, N)
    return viol, worst

# ---------------- extended battery ----------------
def battery_flags(N, eff):
    flags = []
    # B0 exact union bound (PROVED): |B_a|+|B_b|+|B_c| - 2 < N  => no covering.
    # (Each B_w contains 0, so the union loses >= 2 to the triple point.)
    sizes = [len(bad_set(w % N, N)) for w in eff]
    if sum(sizes) - 2 < N:
        flags.append(f'B0:count({"+".join(map(str, sizes))})')
    for j in (2, 3, 4, 5):
        if N % j == 0 and all(w % j != 0 for w in eff):
            flags.append(f'B1:jgon{j}')
    for (x, y) in combinations(eff, 2):
        if x % N == y % N or (x + y) % N == 0:
            flags.append('B2:coincide')
            break
    mods = [w % N for w in eff]
    all_inv = all(m != 0 and gcd(m, N) == 1 for m in mods)
    if all_inv:
        if N not in (7, 11, 13):
            flags.append('B3:nocover-modulus')
        else:
            # B3' proved classification: covering iff class condition holds
            def cls(x):
                r = x % N
                return min(r, N - r)
            m = (N - 1) // 2
            def cprod(x, y):
                return cls(x * y)
            def cdiv(x, y):
                for z in range(1, m + 1):
                    if cprod(z, y) == x:
                        return z
            g = cls(2)
            if N == 7:
                tri = {cls(w) for w in mods}
                if len(tri) < 3:
                    flags.append("B3':cls7")
            elif N == 11:
                def domino(w):
                    X = cdiv(cls(1), cls(w))
                    return {X, cprod(X, g)}
                u = set().union(*[domino(w) for w in mods])
                if len(u) < 5:
                    flags.append("B3':cls11")
            elif N == 13:
                Xs = [cdiv(cls(1), cls(w)) for w in mods]
                def gpow(x, e):
                    r = x
                    for _ in range(e):
                        r = cprod(r, g)
                    return r
                tri = set(Xs)
                tiled = any(tri == {X, gpow(X, 2), gpow(X, 4)}
                            for X in range(1, m + 1))
                if not tiled:
                    flags.append("B3':cls13")
    # B5 subgroup counting (SOUND version): an eff residue of additive order
    # j in {4,5} (j | N, j < N i.e. genuinely non-invertible) has |B_w| = N/j
    # exactly (arc avoids nonzero subgroup elements when j <= 5); if the OTHER
    # two residues are invertible their |B| = 2h+1 exactly; then
    #   N/j + 2(2h+1) - 2 < N   =>  no covering.
    # (The earlier version bounded non-invertible others by 2h+1 -- unsound
    #  when another residue is stuck (|B|=N) or order 2 (|B|=N/2); caught by
    #  the soundness cross-check, 245 failures, all B5-only.)
    mods2 = [w % N for w in eff]
    orders = [N // gcd(m, N) if m else 1 for m in mods2]

    def others_invertible(j):
        return all(m != 0 and gcd(m, N) == 1
                   for m, o in zip(mods2, orders) if o != j)

    if 5 in orders and N % 5 == 0 and 5 < N and others_invertible(5):
        flags.append('B5:order5')
    if (4 in orders and N % 4 == 0 and 4 < N and others_invertible(4)
            and N // 4 + 2 * (2 * h_of(N) + 1) - 2 < N):
        flags.append('B5:order4')
    return flags

def main():
    verify_classifications()
    verify_B5()
    print("\n== 2-runner grid law, N in [3,80] ==")
    allv = []
    wg = (10**9, 1)
    for N in range(3, 81):
        v, worst = two_runner_law(N)
        allv += [(N,) + t for t in v]
        if worst and worst[0] * wg[1] < wg[0] * worst[1]:
            wg = worst
    print(f"  violations: {len(allv)}" + (f" {allv[:5]}" if allv else ""))
    print(f"  worst ratio: {wg[0]}/{wg[1]} = {wg[0]/wg[1]:.4f}")

    print("\n== extended battery on gridonly tier ==")
    tier = Counter()
    gridopen = []
    pair_closed = pair_total = 0
    proved_only = Counter()
    for V in combinations(range(1, VMAX + 1), 4):
        if gcd(*V) != 1:
            continue
        per_pair = []
        for (p, q) in combinations(range(4), 2):
            N = V[p] + V[q]
            others = [V[i] for i in range(4) if i not in (p, q)]
            eff = (V[p], others[0], others[1])
            m3, args = grid_max(N, eff)
            cert = 5 * m3 >= N
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
                                 cert=cert, cls=cls))
        if any(t['cls'] == 'classical' for t in per_pair):
            tier['classical'] += 1
        elif any(t['cls'] == 'semi' for t in per_pair if t['cert']):
            tier['semi'] += 1
        elif any(t['cert'] for t in per_pair):
            tier['gridonly'] += 1
            # battery on certifying pairs
            closed = False
            for t in per_pair:
                if not t['cert']:
                    continue
                pair_total += 1
                fl = battery_flags(t['N'], t['eff'])
                # proved-only: B0, B1, B3', B5 (B2 and B3 are verified/pending)
                proved = [f for f in fl if f.startswith(('B0', 'B1', "B3'", 'B5'))]
                if proved:
                    proved_only[proved[0]] += 1
                if fl:
                    pair_closed += 1
                    closed = True
            if not closed:
                gridopen.append((V, [(t['pq'], t['N'], t['eff'], t['m3'])
                                     for t in per_pair if t['cert']]))
        else:
            tier['fail'] += 1
    print(f"tier split: {dict(tier)}")
    print(f"gridonly certifying pairs closed by extended battery: "
          f"{pair_closed}/{pair_total}")
    print(f"proved-only closures (B1/B3'/B5): {sum(proved_only.values())}"
          f"  {dict(proved_only)}")
    print(f"gridonly sets still OPEN: {len(gridopen)}")
    for (V, certs) in gridopen:
        print(f"  OPEN V={V}  certifying pairs:")
        for c in certs:
            print(f"     pair={c[0]} N={c[1]} eff={c[2]} m3={c[3]}")

if __name__ == "__main__":
    main()
    sum_relation_scan()

# ---- sum-relation micro-experiment (appended) ----
def sum_relation_scan():
    print("\n== sum-relation triples (a,b,a+b): covering moduli ==")
    cov_N = {}
    for a in range(1, 13):
        for b in range(a + 1, 13):
            c = a + b
            for N in range(5, 41):
                eff = (a % N, b % N, c % N)
                if any(w == 0 for w in eff):
                    continue
                full = set(range(N))
                if bad_set(a % N, N) | bad_set(b % N, N) | bad_set(c % N, N) == full:
                    cov_N.setdefault((a, b), []).append(N)
    for key, Ns in sorted(cov_N.items()):
        print(f"  eff=({key[0]},{key[1]},{key[0]+key[1]}): covering at N = {Ns}")
    allN = sorted({n for v in cov_N.values() for n in v})
    print(f"  union of covering moduli for sum-triples: {allN}")
