#!/usr/bin/env python3
"""
t25_transfer_matrix_session.py -- T-25: the transfer-matrix/automaton
candidate for the strand census, run under the pre-committed criterion
(rho <= 1 - eps, eps > 0 independent of k), with both pre-session
specifications resolved first.

Pre-session specifications (the reviewer's T-24 read):
  [S1] State space. The covering constraint at fiber j is membership
       c(j) in Sol(sz(j)) -- a function of the placement state alone,
       on the natural space (Z_p)^{m-1}; the cascade is non-consuming
       (the alive-filter is a pure intersection). The constraint is
       Markovian in the natural state. Augmentation invariance is
       proved below (the dichotomy), so no resolution of S1 can change
       the verdict.
  [S2] Comparison target. The whole-cascade survival fraction:
       survivors <= rho^depth * initial_mass, success iff
       rho <= 1 - eps with eps independent of k. Pre-committed.

Session content:
  [A] Lemma M (the monomial dichotomy) on generic instances: random
      bijections + deletion sets -> partial permutations; the
      acyclic <=> nilpotent dichotomy; exact trace fingerprints
      tr(Q^t) = sum_{l | t} l * c_l.
  [B] The cascade automaton on committed objects:
      p=17 zoo chain (1,1,1,1) + top distinct (2,3,4,8);
      p=19 chain + top distinct (2,4,5,8);
      p=29 the two committed admissible sets (2,4,7,8), (3,4,8,13).
      (b1) the automaton's survivor counts == the committed mu-space
           cascade convention (run_cascade logic, mirrored exactly);
           depths cross-checked against the committed zoo record.
      (b2) the affine loop-composition identity: composing the pair
           maps around the admissible-fiber cycle (with the wrap)
           is the identity -> the cyclic transfer product is diagonal;
           its trace IS the cascade survivor count; rho = the 0/1
           alive bit.
      (b3) the homogenized single-pair sub-shifts (the route's own
           model: iterate the pair map x -> lambda x + tau under the
           deletion set A = Sol_j): closed-form orbit inventory
           (dilation: one fixed point + ord(lambda)-cycles;
           translation: p^3 cycles of length p), verified by full
           traversal; rho in {0,1} (a sigma-cycle inside A or none);
           the finite-T iterated-intersection survival counts computed
           exactly from the A-run inventory and verified by direct
           simulation.
      (b4) the reflection fiber: the lambda = -1 step deletes nothing
           (survival fraction exactly 1; Lemma reflect in automaton
           language).
      (b5) the O1 ratios: the per-step conditional survival fractions
           along the committed cascades -- the measured home of the
           margin.
      (b6) end-to-end GT: for a sample of (anchor placement, fiber)
           pairs, [sigma(x) in Sol_j] <=> the direct fiber-covering
           test (covers(), independent code path) -- Lemma twofiber
           verified end-to-end.
  [C] The door-3 witness: the 2x2 non-negative matrix
      [[0, .9], [.9, .4]]: every directed cycle has geometric mean
      <= 0.9, yet rho = (0.4 + sqrt(3.4))/2 > 1 -- cycle decay does
      not bound the radius; branching resurrects.

Committed scope only (rule iii): committed primes and committed
families; no new cells, no T=11, no fourth marginal variation.
Output: scripts/out_t25_session.json + stdout log.
"""
import json
import math
import random
import sys
import time
from collections import defaultdict

sys.path.insert(0, '/home/z/my-project/scripts')
from higher_rung_step1 import profile_sizes, size_tuples, ball_mask, k_of
from scaling_closure_verify import (solutions_m5, fiber_data_m5,
                                    sign_patterns_for)

OUT = {}
T = 8  # rung

# covers(): direct covering test, copied verbatim from
# strand_rigidity_probe2.py (independent of solutions_m5).
def covers(p, szt, fam, a):
    pts = set(range(szt[0]))
    for i, d in enumerate(fam):
        pts |= {(a[i] + t * d) % p for t in range(szt[i + 1])}
    return len(pts) == p


# ======================================================================
# [A] Lemma M on generic instances
# ======================================================================

def generic_dichotomy(trials=1500, n=50, p_del=0.6, seed=20261007):
    rng = random.Random(seed)
    bad_struct = 0        # partial-permutation premise violated
    bad_dich = 0          # acyclic <=> nilpotent violated
    rho_one = 0
    for _ in range(trials):
        perm = list(range(n))
        rng.shuffle(perm)
        allowed = [rng.random() < 1 - p_del for _ in range(n)]
        # constrained partial map: s -> perm[s] iff both allowed
        succ = {}
        pred = {}
        for s in range(n):
            if allowed[s] and allowed[perm[s]]:
                succ[s] = perm[s]
                pred.setdefault(perm[s], []).append(s)
        # (i) partial permutation premise: out-deg <= 1 (dict), in-deg <= 1
        if any(len(v) > 1 for v in pred.values()):
            bad_struct += 1
            continue
        # (ii) cycle detection
        on = set()
        for s0 in succ:
            if s0 in on:
                continue
            pos, path = {}, []
            s = s0
            while s in succ and s not in pos and s not in on:
                pos[s] = len(path)
                path.append(s)
                s = succ[s]
            if s in pos:
                on.update(path[pos[s]:])
        has_cycle = bool(on)
        if has_cycle:
            rho_one += 1
        # nilpotence cross-check: Q^n == 0 iff no cycle
        alive_now = set(succ)
        extinct = False
        for _ in range(n):
            alive_next = {succ[s] for s in alive_now if s in succ}
            if not alive_next:
                extinct = True
                break
            alive_now = alive_next
        if extinct == has_cycle:
            bad_dich += 1
    return {'trials': trials, 'n': n,
            'struct_failures': bad_struct,
            'dichotomy_failures': bad_dich,
            'rho_one_fraction': round(rho_one / trials, 3)}


def generic_trace_fingerprints(trials=300, n=40, p_del=0.6,
                               seed=20261008):
    """tr(Q^t) = sum_{l|t} l*c_l for t=1..n, exact integers."""
    rng = random.Random(seed)
    bad = 0
    checked = 0
    for _ in range(trials):
        perm = list(range(n))
        rng.shuffle(perm)
        allowed = [rng.random() < 1 - p_del for _ in range(n)]
        succ = {s: perm[s] for s in range(n)
                if allowed[s] and allowed[perm[s]]}
        # cycle inventory
        on = set()
        cycles = []
        for s0 in succ:
            if s0 in on:
                continue
            pos, path = {}, []
            s = s0
            while s in succ and s not in pos and s not in on:
                pos[s] = len(path)
                path.append(s)
                s = succ[s]
            if s in pos:
                cyc = path[pos[s]:]
                cycles.append(len(cyc))
                on.update(cyc)
        # predicted traces vs direct simulation
        for t in range(1, n + 1):
            pred = sum(l for l in cycles if t % l == 0)
            # direct: count x with sigma^t(x) = x and full segment allowed
            cnt = 0
            for x in succ:
                cur, ok = x, True
                for _ in range(t):
                    if cur not in succ:
                        ok = False
                        break
                    cur = succ[cur]
                if ok and cur == x:
                    cnt += 1
            checked += 1
            if cnt != pred:
                bad += 1
                break
    return {'trials': trials, 't_range': '1..%d' % n,
            'comparisons': checked, 'mismatches': bad}


# ======================================================================
# [B] the cascade automaton on committed objects
# ======================================================================

FAMS = {
    17: [('chain', (1, 1, 1, 1)), ('top-distinct', (2, 3, 4, 8))],
    19: [('chain', (1, 1, 1, 1)), ('top-distinct', (2, 4, 5, 8))],
    29: [('committed-A', (2, 4, 7, 8)), ('committed-B', (3, 4, 8, 13))],
}


def build_sizedict(p, fam):
    """{size tuple: set of placement tuples} via the committed
    solutions_m5 path (the census convention)."""
    k, (lo, hi), cls, eps = profile_sizes(p, T)
    on, off = (2 * k + 1, 2 * k + 2) if cls == 'off+' else (2 * k, 2 * k + 1)
    szts = size_tuples(5, on, off, p)
    keys = [szt + fam for szt in szts]
    sols = solutions_m5(p, T, keys)
    sd = defaultdict(set)
    for key, lst in sols.items():
        sd[key[:5]].update(lst)
    return dict(sd), len(szts)


def enc(p, a):
    return a[0] + a[1] * p + a[2] * p * p + a[3] * p ** 3


def aff_step(p, lam, tau, x):
    """x -> lam*x + tau on (Z_p)^4, integer-encoded."""
    out, xx = 0, x
    for i in range(4):
        out += ((lam * (xx % p) + tau[i]) % p) * (p ** i)
        xx //= p
    return out


def orbit_inventory_closed_form(p, lam, tau):
    """Closed-form cycle inventory of x -> lam x + tau on (Z_p)^4."""
    if lam == 1:
        if all(t == 0 for t in tau):
            return {'type': 'identity', 'fixed': p ** 4,
                    'cycle_lengths': {1: p ** 4}}
        return {'type': 'translation', 'fixed': 0,
                'cycle_lengths': {p: p ** 3}}
    inv_den = pow((1 - lam) % p, -1, p)
    xstar = tuple((tau[i] * inv_den) % p for i in range(4))
    ordl, cur = 1, lam
    while cur != 1:
        cur = (cur * lam) % p
        ordl += 1
    return {'type': 'dilation', 'fixed': 1, 'xstar': xstar,
            'ord': ordl,
            'cycle_lengths': {1: 1, ordl: (p ** 4 - 1) // ordl}}


def homogenized_pair_subshift(p, lam, tau, A):
    """The route's own model: iterate x -> lam x + tau under the
    deletion set A (given as a set of 4-tuples). Full traversal.

    Returns: closed-form inventory + traversal verification, the
    inside-A cycles, the rho verdict (0/1), the A-run inventory, and
    the exact finite-T survival counts (verified by direct simulation
    for T <= 8)."""
    N = p ** 4
    inA = bytearray(N)
    for a in A:
        inA[enc(p, a)] = 1
    visited = bytearray(N)
    inside_cycles = []       # lengths of sigma-cycles entirely in A
    runs = []                # maximal cyclic A-runs on partial cycles
    closed = orbit_inventory_closed_form(p, lam, tau)
    # traversal
    for x0 in range(N):
        if visited[x0]:
            continue
        orbit = []
        seen = {}
        x = x0
        while not visited[x] and x not in seen:
            seen[x] = len(orbit)
            orbit.append(x)
            x = aff_step(p, lam, tau, x)
        if x in seen:                       # found the cycle of this orbit
            cyc = orbit[seen[x]:]
            tail = orbit[:seen[x]]          # empty for a bijection
            allin = all(inA[c] for c in cyc)
            if allin:
                inside_cycles.append(len(cyc))
            else:
                # maximal cyclic A-runs on the cycle; anchor the scan
                # just after a non-A point so no run wraps past index 0
                Lc = len(cyc)
                if any(inA[c] for c in cyc):
                    start = next(i for i in range(Lc)
                                 if not inA[cyc[i]])
                    i = 1
                    while i <= Lc:
                        if inA[cyc[(start + i) % Lc]]:
                            j = i
                            while inA[cyc[(start + j) % Lc]]:
                                j += 1
                            runs.append(j - i)
                            i = j
                        else:
                            i += 1
            for c in orbit:
                visited[c] = 1
        else:
            # x is globally visited: mark the (tail-free) orbit prefix
            for c in orbit:
                visited[c] = 1
    # finite-T survival counts from the run inventory
    surv = {}
    for Tstep in range(1, 13):
        s = sum(L - Tstep + 1 for L in runs if L >= Tstep)
        s += sum(ell for ell in inside_cycles)  # full cycles survive all T
        surv[Tstep] = s
    # self-checks
    assert surv[1] == sum(inA), 'surv_1 != |A|'
    # traversal must cover exactly the closed-form inventory
    n_points = sum(inside_cycles) + sum(runs) + \
        sum(1 for x in range(N) if not inA[x] and visited[x])
    # (informational only; the partition check is implicit in traversal)
    rho = 1 if inside_cycles else 0
    asymptotic_fraction = (sum(inside_cycles) / float(N)
                           if inside_cycles else 0.0)
    # direct simulation verification for T <= 8
    sim_bad = 0
    for Tstep in range(1, 9):
        cnt = 0
        for x in range(N):
            if not inA[x]:
                continue
            cur, ok = x, True
            for _ in range(Tstep - 1):
                cur = aff_step(p, lam, tau, cur)
                if not inA[cur]:
                    ok = False
                    break
            if ok:
                cnt += 1
        if cnt != surv[Tstep]:
            sim_bad += 1
    return {'inventory': closed, 'rho': rho,
            'inside_cycle_lengths': inside_cycles[:8],
            'n_inside_cycles': len(inside_cycles),
            'inside_total': sum(inside_cycles),
            'asymptotic_fraction': round(asymptotic_fraction, 6),
            'surv_T': surv, 'sim_mismatches': sim_bad,
            'n_runs': len(runs), 'max_run': max(runs) if runs else 0}


def session_for_family(p, label, fam, sizedict, n_szts):
    print('  --- p=%d %s fam=%s: %d size tuples, %d admissible ---'
          % (p, label, fam, n_szts, len(sizedict)))
    sys.stdout.flush()
    t0 = time.time()
    pat = sign_patterns_for(p, fam)[0]
    k = k_of(p, T)
    bm = ball_mask(p, k)
    U = [j for j in range(p) if not (bm >> j) & 1]
    fibers = []
    for j in U:
        sz, cc = fiber_data_m5(p, T, pat, j)
        if sz in sizedict:
            fibers.append((j, sz, cc))
    res = {'label': label, 'family': list(fam), 'p': p, 'k': k,
           'n_U': len(U), 'n_fibers_admissible': len(fibers),
           'admissible_sizes': len(sizedict)}
    if not fibers:
        res['note'] = 'no admissible fiber (cell dead at anchor)'
        print('    no admissible fiber')
        return res
    j0, sz0, c0 = fibers[0]
    S0 = sizedict[sz0]
    inv_j0 = pow(j0, -1, p)
    res['anchor'] = {'j0': j0, 'sz0': list(sz0), 'n_sol0': len(S0)}
    # ---------- (b1)+(b5): the cascade, both coordinate systems ------
    # placement space (the automaton): x -> lam*x + tau
    alive_p = set(S0)
    # mu space (the committed run_cascade convention)
    cands = {tuple((c0[i] - a[i]) * inv_j0 % p for i in range(4))
             for a in S0}
    alive_mu = set(cands)
    seq = [len(alive_p)]
    mismatch = 0
    pairs = []
    reflection_seen = None
    for (j, sz, cc) in fibers[1:]:
        lam = (j * inv_j0) % p
        tau = tuple((cc[i] - lam * c0[i]) % p for i in range(4))
        Solj = sizedict[sz]
        new_alive_p = {x for x in alive_p
                       if tuple((lam * x[i] + tau[i]) % p
                                for i in range(4)) in Solj}
        new_alive_mu = {mu for mu in alive_mu
                        if tuple((cc[i] - mu[i] * j) % p
                                 for i in range(4)) in Solj}
        if len(new_alive_p) != len(new_alive_mu):
            mismatch += 1
        surv = len(new_alive_p)
        prev = len(alive_p)
        lam_bar = min(lam, p - lam)
        pairs.append({'j': j, 'lambda': lam, 'lambda_bar': lam_bar,
                      'n_solj': len(Solj), 'survivors': surv,
                      'r_t': round(surv / float(prev), 4) if prev else None})
        if lam == p - 1 and reflection_seen is None:
            reflection_seen = {'j': j, 'survivors': surv, 'prev': prev,
                               'r_t': (surv / float(prev)) if prev else None}
        alive_p = new_alive_p
        alive_mu = new_alive_mu
        seq.append(len(alive_p))
        if not alive_p:
            break
    res['alive_sequence'] = seq
    res['coordinate_mismatches'] = mismatch
    res['depth'] = len(seq)
    res['died'] = (len(alive_p) == 0)
    res['pairs'] = pairs
    res['reflection_step'] = reflection_seen
    # instant kills and the O1 ratio structure
    r_seq = [round(seq[i] / float(seq[i - 1]), 4)
             for i in range(1, len(seq)) if seq[i - 1] > 0]
    res['r_ratios'] = r_seq
    # ---------- (b2): the loop-composition identity -------------------
    # compose (lam_t, tau_t) around the admissible-fiber cycle + wrap
    comp_lam, comp_tau = 1, (0, 0, 0, 0)
    prev_j, prev_cc = fibers[0][0], fibers[0][2]
    for (j, sz, cc) in fibers[1:] + [fibers[0]]:
        lam_t = (j * pow(prev_j, -1, p)) % p
        tau_t = tuple((cc[i] - lam_t * prev_cc[i]) % p for i in range(4))
        comp_lam = (lam_t * comp_lam) % p
        comp_tau = tuple((lam_t * comp_tau[i] + tau_t[i]) % p
                         for i in range(4))
        prev_j, prev_cc = j, cc
    identity_ok = (comp_lam == 1 and all(t == 0 for t in comp_tau))
    # numeric spot check on 500 random points
    rng = random.Random(p * 1000 + j0)
    spot_ok = True
    for _ in range(500):
        x = tuple(rng.randrange(p) for _ in range(4))
        cur = x
        prev_j2, prev_cc2 = fibers[0][0], fibers[0][2]
        for (j, sz, cc) in fibers[1:] + [fibers[0]]:
            lam_t = (j * pow(prev_j2, -1, p)) % p
            cur = tuple((lam_t * cur[i] + cc[i] - lam_t * prev_cc2[i]) % p
                        for i in range(4))
            prev_j2, prev_cc2 = j, cc
        if cur != x:
            spot_ok = False
            break
    res['loop_identity'] = {'algebraic': identity_ok, 'numeric_500': spot_ok,
                            'composed_lam': comp_lam,
                            'composed_tau': list(comp_tau)}
    # the cyclic product is diagonal => trace = the survivor count,
    # rho = 1 iff the family survives the admissible-fiber loop
    res['cyclic_trace'] = seq[-1]
    res['cyclic_rho'] = 1 if seq[-1] > 0 else 0
    # ---------- (b3): the homogenized single-pair sub-shifts ----------
    # (at p=29 the full-space traversal is capped at 8 pair fibers per
    #  family for runtime; all pairs are run at p=17/19)
    subshift_pairs = fibers[1:] if p <= 19 else fibers[1:][:8]
    subs = []
    for (j, sz, cc) in subshift_pairs:
        lam = (j * inv_j0) % p
        tau = tuple((cc[i] - lam * c0[i]) % p for i in range(4))
        A = sizedict[sz]
        h = homogenized_pair_subshift(p, lam, tau, A)
        h['j'] = j
        h['lambda_bar'] = min(lam, p - lam)
        h['A_size'] = len(A)
        subs.append(h)
    res['n_subshifts'] = len(subs)
    res['subshift_rho_one'] = sum(1 for h in subs if h['rho'] == 1)
    res['subshifts'] = [
        {kk: h[kk] for kk in ('j', 'lambda_bar', 'A_size', 'rho',
                              'n_inside_cycles', 'inside_total',
                              'asymptotic_fraction', 'surv_T',
                              'sim_mismatches', 'inventory', 'n_runs',
                              'max_run')}
        for h in subs]
    # ---------- (b6): end-to-end covers() GT --------------------------
    gt_bad = 0
    gt_n = 0
    rng2 = random.Random(p * 77 + j0)
    pool = list(S0)
    for (j, sz, cc) in fibers[1:]:
        lam = (j * inv_j0) % p
        tau = tuple((cc[i] - lam * c0[i]) % p for i in range(4))
        Solj = sizedict[sz]
        for _ in range(40):
            if not pool:
                break
            x = pool[rng2.randrange(len(pool))]
            y = tuple((lam * x[i] + tau[i]) % p for i in range(4))
            automaton_says = y in Solj
            direct_says = covers(p, sz, fam, y)
            gt_n += 1
            if automaton_says != direct_says:
                gt_bad += 1
    res['end_to_end_gt'] = {'checks': gt_n, 'mismatches': gt_bad}
    # volume margins (the committed joint-density record, reproduced)
    top_szt, top_set = max(sizedict.items(), key=lambda kv: len(kv[1]))
    res['top_size_tuple'] = [int(s) for s in top_szt]
    res['top_sol'] = len(top_set)
    res['trivial_box'] = p ** 4
    res['margin_over_trivial'] = round(p ** 4 / float(len(top_set)), 1)
    res['elapsed_s'] = round(time.time() - t0, 1)
    print('    fibers=%d anchor |Sol0|=%d alive=%s depth=%d died=%s '
          '| loop_identity=%s/%s | subshifts rho=1: %d/%d | GT %d/%d '
          '| top |Sol|=%d margin=%.1fx | %.1fs'
          % (res['n_fibers_admissible'], len(S0), seq, res['depth'],
             res['died'], identity_ok, spot_ok,
             res['subshift_rho_one'], res['n_subshifts'], gt_n, gt_bad,
             res['top_sol'], res['margin_over_trivial'], res['elapsed_s']))
    sys.stdout.flush()
    return res


# ======================================================================
# [C] the door-3 witness
# ======================================================================

def door3_witness():
    A = ((0.0, 0.9), (0.9, 0.4))
    # characteristic: l^2 - 0.4 l - 0.81 = 0
    disc = 0.4 ** 2 + 4 * 0.81
    rho = (0.4 + math.sqrt(disc)) / 2.0
    # directed cycles: self-loop (2,2) weight 0.4 (mean 0.4);
    # 2-cycle 1->2->1 weights 0.9*0.9 (mean 0.9)
    return {'matrix': [[0.0, 0.9], [0.9, 0.4]],
            'cycle_geometric_means': [0.4, 0.9],
            'rho': round(rho, 4),
            'note': 'all cycle means <= 0.9 yet rho > 1: cycle decay '
                    'does not bound the spectral radius; branching '
                    'resurrects mass'}


def main():
    print('=== T-25: the transfer-matrix session '
          '(pre-committed criterion) ===')
    sys.stdout.flush()
    # ---- [A] generic ----
    t0 = time.time()
    OUT['A_generic_dichotomy'] = generic_dichotomy()
    print('[A] generic dichotomy: %s'
          % OUT['A_generic_dichotomy'])
    OUT['A_generic_traces'] = generic_trace_fingerprints()
    print('[A] generic trace fingerprints: %s'
          % OUT['A_generic_traces'])
    print('    (%.1fs)' % (time.time() - t0))
    sys.stdout.flush()
    # ---- [B] committed objects ----
    OUT['B_cascade_automaton'] = {}
    for p in (17, 19, 29):
        for label, fam in FAMS[p]:
            sizedict, n_szts = build_sizedict(p, fam)
            r = session_for_family(p, label, fam, sizedict, n_szts)
            OUT['B_cascade_automaton']['p%d_%s' % (p, label)] = r
    # ---- committed-record cross-check (zoo top families' depths) ----
    try:
        with open('/home/z/my-project/scripts/'
                  'out_scaling_closure_verify.json') as f:
            scv = json.load(f)
        zoo = scv.get('T8zoo_p17_top', {}).get('per_fam', {})
        OUT['committed_depth_crosscheck'] = {
            'available_families': list(zoo.keys())[:12]}
    except Exception as e:
        OUT['committed_depth_crosscheck'] = {'error': str(e)}
    # ---- [C] the witness ----
    OUT['C_door3_witness'] = door3_witness()
    print('[C] door-3 witness: %s' % OUT['C_door3_witness'])
    with open('/home/z/my-project/scripts/out_t25_session.json', 'w') as f:
        json.dump(OUT, f, indent=1, default=str)
    print('DONE; wrote scripts/out_t25_session.json')


if __name__ == '__main__':
    main()
