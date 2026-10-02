#!/usr/bin/env python3
"""
lrc_theory_experiments.py -- experiments supporting the integrality (TDI /
circulation) and ladder analyses. Uses the committed census (tight lists,
dumps) as test data.

E1  Exact vertex enumeration of P(v) via the difference-constraint network:
    vertices of the z-polytope = potentials of oriented spanning trees of the
    constraint graph (complete digraph on runners + source node 0). Counts
    fractional vertices -> refutes TDI / unimodular-TU-transformability of the
    full system (its data is integral).

E2  Tight-witness structure: for every census witness k, the set of tight
    pair constraints forms a graph -- count edges, acyclicity, and how many
    witnesses are vertices of P(v) (all-box-boundary or pair-supported).

E3  Modular-witness census: rung-1 (v_i != 0 mod n+1) and rung-2 (exists c
    with c*v mod 2n+1 in [2, 2n-1]) statistics against the gap classes,
    for n=4..7 dumps.

E4  Smith normal form of the pair-constraint matrix A = (n+1)W: verify the
    closed form  invariant factors = (n+1)*gcd(v)  (n-1 of them).

Output: out/theory_experiments.json
"""
import itertools
import json
import math
import random
from fractions import Fraction

import numpy as np

OUT = "/home/z/my-project/scripts/out"
SEED = 20261002

DT = np.dtype({
    "names": ["v", "k", "ta", "tb", "num", "den", "feasible", "pad"],
    "formats": ["(8,)u1", "(8,)u1", "<u2", "<u2", "<u2", "<u2", "u1", "(7,)u1"],
    "offsets": [0, 8, 16, 18, 20, 22, 24, 25],
    "itemsize": 32,
})


# ---------------------------------------------------------------- constraints
def arc_lengths(n, v):
    """Constraint graph: nodes 0..n (0 = source, z_0 = 0).
    Returns dict arc (u, w) -> length, meaning z_w <= z_u + len."""
    a = [Fraction(1, (n + 1) * vi) for vi in v]
    arcs = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            arcs[(j + 1, i + 1)] = (n - 1) * a[j]     # z_i <= z_j + (n-1)a_j
    for i in range(n):
        arcs[(0, i + 1)] = Fraction(1) - n * a[i]     # z_i <= 1 - n a_i
        arcs[(i + 1, 0)] = -a[i]                       # z_i >= a_i
    return arcs, a


def z_to_k(n, v, z, a):
    return [v[i] * (z[i + 1] - a[i]) for i in range(n)]


def feasible_z(n, v, z, arcs):
    for (u, w), length in arcs.items():
        if z[w] - z[u] > length:
            return False
    return True


# ---------------------------------------------------------------- E1: vertices
def pruefer_trees(m):
    """All labeled trees on nodes 0..m-1 via Prüfer sequences."""
    for seq in itertools.product(range(m), repeat=m - 2):
        yield seq


def tree_from_pruefer(seq, m):
    deg = [1] * m
    for x in seq:
        deg[x] += 1
    edges = []
    it = iter(range(m))
    leaves = [i for i in range(m) if deg[i] == 1]
    import heapq
    heapq.heapify(leaves)
    seq = list(seq)
    for x in seq:
        leaf = heapq.heappop(leaves)
        edges.append((leaf, x))
        deg[leaf] -= 1
        deg[x] -= 1
        if deg[x] == 1:
            heapq.heappush(leaves, x)
    a = heapq.heappop(leaves)
    b = heapq.heappop(leaves)
    edges.append((a, b))
    return edges


def enumerate_vertices(n, v):
    """All vertices of P(v), exactly, via oriented spanning trees."""
    arcs, a = arc_lengths(n, v)
    m = n + 1
    seen = {}
    for seq in pruefer_trees(m):
        edges = tree_from_pruefer(seq, m)
        # orient each edge both ways
        for orient in range(1 << n):
            # build adjacency: child -> (parent, arc direction)
            # edge t (0-index), bit t: 0 => arc parent->child, 1 => child->parent
            # we need z from root 0 via BFS using tight arcs
            adj = [[] for _ in range(m)]
            for t, (x, y) in enumerate(edges):
                if orient >> t & 1:
                    # arc y->x tight: z_x = z_y + len(y->x)
                    adj[x].append(y)   # x learns from y
                    adj[y].append(x)   # y learns from x
                    # store direction info separately below
            # simpler: recompute with explicit arc choice
            # redo: for each edge pick the directed arc (u,w): z_w = z_u + L
            chosen = []
            for t, (x, y) in enumerate(edges):
                if orient >> t & 1:
                    chosen.append((y, x))   # arc y -> x
                else:
                    chosen.append((x, y))   # arc x -> y
            # valid orientation: every chosen arc must exist in arcs dict
            if any(ch not in arcs for ch in chosen):
                continue
            # propagate z from node 0
            z = {0: Fraction(0)}
            nbr = {i: [] for i in range(m)}
            for (u, w) in chosen:
                nbr[u].append(w)
                nbr[w].append(u)
            stack = [0]
            ok = True
            while stack:
                u = stack.pop()
                for w in nbr[u]:
                    if w not in z:
                        # find arc between them in the tight orientation
                        # if (u,w) is a chosen arc: z_w = z_u + L(u->w)
                        # if (w,u) is a chosen arc: z_w = z_u - L(w->u)
                        if (u, w) in [c for c in chosen]:
                            z[w] = z[u] + arcs[(u, w)]
                        else:
                            z[w] = z[u] - arcs[(w, u)]
                        stack.append(w)
            if len(z) != m:
                continue
            if not feasible_z(n, v, z, arcs):
                continue
            k = z_to_k(n, v, z, a)
            key = tuple(k)
            if key not in seen:
                seen[key] = k
    return list(seen.values())


def e1_vertices(vectors_by_n):
    result = {}
    for n, vecs in vectors_by_n.items():
        n_frac_vectors = 0
        n_int_vertices = 0
        total_vertices = 0
        first_fractional = None
        for v in vecs:
            verts = enumerate_vertices(n, list(v))
            total_vertices += len(verts)
            has_frac = False
            for k in verts:
                if all(float(x).is_integer() for x in k):
                    n_int_vertices += 1
                else:
                    has_frac = True
            if has_frac:
                n_frac_vectors += 1
                if first_fractional is None:
                    for k in verts:
                        if not all(float(x).is_integer() for x in k):
                            first_fractional = {
                                "v": list(v), "vertex_k": [str(x) for x in k]}
                            break
        result[f"n{n}"] = {
            "vectors": len(vecs),
            "total_vertices": total_vertices,
            "integral_vertices": n_int_vertices,
            "vectors_with_fractional_vertex": n_frac_vectors,
            "first_fractional_vertex": first_fractional,
        }
    return result


# ---------------------------------------------------------------- E2: witnesses
def parse_tight_file(path, n):
    out = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        # v=(..) k=(..) t=a/b gap=c/d
        vpart = line.split("v=(")[1].split(")")[0]
        kpart = line.split("k=(")[1].split(")")[0]
        v = [int(x) for x in vpart.split(",")]
        k = [int(x) for x in kpart.split(",")]
        out.append((v, k))
    return out


def tight_pairs(n, v, k):
    """Pair constraints tight at k, plus tight box rows."""
    pairs = []
    box = []
    for i in range(n):
        if k[i] == 0:
            box.append((i, "lo"))
        elif k[i] == v[i] - 1:
            box.append((i, "hi"))
        for j in range(n):
            if i == j:
                continue
            lhs = (n + 1) * (v[j] * k[i] - v[i] * k[j])
            rhs = n * v[i] - v[j]
            if lhs == rhs:
                pairs.append((i, j))
    return pairs, box


def graph_stats(n, pairs):
    # edges as undirected
    ed = set()
    for (i, j) in pairs:
        ed.add((min(i, j), max(i, j)))
    # acyclicity via union-find
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    acyclic = True
    for (i, j) in ed:
        ri, rj = find(i), find(j)
        if ri == rj:
            acyclic = False
        else:
            parent[ri] = rj
    comps = len({find(i) for i in range(n)})
    return {"edges": len(ed), "directed_tight_pairs": len(pairs),
            "acyclic": acyclic, "components": comps}


def e2_witnesses():
    result = {}
    for n in range(2, 8):
        path = f"{OUT}/n{n}_V50_tight.txt"
        import os
        if not os.path.exists(path):
            path = f"{OUT}/n{n}_V38_tight.txt"
        if not os.path.exists(path):
            continue
        wit = parse_tight_file(path, n)
        stats = []
        for v, k in wit:
            pairs, box = tight_pairs(n, v, k)
            gs = graph_stats(n, pairs)
            gs["v"] = v
            gs["box_tight"] = len(box)
            gs["all_box_boundary"] = (len(box) == n)
            stats.append(gs)
        result[f"n{n}"] = {
            "witnesses": len(wit),
            "all_box_boundary": sum(1 for s in stats if s["all_box_boundary"]),
            "tight_pair_edges_distribution": sorted(
                [s["edges"] for s in stats]),
            "all_forests": all(s["acyclic"] for s in stats),
            "max_edges": max((s["edges"] for s in stats), default=0),
        }
    return result


# ---------------------------------------------------------------- E3: modular
def iter_dump(path):
    with open(path, "rb") as f:
        while True:
            recs = np.fromfile(f, dtype=DT, count=1_000_000)
            if recs.size == 0:
                break
            yield recs


def e3_modular(specs):
    result = {}
    for n, prefix, V in specs:
        delta = Fraction(1, n + 1)
        rung2 = Fraction(2, 2 * n + 1)
        q = 2 * n + 1
        total = 0
        tight = 0
        below_rung2 = 0        # gap in (delta, rung2): ladder-counterexamples
        wit_rung2 = 0
        wit_rung2_and_below = 0
        tight_with_wit = 0
        divfree = 0            # no v_i = 0 mod n+1  (rung-1 modular witness)
        tight_and_divfree = 0
        for recs in iter_dump(f"{OUT}/{prefix}_dump.bin"):
            v = recs["v"][:, :n].astype(np.int64)
            num = recs["num"].astype(np.int64)
            den = recs["den"].astype(np.int64)
            g = np.gcd(num, den)
            # gap as float class: tight iff num/g == 1/(n+1) exactly
            is_tight = (num * (n + 1) == den)
            # gap < rung2 (strictly) and not tight:
            # num/den < 2n+1... compare num*(2n+1) < 2*den
            below = (num * (2 * n + 1) < 2 * den) & (~is_tight) & \
                    (num * (n + 1) > den)
            # rung-2 modular witness: exists c
            hasw = np.zeros(recs.size, dtype=bool)
            for c in range(1, q):
                r = (c * v) % q
                ok = np.all((r >= 2) & (r <= q - 2), axis=1)
                hasw |= ok
            divfree_mask = np.all((v % (n + 1)) != 0, axis=1)
            total += recs.size
            tight += int(is_tight.sum())
            below_rung2 += int(below.sum())
            wit_rung2 += int(hasw.sum())
            wit_rung2_and_below += int((hasw & below).sum())
            tight_with_wit += int((hasw & is_tight).sum())
            divfree += int(divfree_mask.sum())
            tight_and_divfree += int((divfree_mask & is_tight).sum())
        result[f"n{n}"] = {
            "V": V, "total": total, "tight": tight,
            "gap_strictly_between_rung1_rung2": below_rung2,
            "rung2_modular_witness_fraction":
                round(wit_rung2 / total, 4),
            "ladder_counterexamples_with_witness": wit_rung2_and_below,
            "tight_with_rung2_witness": tight_with_wit,
            "divfree_fraction": round(divfree / total, 4),
            "tight_divfree": tight_and_divfree,
        }
    return result


# ---------------------------------------------------------------- E4: SNF
def e4_snf(samples):
    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_form
    result = []
    for n, v in samples:
        rows = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                row = [0] * n
                row[i] = (n + 1) * v[j]
                row[j] = -(n + 1) * v[i]
                rows.append(row)
        A = Matrix(rows)
        S = smith_normal_form(A, domain=ZZ)
        inv = [S[i, i] for i in range(min(S.rows, S.cols))]
        g = math.gcd(*[x for x in v]) if n > 1 else v[0]
        predicted = [(n + 1) * g] * (n - 1) + ([0] * (len(inv) - (n - 1)))
        result.append({
            "n": n, "v": list(v),
            "invariant_factors": [int(x) for x in inv if int(x) != 0],
            "predicted": [(n + 1) * g] * (n - 1),
            "match": [int(x) for x in inv if int(x) != 0] ==
                     [(n + 1) * g] * (n - 1),
        })
    return result


# ============================================================================
def main():
    random.seed(SEED)
    res = {}

    # E1 -- exact vertex enumeration
    print("E1: vertex enumeration ...", flush=True)
    vecs2 = [tuple(c) for c in itertools.combinations(range(1, 13), 2)]
    vecs3 = [tuple(c) for c in itertools.combinations(range(1, 10), 3)]
    tight4 = parse_tight_file(f"{OUT}/n4_V50_tight.txt", 4)
    tight5 = parse_tight_file(f"{OUT}/n5_V50_tight.txt", 5)
    vecs4 = [tuple(v) for v, k in tight4] + [(1, 3, 4, 5), (1, 2, 3, 8)]
    vecs5 = [(1, 3, 4, 5, 9), (1, 2, 3, 4, 5), (1, 3, 4, 5, 7)]
    res["E1_vertices"] = e1_vertices({2: vecs2, 3: vecs3, 4: vecs4, 5: vecs5})

    # E2 -- census witness structure
    print("E2: census witnesses ...", flush=True)
    res["E2_witnesses"] = e2_witnesses()

    # E3 -- modular census
    print("E3: modular witness census ...", flush=True)
    res["E3_modular"] = e3_modular([
        (4, "n4_V50", 50), (5, "n5_V50", 50),
        (6, "n6_V50", 50), (7, "n7_V38", 38)])

    # E4 -- Smith normal form
    print("E4: Smith normal form ...", flush=True)
    samples = [
        (2, (1, 2)), (2, (2, 3)), (2, (2, 4)), (2, (6, 10)),
        (3, (1, 2, 3)), (3, (1, 3, 4)), (3, (2, 4, 6)), (3, (6, 10, 15)),
        (4, (1, 2, 3, 4)), (4, (1, 3, 4, 7)), (4, (2, 4, 6, 8)),
        (5, (1, 2, 3, 4, 5)), (5, (1, 3, 4, 5, 9)), (5, (10, 20, 30, 40, 50)),
    ]
    res["E4_snf"] = e4_snf(samples)

    with open(f"{OUT}/theory_experiments.json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1)[:6000])


if __name__ == "__main__":
    main()
