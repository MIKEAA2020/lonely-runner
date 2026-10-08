"""
LR-zonotope Ehrhart / covering-radius library (exact rational arithmetic).

Model (established and hand-verified in prior sessions, see
lonely-runner/audits/glm thought4.txt / glm thought5.txt):

  speeds v = (v_1,...,v_n), positive, gcd 1.  V = R^n / Rv,  pi : R^n -> V.
  Lambda = pi(Z^n)  (rank d = n-1).   Z(v) = pi([0,1]^n).
  Ehrhart:  Ehr_Z(t) = sum_{k=0}^{d} e_k t^k,
            e_k = sum_{S subset [n], |S|=k} gcd(v_{S^c}),  gcd(empty)=0.
  h* from  sum_{t>=0} Ehr(t) z^t = h*(z)/(1-z)^{d+1}.
  lambda(v) = sup_t min_j ||t v_j||_T ;  delta = 1/2 - lambda.
  c-coordinates (requires v_1 = 1):  c_j = x_j - v_j x_1  (j=2..n), Lambda = Z^d.
  V-norm:  ||c||_V = max( i<j ) |phi^{ij}| ,  phi^{1j} = c_j/(1+v_j),
           phi^{ij} = (v_j c_i - v_i c_j)/(v_i + v_j)   (2<=i<j).
  FACT (proved here): ||c||_V <= ||c||_inf  (each form has l1-norm <= 1),
    hence B = unit ball contains [-1,1]^d and the grid certification of rho
    needs no norm-equivalence constant:  rho <= gridmax + h/2.
  delta = dist_V(c*, Lambda),  c* = pi(1/2 * 1) = ((1-v_j)/2)_{j=2..n}
    (a 2-torsion point:  c*_j = 1/2 iff v_j even, 0 iff v_j odd).
  rho = max over torus of D(x) = min_{z in Z^d} ||x - z||_V.
"""
from fractions import Fraction
from itertools import combinations
import math

# ---------------------------------------------------------------- basic helpers

def frac_dist_to_Z(q):
    """||q|| for rational q: distance to nearest integer."""
    fl = math.floor(q)
    cand = [q - fl, fl + 1 - q]
    return min(cand)

def norm_T(v, t):
    return min(frac_dist_to_Z(t * vj) for vj in v)

def lambda_exact(v):
    """Exact loneliness lambda(v) via breakpoint/crossover enumeration."""
    n = len(v)
    cands = {Fraction(0), Fraction(1)}
    for j in range(n):
        for m in range(0, 2 * v[j] + 1):
            cands.add(Fraction(m, 2 * v[j]))
    for i in range(n):
        for j in range(i + 1, n):
            s = v[i] + v[j]
            for m in range(0, s + 1):
                cands.add(Fraction(m, s))
            d0 = abs(v[i] - v[j])
            if d0:
                for m in range(0, d0 + 1):
                    cands.add(Fraction(m, d0))
    best = Fraction(-1)
    best_t = None
    for t in sorted(cands):
        val = min(frac_dist_to_Z(t * vj) for vj in v)
        if val > best:
            best, best_t = val, t
    return best, best_t

# ---------------------------------------------------------------- Ehrhart / h*

def ehr_coeffs(v):
    """e_k = sum_{|S|=k} gcd(v_{S^c});  gcd(empty set) := 0."""
    n = len(v)
    d = n - 1
    e = [0] * (d + 1)
    for mask in range(1 << n):
        k = bin(mask).count("1")
        if k > d:
            continue  # S = full set: gcd(empty) = 0, no contribution
        comp = [v[j] for j in range(n) if not (mask >> j) & 1]
        g = math.gcd(*comp) if comp else 0
        e[k] += g
    return e  # length d+1, e[d] = sum v_j (volume)

def _eulerian_direct(m):
    """A_m(z) = sum_{j=0..m} <m choose j> z^j (Eulerian numbers)."""
    if m == 0:
        return [1]
    E = [[0] * (m + 1) for _ in range(m + 1)]
    # recurrence: <n,k> = (k+1)<n-1,k> + (n-k)<n-1,k-1>
    for n_ in range(1, m + 1):
        for k_ in range(0, n_ + 1):
            if n_ == 1:
                E[n_][k_] = 1 if k_ == 0 else 0
            else:
                a = (k_ + 1) * (E[n_ - 1][k_] if k_ <= n_ - 1 else 0)
                b = (n_ - k_) * (E[n_ - 1][k_ - 1] if k_ - 1 >= 0 else 0)
                E[n_][k_] = a + b
    return [E[m][j] for j in range(m)] if m >= 1 else [1]   # A_m has degree m-1

def polymul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r

def polyadd(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i, x in enumerate(a):
        r[i] += x
    for i, y in enumerate(b):
        r[i] += y
    return r

def polyscale(a, s):
    return [x * s for x in a]

def pow_minus_z(m):
    """(1-z)^m as coefficient list (ascending)."""
    r = [1]
    f = [1, -1]
    for _ in range(m):
        r = polymul(r, f)
    return r

def hstar_from_ehr(e, d):
    """h*(z) = e_0 (1-z)^d + sum_{k>=1} e_k z A_k(z) (1-z)^{d-k},
    where sum_{t>=0} t^k z^t = z A_k(z)/(1-z)^{k+1} (A_k: Eulerian, deg k-1)."""
    h = polyscale(pow_minus_z(d), e[0])
    for k in range(1, min(len(e), d + 1)):
        if e[k] == 0:
            continue
        zk = [0] + _eulerian_direct(k)  # z * A_k(z)  (degree k)
        term = polymul(zk, pow_minus_z(d - k))
        h = polyadd(h, polyscale(term, e[k]))
    while len(h) > 1 and h[-1] == 0:
        h.pop()
    return h

def ehr_poly_eval(e, t):
    return sum(e[k] * t ** k for k in range(len(e)))

def hstar_by_inversion(Ehr_vals, d):
    """h*_k = sum_{j=0..k} (-1)^j C(d+1,j) Ehr(k-j)."""
    from math import comb
    hs = []
    for k in range(d + 1):
        acc = 0
        for j in range(k + 1):
            acc += (-1) ** j * comb(d + 1, j) * Ehr_vals[k - j]
        hs.append(acc)
    while len(hs) > 1 and hs[-1] == 0:
        hs.pop()
    return hs

def ehr_from_hstar(hs, d, t):
    """Ehr(t) = sum_k h*_k C(t + d - k, d)."""
    from math import comb
    acc = 0
    for k, h in enumerate(hs):
        if h:
            acc += h * comb(t + d - k, d)
    return acc

# ---------------------------------------------------------------- DFS counting

def dfs_count(v, t):
    """# of cosets of Zv in Z^n whose Rv-line meets [0,t]^n (t integer >=0)."""
    n = len(v)
    v1 = v[0]
    total = 0
    for u1 in range(v1):
        lo = Fraction(-u1, v[0])
        hi = Fraction(t - u1, v[0])
        total += _dfs_rec(v, 1, lo, hi, t)
    return total

def _dfs_rec(v, j, lo, hi, t):
    if j == len(v):
        return 1 if lo <= hi else 0
    # u_j must allow some s in [lo,hi] with u_j + s v_j in [0,t]
    lo_u = -hi * v[j]           # u_j >= -hi * v_j
    hi_u = t - lo * v[j]        # u_j <= t - lo * v_j
    a = math.ceil(lo_u)
    b = math.floor(hi_u)
    if b < a:
        return 0
    cnt = 0
    for uj in range(a, b + 1):
        nlo = max(lo, Fraction(-uj, v[j]))
        nhi = min(hi, Fraction(t - uj, v[j]))
        if nlo <= nhi:
            cnt += _dfs_rec(v, j + 1, nlo, nhi, t)
    return cnt

# ---------------------------------------------------------------- c-model

def make_forms(v):
    """Linear forms of the V-norm in c-coordinates (v[0] must be 1).
    Returns list of (coeffs, denom): form(c) = (sum coeff_k c_k)/denom,
    with coeff list of ints of length d = n-1 (1- or 2-sparse), denom>0."""
    assert v[0] == 1, "c-model requires v_1 = 1"
    n = len(v)
    d = n - 1
    forms = []
    for j in range(2, n + 1):        # phi^{1j} = c_j / (1+v_j)
        coeffs = [0] * d
        coeffs[j - 2] = 1
        forms.append((tuple(coeffs), 1 + v[j - 1]))
    for i in range(2, n + 1):
        for j in range(i + 1, n + 1):  # phi^{ij} = (v_j c_i - v_i c_j)/(v_i+v_j)
            coeffs = [0] * d
            coeffs[i - 2] = v[j - 1]
            coeffs[j - 2] = -v[i - 1]
            forms.append((tuple(coeffs), v[i - 1] + v[j - 1]))
    return forms

def form_value(form, c):
    coeffs, denom = form
    s = sum(a * x for a, x in zip(coeffs, c))
    return Fraction(s, denom)

def norm_V(c, forms):
    return max(abs(form_value(f, c)) for f in forms)

def cstar(v):
    """c* = pi(1/2 * 1) in c-coords, reduced mod 1 to [0,1)^d."""
    d = len(v) - 1
    return tuple(Fraction(1 - v[j], 2) % 1 for j in range(1, len(v)))

def _round_frac(x):
    fl = math.floor(x)
    if x - fl >= Fraction(1, 2):
        return fl + 1
    return fl

def D_exact(x, v, forms, r_init=None):
    """Exact D(x) = min_{z in Z^d} ||x-z||_V via pruned DFS over z.
    x: tuple of Fractions (length d)."""
    d = len(x)
    beta = [1 + v[j + 1] for j in range(d)]     # exact axis extents of B
    z0 = [_round_frac(xi) for xi in x]
    r = norm_V(tuple(xi - zi for xi, zi in zip(x, z0)), forms)
    # search all z with ||x - z||_V <= r (r shrinks when better z found)
    order = sorted(range(d), key=lambda k: -beta[k])
    pos = {k: i for i, k in enumerate(order)}    # DFS position of each coord
    forms_by_late = {k: [] for k in range(d)}   # forms whose LAST-ASSIGNED coord is k
    for fi, (coeffs, denom) in enumerate(forms):
        supp = [k for k in range(d) if coeffs[k] != 0]
        if not supp:
            continue
        late = max(supp, key=lambda k: pos[k])
        forms_by_late[late].append(fi)

    best_r = r
    assigned = [None] * d

    def rec(idx):
        nonlocal best_r
        if idx == d:
            zz = tuple(assigned)
            val = norm_V(tuple(xi - zi for xi, zi in zip(x, zz)), forms)
            if val < best_r:
                best_r = val
            return
        k = order[idx]
        a = math.ceil(x[k] - best_r * beta[k])
        b = math.floor(x[k] + best_r * beta[k])
        for zk in range(a, b + 1):
            assigned[k] = zk
            ok = True
            for fi in forms_by_late[k]:
                coeffs, denom = forms[fi]
                s = 0
                for kk, cc in enumerate(coeffs):
                    if cc and assigned[kk] is not None:
                        s += cc * (x[kk] - assigned[kk])
                val = abs(Fraction(s, denom))
                # unassigned support coords (at most one, and it is not k)
                rest = [kk for kk, cc in enumerate(coeffs) if cc and assigned[kk] is None]
                if not rest:
                    if val > best_r:
                        ok = False
                else:
                    kk = rest[0]
                    pot = abs(coeffs[kk]) * best_r * beta[kk]
                    if val - pot > best_r:
                        ok = False
            if ok:
                rec(idx + 1)
            assigned[k] = None
        return

    rec(0)
    return best_r

def M_corner(v, forms):
    """M = max over cube corners of ||c||_V  (upper bound for D on the torus)."""
    d = len(v) - 1
    best = Fraction(0)
    for corner in __import__("itertools").product([Fraction(0), Fraction(1)], repeat=d):
        best = max(best, norm_V(corner, forms))
    return best

def enumerate_Z0(v, forms, M):
    """z in Z^d with dist_V(z,[0,1]^d) <= M (necessary per-form filter)."""
    d = len(v) - 1
    beta = [1 + v[j + 1] for j in range(d)]
    # crude box, then filter: for each form, min over box of |form(z-x)| <= M
    out = []
    ranges = []
    for k in range(d):
        a = math.ceil(-M * beta[k] - 1)
        b = math.floor(1 + M * beta[k] + 1)
        ranges.append(range(a, b + 1))
    import itertools as it
    for z in it.product(*ranges):
        ok = True
        for coeffs, denom in forms:
            # min over x in [0,1]^d of |form(z - x)|
            lo = 0
            hi = 0
            for k, cc in enumerate(coeffs):
                if cc:
                    t = cc * (z[k])
                    interval = [min(t - cc * 1, t), max(t - cc * 1, t)]
                    lo += interval[0]
                    hi += interval[1]
            # s = sum cc*(z_k - x_k) ranges over [lo, hi] / denom
            if max(abs(Fraction(lo, denom)), abs(Fraction(hi, denom))) <= M:
                pass  # fine
            elif lo > 0 or hi < 0:
                if min(abs(Fraction(lo, denom)), abs(Fraction(hi, denom))) > M:
                    ok = False
                    break
            # if 0 in [lo,hi]: min |.| = 0 <= M fine
        if ok:
            out.append(z)
    return out

# ---------------------------------------------------------------- float grid

def D_grid_float(v, forms, Z0, K):
    """Vectorized float D on the K^d grid over [0,1)^d.
    Returns (floatmax, argmax_index tuple)."""
    import numpy as np
    d = len(v) - 1
    fdata = []
    for coeffs, denom in forms:
        supp = [k for k in range(d) if coeffs[k] != 0]
        fdata.append((supp, [coeffs[k] / denom for k in supp]))
    grid_axes = [np.arange(K) / K for _ in range(d)]
    # broadcast shape: (K,K,...,K)
    best = np.full((K,) * d, np.inf)
    for z in Z0:
        acc = None
        for supp, cf in fdata:
            term = None
            for k, a in zip(supp, cf):
                sh = [1] * d
                sh[k] = K
                comp = (grid_axes[k].reshape(sh) - z[k]) * a
                term = comp if term is None else term + comp
            term = np.abs(term)
            acc = term if acc is None else np.maximum(acc, term)
        best = np.minimum(best, acc)
    idx = np.unravel_index(np.argmax(best), best.shape)
    return float(best[idx]), tuple(int(i) for i in idx)

def D_float_at(points, v, forms, Z0):
    """float D at given points (list of float tuples)."""
    import numpy as np
    pts = np.array(points, dtype=float)
    out = []
    for p in pts:
        best = np.inf
        for z in Z0:
            m = 0.0
            for coeffs, denom in forms:
                s = sum(a * (p[k] - z[k]) for k, a in enumerate(coeffs) if a)
                m = max(m, abs(s / denom))
                if m >= best:
                    break
            best = min(best, m)
        out.append(float(best))
    return out

# ---------------------------------------------------------------- hyperplanes / vertices

def hyperplanes_for(v, forms, Z0, M, box, include_valleys=True, pair_V_cut=None):
    """All hyperplanes {L.x = c} relevant inside `box` (list of (lo,hi) Fractions):
    valleys {form(x-z)=0}, ridges/exchanges {s form_a(x-z1) = s' form_b(x-z2)}."""
    d = len(v) - 1
    hyp = []  # (L coeffs tuple, c Fraction) -- L includes denominators folded in
    seen = set()

    def interval_of(L):
        lo = hi = Fraction(0)
        for k, a in enumerate(L):
            if a:
                t = a * box[k][0]
                u = a * box[k][1]
                lo += min(t, u)
                hi += max(t, u)
        return lo, hi

    def add(L, c):
        Ls = tuple(L)
        if (Ls, c) in seen:
            return
        # float pre-filter with generous slack (no false drops), then exact
        lo = hi = 0.0
        for k, a in enumerate(Ls):
            if a:
                t = a * float(box[k][0])
                u = a * float(box[k][1])
                lo += min(t, u); hi += max(t, u)
        slack = 1e-9 * (1.0 + abs(lo) + abs(hi) + abs(float(c)))
        if not (lo - slack <= float(c) <= hi + slack):
            return
        loe, hie = interval_of(Ls)
        if loe <= c <= hie:
            seen.add((Ls, c))
            hyp.append((Ls, c))

    if include_valleys:
        for z in Z0:
            for coeffs, denom in forms:
                L = list(coeffs)
                c = sum(coeffs[k] * z[k] for k in range(d))
                add(L, Fraction(c))   # form(x-z)=0 <=> sum coeffs*x = sum coeffs*z
    # pairs
    zs = list(Z0)
    for i1 in range(len(zs)):
        for i2 in range(i1, len(zs)):
            z1, z2 = zs[i1], zs[i2]
            dz = tuple(a - b for a, b in zip(z1, z2))
            if pair_V_cut is not None:
                if norm_V(dz, forms) > pair_V_cut:
                    continue
            for fa in forms:
                for fb in forms:
                    if z1 == z2 and fa == fb:
                        continue
                    for sa in (1, -1):
                        for sb in (1, -1):
                            if z1 == z2 and sa == sb and fa == fb:
                                continue
                            ca, da_ = fa
                            cb, db_ = fb
                            if da_ != db_:
                                # scale to common denom lcm
                                import math as _m
                                g = _m.gcd(da_, db_)
                                l = da_ * db_ // g
                                m1, m2 = l // da_, l // db_
                            else:
                                l = da_
                                m1 = m2 = 1
                            L = [0] * d
                            c = 0
                            for k in range(d):
                                a = sa * ca[k] * m1 - sb * cb[k] * m2
                                if a:
                                    L[k] = a
                                c += (sa * ca[k] * m1) * z1[k] - (sb * cb[k] * m2) * z2[k]
                            add(L, Fraction(c))   # no extra division: l already absorbed in m1,m2

    return hyp

def solve_vertex(hyps, d):
    """Solve d hyperplanes L.x = c exactly (L: plain coefficient list).
    Returns x tuple of Fractions or None."""
    A = []
    b = []
    for L, c in hyps:
        row = [Fraction(a) for a in L]
        if len(row) < d:
            row += [Fraction(0)] * (d - len(row))
        A.append(row)
        b.append(c)
    # Gauss elimination
    for col in range(d):
        piv = None
        for r in range(col, d):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        b[col], b[piv] = b[piv], b[col]
        inv = A[col][col]
        A[col] = [x / inv for x in A[col]]
        b[col] = b[col] / inv
        for r in range(d):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
                b[r] = b[r] - f * b[col]
    return tuple(b)

def _batch_vertices(hyps, box, d, chunk=200000):
    """All d-fold hyperplane intersections inside `box`: batched float solve
    with generous slack, then exact Fraction re-solve of survivors.
    Includes (d-1)-fold intersections with box faces (pad rank with face rows)."""
    import numpy as np
    from itertools import combinations
    m = len(hyps)
    if m < d:
        return
    Lf = np.array([[float(a) for a in L] for L, _ in hyps], dtype=float)
    cf = np.array([float(c) for _, c in hyps], dtype=float)
    blo = np.array([float(b[0]) for b in box], dtype=float)
    bhi = np.array([float(b[1]) for b in box], dtype=float)
    slack = 1e-7
    combos = combinations(range(m), d)
    while True:
        block = []
        try:
            for _ in range(chunk):
                block.append(next(combos))
        except StopIteration:
            pass
        if not block:
            break
        idx = np.array(block, dtype=int)
        A = Lf[idx]                      # (n, d, d)
        b = cf[idx][:, :, None]          # (n, d, 1)
        det = np.linalg.det(A)
        good = np.abs(det) > 1e-9 * np.maximum(1.0, np.abs(A).max(axis=(1, 2)) ** d)
        X = np.full((len(block), d), np.nan)
        if good.any():
            Xg = np.linalg.solve(A[good], b[good])   # (ng, d, 1)
            X[good] = Xg[:, :, 0]
        ok = np.all(np.isfinite(X), axis=1)
        Xok = X[ok]
        inside = np.all((Xok >= blo - slack) & (Xok <= bhi + slack), axis=1)
        # exact re-solve for survivors
        ok_idx = np.where(ok)[0][inside]
        for j in ok_idx:
            combo = block[int(j)]
            x = solve_vertex([hyps[i] for i in combo], d)
            if x is not None and in_box(x, box):
                yield x
        if len(block) < chunk:
            break

def in_box(x, box):
    return all(box[k][0] <= x[k] <= box[k][1] for k in range(len(x)))

# ---------------------------------------------------------------- rho: d=2 exact global

def rho_exact_d2(v, forms):
    """Exact rho for d=2: global arrangement enumeration on the torus box."""
    d = 2
    eps = Fraction(1, 50)
    box = [(Fraction(0) - eps, Fraction(1) + eps)] * d
    M = _box_corner_max(forms, box)
    Z0 = enumerate_Z0(v, forms, M)
    hyps = hyperplanes_for(v, forms, Z0, M, box,
                           include_valleys=True,
                           pair_V_cut=2 * M)
    cands = []
    for h1, h2 in combinations(hyps, 2):
        x = solve_vertex([h1, h2], d)
        if x is not None and in_box(x, box):
            cands.append(tuple(xi % 1 for xi in x))
    # face candidates: hyperplane meets face lines x_k = 0 / 1
    for h in hyps:
        L, c = h
        for faceval in (Fraction(0), Fraction(1)):
            for k in range(d):
                # solve L.x = c with x_k = faceval: 1 unknown
                other = [j for j in range(d) if j != k][0]
                a = L[other]
                if a == 0:
                    continue
                c2 = c - L[k] * faceval
                xo = c2 / a
                x = [None, None]
                x[k] = faceval
                x[other] = xo
                x = tuple(x)
                if in_box(x, box):
                    cands.append(tuple(xi % 1 for xi in x))
    cands += [(Fraction(0), Fraction(0))]
    best = Fraction(-1)
    arg = None
    for x in set(cands):
        val = D_exact(x, v, forms)
        if val > best:
            best, arg = val, x
    return best, arg, len(hyps), len(set(cands))

# ---------------------------------------------------------------- rho: d=3,4 (grid + local exact)

def rho_semid2(v, forms, K, verbose=False):
    """Exact rho for d in {3,4,5}: float grid (computed ONCE) -> certified-away
    + local exact arrangement analysis inside the high boxes.
    Returns dict; rho_exact is not None iff the away-certification succeeded."""
    import numpy as np
    from scipy import ndimage
    d = len(v) - 1
    M = M_corner(v, forms)
    Z0 = enumerate_Z0(v, forms, M)
    grid = _grid_D_numpy(v, forms, Z0, K)          # computed once
    fmax = float(grid.max())
    h = 1.0 / K
    rho_hi_f = fmax + h / 2 + 1e-6

    # torsion points (exact)
    torsion_vals = []
    for bits in range(2 ** d):
        x = tuple(Fraction(1, 2) if (bits >> k) & 1 else Fraction(0)
                  for k in range(d))
        torsion_vals.append((x, D_exact(x, v, forms)))
    rho_lo = max(v_ for _, v_ in torsion_vals)
    rho_lo_pts = [x for x, v_ in torsion_vals if v_ == rho_lo]

    # exact D at top-50 float grid points
    flat = grid.flatten()
    ntop = min(50, flat.size)
    idxs = np.argpartition(flat, -ntop)[-ntop:]
    for fi in idxs:
        idx = np.unravel_index(fi, grid.shape)
        x = tuple(Fraction(int(i), K) for i in idx)
        val = D_exact(x, v, forms)
        if val > rho_lo:
            rho_lo = val
            rho_lo_pts = [x]

    # high-set: grid cells with floatD >= rho_lo - h/2 - 2e-6
    mask = grid >= (float(rho_lo) - h / 2 - 2e-6)
    if not mask.any():
        idx = np.unravel_index(np.argmax(grid), grid.shape)
        mask = np.zeros_like(mask, dtype=bool)
        mask[idx] = True
    structure = np.ones((3,) * d, dtype=bool)      # full neighborhood
    lab, ncomp = ndimage.label(mask, structure=structure)
    boxes = []
    for ci in range(1, ncomp + 1):
        cells = np.argwhere(lab == ci)
        lo_i = cells.min(axis=0)
        hi_i = cells.max(axis=0)
        bnd = []
        for k in range(d):
            lo = Fraction(int(2 * lo_i[k]) - 5, 2 * K)   # lo_i/K - 1/(2K) - 2/K
            hi = Fraction(int(2 * (hi_i[k] + 1)) + 5, 2 * K)
            bnd.append((max(lo, Fraction(-1, 2)), min(hi, Fraction(3, 2))))
        boxes.append(bnd)

    # local exact analysis in each box
    rho_local_max = rho_lo
    local_arg = rho_lo_pts[0] if rho_lo_pts else None
    local_info = []
    radius = rho_hi_f + h
    for box in boxes:
        res = _local_exact(v, forms, Z0, box, radius=radius)
        if res is None:
            local_info.append(None)
            continue
        val, arg, nh, nv = res
        local_info.append((str(val), str(arg), nh, nv))
        if val > rho_local_max:
            rho_local_max = val
            local_arg = arg

    # away-certification: grid max outside dilated boxes (+ h/2) < rho_lo
    inbox = np.zeros(grid.shape, dtype=bool)
    gl = np.arange(K) / K
    for box in boxes:
        sl = []
        for k in range(d):
            i0 = np.searchsorted(gl, float(box[k][0]) - 1e-12, "left")
            i1 = np.searchsorted(gl, float(box[k][1]) + 1e-12, "right")
            sl.append(slice(i0, i1))
        inbox[tuple(sl)] = True
    away = np.where(~inbox, grid, -np.inf)
    away_max = float(away.max()) if np.isfinite(away.max()) else 0.0
    away_ok = bool(away_max + h / 2 + 1e-6 < float(rho_lo))
    return {
        "float_grid_max": fmax, "grid_K": K, "h": h,
        "rho_hi_grid": rho_hi_f,
        "rho_lo_exact": rho_local_max,
        "rho_exact": rho_local_max if away_ok else None,
        "away_certified": away_ok, "away_max": away_max,
        "torsion": [(tuple(str(c) for c in x), str(val)) for x, val in torsion_vals],
        "n_boxes": len(boxes), "boxes": [str(b) for b in boxes],
        "local_info": local_info,
        "argmax": local_arg,
        "Z0_size": len(Z0),
    }

def _grid_D_numpy(v, forms, Z0, K):
    import numpy as np
    d = len(v) - 1
    fdata = []
    for coeffs, denom in forms:
        supp = [k for k in range(d) if coeffs[k] != 0]
        fdata.append((supp, [coeffs[k] / denom for k in supp]))
    grid_axes = [np.arange(K) / K for _ in range(d)]
    best = np.full((K,) * d, np.inf)
    for z in Z0:
        acc = None
        for supp, cf in fdata:
            term = None
            for k, a in zip(supp, cf):
                sh = [1] * d
                sh[k] = K
                comp = (grid_axes[k].reshape(sh) - z[k]) * a
                term = comp if term is None else term + comp
            term = np.abs(term)
            acc = term if acc is None else np.maximum(acc, term)
        best = np.minimum(best, acc)
    return best

def _local_exact(v, forms, Z0, box, radius=None):
    """Exact max of D over `box` via local ridge/exchange arrangement.
    `radius`: relevance radius for z's (default: box-corner max).
    Valleys are excluded only when the box contains no lattice point
    (there D > 0 strictly and argmax-forms cannot flip sign)."""
    import math as _m
    d = len(v) - 1
    M = radius if radius is not None else _box_corner_max(forms, box)
    has_lattice = any(all(box[k][0] <= zk <= box[k][1] for k in range(d))
                      for zk in _box_corners_int(box, d))
    # z's relevant to the box: within M of box
    zloc = []
    for z in Z0:
        ok = True
        for coeffs, denom in forms:
            lo = hi = 0
            for k, cc in enumerate(coeffs):
                if cc:
                    t = cc * z[k]
                    interval = (min(t - cc * box[k][1], t - cc * box[k][0]),
                                max(t - cc * box[k][1], t - cc * box[k][0]))
                    lo += interval[0]
                    hi += interval[1]
            if lo > 0 or hi < 0:
                if min(abs(Fraction(lo, denom)), abs(Fraction(hi, denom))) > M:
                    ok = False
                    break
        if ok:
            zloc.append(z)
    cut = 2 * M
    hyps = hyperplanes_for(v, forms, zloc, M, box, include_valleys=has_lattice,
                           pair_V_cut=cut)
    from itertools import combinations
    cands = list(_batch_vertices(hyps, box, d))
    # box corners
    for corner in _box_corners(box):
        cands.append(corner)
    best = Fraction(-1)
    arg = None
    for x in set(cands):
        val = D_exact(x, v, forms)
        if val > best:
            best, arg = val, x
    return best, arg, len(hyps), len(set(cands))

def _box_corners(box):
    from itertools import product
    d = len(box)
    for c in product(*[b for b in box]):
        yield c

def _box_corner_max(forms, box):
    """max of ||c||_V over box corners (norm is convex -> max at corners)."""
    return max(norm_V(c, forms) for c in _box_corners(box))

def _box_corners_int(box, d):
    import math as _m
    ranges = []
    for k in range(d):
        a = _m.ceil(box[k][0])
        b = _m.floor(box[k][1])
        ranges.append(range(a, b + 1))
    from itertools import product
    for c in product(*ranges):
        yield c
