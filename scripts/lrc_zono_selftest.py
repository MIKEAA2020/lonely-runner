"""Self-test of lrc_zono_lib on hand-verified cases."""
import sys
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
from lrc_zono_lib import *

print("=== Ehrhart formula vs hand values ===")
for v, e_exp in [((1, 2, 6), [1, 4, 9]), ((1, 2, 3), [1, 3, 6]),
                 ((1, 2, 4), [1, 4, 7]), ((2, 3), [1, 5]), ((1, 2), [1, 3])]:
    e = ehr_coeffs(v)
    print(f"v={v}: e={e}  expected={e_exp}  {'OK' if e == e_exp else 'FAIL'}")

print("\n=== h* extraction (two methods) ===")
for v, h_exp in [((1, 2, 6), [1, 11, 6]), ((1, 2, 3), [1, 7, 4]),
                 ((1, 2, 4), [1, 9, 4])]:
    d = len(v) - 1
    e = ehr_coeffs(v)
    h1 = hstar_from_ehr(e, d)
    Ehr_vals = [ehr_poly_eval(e, t) for t in range(d + 2)]
    h2 = hstar_by_inversion(Ehr_vals, d)
    # reconstruction check
    rec = [ehr_from_hstar(h1, d, t) for t in range(d + 2)]
    ok = h1 == h2 == h_exp and rec == Ehr_vals
    print(f"v={v}: h*_eulerian={h1} h*_inversion={h2} expected={h_exp} "
          f"Ehr(0..{d+1})={Ehr_vals} reconstr_ok={rec == Ehr_vals} -> {'OK' if ok else 'FAIL'}")

print("\n=== DFS cross-check ===")
for v in [(1, 2, 6), (1, 2, 3), (1, 2, 4), (2, 3), (1, 2)]:
    e = ehr_coeffs(v)
    for t in (1, 2):
        dfs = dfs_count(v, t)
        poly = ehr_poly_eval(e, t)
        print(f"v={v} t={t}: DFS={dfs} poly={poly} {'OK' if dfs == poly else 'FAIL'}")

print("\n=== lambda / delta ===")
for v, lam_exp in [((1, 2), F(1, 3)), ((2, 3), F(2, 5)), ((1, 3), F(1, 2)),
                   ((1, 4), F(2, 5)), ((1, 2, 3), F(1, 4)), ((1, 2, 6), F(2, 7))]:
    lam, t = lambda_exact(v)
    print(f"v={v}: lambda={lam} (t*={t}) expected>={lam_exp} "
          f"{'OK' if lam >= lam_exp else 'FAIL'}")

print("\n=== c-model: D(c*) == delta? ===")
for v in [(1, 2), (1, 2, 3), (1, 2, 4), (1, 2, 6)]:
    assert v[0] == 1
    lam, _ = lambda_exact(v)
    delta = F(1, 2) - lam
    forms = make_forms(v)
    c = cstar(v)
    dstar = D_exact(c, v, forms)
    print(f"v={v}: lambda={lam} delta={delta} D(c*)={dstar} "
          f"{'OK' if dstar == delta else 'FAIL'}")

print("\n=== Ehr(-1) reciprocity (interior points) ===")
for v in [(1, 2, 6), (1, 2, 3)]:
    e = ehr_coeffs(v)
    print(f"v={v}: Ehr(-1) = {ehr_poly_eval(e, -1)}  (Pick interior hexagon: 6 / 4)")

print("\n=== rho exact d=2: n=2 baseline (1,2) ===")
# d=1: B=[-3,3], rho should be 1/6 = delta(1,2)
forms = make_forms((1, 2))
# manual: rho = 1/6
import numpy as np
Z0 = enumerate_Z0((1, 2), forms, M_corner((1, 2), forms))
print("Z0 for (1,2):", Z0[:10], "... n=", len(Z0))
best = max(D_exact((F(k, 20),), (1, 2), forms) for k in range(20))
print("rho(1,2) ~ grid:", best, " expected 1/6 =", F(1, 6))
