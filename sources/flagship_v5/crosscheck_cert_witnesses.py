#!/usr/bin/env python3
"""Cross-validation of the closure runs' witness-level claims against the
independent linear-form/lattice-search implementation (lrc_zono_lib.py,
the implementation behind this paper's provenance ledger).

Checks (all exact rational arithmetic, no floats in the verdicts):
  1-7. Witness values:
     m=10 torsion center            D = 7/22
     harmonic balance family k=20   D = 16/47
     harmonic balance family k=27   D = 16/47
     harmonic balance family k=35   D = 16/47   (the recorded witness)
     m=9 dyadic search point        D = 90449257437660783/252201579132747776
     m=12 clean witness             D = 17/47
     m=9 clean witness              D = 455/1269
  8-11. Binding structure at the four certified n=5-multiple center
     classes (m=10,20,25,30): the pair (1,m) binds at two translates
     with opposite signs (pair bisection), and D(c*) equals the claimed
     ladder value rho = delta.
Exit code 0 iff all checks pass.
"""
import sys
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F
import lrc_zono_lib as L

def check_value(name, v, c, claimed):
    forms = L.make_forms(v)
    got = L.D_exact(tuple(c), v, forms)
    ok = (got == claimed)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: D(c) = {got}"
          f"{'  == claimed' if ok else '  != claimed ' + str(claimed)}")
    return ok

results = []

# --- witness values ------------------------------------------------------
results.append(check_value("m=10 torsion center 7/22",
    (1, 2, 3, 4, 10),
    [F(1, 2), F(0), F(1, 2), F(1, 2)],
    F(7, 22)))

for k in (20, 27, 35):
    results.append(check_value(f"harmonic balance family c({k}) = 16/47",
        (1, 2, 3, 4, 5),
        [F(2 * k + 160, 235), F(3 * k - 13, 235),
         F(4 * k + 91, 235), F(5 * k, 235)],
        F(16, 47)))

results.append(check_value("m=9 dyadic search witness (deepest on record)",
    (1, 2, 3, 4, 5, 9),
    [F(1099170089337027, 36028797018963968),
     F(8466543868095245, 18014398509481984),
     F(1233377658248327, 9007199254740992),
     F(1848805658902507, 2251799813685248),
     F(8532950000663915, 18014398509481984)],
    F(90449257437660783, 252201579132747776)))

results.append(check_value("m=12 clean witness 17/47",
    (1, 2, 3, 4, 5, 12),
    [F(0), F(27, 47), F(0), F(25, 94), F(22, 47)],
    F(17, 47)))

results.append(check_value("m=9 clean witness 455/1269",
    (1, 2, 3, 4, 5, 9),
    [F(13, 423), F(22, 47), F(58, 423), F(695, 846), F(401, 846)],
    F(455, 1269)))

# --- binding structure at the certified ladder center classes ------------
sys.path.insert(0, '/home/z/my-project/download/paper2_sources_v3')
from rho_closeout import binding_structure

def check_ladder(m, rho):
    v = (1, 2, 3, 4, m)
    cs = L.cstar(v)
    forms = L.make_forms(v)
    d = L.D_exact(cs, v, forms)
    lam, _t = L.lambda_exact(v)
    bs = binding_structure(v, cs)
    pairs = {}
    for b in bs['bindings']:
        key = tuple(b['pair'])
        pairs.setdefault(key, []).append(F(b['value']))
    bisecting = [p for p, vals in pairs.items()
                 if len(vals) >= 2 and any(x < 0 for x in vals)
                 and any(x > 0 for x in vals)]
    bisect_int = 2 * (1 + m) * rho
    ok = (d == rho and F(1, 2) - lam == rho
          and bisect_int.denominator == 1
          and (1, m) in bisecting)
    print(f"[{'PASS' if ok else 'FAIL'}] m={m}: D(c*) = {d} (rho={rho}), "
          f"delta = {F(1,2)-lam}, 2(1+{m})rho = {bisect_int} "
          f"({'integer' if bisect_int.denominator == 1 else 'NOT integer'}), "
          f"pair-bisecting at: {bisecting}")
    return ok

for m, rho in [(10, F(7, 22)), (20, F(13, 42)),
               (25, F(4, 13)), (30, F(19, 62))]:
    results.append(check_ladder(m, rho))

print()
print(f"VERDICT: {sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
