"""Test exact rho (d=2) on n=3 cases — the pivotal deepest-hole test."""
import sys, time
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
from lrc_zono_lib import *

for v in [(1, 2, 3), (1, 2, 6), (1, 2, 4)]:
    t0 = time.time()
    forms = make_forms(v)
    lam, tstar = lambda_exact(v)
    delta = F(1, 2) - lam
    c = cstar(v)
    dstar = D_exact(c, v, forms)
    rho, arg, nhyp, ncand = rho_exact_d2(v, forms)
    thr = F(len(v) - 1, 2 * (len(v) + 1))
    print(f"v={v}: lambda={lam} delta={delta} thr={thr} "
          f"D(c*)={dstar} rho={rho} argmax={arg}")
    print(f"   deepest-hole (rho==delta): {rho == delta}   "
          f"rho-thr margin: {rho - thr}   hyps={nhyp} cands={ncand} "
          f"({time.time()-t0:.1f}s)")
    print(f"   LRC for this v: delta <= thr ? {delta <= thr}")
