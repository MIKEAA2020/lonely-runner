import sys, time
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
from lrc_zono_lib import *

for v in [(1, 2, 3, 4), (1, 2, 3, 8)]:
    t0 = time.time()
    forms = make_forms(v)
    lam, tstar = lambda_exact(v)
    delta = F(1, 2) - lam
    thr = F(len(v) - 1, 2 * (len(v) + 1))
    res = rho_semid2(v, forms, K=48)
    print(f"v={v}: lambda={lam} delta={delta} thr={thr}")
    print(f"   float-grid max={res['float_grid_max']:.7f} (K=48) rho_hi_grid={res['rho_hi_grid']:.7f}")
    print(f"   rho_lo_exact={res['rho_lo_exact']} (= {float(res['rho_lo_exact']):.7f})")
    print(f"   rho_exact={res['rho_exact']} away_certified={res['away_certified']} away_max={res['away_max']:.7f}")
    print(f"   n_boxes={res['n_boxes']} Z0={res['Z0_size']} argmax={res['argmax']}")
    print(f"   deepest-hole: rho==delta -> {res['rho_lo_exact'] == delta}; rho<=thr: {res['rho_lo_exact'] <= thr}")
    print(f"   ({time.time()-t0:.1f}s)", flush=True)
