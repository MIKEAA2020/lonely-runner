import sys, time
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
from lrc_zono_lib import *

for v, K in [((1,2,3,8), 48), ((1,2,3,4,5), 32), ((1,2,3,4,10), 32)]:
    t0 = time.time()
    forms = make_forms(v)
    lam, tstar = lambda_exact(v)
    delta = F(1, 2) - lam
    thr = F(len(v) - 1, 2 * (len(v) + 1))
    res = rho_semid2(v, forms, K=K)
    print(f"v={v}: lambda={lam} delta={delta} thr={thr}", flush=True)
    print(f"   grid(K={K}) max={res['float_grid_max']:.7f} rho_hi={res['rho_hi_grid']:.7f}"
          f" rho_lo_exact={res['rho_lo_exact']} ({float(res['rho_lo_exact']):.7f})", flush=True)
    print(f"   rho_exact={res['rho_exact']} away_ok={res['away_certified']}"
          f" n_boxes={res['n_boxes']} Z0={res['Z0_size']}", flush=True)
    print(f"   argmax={res['argmax']}", flush=True)
    print(f"   deepest-hole: {res['rho_lo_exact'] == delta};  rho<=thr: {res['rho_lo_exact'] <= thr}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
