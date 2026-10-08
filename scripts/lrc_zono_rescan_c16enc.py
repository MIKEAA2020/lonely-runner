"""n=4 m=16: enclosure-level record (exact certification exceeds the session budget)."""
import sys, json, time
sys.path.insert(0, "/home/z/my-project/scripts")
from fractions import Fraction as F
from lrc_zono_lib import lambda_exact, cstar, D_exact
from lrc_zono_rescan import rho_enclosure, _C_save, OUT
import numpy as np  # noqa

v = (1, 2, 3, 16)
t0 = time.time()
lam, _ = lambda_exact(v)
delta = F(1, 2) - lam
thr = F(3, 10)
enc = rho_enclosure(v, K=40)
rho_lo = enc["rho_lo_exact"]
enc.pop("forms")
rec = {"v": list(v), "m": 16, "n": 4, "lambda": str(lam), "delta": str(delta),
       "thr": str(thr), "mode": "enclosure-consistent (exact cert > 10 min, not run)",
       "rho_lo_exact": str(rho_lo), "rho_hi": enc["rho_hi"], "fmax": enc["fmax"],
       "tor_max": str(enc["tor_max"]),
       "witness": [str(t) for t in enc["witness"]],
       "witness_val": str(enc["witness_val"]),
       "deepest_hole_refuted": rho_lo > delta,
       "deepest_hole_certified": False,
       "deepest_hole_consistent": (rho_lo == delta and enc["fmax"] <= float(delta) + 1e-7),
       "rho_leq_thr_lo": rho_lo <= thr, "n_div_m": True,
       "sec": round(time.time() - t0, 1)}
_C_save(rec, "fail")
print(json.dumps(rec, indent=1, default=str))
