import sys, time
sys.path.insert(0, "/home/z/my-project/scripts")
import lrc_zono_lib as L
from fractions import Fraction as F

v = (1, 2, 3, 4)
forms = L.make_forms(v)
t0=time.time(); M = L.M_corner(v, forms); t_M=time.time()-t0
t0=time.time(); Z0 = L.enumerate_Z0(v, forms, M); t_Z=time.time()-t0
t0=time.time(); grid = L._grid_D_numpy(v, forms, Z0, 48); t_g=time.time()-t0
fmax = float(grid.max()); h = 1/48; rho_hi = fmax + h/2 + 1e-6
# torsion
t0=time.time()
tor = []
for bits in range(8):
    x = tuple(F(1,2) if (bits>>k)&1 else F(0) for k in range(3))
    tor.append((x, L.D_exact(x, v, forms)))
t_tor=time.time()-t0
rho_lo = max(v_ for _,v_ in tor)
print(f"M={float(M):.3f} |Z0|={len(Z0)} times: M={t_M:.1f} Z0={t_Z:.1f} grid={t_g:.1f} torsion={t_tor:.1f}", flush=True)
# local exact on the two boxes (reproduce _high_boxes logic quickly)
import numpy as np
from scipy import ndimage
mask = grid >= (float(rho_lo) - h/2 - 2e-6)
lab, ncomp = ndimage.label(mask, structure=np.ones((3,)*3, dtype=bool))
t0=time.time()
boxes = []
for ci in range(1, ncomp+1):
    cells = np.argwhere(lab == ci)
    lo_i = cells.min(axis=0); hi_i = cells.max(axis=0)
    bnd = [(F(int(2*lo_i[k])-5, 96), F(int(2*(hi_i[k]+1))+5, 96)) for k in range(3)]
    boxes.append(bnd)
t_b=time.time()-t0
print(f"n_boxes={ncomp} boxes={boxes} t={t_b:.2f}", flush=True)
radius = rho_hi + h
for b in boxes:
    t0=time.time()
    res = L._local_exact(v, forms, Z0, b, radius=radius)
    print(f"local box -> val={res[0]} arg={res[1]} nhyp={res[2]} ncand={res[3]} t={time.time()-t0:.1f}s", flush=True)
