import sys, time
sys.path.insert(0, "/home/z/my-project/scripts")
import lrc_zono_lib as L
from fractions import Fraction as F
import numpy as np

v = (1, 2, 3, 4, 5)
forms = L.make_forms(v)
t0=time.time(); M = L.M_corner(v, forms); t_M=time.time()-t0
t0=time.time(); Z0 = L.enumerate_Z0(v, forms, M); t_Z=time.time()-t0
t0=time.time(); grid = L._grid_D_numpy(v, forms, Z0, 32); t_g=time.time()-t0
fmax = float(grid.max()); h=1/32; rho_hi = fmax + h/2 + 1e-6
print(f"M={float(M):.3f} |Z0|={len(Z0)} tM={t_M:.1f} tZ={t_Z:.1f} tgrid={t_g:.1f} fmax={fmax:.6f}", flush=True)
t0=time.time()
tor = []
for bits in range(16):
    x = tuple(F(1,2) if (bits>>k)&1 else F(0) for k in range(4))
    tor.append((x, L.D_exact(x, v, forms)))
rho_lo = max(v_ for _,v_ in tor)
print(f"torsion max = {rho_lo} ({float(rho_lo):.6f}) t={time.time()-t0:.1f}s", flush=True)
from scipy import ndimage
mask = grid >= (float(rho_lo) - h/2 - 2e-6)
lab, ncomp = ndimage.label(mask, structure=np.ones((3,)*4, dtype=bool))
boxes = []
for ci in range(1, ncomp+1):
    cells = np.argwhere(lab == ci)
    lo_i = cells.min(axis=0); hi_i = cells.max(axis=0)
    bnd = [(F(int(2*lo_i[k])-5, 64), F(int(2*(hi_i[k]+1))+5, 64)) for k in range(4)]
    boxes.append(bnd)
print(f"n_boxes={ncomp} sizes={[ [round(float(a),3),round(float(b),3)] for a,b in boxes[0] ]}", flush=True)
radius = rho_hi + h
# zloc count check
box = boxes[0]
t0=time.time()
zloc = []
for z in Z0:
    ok = True
    for coeffs, denom in forms:
        lo = hi = 0
        for k, cc in enumerate(coeffs):
            if cc:
                t = cc * z[k]
                interval = (min(t - cc*box[k][1], t - cc*box[k][0]), max(t - cc*box[k][1], t - cc*box[k][0]))
                lo += interval[0]; hi += interval[1]
        if lo > 0 or hi < 0:
            if min(abs(F(lo, denom)), abs(F(hi, denom))) > radius:
                ok = False; break
    if ok: zloc.append(z)
print(f"zloc(|radius={float(radius):.3f}|)={len(zloc)} t={time.time()-t0:.1f}s", flush=True)
t0=time.time()
res = L._local_exact(v, forms, Z0, box, radius=radius)
print(f"local: val={res[0]} arg={res[1]} nhyp={res[2]} ncand={res[3]} t={time.time()-t0:.1f}s", flush=True)
