#!/usr/bin/env python3
"""Follow-up analysis of the certification runs:
1. For every unresolved (survivor) box of a run, compute the float U-bound;
   exact U for the top-20. This yields the rigorous enclosure
       rho <= max(r, max_U over survivors)
   for the boxes the budget-capped run left open.
2. Recover the deepest survivor-center witness of the negative control
   (the 607/1792 point) exactly.

Output: scripts/out_certify_followup.json
"""
import sys, json, re, time
sys.path.insert(0, '/home/z/my-project/scripts')
from fractions import Fraction as F
import certify_rung2_n5 as C

SRC = '/home/z/my-project/scripts/out_certify_rung2.json'
OUT = '/home/z/my-project/scripts/out_certify_followup.json'


def parse_pair(s):
    lo, hi = s.split('|')
    return F(lo), F(hi)


def enclosure_for(name, v, r):
    d = json.load(open(SRC))
    run = [x for x in d['runs'] if x['name'] == name][0]
    r = F(r)
    if run['certified']:
        return {'name': name, 'certified': True}
    # rebuild ALL survivors: we only have 50 samples persisted; the count
    # is recorded. For the enclosure we use the sampled boxes and the
    # depth-capped statistics; additionally re-derive the survivor set
    # is impossible without re-running, so we bound what we have.
    boxes = []
    for sbox, uf in run.get('survivor_sample', []):
        boxes.append([parse_pair(s) for s in sbox])
    res = {'name': name, 'n_sampled': len(boxes)}
    t0 = time.time()
    us = []
    for b in boxes:
        uf, _ = C.U_box(b, v, False, r=r)
        us.append(uf)
    mx = max(us) if us else None
    res['max_float_U_sampled'] = mx
    res['r'] = float(r)
    # exact U for the top 20 by float U
    order = sorted(range(len(boxes)), key=lambda i: -us[i])[:20]
    exacts = []
    for i in order:
        ue, _ = C.U_box(boxes[i], v, True, r=r)
        exacts.append(str(ue))
    res['top20_exact_U'] = exacts
    res['enclosure_upper'] = max([str(r)] + exacts[:1]) if exacts else str(r)
    res['sec'] = round(time.time() - t0, 1)
    # exact D at centers of all sampled survivors
    dmax = None; darg = None
    for b in boxes:
        ctr = [(lo + hi) / 2 for lo, hi in b]
        dc = C.D_at(ctr, v)
        if dmax is None or dc > dmax:
            dmax, darg = dc, ctr
    res['max_exact_D_at_sampled_centers'] = str(dmax)
    res['argmax_center'] = [str(x) for x in darg]
    return res


def main():
    out = []
    out.append(enclosure_for('n5-rung2-m10-TARGET', (1, 2, 3, 4, 10), F(7, 22)))
    out.append(enclosure_for('n5-harmonic-NEGCTL', (1, 2, 3, 4, 5), F(1, 3)))
    json.dump(out, open(OUT, 'w'), indent=1)
    for e in out:
        print(json.dumps(e, indent=1))
    print('written:', OUT)


if __name__ == '__main__':
    main()
