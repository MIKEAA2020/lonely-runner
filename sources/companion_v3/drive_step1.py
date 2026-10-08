#!/usr/bin/env python3
"""Resumable driver: runs the step-1 censuses one prime at a time,
flushing to out_higher_rung_step1.json after each; skips completed work."""
import json, sys, time
sys.path.insert(0, '/home/z/my-project/scripts')
import higher_rung_step1 as H

PATH = '/home/z/my-project/scripts/out_higher_rung_step1.json'
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 420.0

def load():
    try:
        return json.load(open(PATH))
    except Exception:
        return {}

OUT = load()
if 'machinery_T6_control' not in OUT:
    OUT['machinery_T6_control'] = H.machinery_control()
    OUT['footprint_control'] = H.footprint_control()

T7 = [29, 43, 71, 23, 37, 79, 17, 31, 59, 11, 53, 19, 47, 61, 13, 41, 83]
T8 = [17, 41, 73, 19, 43, 29, 31, 47]

def done(res_key, p):
    r = OUT.get(res_key, {})
    return ('p%d' % p) in r

t0 = time.time()
for p in T7:
    if done('T7_census', p):
        continue
    if time.time() - t0 > BUDGET:
        print('BUDGET OUT; resume later'); break
    if 'T7_census' not in OUT: OUT['T7_census'] = {}
    fams, meta = H.census_m4(p, 7)
    st = H.family_stats(fams)
    pat, patd = st['pat'], st['pat_distinct']
    OUT['T7_census']['p%d' % p] = {
        'eps': meta['eps'], 'k': meta['k'], 'cls': meta['cls'],
        'n_families': len(pat), 'n_configs': sum(pat.values()),
        'n_distinct_families': len(patd),
        'n_distinct_configs': sum(patd.values()),
        'top_families': [[list(k), v] for k, v in pat.most_common(15)],
        'max_proj_per_family': {str(k): dict(v) for k, v in
                                st['maxproj'].items()},
        'families': {str(k): v for k, v in sorted(pat.items())},
        'famsize_detail': {str(k): {'n': v['n'],
                                    **{nm: s for nm, s in v.items()
                                       if nm != 'n'}} for k, v in
                           fams.items()},
    }
    print('p=%3d eps=%d k=%d: %4d fams %6d cfgs; distinct %d; maxproj %s' %
          (p, meta['eps'], meta['k'], len(pat), sum(pat.values()),
           len(patd), {nm: max(mp[nm] for mp in st['maxproj'].values()
                               if nm in mp)
                       for nm in ('a2', 'd32', 'd43', 'd42')}))
    sys.stdout.flush()
    json.dump(OUT, open(PATH, 'w'), indent=1)

for p in T8:
    if done('T8critical_census', p):
        continue
    if time.time() - t0 > BUDGET:
        print('BUDGET OUT; resume later'); break
    if 'T8critical_census' not in OUT: OUT['T8critical_census'] = {}
    fams, meta = H.census_m4(p, 8)
    st = H.family_stats(fams)
    pat, patd = st['pat'], st['pat_distinct']
    OUT['T8critical_census']['p%d' % p] = {
        'eps': meta['eps'], 'k': meta['k'], 'cls': meta['cls'],
        'n_families': len(pat), 'n_configs': sum(pat.values()),
        'n_distinct_families': len(patd),
        'n_distinct_configs': sum(patd.values()),
        'top_families': [[list(k), v] for k, v in pat.most_common(15)],
        'max_proj_per_family': {str(k): dict(v) for k, v in
                                st['maxproj'].items()},
        'families': {str(k): v for k, v in sorted(pat.items())},
        'famsize_detail': {str(k): {'n': v['n'],
                                    **{nm: s for nm, s in v.items()
                                       if nm != 'n'}} for k, v in
                           fams.items()},
    }
    print('T8c p=%3d eps=%d k=%d: %4d fams %6d cfgs; distinct %d' %
          (p, meta['eps'], meta['k'], len(pat), sum(pat.values()),
           len(patd)))
    sys.stdout.flush()
    json.dump(OUT, open(PATH, 'w'), indent=1)

if all(done('T7_census', p) for p in T7) and \
   all(done('T8critical_census', p) for p in T8):
    if 'T8_m5_census' not in OUT:
        # m5 census is expensive; do p=17 only in this driver
        res = {}
        counts, meta = H.census_m5(17, 8)
        pat = {}
        for (s1, d2, d3, d4, d5), n in counts.items():
            pat.setdefault((d2, d3, d4, d5), 0)
            pat[(d2, d3, d4, d5)] += n
        patd = {k: v for k, v in pat.items()
                if len({k[0], k[1], k[2], k[3]}) == 4 and 1 not in k}
        res['p17'] = {'eps': meta['eps'], 'k': meta['k'],
                      'n_families': len(pat),
                      'n_configs': sum(pat.values()),
                      'n_distinct_families': len(patd),
                      'n_distinct_configs': sum(patd.values()),
                      'top_families': [[list(k), v] for k, v in
                                       sorted(pat.items(),
                                              key=lambda kv: -kv[1])[:15]],
                      'families': {str(k): v for k, v in
                                   sorted(pat.items())}}
        OUT['T8_m5_census'] = res
        print('m5 p=17: %d fams %d cfgs; distinct %d fams' %
              (len(pat), sum(pat.values()), len(patd)))
        json.dump(OUT, open(PATH, 'w'), indent=1)
    print('ALL STEP1 WORK COMPLETE')
print('driver exit; elapsed %.0fs' % (time.time() - t0))
