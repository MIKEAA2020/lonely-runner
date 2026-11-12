#!/usr/bin/env python3
"""Resumable step-2 driver."""
import json, sys, time
sys.path.insert(0, '/home/z/my-project/scripts')
import higher_rung_step2 as S

PATH = '/home/z/my-project/scripts/out_higher_rung_step2.json'
try:
    OUT = json.load(open(PATH))
except Exception:
    OUT = {}

t0 = time.time()
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 480

if 'drift_law' not in OUT:
    OUT['drift_law'] = S.drift_law_check()
    json.dump(OUT, open(PATH, 'w'), indent=1)

for p in (5, 7, 11):
    k = 'gt_bf_p%d' % p
    if k not in OUT:
        OUT[k] = S.ground_truth_bruteforce(p, 7)
        json.dump(OUT, open(PATH, 'w'), indent=1)

for p in (5, 11, 13, 17, 19, 23, 29, 31, 37):
    k = 'gt_p%d' % p
    if k not in OUT:
        if time.time() - t0 > BUDGET:
            print('BUDGET OUT'); break
        OUT[k] = S.ground_truth_census_restricted(p, 7)
        json.dump(OUT, open(PATH, 'w'), indent=1)

for p in (17, 29, 43):
    k = 'exit_p%d' % p
    if k not in OUT:
        if time.time() - t0 > BUDGET:
            print('BUDGET OUT'); break
        OUT[k] = S.exit_time_top(p, 7)
        json.dump(OUT, open(PATH, 'w'), indent=1)

print('driver2 exit; %.0fs' % (time.time() - t0))
