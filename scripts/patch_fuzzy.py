#!/usr/bin/env python3
"""patch_fuzzy.py — whitespace-robust patcher.
Each patch: (file, old, new). Matching ignores the exact positions of
line breaks (compares on whitespace-normalized text) but replaces the
full original byte range with the new text verbatim."""
import sys, os, re
import os as _os
D = _os.environ.get('PATCH_DIR', '/home/z/my-project/download/paper2_sources_v2')

def norm(s):
    return re.sub(r'\s+', ' ', s)

def apply(path, old, new):
    s = open(path).read()
    n_old, n_new = norm(old), norm(new)
    if n_old == n_new:
        return 'NO-OP (identical after normalize)'
    if old in s:
        open(path, 'w').write(s.replace(old, new, 1))
        return 'OK (exact)'
    # fuzzy: find the span whose normalized form matches
    target = n_old
    # build char-offset map of normalized text
    out, spans = [], []
    i = 0
    for m in re.finditer(r'\s+', s):
        out.append(s[i:m.start()]); spans.append((len(''.join(out)) + len(out)*0,))
        i = m.end()
    # simpler: incremental scan
    npos = []
    nstr = []
    j = 0
    last = 0
    posmap = []  # (index_in_norm, index_in_s)
    for m in re.finditer(r'\s+', s):
        chunk = s[last:m.start()]
        for k, ch in enumerate(chunk):
            posmap.append((j, last + k)); j += 1
        # one space for the whitespace run
        posmap.append((j, m.start())); j += 1
        last = m.end()
    chunk = s[last:]
    for k, ch in enumerate(chunk):
        posmap.append((j, last + k)); j += 1
    normed = ''.join(s[p[1]] if False else (' ' if False else '') for p in posmap)
    # reconstruct normalized string directly
    normed = re.sub(r'\s+', ' ', s)
    idx = normed.find(target)
    if idx < 0:
        return 'MISS'
    # map normalized indices back to raw indices
    # posmap: list of (norm_idx, raw_idx) built above
    back = {}
    for n_i, r_i in posmap:
        back.setdefault(n_i, r_i)
    start = back.get(idx)
    end_n = idx + len(target) - 1
    end = back.get(end_n)
    if start is None or end is None:
        return 'MAPFAIL'
    end_raw = end + 1
    # verify
    if norm(s[start:end_raw]) != target:
        return 'VERIFYFAIL'
    open(path, 'w').write(s[:start] + new + s[end_raw:])
    return 'OK (fuzzy)'

if __name__ == '__main__':
    import importlib.util
    spec = importlib.util.spec_from_file_location('patches', sys.argv[1])
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for f, old, new in mod.P:
        print(apply(os.path.join(D, f), old, new), '|', f, '|', old[:55].replace('\n', ' '))
