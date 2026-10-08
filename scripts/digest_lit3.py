#!/usr/bin/env python3
"""digest_lit3.py — digest round-2 lit search (r01-r16) into readable form."""
import json, glob, os

out = []
for f in sorted(glob.glob('/home/z/my-project/scripts/lit/r*.json')):
    key = os.path.basename(f).replace('.json', '')
    try:
        data = json.load(open(f))
    except Exception as e:
        out.append(f"### {key}  (ERROR: {e})")
        continue
    # data may be a list or wrapped
    items = data if isinstance(data, list) else data.get('results', data.get('data', []))
    out.append(f"### {key}  ({len(items)} results)")
    for it in items[:8]:
        name = (it.get('name') or '')[:110]
        url = it.get('url') or ''
        snip = (it.get('snippet') or '').replace('\n', ' ')[:260]
        date = it.get('date') or ''
        out.append(f"- [{name}] ({url}) {date}\n  {snip}")
    out.append("")

open('/home/z/my-project/scripts/lit/digest2.txt', 'w').write('\n'.join(out))
print('\n'.join(out[:80]))
print("...")
print(f"[full digest: scripts/lit/digest2.txt, {len(out)} lines]")
