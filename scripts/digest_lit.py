#!/usr/bin/env python3
"""digest_lit.py — compact digest of the LRC literature search JSONs."""
import json, glob, os

OUT = "/home/z/my-project/scripts/lit"
files = sorted(glob.glob(os.path.join(OUT, "q*.json")))
lines = []
for f in files:
    try:
        data = json.load(open(f))
    except Exception as e:
        lines.append(f"### {os.path.basename(f)} — PARSE ERROR {e}")
        continue
    # data may be a list or a dict wrapper
    items = data if isinstance(data, list) else data.get("results", data.get("data", []))
    lines.append(f"### {os.path.basename(f)}  ({len(items)} results)")
    for it in items:
        if not isinstance(it, dict):
            continue
        name = (it.get("name") or "").strip()[:110]
        url = (it.get("url") or "").strip()[:120]
        host = (it.get("host_name") or "").strip()[:40]
        date = (it.get("date") or "").strip()[:12]
        snip = (it.get("snippet") or "").replace("\n", " ").strip()[:320]
        lines.append(f"- [{name}] ({host}, {date})\n  {url}\n  {snip}")
    lines.append("")
open("/home/z/my-project/scripts/lit/digest.txt", "w").write("\n".join(lines))
print(f"wrote digest: {len(lines)} lines from {len(files)} files")
