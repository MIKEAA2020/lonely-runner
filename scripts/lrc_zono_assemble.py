"""Assemble all results into one JSON + markdown report."""
import sys, json, math
sys.path.insert(0, "/home/z/my-project/scripts")

def spearman(xs, ys):
    n = len(xs)
    if n < 3: return None
    def ranks(a):
        order = sorted(range(n), key=lambda i: a[i]); r = [0.0]*n; i = 0
        while i < n:
            j = i
            while j+1 < n and a[order[j+1]] == a[order[i]]: j += 1
            avg = (i+j)/2 + 1
            for k in range(i, j+1): r[order[k]] = avg
            i = j+1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx)/n, sum(ry)/n
    num = sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    den = math.sqrt(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))
    return num/den if den else None

out1 = json.load(open("out_lrc_zono.json"))
out2b = json.load(open("out_lrc_zono2b.json"))
out3 = json.load(open("out_lrc_zono3.json"))

# n=4 exact core results (from the validated runs)
n4_core = [
 {"v": [1,2,3,4], "family": "harmonic", "lambda": "1/5", "delta": "3/10",
  "thr": "3/10", "rho": "3/10", "rho_f": 0.3, "method": "exact (d=3 grid+local, away-certified)",
  "deepest_hole": True, "rho_leq_thr": True, "argmax": "torsion c*=(1/2,0,1/2) (+ ridge)"},
 {"v": [1,2,3,8], "family": "rung2", "lambda": "2/9", "delta": "5/18",
  "thr": "3/10", "rho": "5/18", "rho_f": 0.2777778, "method": "exact (d=3 grid+local, away-certified)",
  "deepest_hole": True, "rho_leq_thr": True, "argmax": "torsion c*=(1/2,0,1/2)"},
 {"v": [1,2,3,4], "family": "sweep m=4", "lambda": "1/5", "delta": "3/10",
  "thr": "3/10", "rho": "3/10", "rho_f": 0.3, "method": "exact",
  "deepest_hole": True, "rho_leq_thr": True},
]

# n=4 correlation from sweep
sw4 = out2b["rho_n4_sweep_enclosure"]
sw4h = out2b["hstar_n4_sweep"]
sp4 = spearman([r["strictness_f"] for r in sw4h], [r["delta_margin_f"] for r in sw4h])

final = {
 "hstar_table": out1["hstar_table"],
 "rho_n3_core": out1["rho_n3_core"],
 "rho_n3_sweep": out1["rho_n3_sweep"],
 "hstar_n3_sweep": out1["hstar_n3_sweep"],
 "correlation_n3": out1["correlation_n3"],
 "rho_n4_core": n4_core,
 "rho_n4_sweep_enclosure": sw4,
 "hstar_n4_sweep": sw4h,
 "correlation_n4": {"strictness": [r["strictness_f"] for r in sw4h],
                    "delta_margin": [r["delta_margin_f"] for r in sw4h],
                    "spearman": sp4},
 "rho_n5_n6": out3,
}
with open("out_lrc_zono_final.json", "w") as f:
    json.dump(final, f, indent=1, default=str)
print("spearman n4:", sp4)
print("assembled -> out_lrc_zono_final.json")

# quick digest for the report
print("\n--- DIGEST ---")
for r in out1["hstar_table"]:
    print(f"{r['family']:8s} n={r['n']}: h*={r['hstar_str']} RR={r['real_rooted']} "
          f"LC={r['log_concave']} strict={r['strictness']}")
print("n3 covering holds everywhere:", all(r["rho_leq_thr"] for r in out1["rho_n3_sweep"]))
print("n3 deepest iff 3|m:", [(r["v"][2], r["deepest_hole"]) for r in out1["rho_n3_sweep"]])
print("n4 covering holds everywhere (sweep+core):",
      all(r["rho_leq_thr"] for r in sw4))
print("n4 deepest pattern:", [(r["m"], r["deepest_hole"]) for r in sw4])
for r in out3:
    print(f"n={len(r['v'])} {r['v']}: rho_lo={r['rho_lo_exact']} thr={r['thr']} "
          f"rho>thr={r['rho_gt_thr']} deepest={r['deepest_hole_possible']}")
