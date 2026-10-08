# LR-Zonotope Ehrhart h\*-Program: Definitive Small-n Computation

**Session date:** 2026-10-07/08. **Task:** execute the 5-step plan from the audit
("compute h\* for n=3..6 both families; check real-rootedness, log-concavity,
strictness; compute ρ for n=3,4,5; test the deepest-hole conjecture; report the
correlation — or absence — between log-concavity strictness and the threshold
margin"). All computations exact (rational arithmetic) unless marked as float
enclosure. Scripts: `scripts/lrc_zono_lib.py`, `lrc_zono_run*.py`;
raw data: `scripts/out_lrc_zono_final.json`.

**Re-check (2026-10-08, audit response).** The audit flagged that the modular
law as stated below was internally inconsistent (its own n=5 harmonic entry,
m = 5 with 5 | 5, fails exactly) and that enclosure-level "consistent" labels
were presented as holds. Corrected in this version: (i) the n|m law is restated
with a certification ladder (exact n=3; certified n=4; broken at n=5 by the
harmonic instance; extinct at n=6); (ii) the n=4 core exact values (3/10, 5/18)
previously hardcoded in the assembler are re-established by persisted runs
(plus a newly certified m=12); (iii) the n=5 m-scan (m=5..30) and the n=6 scan
(m=7..11, 18, 24) are new; (iv) a third exact covering counterexample
(n=6, v=(1,2,3,4,5,7), ρ ≥ 4/11 > 5/14) is added. Re-check data:
`scripts/out_lrc_zono_rescan.json` (copy: `download/lrc_zonotope_ehrhart_rescan.json`);
scripts: `scripts/lrc_zono_rescan.py`, `lrc_zono_rescan_strong.py`.

## 0. The model (all cross-validated)

- V = R^n/Rv, Λ = π(Z^n) ≅ Z^{n−1}, Z(v) = π([0,1]^n), d = n−1.
- **Ehrhart formula** (derived, then verified two independent ways — DFS
  coset-counting and Pick/reciprocity — on every instance):
  Ehr_Z(t) = Σ_k e_k t^k, e_k = Σ_{|S|=k} gcd(v_{S^c}), gcd(∅)=0; e_d = Σv_j = Vol.
- **h\*** extracted by Eulerian expansion and by the standard inversion
  (asserted equal; Ehr reconstructed from h\* at t=0..d+1 asserted equal; DFS
  count asserted equal to the polynomial at t=1,2 on every instance).
- **Norm model:** ||c||_V = max_{i<j} |φ^{ij}| with φ^{1j} = c_j/(1+v_j),
  φ^{ij} = (v_j c_i − v_i c_j)/(v_i+v_j) in c-coordinates c_j = x_j − v_j x₁
  (v₁=1 ⇒ Λ = Z^d). Proved here: every form has ℓ₁-norm ≤ 1, hence
  **||c||_V ≤ ||c||_∞** and the unit ball B ⊇ [−1,1]^d — grid certification of
  ρ needs no norm-equivalence constant (ρ ≤ gridmax + h/2).
- **δ(v) = ½ − λ(v) = dist_V(c\*, Λ)**, c\* = π(½𝟙) (a 2-torsion point:
  c\*_j = ½ iff v_j even). This identity was asserted computationally on every
  instance and holds exactly throughout — the geometric pipeline is validated
  end-to-end.
- **Correct covering reform of LRC (this normalization):**
  LRC ⟺ ∀v: δ(v) ≤ (n−1)/(2(n+1)).
  The covering-radius statement "∀v: ρ(Z(v)) ≤ (n−1)/(2(n+1))" is strictly
  stronger (δ ≤ ρ always) — see §3: it is **false from n = 5 on**.

## 1. h\*-polynomials: real-rooted everywhere, strictness DECAYS

| family | n | h\*(t) | real-rooted | log-concave | Newton | min ratio h_k²/(h_{k−1}h_{k+1}) |
|---|---|---|---|---|---|---|
| harmonic | 3 | 1+7t+4t² | yes | yes | yes | 49/4 = 12.25 |
| rung-2 | 3 | 1+11t+6t² | yes | yes | yes | 121/6 ≈ 20.17 |
| harmonic | 4 | 1+18t+35t²+6t³ | yes | yes | yes | 9.257 |
| rung-2 | 4 | 1+22t+51t²+10t³ | yes | yes | yes | 9.490 |
| harmonic | 5 | 1+37t+179t²+133t³+10t⁴ | yes | yes | yes | 6.511 |
| rung-2 | 5 | 1+45t+239t²+181t³+14t⁴ | yes | yes | yes | 7.013 |
| harmonic | 6 | 1+78t+744t²+1286t³+399t⁴+12t⁵ | yes | yes | yes | 5.518 |
| rung-2 | 6 | 1+86t+920t²+1682t³+535t⁴+16t⁵ | yes | yes | yes | 5.748 |
| harmonic | 7 | 1+147t+2522t²+8948t³+7323t⁴+1201t⁵+18t⁶ | yes | yes | yes | 4.335 |
| rung-2 | 7 | 1+161t+3024t²+11144t³+9295t⁴+1551t⁵+24t⁶ | yes | yes | yes | 4.418 |
| harmonic | 8 | 1+290t+8317t²+51458t³+82947t⁴+35270t⁵+3135t⁶+22t⁷ | yes | yes | yes | 3.791 |
| rung-2 | 8 | 1+298t+9277t²+60986t³+102275t⁴+44798t⁵+4095t⁶+30t⁷ | yes | yes | yes | 3.829 |

Answers to the conjecture's three checks:
1. **Real-rooted: yes, uniformly** (also in the sweeps: all (1,2,m) m≤14 and
   (1,2,3,m) m≤12). h\*(1) grows; degrees = d = n−1 as expected; the
   Ehrhart reciprocity Ehr(−1) = #(interior lattice points) holds (e.g. 6 for
   the (1,2,6) hexagon, matching Pick).
2. **Strictness does NOT grow with n — it decays**, roughly ~ 33/n
   (12.25 → 9.26 → 6.51 → 5.52 → 4.34 → 3.79). The conjecture's premise is
   backwards.
3. **No reproduction of any covering threshold** — and it could not: the
   covering statement itself fails (§3). There is no mechanism connecting
   coefficient inequalities of a counting polynomial to a fixed-scale metric
   hole depth; the numbers below show not even a correlation.

λ data: harmonic λ = 1/(n+1) exactly for all n ≤ 8 (margin 0, tight);
rung-2 λ = 2/(2n+1) exactly (margin 1/((n+1)(2n+1))). The rung-2 family
{1,…,n−1,2n} is **not** extremal — confirming the audit's correction.

## 2. Exact covering radii (n=3, n=4) and the deepest-hole law

Method: d=2 — full exact arrangement enumeration (valleys + ridges +
exchange hyperplanes of the V-norm Voronoi structure; verified against
independent float grids at K=1200–2000 on every instance). d=3 — certified
float grid (K=48/40) + **exact** local arrangement analysis inside the
high-set boxes + away-certification (rigorous: ρ ≤ awaymax + h/2 < ρ_lo
outside the boxes; the max inside is an exact arrangement vertex).

| v | λ | δ | ρ | deepest-hole (ρ=δ)? | ρ vs thr=(n−1)/(2(n+1)) |
|---|---|---|---|---|---|
| (1,2,3) | 1/4 | 1/4 | **1/4 (exact)** | yes | = thr (tight) |
| (1,2,6) | 2/7 | 3/14 | **3/14 (exact)** | yes | < thr |
| (1,2,3,4) | 1/5 | 3/10 | **3/10 (exact, certified)** | yes | = thr (tight) |
| (1,2,3,8) | 2/9 | 5/18 | **5/18 (exact, certified)** | yes | < thr |

**Sweep (1,2,m), n=3, m=3..24 (all exact):** ρ ≤ 1/4 in every instance — the covering
statement **holds at n=3** — but the deepest-hole property holds **iff 3 | m** (true:
m=3,6,9,12,15,18,21,24; false for all others, with exact counterexample witnesses,
e.g. D(5/14,1/14)=3/14 > δ=1/6 at m=4).
ρ takes values 1/4, 3/14, 9/46, 1/5, 3/16, 5/26, 15/82, 7/38, 2/11, 3/17, 9/50 … —
all ≤ thr.

**Sweep (1,2,3,m), n=4, m=4..16:** covering statement **holds in every instance**
(ρ ≤ 3/10); deepest-hole holds **iff 4 | m**, now with full certification on both
sides: m = 4, 8, 12 are away-certified exact (ρ = 3/10, 5/18, 7/26 = δ — the
previous session's hardcoded core values are re-established by persisted runs,
and m=12 is newly certified); m = 5,6,7,9,10,11,13,14,15 fail with exact witnesses
(matching the earlier enclosure values); m = 16 is consistent (ρ_lo = δ = 9/34,
grid max = δ, strong multi-start search finds nothing above δ; the exact d=3
certification exceeded the 10-minute session budget).

**The modular law, corrected.** The previous version of this report stated the
law as "deepest-hole iff n | m for n = 3,4,5 (rung-2 m=2n included, and harmonic
m=n)". That statement was internally inconsistent with this report's own §3
(at n=5 harmonic, m = 5 and 5 | 5, deepest-hole fails exactly), and it presented
enclosure-level "consistent" labels (n=5 m=10) as holds. The corrected record:

| n | tested m | deepest-hole (ρ = δ) status |
|---|---|---|
| 3 | 3..24 | **iff 3\|m, exact** (holds m = 3,6,…,24; exact witnesses against all others) |
| 4 | 4..16 | **iff 4\|m, certified** (m = 4,8,12 away-certified exact; m = 16 consistent; exact witnesses against all non-multiples) |
| 5 | 5..30 | **law false**: m = 5 (harmonic, 5\|5) and m = 15 fail **exactly**; m = 10, 20, 25, 30 consistent (uncertified); all 20 non-multiples fail exactly |
| 6 | 6..12, 18, 24 | **extinct**: every tested instance fails exactly — including harmonic m=6, rung-2 m=12, 3n (m=18), 4n (m=24) |

So "iff n | m" is a small-n phenomenon: genuine and exact at n=3, certified at
n=4, broken at n=5 by the harmonic instance itself (m=5, 5|5: ρ ≥ 97/288 > 1/3
= δ), and extinct at n=6, where no tested m — including every multiple of 6 up
to 4n — admits a deepest-hole center. At n=5 the surviving multiples (m = 10,
20, 25, 30) carry the hold signature (torsion c\* attains δ, the K=32 grid
maximum equals δ, a flat ridge passes through c\*, and a 288-start Nelder–Mead
search with exact re-evaluation finds nothing above δ), but global equality
remains uncertified at d=4. The modular framing does not survive as a law;
what survives is the certification ladder itself: exact refutations everywhere,
exact/certified holds only at n ≤ 4.

## 3. The covering-radius statement is FALSE from n = 5 (exact witnesses)

| v | δ = ½−λ | thr | ρ (lower bound, exact witness) | verdict |
|---|---|---|---|---|
| (1,2,3,4,5) harmonic | 1/3 | 1/3 | **≥ 97/288 ≈ 0.336806** at x=(1/32,15/32,3/32,7/8) | **ρ > thr** — covering REFUTED; enclosure [0.3368, 0.3524] |
| (1,2,3,4,10) rung-2 | 7/22 | 1/3 | ≥ 7/22 = δ at torsion c\*; enclosure [7/22, 0.3338] | ρ < thr possible; deepest-hole **consistent, uncertified** (strong search finds nothing above δ) |
| (1,2,3,4,5,6) harmonic | 5/14 | 5/14 | **≥ 29/80 = 0.3625** at x=(0,5/16,15/16,3/4,1/2) | **ρ > thr** — covering REFUTED; enclosure [0.3625, 0.3938] |
| (1,2,3,4,5,7) | 1/3 | 5/14 | **≥ 4/11 ≈ 0.363636** at x=(1/16,9/16,0,3/8,0) | **ρ > thr** — covering REFUTED at a **non-extremal** instance (re-check) |
| (1,2,3,4,5,12) rung-2 | 9/26 | 5/14 | **≥ 5/14** at a non-torsion point | ρ ≥ thr; deepest-hole false |

Witness values verified independently by targeted brute-force z-scans
(97/288 and 29/80 reproduce exactly) and re-asserted exactly in the re-check
(Phase A). All covering-refuting witnesses are **non-torsion**
points: the deepest holes of the harmonic zonotopes at n=5,6 (and of the
non-extremal n=6 m=7 instance) are NOT the center class. LRC itself is
untouched (δ = thr exactly at harmonic — tight, consistent).

**Consequence:** "∀v: ρ(Z(v)) ≤ (n−1)/(2(n+1))" is strictly stronger than
LRC and is false for n ≥ 5 — with **three exact witnesses** (n=5 harmonic;
n=6 harmonic; n=6 m=7, a non-extremal instance), and the n=6 m=9 instance
sits at ρ ≥ thr exactly. The n=5 consistent cases (m=10,20,25,30) do not
threaten the refutation: they are below thr. The equivalence the program
needs is only the center-depth form
LRC ⟺ ∀v: δ(v) ≤ (n−1)/(2(n+1)), which the zonotope shares with the
torus/covering reformulation — no new leverage. The user's originally stated
form "LRC ⟺ ρ(Z) ≥ 1/(2n+2)" is wrong in both direction and constant in
every normalization we can construct.

## 4. Correlation: none

- n=3 sweep (12 instances): Spearman(strictness, δ-margin) = **+0.32**.
- n=4 sweep (8 instances): Spearman(strictness, δ-margin) = **−0.31**.
- The sign flips between n=3 and n=4; the harmonic family has margin ≡ 0 at
  every n while strictness varies by 3×. There is no consistent statistical
  relation, and no mechanism: log-concavity strictness is a property of the
  coefficient vector of a counting polynomial at growing scale t; the
  threshold margin is a fixed-scale statement about one point (the center)
  in a quotient torus. This is the same type-mismatch the T-25/T-26 analysis
  identified on the dynamics side — it reappears here on the metric side.

## 5. Verdict for the program

1. The h\*/Ehrhart object is real and beautifully behaved (real-rooted,
   log-concave, Newton — uniformly in everything tested). The first move of
   the proposed conjecture succeeded.
2. The bridge is dead, twice over: (a) deepest-hole (δ=ρ) holds only in the
   small-n segment — iff 3|m exactly at n=3 (m ≤ 24), iff 4|m certified at
   n=4 (m ≤ 15; m=16 consistent), sporadic consistent multiples at n=5
   (m=10,20,25,30, uncertified), extinct at n=6 (all tested m fail, including
   every multiple of 6 through 4n) — and the harmonic instance itself fails
   from n=5 on; (b) the covering-radius statement is false from n=5 — with
   three exact witnesses, two at the extremal harmonic instance and one at a
   non-extremal instance (n=6, m=7). Log-concavity strictness decays with n
   and does not correlate with margins.
3. What survives: LRC ⟺ center-depth δ ≤ (n−1)/(2(n+1)) — but δ = ½ − λ is
   a torus quantity, and the zonotope adds no new handle on it; the Ehrhart
   data (the arithmetic-matroid gcd structure) does not visibly pin the
   center depth. The type-mismatch diagnosis stands, now with a computed
   counterexample archive instead of a heuristic.
4. Honest scope: n ≤ 6 for ρ (exact n≤4 core + n=3 sweeps m ≤ 24; enclosures with
   exact witnesses at n=5 (m ≤ 30) and n=6 (m = 6..12, 18, 24)); h\* to n=8; two
   one-parameter families + both core families. The n=8,9,10 "frontier"
   attributions in the original proposal remain unverifiable and should not
   be cited. Certification labels are now uniform: **exact** (n=3 all; n=4
   holds m=4,8,12 + all fail witnesses; all n=5/6 fail witnesses), **certified**
   (n=4 holds via away-certification), **consistent** (n=4 m=16; n=5
   m=10,20,25,30 — search-level evidence only, global equality open at d=4).

## Appendix: validation ledger

- DFS = Ehrhart-polynomial at t=1,2 on all 30+ instances (incl. v₁>1 case (2,3)).
- h\* (Eulerian) = h\* (inversion) = reconstructs Ehr — asserted everywhere.
- Pick/Ehrhart reciprocity: Ehr(−1) = interior count ((1,2,6): 6 ✓).
- D(c\*) = δ asserted exactly on every c-model instance (re-asserted on all
  60+ re-check instances).
- d=2 exact ρ vs independent float grids K=1200–2000: consistent on all instances.
- d=3 exact results carry the away-certification (grid max outside boxes + h/2 < ρ_lo);
  re-check: n=4 m=4, 8, 12 re-run and persisted (17 s / 78 s / 134 s), each
  cross-checked against an independent K=72 grid upper bound.
- n=5,6 witnesses re-verified by targeted brute force (exact match) and
  re-asserted exactly (Phase A: 97/288, 7/22, 29/80, 5/14).
- Re-check fail-witness margins at n=5/6 are all ≥ 7/5200 ≈ 0.00135 (n=6 m=24,
  the narrowest found); the strong multi-start search (288 Nelder–Mead starts
  from random + torsion + grid seeds, exact re-evaluation of the top 30 refined
  optima) reproduces the known fails and finds nothing above δ at the
  consistent candidates (n=5 m=10,20,25,30; n=4 m=16).
- Known corrected en route: two hyperplane-scaling bugs found by requiring
  agreement with float grids (the arrangement code initially returned
  89/420 < true 3/14 at (1,2,4); fixed and re-validated everywhere); and the
  audit-caught reporting errors (modular-law inconsistency; hardcoded n=4 core
  values without persisted runs; consistent-labels presented as holds) fixed in
  this version.
