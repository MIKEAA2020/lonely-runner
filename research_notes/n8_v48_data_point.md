# The V=48 Data Point: Both Decay Series Completed

**Program:** LRC residue-language reformulation (pair-sum finitization).
**Charge (reviewer, optional):** "If compute allows, running the augmented
battery at V=48 would confirm the extrapolation. If augmented coverage at
V=48 is still above 80%, the decay is slower than extrapolated and the
paper's framing can be more confident. If it drops below 78%, the decay is
on the extrapolated track, and the paper's framing should be cautious."
**Committed script:** `scripts/N8_v48_extension.py`, log
`scripts/n8_run.log` (538.8s, deterministic). Battery byte-identical to
committed N6 (same: P4 at {7,13,17}) and N7 (augmented: P4 at
{7,13,17,19,37}); the only change is range — tables/capability/sigcap for
all N ≤ 95 (pair sums reach 2·48−1 = 95).

---

## 1. Controls (all exact)

- n = 4 (v ≤ 16): 1651/1661/1744 of 1745, OPEN 1, open set (3,5,8,13).
- Same battery: V=16 4179/4193/4311; V=24 35778/36362/41656 (V-tier
  4878+416); V=32 156461/161596/196751; V=40 480454/500280/641166 — every
  committed n6_run.log number.
- Augmented battery: V=24 38611/38924/41656; V=32 173337/176376/196751;
  V=40 526607/540415/641166 — every committed n7_run.log number; the
  unflagged-allinv-at-19/37 self-check is 0 at every range.
- Capability: n=4 {7,11,13}; n=5 N≤47 {7,13,17,19,37} with counts
  20/68/16/64/32; capable moduli in 48..79 = [] (committed); sigcap map
  for 48..79 exact vs committed.
- Soundness (new range): Lemma S5 for all N ≤ 95 (0 violations); P4
  predicate exactness at {7,13,17,19,37} (0 mismatches); flag⟹cert
  sweeps — 80..95 random 0/64,000, N=81 exhaustive 0/1,929,501, N=95
  exhaustive 0/3,612,280.

## 2. The measurement

| V | sets | same ports | same +D3 | aug ports | aug +D3 | OPEN |
|---|---|---|---|---|---|---|
| 16 | 4,311 | 96.94% | 97.26% | — | — | 0 |
| 24 | 41,656 | 85.89% | 87.29% | 92.69% | 93.44% | 0 |
| 32 | 196,751 | 79.52% | 82.13% | 88.10% | 89.64% | 0 |
| 40 | 641,166 | 74.93% | 78.03% | 82.13% | 84.29% | 0 |
| 48 | 1,665,356 | **72.10%** | **75.25%** | **77.79%** | **80.21%** | 0 |

**Decision (reviewer's rule): augmented +D3 at V=48 = 80.21% > 80% — the
decay is slower than extrapolated; the paper's framing can be confident.**
Honest margins: it cleared the bar by 0.21 points; augmented ports is
77.79% (just below 78); the same-battery series continues its decelerating
decay (−6.37, −4.59, −2.83 pts per 8 V-units on ports). The augmented +D3
decay per 8 V-units is 3.80 → 5.36 → **4.07** — the 40→48 interval
decelerated back below the 32→40 interval, contradicting the "mildly
accelerating" extrapolation. Under the last interval's rate the 60%
crossing moves from V ≈ 70 to V ≈ 88; the series is not fit by a single
trend, and the paper should say so rather than fit one.

- V-tier (augmented) share: 6.6% → 10.4% → 15.7% → 19.8% (322,105
  V1-closed + 7,394 V2-only at V=48). Growth continues, slower than
  linear.
- Residual decomposition at V=48 (augmented): 682,490 unflagged
  certifying pairs = 532,700 kernel at N ≤ 47 + 137,484 kernel at
  composites 48–79 + 12,306 kernel at composites 80–95. **The residual is
  still purely the kernel world** (mixed order-signatures at composite
  moduli); zero unflagged all-invertible pairs anywhere.
- Ground truth: (TAU-5) holds on all 1,665,356 sets at v ≤ 48 (every set
  has a certifying pair); the V2 empirical law (every set has a
  certifying pair at a never-covering (N, order-signature)) survives the
  386× extension of the original matched-range base.

## 3. The rigid zone does not grow in 80–95 either

No capable modulus exists anywhere in 80..95. New odd composites, with the
reachability structure (T = 6):

| N | h | non-unit classes | reachable | slots/set | capable | # coverings |
|---|---|---|---|---|---|---|
| 81 | 13 | 13 | 12 | 9+4 | no | 0 |
| 85 | 14 | 10 | 8 | 12+2 | no | 0 |
| 87 | 14 | 15 | 14 | 10+4 | no | 0 |
| 91 | 15 | 9 | **9** | 12+3 | no | 0 |
| 93 | 15 | 16 | 15 | 10+5 | no | 0 |
| 95 | 15 | 11 | 9 | 12+3 | no | 0 |

Primes 83 and 89: safe (inside the committed (47,150] scan; re-derived
and set-based cross-checked as controls). All counts cross-checked by the
independent set-based scanner (exact agreement at all eight moduli).

**The open composite micro-cases are now {49, 77, 91} = {7², 7·11, 7·13}**
— the moduli where ALL non-unit classes are reachable by unit bad sets,
no covering exists, and the obstruction is genuinely multiplicative. The
pattern: the three honest micro-cases are exactly the products of 7 with
the n = 4 rigid-zone primes {7, 11, 13}. Natural predicted next members
(7·17 = 119, 7·19 = 133, 7·37 = 259) lie outside the scanned composite
range (≤ 95). This upgrades the paper's composite open-problem entry from
"micro-cases 49/77" to "the 7·{7,11,13} family, verified safe to 95".

Verified-range statement for the paper: *within the verified range —
primes ≤ 150, composites ≤ 95 — the n = 5 unit-covering-capable moduli
are exactly {7, 13, 17, 19, 37}; beyond that range the question is open.*

## 4. Verdict

The reviewer's optional data point was taken (compute allowed: 538.8s).
Result: **augmented +D3 at V=48 = 80.21%, above the 80% bar — confident
framing for the paper**, with the razor-thin margin and the non-monotonic
decay rate reported as findings, not smoothed away. The rigid zone is
empirically constant over the entire verified range (primes ≤ 150,
composites ≤ 95). The residual remains one named class (the kernel
world). The paper is now written with this as the final data point.
