# Assessment of the Pair-Collapse Reduction: Novelty, Generalization, Recovery

**Session:** LRC residue-language program. **Charge (reviewer):** before any
plain-regime chase, settle three questions — (1) is the Reduction Theorem /
(TAU) novel? (2) does the method generalize to n = 5? (3) can the classical
n = 4 proof be recovered through it? This note reports all three, plus one
result that appeared during the diagnostics and changes the framing: the
reduction is **lossless** — an assembled exact equality
M(V) = max over pairs of tau, making (TAU-n) *equivalent* to LRC-n.

**Committed scripts:** `scripts/lit_search_lrc.sh`, `scripts/lit_search_lrc2.sh`
(+ `lit/` raw JSON, `lit/digest.txt`), `scripts/E2_generalize_n5.py`,
`scripts/E2b_wider_base.py`, `scripts/E3_recovery_diagnostic.py`,
`scripts/E4_cross_check.py`, `scripts/E5_tight_check.py`,
`scripts/E6_two_sided_binding.py`.

---

## 1. Step 1 — Literature verdict

Searched (22 queries): complementary/pair-sum speeds, mirror reductions,
finite cyclic formulations, Cusick view obstruction, BHK, Barajas–Serra,
tight sets, spectra, small-case proofs, induction routes, recent arXiv.

**(a) The mirror technique is folklore.** The pair-collapse identity
(v_q ≡ −v_p on Z_{v_p+v_q}) is elementary, and the known breakpoint fact —
Tao's 2015 blog remark states it in passing: *the optimal time may be assumed
rational, with denominator being the sum or difference of two speeds* — is
exactly the lattice structure the reduction exploits. Small-case proofs have
long used pair-sum type times. The technique must be treated as known.

**(b) The formulation appears new as stated.** No source found states the
constructive per-pair bound M(V) ≥ tau(p,q) with the (n−1)-effective-runner
reading on Z_N, the (TAU-n) statement, or the empirical pair-selection laws
(max-speed pairs suffice; which pair carries the certificate). The closest
active lines — tight-set classification (Kravitz; Giri–Kravitz; Tin Fan's
spectrum amendments), Tao's compactness reduction (2018), the
Rosenfeld/Trakulthongchai/Sungkawichai finite-checking program — study
different statements. Verdict: **technique known, formulation new; treat the
reduction as a reformulation whose value is structural, and cite the
folklore honestly.**

**(c) The frontier has moved — framing update.** LRC is now proved through
*ten* runners (Rosenfeld 2025, arXiv:2509.14111 / 2511.22427;
Trakulthongchai, EJC 2026), with computer-assisted claims to thirteen–fifteen
(arXiv:2609.02604; Quanta, Mar 2026). Our session summary's "known ≤ 7" was
stale. Consequence: an n = 4 result has no rung value; the program's only
possible value is the lossless finite reformulation and its obstruction
geometry. (Also noted: the *shifted* LRC was disproved 2025–26; the
loneliness-spectrum conjecture was disproved at n = 4 by Tin Fan — the
"empirical laws in this family can fail" reminder is live.)

## 2. Step 2 — Generalization (the decisive experiment)

(TAU-n): some pair has (n+1)·tau_{n−1}(v_p+v_q) ≥ v_p+v_q. Exact integer
arithmetic, gcd-1 sets, conventions identical to the committed k = 3 scripts.

| n | target | bound | sets | any pair | max-speed pair | top pair |
|---|--------|-------|------|----------|----------------|----------|
| 3 | 1/4 | v ≤ 30 | 3,472 | **3,472/3,472** | 3,472/3,472 | 3,471/3,472 |
| 4 | 1/5 | v ≤ 16 | 1,745 | **1,745/1,745** | 1,745/1,745 | 1,737/1,745 |
| 5 | 1/6 | v ≤ 24 | 41,656 | **41,656/41,656** | 41,655/41,656 | 41,480/41,656 |
| 6 | 1/7 | v ≤ 18 | 18,479 | **18,479/18,479** | 18,477/18,479 | 18,152/18,479 |

Zero failures anywhere: **the method generalizes; it is a tool, not an n = 4
coincidence.** Additional structure:

- The regular sets are exactly tight at every n — (1,2,3), (1,2,3,4),
  (1,2,3,4,5), (1,…,6) all have best-pair ratio exactly 1 — plus non-regular
  tight sets (1,3,4,7) at n = 4 and (1,3,4,5,9) at n = 5, matching the
  tight-set literature's examples.
- The **pair-selection rule is robust but not rigid**: exactly one n = 5 set
  in v ≤ 24 is certified by no max-speed pair — V = (2,6,8,10,11), certified
  by pairs (2,10) and (8,10) — hand-verified (pair (2,10): N = 12, k = 5,
  residuals (2,6,4,5), 6·2 = 12). Two such sets at n = 6. The (min,max) pair
  degrades fast (98.0% at n = 5, 92.5% at n = 6): selection is genuinely
  set-dependent beyond n = 4.
- **Iteration is available but not universal**: 57% of n = 5 sets contain two
  disjoint pairs of equal sum (both mirrored on a common grid); that double
  collapse certifies 96% of them. The rest certify through plain 4-effective
  grid values — the collapse-to-3 induction does *not* close n = 5 on its own.
- The j-gon (regular-time) share of certificates *shrinks* with n (29% at
  n = 5, 35% at n = 6, "plain" 73% at n = 5 v ≤ 24): case-chasing does not
  scale; only a structural proof of the plain regime would.
- Cross-validation (E4): for 150 random n = 5 and 60 random n = 6 sets, M(V)
  recomputed by exact unrestricted breakpoint arithmetic — Theorem 1 held
  with zero violations, and, unexpectedly, **equality in all 210/210 cases**
  (see §4).

## 3. Step 3 — Recovery diagnostic: the classical proof does NOT re-emerge

For every n = 4 set and pair, the unrestricted 3-effective optimum mu was
computed exactly (breakpoint enumeration over d ∈ {2a, 2b, 2c, a±b, a±c,
b±c}), and compared with the grid certificate. Tiers:

- **classical** (grid captures the unrestricted optimum, mu ≥ 1/4):
  628/1745 sets = **36%**;
- **semi** (certificate at an unrestricted breakpoint — mirror/pole/meeting
  structure — with value in [1/5, 1/4)): 1,047/1745 = **60%**;
- **gridonly** (certificate only at non-breakpoint grid times):
  70/1745 = **4%**.

Sanity: 0 violations of grid ≤ unrestricted; 0 violations of mu ≥ 1/4
(classical k = 3 re-confirmed inside the harness). Hand-verified tier
assignments: V = (1,2,3,4) is semi (its certificate t = 1/5 is the mirror
breakpoint of the *other* pair (2,3) — the pentagon structure, i.e., how the
classical regular-set proof actually works); V = (1,3,4,5) is gridonly (both
certifying times k = 4, 5 on Z_9 are non-breakpoints); V = (1,2,3,16) pair
(1,16) is gridonly (only k = 4, 13 on Z_17, neither a breakpoint).

**Verdict: the classical k = 4 theorem is not recovered from classical k = 3
content.** 64% of sets need grid-specific value laws. The reduction is a
genuinely different route — not a repackaging of the classical case analysis.

## 4. The new result: the reformulation is LOSSLESS

E4's equality finding (M = best-pair tau, 210/210) was confirmed on targeted
sets (E5: five tight/exceptional sets, all exact) and explained:

**Theorem (exact grid reformulation).** For every speed set V,
    M(V) = max over pairs (p,q) of tau_{n−1}(v_p+v_q)/(v_p+v_q).

*Assembly (n = 4, from last turn's proved pieces):* (≥) is Theorem 1. For
(≤), let t* = x/B (reduced) be a global argmax with m/B = M(V) < 1/2.
Two-sided binding (proved last turn) gives binders u, w with residues m and
B−m; the binder pair-sum identity gives B | v_u + v_w; hence t* lies on the
pair-sum lattice of (u, w), where the pair mirrors, so the (u,w)-grid
attains M(V): tau(u,w) ≥ M(V). If M(V) = 1/2, every pair's lattice contains
the pole time. ∎ (Hypotheses of the two-sided-binding proof re-checked; its
continuity argument concerns the active set at the argmax and appears
dimension-free.)

**n-genericity of the missing ingredient (E6):** two-sided binding at *every*
exact argmax breakpoint held in 480/480 random sets at n = 4, 5, 6, with
zero single-binder argmaxes. The equality is thus empirically a law of the
problem, not an n = 4 artifact.

**Consequences.**
1. **(TAU-n) ⟺ LRC-n.** The reduction is not a relaxation or a translation
   of the classical proof — it is a lossless re-encoding of LRC into one
   finite cyclic problem per pair: "some pair-sum grid Z_{v_p+v_q} carries an
   (n−1)-effective configuration at ratio ≥ 1/(n+1)". The reviewer's "LRC in
   disguise" worry is confirmed in the strongest, benign sense: the disguise
   is exact, and the plain regime *is* the whole difficulty of LRC-n, now in
   finite form.
2. The equality also upgrades the empirical laws: "max-speed pairs suffice"
   now reads "the binding pair at the argmax can be chosen to contain the
   maximum speed" — a structural statement about argmaxes.
3. Novelty positioning sharpens: the *equality* is adjacent to the folklore
   "optimal denominator = sum or difference of two speeds" (Tao's remark),
   but the specific statement — a *pair-sum* lattice always suffices, via
   the two-sided binder pair — was not found stated anywhere.

## 5. Updated honest status

| Claim | Status |
|---|---|
| Pair-collapse identity / Theorem 1 (M ≥ tau, every pair) | Proved (elementary; technique folklore). |
| **M(V) = max_pair tau (lossless reformulation)** | **Assembled at n = 4 from proved pieces; two-sided binding verified n-generically 480/480; needs the continuity argument re-written for general n.** |
| (TAU-n) ⟺ LRC-n | Follows from the equality. |
| (TAU-n) empirics | 0 failures: 3,472 + 1,745 + 41,656 + 18,479 sets (n = 3..6). |
| Max-speed-pair selection rule | Holds 99.998% (n = 5, v ≤ 24); exactly 1 + 2 known counterexamples; set-dependent in general. |
| j-gon / pentagon cases of (TAU) | Proved (last turn). Plain share grows with n (73% at n = 5). |
| Recovery of classical k = 4 from k = 3 content | **No** — 36% / 60% / 4% tier split. |
| Plain regime | = LRC-n in finite form; the single remaining content, at every n. |

## 6. Verdict and recommendation

Per the reviewer's decision rule: the reduction generalizes (losslessly) and
its formulation is new-as-stated, so this is paper-shaped — but with the
field frontier at 10–15 runners, only as a **methodology/structure**
contribution: the lossless finite reformulation (equality theorem, with the
n-generic two-sided-binding proof written out), the pair-selection laws, and
the tiered obstruction geometry at the base rung.

Recommended next session (in order):
1. **Write the equality theorem properly**: re-prove two-sided binding for
   general n (the one analytic gap; 480/480 empirical), assemble, state
   (TAU-n) ⟺ LRC-n. This is the citable core.
2. Then, and only then, the n = 4 plain-regime chase — now understood as a
   *validation exercise* for the finite method on a known case, not a path
   to new rungs; abandon it if it needs classical machinery (the 4%
   gridonly core is the honest stress test).
3. Do not chase n ≥ 5 by cases (73% plain share). If the plain regime
   resists structure, pivot to R2 with a clean conscience — the lossless
   reformulation is itself the deliverable of this thread.
