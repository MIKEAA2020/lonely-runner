# The Top-Pair Rescue Lemma at k = 3: Proof Attempt in the Residue Language

**Session:** LRC residue-language program, rung k = 3 (n = 4).
**Directive:** focused proof attempt on the top-pair rescue lemma, in isolation,
in the residue language, without importing the classical k = 3 proof. Success
= real theorem; failure with a clear obstruction = pivot to R2.
**Status of this note:** the witness-level lemma is *reduced* to a pure speed
statement (Theorem 1 + Corollary 1, proved); the pentagon case of that
statement is proved (Lemmas J and P); the remaining "plain regime" is a
finite arithmetic chase, empirically verified 1745/1745, not yet closed.

---

## 0. Vocabulary (residue language)

For a speed set V = {v_1, ..., v_n}, a grid B >= 1 and a grid time x:
- residues r_i = v_i x mod B; residuals c_i = ||r_i||_B = min(r_i, B - r_i);
  m = min_i c_i. The witness (x, B) is 1/5-good iff 5m >= B.
- M(V) = sup over real t of min_i ||v_i t|| (max loneliness). LRC at n = 4:
  M(V) >= 1/5.
- Pair-sum grid of a pair (p, q): N = v_p + v_q. Its lattice times k/N.
- m3*(N; w_1, w_2, w_3) := max over k in Z_N of min(||w_1 k||_N, ||w_2 k||_N,
  ||w_3 k||_N) — the per-grid 3-runner maximum.

## 1. The Reduction Theorem (the attempt's main new result)

**Lemma 1 (pair-collapse identity).** On the grid N = v_p + v_q one has
v_q = -v_p (mod N), hence ||v_q k||_N = ||v_p k||_N for every k.

*Proof.* v_q ≡ N - v_p ≡ -v_p (mod N), and the residual ||·||_N is even. ∎

**Theorem 1 (pair-collapse bound).** For every speed set V and every pair
(p, q) with N = v_p + v_q, letting (w_1, ..., w_{n-2}) be the speeds of the
other runners,

    M(V) >= (1/N) · max_k min( ||v_p k||_N, ||w_1 k||_N, ..., ||w_{n-2}||_N ).

In particular at n = 4: M(V) >= max over the 6 pairs of
tau(p,q) := m3*(v_p + v_q; v_p, v_y, v_z) / (v_p + v_q).

*Proof.* For the k achieving the inner maximum, evaluate the loneliness at the
time t = k/N: the n-runner minimum equals the (n-1)-effective minimum by
Lemma 1, i.e. the value m3*/N. M(V) is a supremum over all times. ∎

Constructive: the bound is exhibited by explicit lattice times. Verified
exactly (breakpoint arithmetic, 71/71 random sets; the earlier apparent
violations were a sampling artifact — a coarse grid misses denominators).
The bound is sometimes tight: V = (1, 2, 3, 5) has M = tau = 1/4; the
regular sets V = (1, 2, 3, 4), (1, 3, 4, 7) have M = tau = 1/5 exactly.

**Corollary 1 (k = 3 reduction).** LRC at n = 4 follows from the pure speed
statement

    (TAU)  for every 4-speed set V with gcd(V) = 1 there is a pair (p, q)
          with 5 · m3*(v_p + v_q; v_p, v_y, v_z) >= v_p + v_q.

Moreover it suffices that (TAU) hold for the three pairs containing the
maximum speed (verified: max-speed-pair form holds 1745/1745 for v <= 16).

*Proof.* gcd(V) = 1 is WLOG (scaling all speeds by c scales M and every tau
identically). If (TAU) holds, Theorem 1 gives M(V) >= 1/5. ∎

**Remark (what dissolves).** The witness-level machinery — deficit law,
covering window, two-sided binding, binder pair-sum law, orbit and divisor
lemmas — is *not needed* for this route: the contradiction form of the
classical argument (assume M < 1/5 at the global argmax, derive a beating
time) has been replaced by a direct constructive bound. The empirical law
"top-pair rescue 116/116 at n = 4" is explained by (TAU): at every failing
witness some pair-sum lattice time beats the witness, because some pair-sum
grid carries a 3-effective configuration at ratio >= 1/5, which dominates
every witness with m/B < 1/5 (rescue ⟺ m/B < tau(p,q)).

## 2. Proved cases of (TAU)

**Lemma J (j-gon times).** Let j >= 2, j | N, and let the effective speeds
satisfy w ≢ 0 (mod j) for all three. Then k = N/j achieves every effective
residual >= N/j; hence tau >= 1/j.

*Proof.* w·(N/j) mod N = (N/j)·(w mod j) since (N/j)·j = N ≡ 0 (mod N);
the residual is (N/j)·min(w mod j, j - w mod j) >= N/j. ∎

Useful j for the 1/5 target: j ∈ {2, 3, 4, 5}. (j = 2 recovers the classical
"all speeds odd ⟹ t = 1/2" fact; j = 5 is the pentagon.)

**Lemma P (pentagon case of TAU).** If 5 ∤ v_i for all four speeds and the
residues v_i mod 5 contain an opposite pair (r, -r), then (TAU) holds.

*Proof.* Take the opposite pair (p, q): v_p + v_q ≡ 0 (mod 5), so 5 | N, and
all effective speeds (pair member + the other two) are ≢ 0 mod 5; Lemma J
with j = 5 gives tau >= 1/5. ∎

Caveat (honest): the pigeonhole gives only same-±-class pairs; a residue
multiset contained in a "half-set" {1,2}, {1,3}, {4,2}, {4,3} has no opposite
pair and the pentagon route is unavailable — e.g. V = (1, 6, 11, 16) ≡
(1, 1, 1, 1) mod 5. Those cases are carried by the plain regime (§3).

**Lemma S (size facts for max-speed pairs).** For the pair (v_4, v_x):
(a) no effective speed is ≡ 0 mod N (stuck is impossible: every other speed
is < v_4 < N, and v_4, v_x < N);
(b) for the top-speed pair (v_3, v_4) no effective speed equals N/2, since
(v_3 + v_4)/2 lies strictly between v_3 and v_4;
(c) if m N / B < 1 (small grids relative to the witness), the rescue needs
only a pole-free time, which exists unless all three effective speeds are
stuck — impossible by (a).

**Lemma M (moat counting, sufficient).** If N(B - 6m) > B·G, where
G = sum of gcd(effective speed, N), the pair rescues the witness. (Three
moats of radius R-1 cannot cover Z_N.) Covers deep bridges; provably
insufficient in the shallow regime (B - 6m = δ - m <= 0).

**Lemma 2 (two-runner base, inclusion–exclusion).** Speeds α, β on grid N:
if N >= 4R - 3 + (g_α + g_β) then some k has both residuals >= R. (Safe-set
sum >= N + 1.) Part of the descent toolkit; not on the critical path of
Corollary 1.

## 3. The plain regime (open)

(TAU) remains open when the j-gon routes are blocked:
- (B) all speeds ≢ 0 mod 5 but residues lie in a half-set (no opposite pair);
- (C) some speed ≡ 0 mod 5.

Empirics (committed scripts, v <= 16, 1745 sets): (TAU) holds always; the
max-speed-pair form holds always; the specific top-speed pair fails 8/1745.
Obstruction taxonomy of *failing* pairs (witness-level data, 824 two-sided
points): stuck (N | v_y), half (v_y = N/2), third (v_y ∈ {N/3, 2N/3}),
5blocked (5 | N with an effective speed ≡ 0 mod 5), covering. The covering
class is confined, in-range, to N ∈ {7, 11, 13} with all three effective
speeds coprime to N (G = 3) — the "rigid zone". Winning times are ~2/3
(near-)pentagon-aligned.

Remaining obligation, precisely: for every V in regimes (B)/(C), exhibit a
max-speed pair (v_4, v_x) and prove m3*(v_4 + v_x) >= ceil((v_4 + v_x)/5)
in the residue language. The obstructions to rule out are exactly the
taxonomy above; stuck and half are already excluded by Lemma S for the right
pairs; what remains is the finite chase over third/5blocked/covering for the
three max-speed pairs simultaneously — a bounded divisibility problem.

## 4. Verdict

- The reviewer's question — "is the top-pair rescue lemma provable in the
  residue language alone?" — has changed shape: the lemma is now the pure
  speed statement (TAU) behind a proved reduction (Theorem 1). The
  witness-level empirics are explained, not just confirmed.
- Translation-detector audit: Lemma J/P (regular times) *is* the classical
  easy direction rediscovered in-language (flagged, not hidden). The novel
  content is Theorem 1 (pair-collapse constructive bound, tight at regular
  sets) and the plain-regime chase, which has no classical counterpart we
  know of. The residual risk that the plain chase secretly needs the
  classical regularity induction is real and is the thing to watch.
- Recommendation: do not pivot to R2 yet. The gap is now a finite arithmetic
  chase with a mapped obstruction taxonomy — qualitatively different from
  the open-ended lemma the turn started with. One more focused session on
  regimes (B)/(C) is warranted; R2's group dictionary may even assist
  (unit arithmetic |moat ∩ μ·moat| is character-theoretic in flavor).

## 5. Committed artifacts

- scripts/verify_top_pair_rescue_k3.py (E0–E7: witness-level verification)
- scripts/classify_rescue_k3.py (T1–T6: taxonomy; fixed a bitwise-OR bug
  caught by hand-verification)
- scripts/verify_tau_statement_k3.py (TAU verification, witness-free)
- this note: download/k3_top_pair_rescue_proof_attempt.md

Key numbers: (TAU) 1745/1745; max-speed-pair form 1745/1745; top-speed pair
1737/1745; witness-level mixed-pair rescue 824/824; pentagon lemma 586/586;
exact Reduction cross-check 71/71; covering class confined to N ∈ {7, 11, 13}.
