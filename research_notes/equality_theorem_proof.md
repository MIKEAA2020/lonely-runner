# The Pair-Sum Equality Theorem: Complete Proof and Novelty Status

**Program:** LRC residue-language reformulation. **This note:** the reviewer's
directive 2 — write out the equality theorem cleanly for general n, with the
case analysis (opposite-slope binding vs. peak binding) in place of the
empirical two-sided-binding ingredient (E6, 480/480). Companion verification:
`scripts/V1_equality_fresh_check.py` (fresh, non-circular). Novelty: §5.

---

## 0. Conventions and statement

Let k >= 2 and let V = (v_1, ..., v_k) be positive integers (not necessarily
distinct — the argument never uses distinctness; LRC adds it). Work on the
circle T = R/Z. For t in T write

  f_p(t) = || v_p t ||  (distance of v_p t to the nearest integer, in [0, 1/2]),
  g(t)   = min_{1<=p<=k} f_p(t),
  M(V)   = max_{t in T} g(t).

g is continuous (min of finitely many continuous functions) on the compact
circle, so the maximum is attained. Dividing all speeds by d = gcd(V) does not
change the problem (t -> dt is a bijection of T and f_p(t) = f_{p/d}(dt)), so
we assume gcd(V) = 1 throughout. Note M(V) > 0: at t = 1/(2 v_max) every f_p
equals v_p/(2 v_max) <= 1/2, so g >= v_min/(2 v_max) > 0.

**Pairs and grids.** For an unordered pair {u, w}, u != w, of indices, let

  N = N_{uw} := v_u + v_w  (>= 3),  the *pair sum*.

The *pair-sum lattice* of {u,w} is the set of times { k/N : 1 <= k <= N-1 }.

**The reduced problem on a pair grid.** For a pair {u, w} define

  tau(u,w) := max_{1 <= k <= N-1}  min_{p != w}  || v_p k / N ||.

This is the (k-1)-runner maximum-loneliness of the *effective speeds*
(v_p mod N)_{p != w} on the cyclic group Z_N — the runner w has been deleted,
and by the mirror identity below the remaining runner u carries its partner.

**THEOREM (pair-sum equality; lossless reformulation).** For all V as above,

  M(V) = max over pairs {u,w} of tau(u,w).

Moreover, the maximum defining M(V) is attained at a time t* = k*/(v_u + v_w)
lying on the pair-sum lattice of some pair; and if M(V) < 1/2 the binders at
t* can be chosen two-sided (one at residue m, one at residue -m).

**COROLLARY.** LRC for k runners (target 1/(k+1)) is equivalent to the finite
statement (TAU-k): for every gcd-1 speed set V there is a pair {u,w} with

  (k+1) · max_{1<=k<=N-1} min_{p != w} || v_p k/N ||  >=  v_u + v_w .

---

## 1. Lemma 1 (mirror / pair collapse)

**Lemma 1.** Let N = v_u + v_w. For every integer k,
  (i)  v_u k/N + v_w k/N = k in Z, hence v_u k ≡ -v_w k (mod N) and
       || v_u k/N || = || v_w k/N ||;
  (ii) min_{p in V} || v_p k/N || = min_{p != w} || v_p k/N ||.

*Proof.* (i) is immediate: the two positions sum to the integer k, so they are
reflections of each other about 0 (and distances to 0 coincide). For (ii): the
left side is the minimum over a superset, so it is <= the right side. For >= :
the left minimum equals min( min_{p not in {u,w}} f_p , f_u , f_w ), and since
f_u = f_w by (i) this equals min( min_{p not in {u,w}} f_p , f_u ), which is
exactly the right side (the set {p != w} consists of {p not in {u,w}} plus u).
∎

Part (ii) says: **at every time of the pair-sum lattice, the full k-runner
minimum and the (k-1)-runner reduced minimum coincide.** Deleting w costs
nothing exactly at the times when w is mirrored by u.

---

## 2. Lemma 2 (attainment at a pair-sum time)

**Lemma 2.** There exists a pair {u,w} and an integer k* with 1 <= k* <= N-1
such that g(k*/N) = M(V). Specifically, let t* be any global argmax (g(t*) =
M(V) =: m) and define the binders B = { p : f_p(t*) = m }. Then:

  (a) if 0 < m < 1/2: some binder is RISING (v_p t* mod 1 in (0,1/2)) and
      some binder is FALLING (v_p t* mod 1 in (1/2,1)); picking u rising and
      w falling, (v_u + v_w) t* is an integer, and t* = k*/N for that pair;
  (b) if m = 1/2: all speeds are odd, t* = 1/2, and every pair has even N and
      k* = N/2;
  (c) m = 0 cannot occur.

*Proof.* (c): M(V) > 0 as noted in §0.

Setup for (a): assume 0 < m < 1/2. Every binder p in B has residue
v_p t* mod 1 equal to m (call p *rising*: on the segment (0,1/2) the position
x = v_p t mod 1 satisfies f_p = x, which is increasing at slope +v_p) or to
1-m (call p *falling*: x in (1/2,1), f_p = 1-x, slope -v_p). No binder is at
a peak (f = 1/2) or a valley (f = 0) because m is neither 1/2 nor 0.

**Claim: not all binders rise; not all binders fall.**
Suppose all binders rise (the other case is symmetric). For p not in B let
delta_p = f_p(t*) - m > 0. All f_p are v_p-Lipschitz. Choose s > 0 with

  s < min( min_{p not in B} delta_p / v_p ,  min_{p in B} (1/2 - m)/(2 v_p) ).

For binders (all rising): for 0 < s' <= s the position stays in (0,1/2), so
f_p(t* + s') = m + v_p s' > m. For non-binders: f_p(t* + s') >= f_p(t*) -
v_p s' = m + (delta_p - v_p s') > m by the choice of s. Hence
g(t* + s') > m = M(V), contradicting maximality of M(V). If all binders fall,
perturb backwards (s' < 0) — symmetric. This proves the claim, so (a) holds
with u rising, w falling in B: v_u t* ≡ m and v_w t* ≡ -m (mod 1), hence
(v_u + v_w) t* ≡ m - m ≡ 0 (mod 1). So k* := (v_u + v_w) t* is an integer;
since t* in T \ {0} (g(0) = 0 < m) we may represent t* in (0,1), whence
0 < k* < N, i.e. 1 <= k* <= N-1. By construction g(t*) = M(V).

(b): if m = 1/2 then every runner satisfies f_p(t*) = 1/2 (the minimum is 1/2
and no f_p exceeds 1/2), i.e. v_p t* ≡ 1/2 (mod 1) for all p. Write t* = a/b
in lowest terms. Then (v_p - v_1) a/b in Z for all p, so b | v_p - v_1 (as
gcd(a,b) = 1); and 2 v_1 a = b(2j+1) for some integer j, so b | 2 v_1.
If b is odd, b | v_1, hence b | v_p for all p, so b | gcd(V) = 1 and b = 1 —
but then t* = 0, excluded. If b = 2c is even, then c | v_1 and
v_p ≡ v_1 (mod 2c) implies v_p ≡ 0 (mod c) for all p, so c | gcd(V) = 1,
c = 1, b = 2. So t* = 1/2, and v_1/2 ≡ 1/2 (mod 1) forces v_1 odd; then
v_p ≡ v_1 (mod 2) makes every speed odd. Conversely (used in §3) if all
speeds are odd then t = 1/2 gives f_p = 1/2 for every p, so M(V) = 1/2.
Finally, every pair sum N = v_u + v_w of two odd speeds is even; distinctness
gives v_u + v_w >= 1 + 3 = 4, so k* = N/2 in {2, ..., N-1}; wait — even
without distinctness, N >= 1+1 = 2 and k* = N/2 in {1, ..., N-1} unless
N = 2, i.e. v_u = v_w = 1, in which case take any other pair if k >= 3 (some
pair of distinct indices with N >= 4 exists when k >= 2 and not all speeds
equal 1; if all speeds equal 1 the normalization gcd = 1 already holds and
M = 1/2 with t* = 1/2 = k*/N for k* = 1, N = 2). In every case t* = 1/2 =
(N/2)/N is a pair-lattice time. ∎

*Remark (what forced two-sidedness means).* The perturbation argument is the
rigorous form of the empirical E6 observation (two-sided binding at every
exact argmax, 480/480). It also shows the folklore's "difference
denominators" |v_u - v_w| never have to carry the optimum: an argmax time
always lies on a SUM lattice, because same-sided binder pairs (which produce
difference relations) can never be the only binders. V1's check C5 confirms:
in 885/885 sets, every unrestricted argmax lies on a pair-sum lattice.

---

## 3. Proof of the Theorem

**(<= direction: M >= tau for every pair — the constructive bound.)**
For every pair {u,w} and every k in {1,...,N-1}, Lemma 1(ii) gives
min_{p != w} ||v_p k/N|| = g(k/N) <= M(V). Taking the max over k:
tau(u,w) <= M(V). Hence max_pair tau <= M(V).

**(=> direction: some pair attains M.)**
By Lemma 2 there is a pair {u,w} with k* in {1,...,N-1} and
g(k*/N) = M(V). By Lemma 1(ii) again,

  min_{p != w} || v_p k*/N || = g(k*/N) = M(V),

so tau(u,w) >= M(V). (In case (b) — all speeds odd — this works for EVERY
pair, with k* = N/2: all runners sit at 1/2, the mirror is trivially exact,
and the reduced minimum is 1/2 = M(V).)

Combining: max_pair tau(u,w) = M(V), with equality attained at the pair {u,w}
produced by Lemma 2. ∎

**Corollary (equivalence with LRC).** LRC-k asserts M(V) >= 1/(k+1) for all
gcd-1 V. By the Theorem this is equivalent to: some pair {u,w} has
tau(u,w) >= 1/(k+1), i.e. (k+1) · max_k min_{p != w} ||v_p k/N|| >= N with
N = v_u + v_w. That is (TAU-k). ∎

---

## 4. Verification record (fresh, non-circular)

`scripts/V1_equality_fresh_check.py` computes M two independent ways:
M_unrestricted from the folklore breakpoint set
D = {2 v_p} U {v_p + v_q} U {|v_p - v_q|} (exhaustive for vertices of the PL
lower envelope — no pair-sum assumption), and M_pairsum from pair-sum
lattices only. Exact integer arithmetic throughout.

| corpus | C1 equality | C2 argmax on sum lattice | C3 all-odd <=> M=1/2 | C4 two-sided at sum argmax |
|---|---|---|---|---|
| n=4, v<=12, gcd 1 (479 sets) | 479/479 | 479/479 | 479/479 | 464/464 |
| n=5, v<=20, random 400 (393 gcd-1) | 393/393 | 393/393 | 393/393 | 389/389 |
| 13 targeted edge sets (n=2, all-odd, tight, exception (2,6,8,10,11)) | 13/13 | 13/13 | 13/13 | 9/9 |

C5 (structural sharpening): **0 sets in 885 have a difference-only argmax** —
every unrestricted argmax lies on a pair-sum lattice, exactly as Lemma 2
predicts (the folklore's "sum or difference" is sharpened to "sum").
Combined with last turn's E4/E5 (215/215 + tight sets at n=5,6) and E6
(480/480), the equality is verified on >1,100 distinct sets by three
independent computations.

---

## 5. Novelty status after literature round 2

(36 web queries total across both rounds + arXiv API enumeration of the full
2024–2026 ring; raw JSON in `scripts/lit/`, digests in `lit/digest.txt`,
`lit/digest2.txt`, `lit/arxiv_recent.xml`.)

**What is known (must be cited as folklore):**
- The analytic core of Lemma 2 — "the maximum is attained at an intersection
  of two sine curves (or a peak), so the optimal time is rational with
  denominator the sum or difference of two speeds" — appears verbatim in a
  comment on Tao's 13 May 2015 blog post ("A remark on the lonely runner
  conjecture"). The technique must be treated as known.
- Exact values ML = p/q and their denominators are an active topic: Kravitz's
  conjecture (q = np+1 when ML < 1/n), Fan–Sun counterexamples and amendment
  (q = np+k), Cordella (arXiv:2609.03444) proving the six-speed case with
  ML = (P-1)/(6P), P ≡ ±1 mod 6. The community computes exact optima on
  infinite families; the breakpoint/lattice structure is certainly in hand
  there, at least implicitly.
- Beck–Everett (arXiv:2609.06259, "Lonely Runner Relations"): counterexamples
  and tight instances lie on hyperplanes m·n = 0 with ||m||_1 <= min(2k+3,
  ((k+1)/(k-1)) flt(k)) — bounded-weight integer relations among speeds,
  via Fourier analysis. Related in spirit (weight-2 relations are pair sums /
  differences) but a different statement: it constrains the speed vector in
  parameter space, not the optimal time.

**What was NOT found anywhere in the reachable literature:**
- The mirror/collapse reading: on Z_{v_u+v_q} the pair is reflected at EVERY
  lattice time, so the k-runner problem on the pair grid IS a (k-1)-runner
  problem (Lemma 1(ii)) — no source states this.
- The exact identity M(V) = max_pair tau (the equality, with the forced
  SUM-only lattice and the all-odd pole case isolated).
- The (TAU-k) equivalent finite statement and the pair-selection laws.

**Honest verdict (unchanged from last turn, now better evidenced).** The
theorem is correct and its proof is elementary; its analytic core is folklore
(blog-comment grade); the specific formulation — mirror grids as (k-1)-runner
problems, the exact equality, (TAU-k) — is new as stated, but with the field
actively computing exact optima (Cordella 2026 computes ML exactly on infinite
families), the absence of the statement from reachable sources is *moderate*
evidence of novelty, not strong. The correct positioning: a reformulation
paper whose citable core is the lossless equivalence, with the folklore
honestly cited, and whose empirical contribution is the pair-selection
structure and the obstruction geometry. The reviewer's caution — "probably a
known result in a different language" — remains the right prior.

**Frontier note (framing).** LRC is now proved through 10 runners
(Rosenfeld 2025; Trakulthongchai, EJC 2026, "Nine and ten lonely runners"),
computer-assisted through 14–15 (Allikvere, arXiv:2609.02604, with public
certificates; Malikiosis–Santos–Schymura finite-checking bound
binom(n+1,2)^{n-1}). An n=4 theorem has no rung value; the equality theorem
matters only as the entry point to the finite cyclic obstruction geometry.
