r"""Companion v2 part 2: sec1 abstract, secscale, sec2, sec7, sec9, bib."""
P = []

# ---- sec1 abstract: three hypotheses; megasentence split (gemini)
P.append(('sec1.tex',
r"""The paper's culminating result is stated honestly as a \emph{conditional}
theorem (Theorem~\ref{thm:SC}): under two named open hypotheses---the
distinct-family volume law H1g and the dilate-sparsity pair H2m/H2r---the
fiber cascade closes the covering obstruction at the frontier cells
$(1K,T{-}3U)$ of every rung $T\ge 8$ at $N=p^2$ within a depth that is
$T$-flat in practice and never exceeds $T+3$ in the worst case; the
reduction chain and every other input (normal form, reflection fiber,
chain-class theorem, renormalization lemma, lossless equivalence) are
proved or machine-validated, and each constant's provenance---derived,
measured, or fitted, with the data it is fit to---is disclosed in a
six-entry hypothesis ledger.
That boundary between the proved and the hypothesized has a structural
type, verified statement-by-statement in \S\ref{sec:scale-status}: every
proved result in this program is a closure, invariance, exhaustiveness,
or machine-verification statement, while every remaining obligation is a
uniform growth-in-$k$ margin over a trivial or measured baseline---a
statement type no turn of the program has proved. This type mismatch is
why the paper stops where it stops. The one candidate instrument
shaped to produce that missing type---a transfer-matrix or automaton
model of the strand census---was then executed as a proof attempt under
a pre-committed success criterion and closed by theorem
(Lemma~\ref{lem:monomial}): the cascade's constrained transfer
matrices are partial permutations with spectral radius exactly $0$ or
$1$ at every scale, the cyclic product is diagonal with trace equal to
the survivor count itself, and the missing margin is provably external
to the lossless dynamics---the reformulation is bounded by its own
fidelity, the losslessness that makes it a faithful rewrite being the
property that withholds the decay a margin would require.""",
r"""The paper's culminating result is stated as a \emph{conditional}
theorem (Theorem~\ref{thm:SC}): under three named open hypotheses---the
multi-fiber survivor inequality Hq (quantified in this revision; it
subsumes the former dilate-sparsity input H2m as its one-step case),
the distinct-family volume law H1g with its uniform margin, and the
resonance budget H2r---the fiber cascade closes the covering
obstruction at the frontier cells $(1K,T{-}3U)$ of every rung $T\ge 8$
at $N=p^2$ within a depth that is $T$-flat in practice and never
exceeds $T+3$ in the worst case. The reduction chain and every other
input (normal form, reflection fiber, chain-class theorem,
renormalization lemma, lossless equivalence) are proved or
machine-validated, and each constant's provenance---derived, measured,
or fitted, with the data it is fit to---is disclosed in the hypothesis
ledger.
That boundary between the proved and the hypothesized has a structural
type, verified statement-by-statement in \S\ref{sec:scale-status}: every
proved result in this program is a closure, invariance, exhaustiveness,
or machine-verification statement, while every remaining obligation is a
uniform growth-in-$k$ margin over a trivial or measured baseline---a
statement type no turn of the program has proved. This type mismatch is
why the paper stops where it stops. The one candidate instrument
shaped to produce that missing type was a transfer-matrix or automaton
model of the strand census. It was executed as a proof attempt under a
pre-committed success criterion and closed by theorem
(Lemma~\ref{lem:monomial}). The cascade's constrained transfer
matrices are partial permutations with a spectral radius of exactly $0$
or $1$ at every scale. The cyclic product is diagonal, with trace equal
to the survivor count itself, and the missing margin is provably
external to the lossless dynamics. The reformulation is bounded by its
own fidelity: the losslessness that makes it a faithful rewrite is the
property that withholds the decay a margin would require."""))

# ---- sec1: "closed by theorem" narrowed + revision-history sentence
P.append(('sec1.tex',
r"""matrices of the lossless cascade are partial permutations
(Lemma~\ref{lem:monomial}), every spectral quantity in the framework
is exactly $0$ or $1$, the cyclic product is diagonal---its trace is
the cascade's own survivor count---and the survival-fraction deficit
the criterion demanded is, in this language, precisely the pair bound
H2m and the volume law H1g restated. The instrument is the census in
different notation; it consumes the hypotheses it was commissioned to
prove. That is the mechanism behind the type mismatch: losslessness
forbids generated decay, and the reformulation is bounded by its own
fidelity.""",
r"""matrices of the lossless cascade are partial permutations
(Lemma~\ref{lem:monomial}), every spectral quantity in the framework
is exactly $0$ or $1$, the cyclic product is diagonal---its trace is
the cascade's own survivor count---and the survival-fraction deficit
the criterion demanded is, in this language, precisely the survivor
inequality Hq and the volume law H1g restated. The instrument is the
census in different notation; it consumes the hypotheses it was
commissioned to prove. That is the mechanism behind the type mismatch:
losslessness forbids generated decay, and the reformulation is bounded
by its own fidelity. (The closure theorem covers automata built on the
cascade's bijective drift---the lossless state space; instruments that
do not inherit the drift are outside its scope, as
\S\ref{sec:scale-status} records.)"""))

# ---- sec1 notation: h added; exact-integer claim fixed
P.append(('sec1.tex',
r"""Throughout, $n\ge 2$ is the number of speeds, $T=n+1$ is the loneliness
target denominator, and $V$ has $\gcd(V)=1$ (dividing by the gcd does
not change $M$). For a modulus $N$ we write $\wn{x}=\min(x\bmod N,\;
N-(x\bmod N))$ and $\ord_N(w)=N/\gcd(w,N)$ for the additive order of $w$
in $\Z_N$. The letter $k$ always denotes a lattice time index. All
computations reported in this paper are exact integer computations; their
verification record is the subject of \S\ref{sec:sound}.""",
r"""Throughout, $n\ge 2$ is the number of speeds, $T=n+1$ is the loneliness
target denominator, and $V$ has $\gcd(V)=1$ (dividing by the gcd does
not change $M$). For a modulus $N$ we write $\wn{x}=\min(x\bmod N,\;
N-(x\bmod N))$, $\ord_N(w)=N/\gcd(w,N)$ for the additive order of $w$
in $\Z_N$, and $h=\floor{(N-1)/T}$ for the covering threshold index
of \S\ref{sec:cover}. The letter $k$ always denotes a lattice time
index. Computations are exact integer/rational arithmetic wherever a
theorem or table entry asserts an exact value; the Wilson intervals,
the fitted $\Lambda_0$, the projected pools of Table~\ref{tab:pool},
and displayed asymptotics are statistical or asymptotic estimates and
are labeled as such in the free-parameter disclosure. The verification
record is the subject of \S\ref{sec:sound}."""))

# ---- secscale: prop:depth relabeled as the independent-extension model
P.append(('secscale.tex',
r"""\begin{proposition}[comparison; under H0, H1g, H2m, H2r]
\label{prop:depth}
In the independent-extension composition the candidate count after $i$
processed fibers is
$C_i\approx\bar\gamma^{\,i}\,|\mathrm{Sol}|^{i+1}/p^{(m-1)i}$, and the
family dies within""",
r"""\begin{proposition}[comparison; the independent-extension model,
evaluated under the calibrated constants]
\label{prop:depth}
\emph{Model statement, not a cascade theorem}: in the
independent-extension composition---the idealization in which the
per-fiber cuts act on independent copies rather than on the actual
history-dependent survivor sets---the candidate count after $i$
processed fibers is
$C_i\approx\bar\gamma^{\,i}\,|\mathrm{Sol}|^{i+1}/p^{(m-1)i}$. The
actual-cascade statement is Theorem~\ref{thm:SC}'s clause (i) with
Proposition~\ref{prop:Hq}'s induction, conditional on the quantified
ledger; the model below is its calibration. In the model, the family
dies within"""))

# ---- secscale: honest-tone reduction (load-bearing kept)
P.append(('secscale.tex',
r"""The calibration is honest about its provenance: the per-fiber cut""",
r"""The calibration carries its provenance explicitly: the per-fiber cut"""))

# ---- secscale: the chain-volume clarification
P.append(('secscale.tex',
r"""which is exactly H1g's content.""",
r"""which is exactly H1g's content; the $k$-growing $s_{\max}$ factors of
the proved chain-class bound (Theorem~\ref{thm:chaindilate}) are loose
by orders of magnitude, and the $\approx\mathrm{ov}^{m-1}$ shape with
an $O(1)$ constant remains empirical---H1g's open content---not a
consequence of the chain theorem."""))
