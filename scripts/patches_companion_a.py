r"""Companion v2 — the Hq repair of Theorem SC (secsc.tex) + abstract (sec1)."""
P = []

# ---- 1. Theorem SC statement: quantified hypotheses, epsilon range, Hq
P.append(('secsc.tex',
r"""\begin{theorem}[scaling closure; conditional]\label{thm:SC}
Assume the hypotheses $\mathrm{NF}$, $\mathrm{H0}$, $\mathrm{H1g}$,
$\mathrm{H2m}$, $\mathrm{H2r}$ of Table~\ref{tab:ledger} (statuses as in
the table), together with the committed equivalence $\mathrm{EQ}$ of
the same table. Let $T\ge 8$, let $p=Tk+\epsilon$ be prime with $k\ge 2$
(so $p\ge p_0(T)$, the least such prime; $p_0(8)=17$), and let the
frontier cell $(1K,T{-}3U)$ at $N=p^2$ be given, with $m=T-3$ unit
arcs, the two-value size window of Lemma~\ref{lem:rungT}, and overlap
budget $\mathrm{ov}=\Sigma s-p\approx(T-6)k+O(1)$. Then:
\smallskip
\emph{(i)} every census-admissible family in the cell dies in the fiber
cascade within
\begin{multline*}
C(T,k)\;=\;
\frac{\log\!\bigl(c_1\,\mathrm{ov}^{\,m-1}\bigr)}
     {(m-1)\,\log\!\bigl(T/(T-6)\bigr)-\log c_1}\\
+\;\frac{2\Lambda_0\,(T-2)}{T}\;+\;1\;+\;O(1)
\end{multline*}
fibers---numerically a $T$-flat plateau $\approx 9$, $9.2$, $10.1$ at
$T=8,9,10$ (at $k=2,3,3$)---with the never-exceeded worst case
\[
(T-4)+N_{\mathrm{res}}+O(1)\;\le\;T+3,
\qquad
N_{\mathrm{res}}=1+\tfrac{2\Lambda_0(T-2)}{T};
\]
\smallskip
\emph{(ii)} since $C(T,k)\ll|U|=(T-2)k$ for every $T\ge 8$, $k\ge 2$,
with a margin growing linearly in $k$, the cascade exhausts the fiber
budget: no census-admissible family admits a covering placement valid on
all of $U$, and hence
\smallskip
\emph{(iii)} the cell $(1K,T{-}3U)$ at $p^2$ is empty; and by the
lossless reformulation equivalence, $\mathrm{LRC}$ holds at that rung.
The measured range of the calibration is $T=8,9,10$ at $p\le 73$, with
committed ground-truth cells at $p\le 61$ for $T\in\{7,8\}$; the
statement at primes beyond the calibration is exactly what H1g and H2m
purchase.
\end{theorem}""",
r"""\begin{theorem}[scaling closure; conditional]\label{thm:SC}
Assume the hypotheses $\mathrm{NF}$, $\mathrm{Hq}$, $\mathrm{H1g}$,
$\mathrm{H2r}$ of Table~\ref{tab:ledger} (statuses as in the table;
$\mathrm{Hq}$ subsumes the former pairwise input $\mathrm{H2m}$ as its
one-step case), together with the committed equivalence $\mathrm{EQ}$
of the same table. Let $T\ge 8$, let $p=Tk+\epsilon$ be prime with
$k\ge 2$ and $0\le\epsilon<T$ (so $k=\floor{p/T}$ is canonical; $p\ge
p_0(T)$, the least such prime, $p_0(8)=17$), and let the frontier cell
$(1K,T{-}3U)$ at $N=p^2$ be given, with $m=T-3$ unit arcs, the
two-value size window of Lemma~\ref{lem:rungT}, and overlap budget
$\mathrm{ov}=\Sigma s-p=(T-6)k+O(1)$. Then:
\smallskip
\emph{(i)} every census-admissible family in the cell dies in the fiber
cascade within
\begin{multline*}
C(T,k)\;=\;
\frac{\log\!\bigl(c_1\,\mathrm{ov}^{\,m-1}\bigr)}
     {(m-1)\,\log\!\bigl(p/\mathrm{ov}\bigr)-\log(c_1\Gamma)}
+\;\frac{2\Lambda_0\,(T-2)}{T}\;+\;1\;+\;O(1/k)
\end{multline*}
fibers, where $\Gamma$ is the uniform constant of $\mathrm{Hq}$; the
denominator is uniformly positive under $\mathrm{H1g}$'s margin clause,
and its $k\to\infty$ form $(m-1)\log(T/(T-6))-\log(c_1\Gamma)$ gives
the $T$-flat plateau $\approx 9$, $9.2$, $10.1$ at $T=8,9,10$ (at
$k=2,3,3$)---with the never-exceeded worst case
\[
(T-4)+N_{\mathrm{res}}+O(1)\;\le\;T+3,
\qquad
N_{\mathrm{res}}=1+\tfrac{2\Lambda_0(T-2)}{T};
\]
\smallskip
\emph{(ii)} since $C(T,k)\ll|U|=(T-2)k$ for every $T\ge 8$, $k\ge 2$,
with a surplus growing linearly in $k$ (conditional on the ledger: the
depth bound is what (i) provides; the surplus language refers to the
fiber budget, not to an independently proved decay rate), the cascade
exhausts the fiber budget: no census-admissible family admits a
covering placement valid on all of $U$, and hence
\smallskip
\emph{(iii)} the cell $(1K,T{-}3U)$ at $p^2$ is empty; and by the
lossless reformulation equivalence, $\mathrm{LRC}$ holds at that rung.
The measured range of the calibration is $T=8,9,10$ at $p\le 73$, with
committed ground-truth cells at $p\le 61$ for $T\in\{7,8\}$; the
statement at primes beyond the calibration is exactly what the
hypotheses purchase.
\end{theorem}

\begin{proposition}[the Hq-induction; the derivation step]\label{prop:Hq}
Let the frontier cell and notation be as in Theorem~\ref{thm:SC}, and
assume $\mathrm{Hq}$ with constant $\Gamma$ and the volume law
$\mathrm{H1g}$ with $|\mathrm{Sol}(f,\mathrm{sz})|\le
c_1\,\mathrm{ov}^{\,m-1}$ for every admissible family and size tuple.
Write $S_i$ for the survivor set after $i$ processed fibers,
$A_i=\lambda_i^{-1}(\mathrm{Sol}_i-\tau_i)$ for the $i$-th fiber's
dilate constraint, $d=m-1$, $B=c_1\,\mathrm{ov}^{\,d}$, and
$r=\Gamma B/p^{\,d}$. Then $|S_0|\le B$, each step contracts,
$|S_i|\le r\,|S_{i-1}|$, and consequently
\[
|S_q|\;\le\;B\,r^{\,q},
\qquad\text{so the cascade is empty whenever}\qquad
q\;>\;\frac{\log B}{(m-1)\,\log(p/\mathrm{ov})-\log(c_1\Gamma)}.
\]
The denominator is uniformly positive iff
$c_1\Gamma<\bigl(p/\mathrm{ov}\bigr)^{m-1}$ with a margin uniform in
$(T,k)$; in the $k\to\infty$ limit $p/\mathrm{ov}\to T/(T-6)$, and the
uniform condition is the threshold clause of $\mathrm{H1g}$.
\end{proposition}

\begin{proof}
$S_0$ is the family's initial placement set, bounded by the volume law.
For each $i$, $S_i=S_{i-1}\cap A_i$ with $A_i$ the $i$-th fiber's
constraint set, so $\mathrm{Hq}$ applied to the admissible history
$S_{i-1}$ gives $|S_i|\le\Gamma\,|S_{i-1}|\,|\mathrm{Sol}_i|/p^{d}\le
r\,|S_{i-1}|$ using $|\mathrm{Sol}_i|\le B$. Induction gives
$|S_q|\le B\,r^{q}$, and $|S_q|<1$ forces emptiness. The resonance tail
(the weak-cut fibers with $\bar\lambda\le\Lambda_0$, each of which may
fail to contract) adds $N_{\mathrm{res}}=1+2\Lambda_0(T-2)/T$ fibers
before the contraction count begins, which is the tail of $C(T,k)$.
\end{proof}

\paragraph{What changed in this revision, and why.}
The theorem's predecessor assumed an unquantified ``sequential
composition'' input (the former $\mathrm{H0}$ row: $C_i\approx
\bar\gamma^{\,i}|\mathrm{Sol}|^{i+1}/p^{(m-1)i}$) and a \emph{pairwise}
dilate-sparsity input (the former $\mathrm{H2m}$). As an audit of this
paper observed, a recurrence for the \emph{independent-extension model}
and a bound on \emph{unconditional} pair intersections do not compose
into a bound on the \emph{iterated, history-dependent} intersections the
actual cascade performs: pairwise sparsity does not control
higher-order correlations, and three sets can have acceptable pairwise
overlaps while sharing a structured common core. The repair is minimal
and explicit: the multi-fiber survivor inequality $\mathrm{Hq}$ (above)
replaces both, with $\mathrm{H2m}$ retained as its one-step case and
$\mathrm{H0}$'s measurements retained as its calibration support. The
derivation from the ledger is then exactly
Proposition~\ref{prop:Hq}'s induction---complete, and conditional on
named open inputs only."""))

# ---- 2. ledger: replace H0 row with Hq; H1g sharpened; H2m base case; H2r OPEN
P.append(('secsc.tex',
r"""$\mathrm{H0}$ &
Sequential composition: the per-fiber cuts compose,
$C_i\approx\bar\gamma^{\,i}\,|\mathrm{Sol}|^{i+1}/p^{(m-1)i}$; no
tangent-drift trajectory stays inside the solution family for many
fibers. &
measured benign in range &
Depths $\le 4$ at the committed cells ($764$-pair record); $\le 11$
across the three calibration rungs; exit time $\le 2$ exhaustive over
admissible families at $T=7$ and $T=8$-critical. Route: tangency
codimension $\ge 1$. \\""",
r"""$\mathrm{Hq}$ &
Multi-fiber survivor inequality (this revision; replaces the former
unquantified H0 and subsumes H2m as its one-step case): for every
admissible cascade history $S_{i-1}$ and next fiber $i$,
$|S_{i-1}\cap A_i|\le\Gamma\,|S_{i-1}|\,|\mathrm{Sol}_i|/p^{m-1}$,
where $A_i=\lambda_i^{-1}(\mathrm{Sol}_i-\tau_i)$, with $\Gamma$
uniform over families, size tuples, dilations, translations, and prior
surviving fibers. &
OPEN (measured benign in range) &
Depths $\le 4$ at the committed cells ($764$-pair record); $\le 11$
across the three calibration rungs; exit time $\le 2$ exhaustive over
admissible families at $T=7$ and $T=8$-critical; the $i=1$ case is the
measured H2m. Route: tangency codimension $\ge 1$; the pairwise
measurements are the base case, the iterated conditional bound is the
open content. \\"""))

P.append(('secsc.tex',
r"""$\mathrm{H1g}$ &
Distinct-family volume law:
$|\mathrm{Sol}(f,\mathrm{sz})|\le c_1\,\mathrm{ov}^{\,m-1}$ for every
$\pm$-distinct family $f$ and admissible size tuple, with
$c_1=O(1)$ uniform in $(T,k,f)$; equivalently
(Proposition~\ref{prop:threshold}), $|\mathrm{Sol}|$ beats the trivial
bound $p^{\,m-1}$ by a factor growing in $k$, over the threshold
$c_1=\bigl(T/(T-6)\bigr)^{\,m-1}$ ($256/243/244.1$ at $T=8/9/10$). &
measured, margin growing in $k$; chain class proved; three attack
routes closed (the automaton route by theorem) &""",
r"""$\mathrm{H1g}$ &
Distinct-family volume law with uniform margin (sharpened in this
revision): $|\mathrm{Sol}(f,\mathrm{sz})|\le
c_1\,\mathrm{ov}^{\,m-1}$ for every $\pm$-distinct family $f$ and
admissible size tuple, with $c_1\le\bigl(T/(T-6)\bigr)^{\,m-1}
e^{-\eta}$ for some $\eta>0$ uniform in $(T,k,f,\mathrm{sz})$---i.e.,
the constant beats the threshold
$c_1=\bigl(T/(T-6)\bigr)^{\,m-1}$ ($256/243/244.1$ at $T=8/9/10$,
where the depth denominator of Proposition~\ref{prop:Hq} vanishes) by
a fixed multiplicative margin. The earlier phrasing ``$c_1=O(1)$
uniform, equivalently beats the trivial bound by a factor growing in
$k$'' is retained as calibration language, not as an equivalent
hypothesis. &
OPEN (measured, margin growing in $k$; chain class proved; three
attack routes closed, the automaton route by theorem) &"""))

P.append(('secsc.tex',
r"""$\mathrm{H2m}$ &
Mid-region dilate-sparsity:
$|\mathrm{Sol}_0\cap\lambda^{-1}(\mathrm{Sol}_j-\tau)|\le
\gamma_0\,|\mathrm{Sol}_0|\,|\mathrm{Sol}_j|/p^{m-1}$ for every fiber
pair with $\bar\lambda>\Lambda_0$, $\gamma_0=O(1)$ uniform. &
measured; proved for the chain class (loose constants) &""",
r"""$\mathrm{H2m}$ &
Mid-region dilate-sparsity---now the $i=1$ (one-step, unconditional)
case of $\mathrm{Hq}$, retained for its stronger measured record:
$|\mathrm{Sol}_0\cap\lambda^{-1}(\mathrm{Sol}_j-\tau)|\le
\gamma_0\,|\mathrm{Sol}_0|\,|\mathrm{Sol}_j|/p^{m-1}$ for every fiber
pair with $\bar\lambda>\Lambda_0$, $\gamma_0=O(1)$ uniform. &
OPEN (measured; proved for the chain class with loose constants) &"""))

P.append(('secsc.tex',
r"""$\mathrm{H2r}$ &
Resonance budget: the weak-cut fibers ($\bar\lambda\le\Lambda_0$ or
small-rational $\lambda$) number
$\le 2\Lambda_0(T-2)/T+1$, $k$-free; the $\lambda=-1$ fiber is never a
cut. &
derived given $\Lambda_0$; the $\lambda=-1$ clause proved &""",
r"""$\mathrm{H2r}$ &
Resonance budget, now stated as an open assumption rather than
``derived'': the weak-cut fibers---those whose dilation ratio
$\lambda_i$ has small rational height, $\bar\lambda\le\Lambda_0$ in
the height function of \S\ref{sec:scale-model}---number $\le
2\Lambda_0(T-2)/T+1$, $k$-free; the $\lambda=-1$ fiber is never a cut
(proved, Lemma~\ref{lem:reflect}). The count's formula shape is
derived from the reflection lemma and a Gauss-circle-type heuristic;
$\Lambda_0\approx 4$ is data-chosen, not proved. &
OPEN (measured support; formula shape derived; the $\lambda=-1$ clause
proved) &"""))

# ---- 3. ledger preamble + caption: statuses as conditional inputs
P.append(('secsc.tex',
r"""\paragraph{The hypothesis ledger.}
Each row of Table~\ref{tab:ledger} names one input, its precise form,
its status, and its evidence or attack route. Two of the six (H1g and
the H2 pair) are the substantive open content; the others are proved,
machine-validated structure (NF), a committed theorem (EQ), or measured
benign in the tested range with a named route (H0). The statuses are
stated once here and are not restated with each citation.""",
r"""\paragraph{The hypothesis ledger.}
Each row of Table~\ref{tab:ledger} names one input, its precise form,
its status, and its evidence or attack route. The rows are
\emph{conditional inputs}, not evidence toward the theorem's validity:
three are open hypotheses with measured support ($\mathrm{Hq}$,
$\mathrm{H1g}$, $\mathrm{H2r}$, the first sharpened and quantified in
this revision); one is the measured one-step base case
($\mathrm{H2m}$); the others are proved, machine-validated structure
($\mathrm{NF}$) and a committed theorem ($\mathrm{EQ}$). The statuses
are stated once here and are not restated with each citation."""))

P.append(('secsc.tex',
r"""\caption{The hypothesis ledger of Theorem~\ref{thm:SC}. Statuses:
\emph{proved} (with the verification record), \emph{measured} (with the
tested range), or \emph{derived} from a disclosed constant. The two
substantive open inputs are H1g and the H2 pair; everything else is
proved, committed, or measured benign in range.}\label{tab:ledger}\\""",
r"""\caption{The hypothesis ledger of Theorem~\ref{thm:SC}. Statuses:
\emph{OPEN} (a conditional input with measured support), \emph{proved}
(with the verification record), \emph{measured} (base case with tested
range), or \emph{committed}. The substantive open inputs are Hq, H1g,
and H2r; every row is an input the theorem assumes, not evidence that
the theorem holds.}\label{tab:ledger}\\"""))

# ---- 4. "What is conditional" paragraph: the reduction chain now cites prop:Hq
P.append(('secsc.tex',
r"""\paragraph{What is conditional and what is not.}
The reduction chain \emph{(i)}$\to$\emph{(ii)}$\to$\emph{(iii)} is
proved from committed machinery: the two-fiber reduction
(Lemma~\ref{lem:twofiber}: the cascade step \emph{is} a
dilate-intersection of census solution sets), the drift law
(Lemma~\ref{lem:rungT}, machine-verified $3{,}528/3{,}528$), the
$T$-free fiber machinery, and the lossless equivalence of
Corollary~\ref{cor:tau}. The arithmetical content---that the per-fiber
cut beats the pool within $C(T,k)$ fibers---is proved \emph{given} the
ledger's inputs.""",
r"""\paragraph{What is conditional and what is not.}
The reduction chain \emph{(i)}$\to$\emph{(ii)}$\to$\emph{(iii)} is
proved from committed machinery: the two-fiber reduction
(Lemma~\ref{lem:twofiber}: the cascade step \emph{is} a
dilate-intersection of census solution sets), the drift law
(Lemma~\ref{lem:rungT}, machine-verified $3{,}528/3{,}528$), the
$T$-free fiber machinery, and the lossless equivalence of
Corollary~\ref{cor:tau}. The arithmetical content---that the per-fiber
cut beats the pool within $C(T,k)$ fibers---is proved \emph{given} the
ledger's inputs, by the induction of Proposition~\ref{prop:Hq} (in
this revision an explicit theorem-level step; its predecessor leaned
on the independent-extension recurrence without the quantified Hq)."""))

# ---- 5. scope paragraph (p2-side/frontier/rigid-zone levels)
P.append(('secsc.tex',
r"""\paragraph{Where the derivation lives.}
Section~\ref{sec:scale} derives the cascade model and proves the proved
pieces (Lemmas~\ref{lem:twofiber}, \ref{lem:pool}, \ref{lem:spacing},
\ref{lem:reflect}, \ref{lem:renorm}; Theorem~\ref{thm:chaindilate});""",
r"""\paragraph{Scope levels, fixed once.}
Three different closures circulate in this paper and must not be
conflated: \emph{cells} (the $T=7$ and $T=8$-critical cells are closed
within the tested primes, $p\le 61$, by exact census plus the shear
step); \emph{rungs} (the frontier cells of every $T\ge 8$ are closed
only conditionally, by Theorem~\ref{thm:SC} under its ledger); and
\emph{rigid-zone finiteness} (open; the five-prime zone of
\S\ref{sec:rigid} is an empirical finding of the verified range).
Statements below use the matching noun.

\paragraph{Where the derivation lives.}
Section~\ref{sec:scale} derives the cascade model and proves the proved
pieces (Lemmas~\ref{lem:twofiber}, \ref{lem:pool}, \ref{lem:spacing},
\ref{lem:reflect}, \ref{lem:renorm}; Theorem~\ref{thm:chaindilate});"""))
