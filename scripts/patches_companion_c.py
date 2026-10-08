r"""Companion v2 part 3: H1g statement in secscale, sec2 (attainment, half case,
novelty), sec7, sec9 placement table, bib."""
P = []

# ---- secscale: H1g statement sharpened (mirrors the ledger row)
P.append(('secscale.tex',
r"""H1g is the volume law
$|\mathrm{Sol}(f,\mathrm{sz})|\le c_1\,\mathrm{ov}^{\,m-1}$ for the
$\pm$-distinct families.""",
r"""H1g is the volume law with uniform margin,
$|\mathrm{Sol}(f,\mathrm{sz})|\le c_1\,\mathrm{ov}^{\,m-1}$ for the
$\pm$-distinct families, with $c_1\le\bigl(T/(T-6)\bigr)^{m-1}e^{-\eta}$
for an $\eta>0$ uniform in $(T,k,f,\mathrm{sz})$---the margin clause
that keeps Proposition~\ref{prop:Hq}'s denominator positive. The
bare-$c_1$ form below is the measured content; the margin clause is
the open content."""))

# ---- sec2: attainment proof expanded (audit point 10)
P.append(('sec2.tex',
r"""\emph{Claim: not all binders rise, and not all binders fall.} Suppose all
binders rise; the falling case is symmetric. For non-binders let
$\delta_p=f_p(t^\ast)-m>0$. All $f_p$ are $v_p$-Lipschitz. Choose
$s>0$ with
\[
s<\min\Bigl(\;\min_{p\notin B}\delta_p/v_p,\;\;
                \min_{p\in B}\tfrac{1}{2}(1/2-m)/v_p\Bigr).
\]
For binders (all rising), $0<s'\le s$ keeps the position inside
$(0,\tfrac12)$ where $f_p$ is linear, so
$f_p(t^\ast+s')=m+v_p s'>m$. For non-binders,
$f_p(t^\ast+s')\ge f_p(t^\ast)-v_p s'=m+(\delta_p-v_p s')>m$ by the
choice of $s$. Hence $g(t^\ast+s')>m=M(V)$, contradicting maximality. If
all binders fall, perturb backwards ($s'<0$); the symmetric argument gives
the same contradiction. This proves the claim.""",
r"""\emph{Claim: not all binders rise, and not all binders fall.} Suppose all
binders rise; the falling case is symmetric. For non-binders let
$\delta_p=f_p(t^\ast)-m>0$. All $f_p$ are $v_p$-Lipschitz. Choose
$s>0$ with
\[
s<\min\Bigl(\;\min_{p\notin B}\delta_p/v_p,\;\;
                \min_{p\in B}\tfrac{1}{2}(1/2-m)/v_p\Bigr).
\]
Two remarks make the estimates airtight. First, the non-binder bound
uses only the \emph{global} Lipschitz property
$f_p(t\pm s)\ge f_p(t)-v_p|s|$, which holds across the cusps of
$f_p=\wnc{v_pt}$---a non-binder crossing a cusp during the
perturbation is immaterial, because the one-sided Lipschitz lower
bound is cusp-free. Second, the binder linearity window: a rising
binder's position $x_p=v_pt^\ast\bmod 1$ equals $m\in(0,\tfrac12)$,
and $f_p$ is linear on $(0,\tfrac12)$; the window survives the shift
iff $x_p+v_ps'<\tfrac12$, which the $s$-bound guarantees per binder
(v_j s'\le\tfrac12\,(\tfrac12-m)<\tfrac12-m$)---this is also where
binders with \emph{different} slopes are handled: each binder gets its
own window from the same minimum.
For binders (all rising), $0<s'\le s$ keeps the position inside
$(0,\tfrac12)$ where $f_p$ is linear, so
$f_p(t^\ast+s')=m+v_p s'>m$. For non-binders,
$f_p(t^\ast+s')\ge f_p(t^\ast)-v_p s'=m+(\delta_p-v_p s')>m$ by the
choice of $s$. Hence $g(t^\ast+s')>m=M(V)$, contradicting maximality. If
all binders fall, perturb backwards ($s'<0$); the symmetric argument gives
the same contradiction. This proves the claim."""))

# ---- sec2: split the m=1/2 case into its own lemma
P.append(('sec2.tex',
r"""\item[(b)] if $m=\tfrac12$: all speeds are odd, $t^\ast=\tfrac12$, and
  every pair has even $N$ with $k^\ast=N/2$;
\item[(c)] $m=0$ cannot occur.
\end{enumerate}
\end{lemma}""",
r"""\item[(b)] if $m=\tfrac12$, the half case of
Lemma~\ref{lem:halfcase} applies: all speeds are odd, $t^\ast=\tfrac12$,
and every pair has even $N$ with $k^\ast=N/2$;
\item[(c)] $m=0$ cannot occur.
\end{enumerate}
\end{lemma}

\begin{lemma}[the half case]\label{lem:halfcase}
If $M(V)=\tfrac12$ then every speed is odd and $t^\ast=\tfrac12$
attains it; every pair sum is even, so $t^\ast=(N/2)/N$ is a
pair-lattice time for any pair with $N\ge 4$ (such a pair exists for
$n\ge2$ unless all speeds equal $1$; in that degenerate case
$M=\tfrac12=t^\ast$ with $N=2$, $k^\ast=1$).
\end{lemma}"""))

# ---- sec2: move the (b) proof into the new lemma (adjust proof text)
P.append(('sec2.tex',
r"""(b) If $m=\tfrac12$ then \emph{every} runner binds (no $f_p$ exceeds
$\tfrac12$), i.e.\ $v_p t^\ast\equiv \tfrac12 \pmod 1$ for all $p$. Write
$t^\ast=a/b$ in lowest terms. Then $b\mid v_p-v_1$ for all $p$ (as
$\gcd(a,b)=1$), and $2v_1 a=b(2j+1)$ for some integer $j$, so $b\mid
2v_1$. If $b$ were odd, $b\mid v_1$, hence $b\mid v_p$ for all $p$, hence
$b\mid\gcd(V)=1$ and $b=1$---but then $t^\ast=0$, excluded. So $b=2c$ is
even; then $c\mid v_1$ and $v_p\equiv v_1\pmod{2c}$ forces
$v_p\equiv 0\pmod c$ for all $p$, so $c\mid\gcd(V)=1$, $c=1$, $b=2$:
$t^\ast=\tfrac12$. Then $v_1/2\equiv\tfrac12\pmod 1$ forces $v_1$ odd,
and $v_p\equiv v_1\pmod 2$ makes every speed odd. Conversely, if all
speeds are odd then $t=\tfrac12$ has $f_p=\tfrac12$ for every $p$. Every
pair sum of two odd speeds is even, so $t^\ast=\tfrac12=(N/2)/N$ is a
pair-lattice time for any pair with $N\ge 4$ (such a pair exists for
$n\ge2$ unless all speeds equal $1$; in that degenerate case
$M=\tfrac12=t^\ast$ with $N=2$, $k^\ast=1$).
\end{proof}""",
r"""\end{proof}

\begin{proof}[Proof of Lemma~\ref{lem:halfcase}]
If $m=\tfrac12$ then \emph{every} runner binds (no $f_p$ exceeds
$\tfrac12$), i.e.\ $v_p t^\ast\equiv \tfrac12 \pmod 1$ for all $p$. Write
$t^\ast=a/b$ in lowest terms. Then $b\mid v_p-v_1$ for all $p$ (as
$\gcd(a,b)=1$), and $2v_1 a=b(2j+1)$ for some integer $j$, so $b\mid
2v_1$. If $b$ were odd, $b\mid v_1$, hence $b\mid v_p$ for all $p$, hence
$b\mid\gcd(V)=1$ and $b=1$---but then $t^\ast=0$, excluded. So $b=2c$ is
even; then $c\mid v_1$ and $v_p\equiv v_1\pmod{2c}$ forces
$v_p\equiv 0\pmod c$ for all $p$, so $c\mid\gcd(V)=1$, $c=1$, $b=2$:
$t^\ast=\tfrac12$. Then $v_1/2\equiv\tfrac12\pmod 1$ forces $v_1$ odd,
and $v_p\equiv v_1\pmod 2$ makes every speed odd. Conversely, if all
speeds are odd then $t=\tfrac12$ has $f_p=\tfrac12$ for every $p$,
which proves attainment.
\end{proof}"""))

# ---- sec2: novelty wording
P.append(('sec2.tex',
r"""new as stated, with moderate evidence of novelty.""",
r"""not located in the searched literature (two rounds, $56$ queries, the
arXiv ring enumerated for 2024--2026); the correct prior, given how
actively the spectrum line computes exact optima \cite{Cordella}, is
that pieces may exist in another language."""))

# ---- sec2: distinctness remark
P.append(('sec2.tex',
r"""\begin{remark}[what forced two-sidedness means]\label{rem:twosided}""",
r"""\begin{remark}[which speed convention is proved where]
The equality theorem and both lemmas above are proved for possibly
non-distinct speeds; nothing in this section uses distinctness. The
conjecture itself is stated for distinct speeds (WLOG after collapsing
equal speeds, which can only decrease $M$), and the covering framework
downstream assumes distinct speeds where it says so. The half case
above is the one place where non-distinctness changes the picture
(all speeds equal $1$), and it is handled inside
Lemma~\ref{lem:halfcase}.
\end{remark}

\begin{remark}[what forced two-sidedness means]\label{rem:twosided}"""))
