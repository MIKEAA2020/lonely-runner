#!/usr/bin/env python3
"""patch_flagship_v2.py — apply the v2 audit revisions to paper2_sources_v2.
Each patch is (file, old, new); all must match exactly. Run from anywhere."""
import sys, os
D = '/home/z/my-project/download/paper2_sources_v2'

P = []
def patch(f, old, new):
    P.append((f, old, new))

# ---------------- sec2: quotient space rename, volume, three-ways, gcd clause
patch('sec2.tex',
"""\\item \\textbf{Quotient-torus form.} Let $V=\\R^n/\\R v$, let
$\\pi\\colon\\R^n\\to V$ be the projection, $\\Lambda=\\pi(\\Z^n)\\cong\\Z^{n-1}$,
and let $\\|\\cdot\\|_V$ be the quotient norm
$\\|c\\|_V=\\operatorname{dist}_\\infty(x,\\R v)$ for any lift $x$ of $c$.
Then $\\delta(V)=\\operatorname{dist}_V(c^*,\\Lambda)$ where
$c^*=\\pi(\\tfrac12\\mathbf{1})$ is the \\emph{center class}, a $2$-torsion
point whose coordinates are $\\tfrac12$ on even speeds and $0$ on odd
speeds. The $\\delta$-form reads: the center class has quotient-norm
depth at most $\\mathrm{thr}_n$ in every $V$.
\\item \\textbf{Covering-radius form ($\\rho$-form).}
Let $\\rho(V)=\\max_{c\\in V}\\operatorname{dist}_V(c,\\Lambda)$ be the
covering radius of $\\Lambda$ in the quotient norm. The $\\rho$-form of the
conjecture asserts $\\rho(V)\\le\\mathrm{thr}_n$ for all $V$.""",
"""\\item \\textbf{Quotient-torus form.} Let $\\mathcal{X}_v=\\R^n/\\R v$ be
the quotient vector space (the calligraphic symbol avoids colliding with
the speed vector $V$), let $\\pi\\colon\\R^n\\to\\mathcal{X}_v$ be the
projection, $\\Lambda=\\pi(\\Z^n)\\cong\\Z^{n-1}$, and let
$\\|\\cdot\\|_{\\mathcal{X}_v}$ be the quotient norm
$\\|c\\|_{\\mathcal{X}_v}=\\operatorname{dist}_\\infty(x,\\R v)$ for any lift
$x$ of $c$. Then $\\delta(V)=\\operatorname{dist}_{\\mathcal{X}_v}(c^*,\\Lambda)$
where $c^*=\\pi(\\tfrac12\\mathbf{1})$ is the \\emph{center class}, a
$2$-torsion point whose coordinates are $\\tfrac12$ on even speeds and $0$
on odd speeds. The $\\delta$-form reads: the center class has quotient-norm
depth at most $\\mathrm{thr}_n$ in every $\\mathcal{X}_v$.
\\item \\textbf{Covering-radius form ($\\rho$-form).}
Let $\\rho(V)=\\max_{c\\in\\mathcal{X}_v}\\operatorname{dist}_{\\mathcal{X}_v}(c,
\\Lambda)$ be the covering radius of $\\Lambda$ in the quotient norm. The
$\\rho$-form of the conjecture asserts $\\rho(V)\\le\\mathrm{thr}_n$ for all
$V$.""")

patch('sec2.tex',
"""each of $\\ell_1$-norm at most $1$; hence
$\\|c\\|_V\\le\\|c\\|_\\infty$ and Lipschitz certification needs no
norm-equivalence constant.""",
"""each of $\\ell_1$-norm at most $1$; hence
$\\|c\\|_{\\mathcal{X}_v}\\le\\|c\\|_\\infty$ and Lipschitz certification
needs no norm-equivalence constant.""")

patch('sec2.tex',
"""The \\emph{LR-zonotope} of $V$ is $Z(V)=\\pi([0,1]^n)\\subset V$; it is an
$(n{-}1)$-dimensional lattice zonotope with respect to $\\Lambda$, with
normalized volume $\\sum_j v_j$. Its lattice-point geometry encodes the""",
"""The \\emph{LR-zonotope} of $V$ is $Z(V)=\\pi([0,1]^n)\\subset
\\mathcal{X}_v$; it is an $(n{-}1)$-dimensional lattice zonotope with
respect to $\\Lambda$, with relative volume $\\sum_j v_j$ (normalized
volume $d!\\,\\sum_j v_j$, so that $h^*(1)=d!\\,\\sum_j v_j$)
\\cite{Shephard,BeckRobins}. Its lattice-point geometry encodes the""")

patch('sec2.tex',
"""with $\\gcd$ of the empty set $0$, and the $h^*$-polynomial is obtained
from the Eulerian expansion""",
"""and the $h^*$-polynomial is obtained
from the Eulerian expansion""")

patch('sec2.tex',
"""$h^*(z)=\\sum_k e_k\\,(1-z)^{d-k}A_k(z)$ with $A_k$ the Eulerian
polynomials. Formula \\eqref{eq:ehr} was derived in the program's records
and verified in two independent ways on every instance used here (a DFS
coset count and the standard Ehrhart inversion;
Appendix~\\ref{sec:prov}).""",
"""$h^*(z)=\\sum_k e_k\\,(1-z)^{d-k}A_k(z)$ with $A_k$ the Eulerian
polynomials \\cite{Frobenius}. Formula \\eqref{eq:ehr} was derived in the
program's records and verified in three independent ways on every
instance used here (the Eulerian expansion, the standard Ehrhart
inversion, and a depth-first coset count at $t=1,2$;
Appendix~\\ref{sec:prov}).""")

patch('sec2.tex',
"""two statements coincide on all tested instances for $n\\le 4$---which is
why the identification was easy to believe---and separate from $n=5$ on,
with exact witnesses.""",
"""two statements coincide on all tested instances for $n\\le 4$---more
precisely, the $\\rho$-form \\emph{holds} on every instance tested at
$n\\le 4$, certified through $m=24$ at $n=3$ and $m=16$ at $n=4$, while
the center-is-deepest property holds exactly at the multiples---which is
why the identification was easy to believe---and separate from $n=5$ on,
with exact witnesses.""")

# ---------------- sec3: mechanism observations, instance 3 m-def, instance 4 wording
patch('sec3.tex',
"""\\item \\textbf{Losslessness is as-hard-ness.} A lossless reformulation
carries the full difficulty of $\\lrc_n$ by construction. No instrument
that respects the equivalence---that is, no instrument whose outputs are
invariant under the lossless rewriting---can distinguish the
reformulation from the original.""",
"""\\item \\textbf{Losslessness transfers truth exactly.} A proof of the
reformulated statement is ipso facto a proof of $\\lrc_n$, and conversely.
A rewrite may still change tractability---the program's own finitization
is an equivalent rewrite that made the problem computable instance by
instance---but it cannot discard the conjecture's content: no lossless
rewrite can be simultaneously strictly stronger than $\\lrc_n$ and true
at the extremal instances.""")

patch('sec3.tex',
"""\\item \\textbf{The strengthening is false where it matters.} The
toolkit-relevant strengthenings are strictly stronger than $\\lrc$, and at
the extremal instances (the harmonic and rung-2 families, the
reflection-symmetric cascade families, the tight covering
configurations) the extra margin they demand does not exist. When the
demand is a definite statement, it is false, and small $n$ exhibits the
falsity with exact witnesses.""",
"""\\item \\textbf{The strengthening is false where it matters.} The
toolkit-relevant strengthenings are strictly stronger than $\\lrc$ (they
imply it; the converse fails), and at the extremal instances (the
harmonic and rung-2 families, the reflection-symmetric cascade families,
the tight covering configurations) the extra margin they demand does not
exist. When the demand is a definite statement, it is false, and small
$n$ exhibits the falsity with exact witnesses.""")

patch('sec3.tex',
"""\\item \\textbf{The instruments land on the closed side.} The instruments a
lossless reformulation naturally supports---exact counts,
classifications, invariances, exhaustiveness checks---provably compose
into statements that close everything \\emph{except} the missing margin.
In the setting of Instance 1 this is a theorem-level audit finding; in
Instance 2 it is a theorem; in Instances 3 and 4 it is measured, with the
saturation and decay records below.
\\end{enumerate}""",
"""\\item \\textbf{The instruments land on the closed side.} The instruments a
lossless reformulation naturally supports---exact counts,
classifications, invariances, exhaustiveness checks---compose into
statements that close everything \\emph{except} the missing margin: by
theorem in Instance 2 (Lemma~\\ref{lem:monomial}), by audit in Instance 1
(the inventory classification), and by measurement in Instances 3 and 4
(the saturation and decay records below).
\\end{enumerate}""")

patch('sec3.tex',
"""The four subsections that follow develop one instance each. They are
independent: different reformulations, different hoped-for toolkits,
different witnesses. What they share is the mechanism, and we state it
once, here, as the paper's organizing claim: \\emph{in each setting, the
falsity of the strengthening is not an accident of the example but a
consequence of the losslessness that made the reformulation faithful.}""",
"""The four subsections that follow develop one instance each. Three of
them sit on the same pair-sum cascade (an inventory audit, an instrument
closure, a rigidity falsification) and the fourth is independent of it
(the zonotope); what is common to all four is the mechanism, and we
state it once, here, as the paper's organizing claim: \\emph{in each
setting, the falsity or unavailability of the strengthening is not an
accident of the example but a consequence of the losslessness that made
the reformulation faithful.}""")

patch('sec3.tex',
"""Under this reduction $\\lrc_n$
becomes the statement that for every $V$ some pair's bad sets
$B_w=\\{k\\in\\Z_N:(n{+}1)\\wn{wk}<N\\}$ fail to cover $\\Z_N$: a finite
covering problem, exact, and apparently made for the covering toolkit.""",
"""Under this reduction $\\lrc_n$
becomes the statement that for every $V$ some pair's bad sets
$B_p=\\{k\\in\\Z_N:(n{+}1)\\wn{v_pk}<N\\}$, $p\\notin\\{u,w\\}$, fail to
cover $\\Z_N$: a finite covering problem, exact, and apparently made for
the covering toolkit.""")

patch('sec3.tex',
"""Finally, the trivial bound sits exactly at the cascade threshold. The
baseline $|Sol|\\le p^{m-1}$ realizes the constant
$c_1=(T/(T-6)+o(1))^{m-1}$---numerically $256$, $243$, $244$ at $T=8,9,10$
---which is exactly the value at which the depth denominator of the
conditional scaling theorem vanishes. The hypothesis H1g's content is to
beat the trivial bound by a factor growing in $k$; the measured factors
($95\\times$, $1{,}207\\times$, $3{,}803\\times$) grow, but the proved ones
do not exist, and the sharp-constant question and the margin question
are one criterion stated twice.""",
"""Finally, the trivial bound sits exactly at the cascade threshold. The
baseline $|Sol|\\le p^{m-1}$ realizes the constant
$c_1=(T/(T-6)+o(1))^{m-1}$ (here and below $m=T-3$ is the arc count of
the frontier cell, so $m-1=T-4$)---numerically $256$, $243$, $244.14$ at
$T=8,9,10$---which is exactly the value at which the depth denominator
of the conditional scaling theorem vanishes. The hypothesis H1g's
content is to beat the trivial bound by a factor growing in $k$; the
measured factors ($95\\times$, $1{,}207\\times$, $3{,}803\\times$ at
$k=2,3,4$) grow, but the first of them starts \\emph{below} the
threshold constant ($95<244$), which is precisely why the margin must
grow in $k$ and why $k=2$ alone cannot certify anything; the proved
factors do not exist, and the sharp-constant question and the margin
question are one criterion stated twice.""")

patch('sec3.tex',
"""It is false, with exact witnesses, from $n=5$. At the extremal harmonic
instance $v=(1,2,3,4,5)$, the rational point
$c=(\\tfrac1{32},\\tfrac{15}{32},\\tfrac3{32},\\tfrac78)$ has
$D(c)=97/288>\\tfrac13=\\mathrm{thr}_5$: the deepest hole is not the
center, and the covering-radius statement fails at the very instance
where the $\\delta$-form is tight. At $n=6$ the harmonic witness
$c=(0,\\tfrac5{16},\\tfrac{15}{16},\\tfrac34,\\tfrac12)$ gives
$D=29/80>\\tfrac5{14}=\\mathrm{thr}_6$, and the non-extremal instance
$v=(1,2,3,4,5,7)$ carries $c=(\\tfrac1{16},\\tfrac9{16},0,\\tfrac38,0)$
with $D=4/11>\\mathrm{thr}_6$: the falsity is not an extremal-curiosity
confined to the tight family. The rung-2 instance
$v=(1,2,3,4,5,12)$ has a witness at exactly $5/14=\\mathrm{thr}_6$. All
witnesses are non-torsion points---the deepest holes are not the center
classes---and all are single exact evaluations of \\eqref{eq:runnerform},
checkable by hand.""",
"""It is false, with exact witnesses, from $n=5$. At the extremal harmonic
instance $v=(1,2,3,4,5)$, the rational point
$c=(\\tfrac1{32},\\tfrac{15}{32},\\tfrac3{32},\\tfrac78)$ has
$D(c)=97/288>\\tfrac13=\\mathrm{thr}_5$ (the branch-and-bound follow-up
found $607/1792$ at
$c=(\\tfrac{11}{256},\\tfrac{163}{256},\\tfrac{15}{256},
\\tfrac{79}{256})$, and the audit of this paper found
$4183/12288$ at
$c=(\\tfrac{465}{4096},\\tfrac{605}{1024},\\tfrac{1019}{4096},
\\tfrac{315}{4096})$): the deepest hole is not the center, and the
covering-radius statement fails at the very instance where the
$\\delta$-form is tight. At $n=6$ the harmonic witness
$c=(0,\\tfrac5{16},\\tfrac{15}{16},\\tfrac34,\\tfrac12)$ gives
$D=29/80>\\tfrac5{14}=\\mathrm{thr}_6$, and the non-extremal instance
$v=(1,2,3,4,5,7)$---whose center depth is the \\emph{computed}
$\\delta=1/3$, strictly below $\\mathrm{thr}_6$ with slack---carries
$c=(\\tfrac1{16},\\tfrac9{16},0,\\tfrac38,0)$
with $D=4/11>\\mathrm{thr}_6$: the falsity is not an extremal-curiosity
confined to the tight family, and the covering-radius failure there is
not an artifact of center-class tightness. The rung-2 instance
$v=(1,2,3,4,5,12)$ carries $c=(\\tfrac1{16},\\tfrac{11}{16},\\tfrac18,
\\tfrac7{16},\\tfrac78)$ at exactly $D=5/14=\\mathrm{thr}_6$: a
zero-slack boundary case that already separates the deepest hole from
the center ($\\delta=9/26<5/14$) without strictly exceeding the
threshold. All witnesses are non-torsion points---the deepest holes are
not the center classes---and all are single exact evaluations of
\\eqref{eq:runnerform}, checkable by hand (Appendix~\\ref{sec:self}).""")

patch('sec3.tex',
"""For $n\\le4$ the two statements coincide on every instance tested: at
$n=3$ the deepest-hole property (that the maximum of $D$ is attained at
the center class) holds exactly when $3\\mid m$ in the family
$(1,2,m)$ for $m\\le24$, and at $n=4$ it holds when
$4\\mid m$ for $m\\le16$ (away-certified at $m=4,8,12$;
branch-and-bound certified at $m=16$, \\S\\ref{sec:bb}), with exact
witnesses against every other $m$.
The coincidence of the $\\delta$- and $\\rho$-forms in the small cases is
precisely why the stronger form was tempting; the $n=5$ witness is
precisely why it is wrong.""",
"""For $n\\le4$ the $\\rho$-form holds on every instance tested: at $n=3$
the deepest-hole property (that the maximum of $D$ is attained at the
center class) holds exactly when $3\\mid m$ in the family
$(1,2,m)$ for $m\\le24$, and at $n=4$ it holds when $4\\mid m$ for
$m\\le16$ (away-certified at $m=4,8,12$; branch-and-bound certified at
$m=16$, \\S\\ref{sec:bb}), with exact witnesses against every other $m$.
The two forms have the same truth value throughout the tested small
range---which is precisely why the stronger form was tempting; the
$n=5$ witness is precisely why it is wrong.""")

patch('sec3.tex',
"""$h^*$-polynomials are real-rooted, log-concave, and Newton, uniformly,
in every instance tested through $n=10$ in both core families
(Appendix~\\ref{sec:hstar}). But the log-concavity \\emph{strictness}
decays like $\\Theta(1/n)$ ($n\\cdot\\text{strictness}$ drifting from
$\\approx37$ to $\\approx30$ over $n=3,\\dots,10$), so there is no growing
quantitative margin to spend; and across the instance sweeps the
strictness has no measurable correlation with the covering margins
(Spearman $+0.32$ at $n=3$, $-0.31$ at $n=4$; the sign flips).""",
"""The $h^*$-polynomials are real-rooted, log-concave, and Newton, uniformly,
in every instance tested through $n=10$ in both core families
(Appendix~\\ref{sec:hstar}). But the log-concavity \\emph{strictness} $S$
decays roughly like $30/n$ in the measured range ($nS$ drifting from
$\\approx37$ to $\\approx30$ in the harmonic family over
$n=3,\\dots,10$, and from $60.5$ to $\\approx30$ in the rung-2 family),
and Newton's inequalities force $S$ above a floor $\\approx1+4/d$ that
itself declines, so much of the decay is generic and the $1/n$ behavior
cannot persist indefinitely---there is no growing quantitative margin to
spend; and across the instance sweeps the strictness has no
statistically distinguishable correlation with the covering margins at
the achievable sample sizes (Spearman $+0.32$, $p=0.32$, $N=12$ at
$n=3$; $-0.31$, $p=0.45$, $N=8$ at $n=4$; the sign flip is consistent
with noise; these arrays are too small to establish correlation or its
absence).""")

patch('sec3.tex',
"""A
coefficient inequality of a counting polynomial at growing scale $t$ is
the wrong type to pin a fixed-scale metric hole depth; the computation
says so twice, once by decay and once by correlation, and the
covering-radius statement it was meant to feed is false at the extremal
instances anyway. There is no mechanism, and none is needed: the bridge
is dead on both ends.""",
"""A
coefficient inequality of a counting polynomial at growing scale $t$ is
the wrong type to pin a fixed-scale metric hole depth; the computation
says so twice, once by decay and once by the absence of a measurable
signal, and the covering-radius statement it was meant to feed is false
at the extremal instances anyway. There is no mechanism, and none is
needed: the bridge is dead on both ends.""")

patch('sec3.tex',
"""\\caption{The four instances of type mismatch. ``Margin type'' is the
statement class the toolkit needed; ``status'' is its verified fate.}""",
"""\\caption{The four instances of type mismatch. ``Margin type'' is the
statement class the toolkit needed; ``status'' is its verified fate,
with the kind of verification distinguished: \\emph{false} (an exact
witness falsifies a definite statement), \\emph{closed by theorem} (a
theorem shows the instrument cannot produce the output), and
\\emph{zero in inventory} (an audit finding: no statement of the type
exists among the proved results). The first three rows sit on the same
pair-sum cascade; the fourth is the independent zonotope setting.}""")

patch('sec3.tex',
"""Pair-sum covering & growth-in-$k$ margin & no proved instance in full
inventory; $k$-cancellation is a $k$-invariance & audit record
($34$ statements) \\\\
Transfer cascade & spectral gap $1-\\epsilon$ & false by Lemma~M;
survival $=1$ at reflection families & six families, exact \\\\
Character rigidity & per-character profile laws & V$'$/V$''$ false;
envelopes saturated; trivial bound at threshold & exhaustive orbits;
$0/32$ \\\\
LR-zonotope & $\\rho\\le\\mathrm{thr}_n$ (covering) & false from $n=5$;
$h^*$ strictness decays, no correlation & $97/288$, $29/80$, $4/11$
(exact) \\\\""",
"""Pair-sum covering & growth-in-$k$ margin & zero in inventory (audit;
$34$ statements); $k$-cancellation is a $k$-invariance & audit record
($34$ statements) \\\\
Transfer cascade & spectral gap $1-\\epsilon$ & closed by theorem
(Lemma~M; survival $=1$ at reflection families) & six families, exact \\\\
Character rigidity & per-character profile laws & false (V$'$/V$''$
witnesses); envelopes saturated; trivial bound at threshold &
exhaustive orbits; $0/32$ \\\\
LR-zonotope & $\\rho\\le\\mathrm{thr}_n$ (covering) & false from $n=5$;
$h^*$ strictness decays, no measurable signal & $607/1792$,
$4183/12288$, $29/80$, $4/11$ (exact) \\\\""")

# run and report
ok = fail = 0
for f, old, new in P:
    path = os.path.join(D, f)
    s = open(path).read()
    if old in s:
        open(path, 'w').write(s.replace(old, new, 1))
        ok += 1
    else:
        fail += 1
        print('FAILED:', f, '|', old[:70].replace('\n', '⏎'))
print(f'patched {ok}, failed {fail}')
sys.exit(1 if fail else 0)
