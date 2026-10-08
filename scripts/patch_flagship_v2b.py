#!/usr/bin/env python3
"""patch_flagship_v2b.py — part 2: sec4, sec5, sec6, secA, secB, main.tex."""
import sys, os
D = '/home/z/my-project/download/paper2_sources_v2'
P = []
def patch(f, old, new): P.append((f, old, new))

# ---------------- sec4: classical references fix, tau explicit
patch('sec4.tex',
"""a chain-class dilate theorem, and a renormalization lemma;
the Mirsky--Newman and de Bruijn--Schoenberg covering-system identities
\\cite{MirskyNewman,deBruijnSchoenberg} play the classical supporting
roles.""",
"""a chain-class dilate theorem, and a renormalization lemma;
the Mirsky--Newman covering-system theorem and the
R\\'edei--de Bruijn--Schoenberg structure theorem for vanishing sums of
roots of unity \\cite{MirskyNewman,LamLeung} play the classical
supporting roles.""")

patch('sec4.tex',
"""\\begin{theorem}[Sum-only equality; Theorem 1.1 of the companion
monograph]\\label{thm:equality}
Let $V=(v_1,\\dots,v_n)$ have distinct positive entries with
$\\gcd(V)=1$ and $0<M(V)<\\tfrac12$. For each unordered pair $\\{u,w\\}$
let $\\tau(u,w)$ denote the optimal loneliness of the
$(n-1)$-runner instance obtained by restricting to the pair-sum lattice
of denominator $N=v_u+v_w$. Then""",
"""\\begin{theorem}[Sum-only equality; Theorem 1.1 of the companion
monograph]\\label{thm:equality}
Let $V=(v_1,\\dots,v_n)$ have distinct positive entries with
$\\gcd(V)=1$ and $0<M(V)<\\tfrac12$. For each unordered pair $\\{u,w\\}$
let $\\tau(u,w)$ denote the optimal loneliness of the $(n-1)$-runner
instance obtained by restricting to the pair-sum lattice of denominator
$N=v_u+v_w$---explicitly (Appendix~\\ref{sec:self}),
\\[
\\tau(u,w)\\;=\\;\\max_{k\\in\\Z_{N}}\\;\\min_{p\\in[n]\\setminus\\{w\\}}
\\wnc{\\tfrac{v_pk}{N}}
\\]
(runners $u$ and $w$ coincide on the lattice, since
$\\wnc{v_uk/N}=\\wnc{v_wk/N}$ there). Then""")

# ---------------- sec5: witness table, three-families sentence, ladder, U(C), naming
patch('sec5.tex',
"""Table~\\ref{tab:refute} lists the witnesses. Each witness $c$ is a
rational point of the quotient torus in $c$-coordinates; the claimed
value is $D(c)$ from \\eqref{eq:runnerform}, computed exactly and
independently by the linear-form/lattice-search implementation of the
program's records. The harmonic instances are extremal for the
$\\delta$-form ($\\delta=\\mathrm{thr}_n$ exactly), so these witnesses kill
the $\\rho$-form at the very instances where the surviving equivalent form
is tight; the $n=6$, $m=7$ witness kills it at a generic, non-extremal
instance; and the rung-2 instance at $m=12$ sits exactly at the
threshold. The witnesses are non-torsion in every case: the deepest
holes are not the center classes.""",
"""Table~\\ref{tab:refute} lists the witnesses. Each witness $c$ is a
rational point of the quotient torus in $c$-coordinates; the claimed
value is $D(c)$ from \\eqref{eq:runnerform}, computed exactly and
independently by the linear-form/lattice-search implementation of the
program's records. The harmonic instances are extremal for the
$\\delta$-form ($\\delta=\\mathrm{thr}_n$ exactly), so these witnesses kill
the $\\rho$-form at the very instances where the surviving equivalent form
is tight; the $n=6$, $m=7$ witness kills it at a generic, non-extremal
instance. The rung-2 instance at $m=12$ sits exactly at the threshold:
its row is a \\emph{boundary} row (zero slack, not a strict refutation),
recorded because it already separates the deepest hole from the center
($\\delta=9/26<5/14$). Three boundary rows at $n=5$ ($m=6,7,8$, witness
value exactly $1/3=\\mathrm{thr}_5$ against $\\delta=3/10$) and one at
$n=6$ ($m=9$) are recorded in the same spirit in the program data. The
witnesses are non-torsion in every case: the deepest holes are not the
center classes.""")

patch('sec5.tex',
"""\\caption{Exact witnesses against the covering-radius form
$\\rho(V)\\le\\mathrm{thr}_n$. All values are exact rational lower bounds
for $\\rho$; enclosures from the grid searches are recorded in the
program's data. Witness points are in $c$-coordinates
(\\S\\ref{sec:prelim}).}
\\label{tab:refute}
\\begin{tabularx}{\\textwidth}{@{} l l l >{\\raggedright\\arraybackslash}X l @{}}
\\toprule
$V$ & $\\delta=\\tfrac12-\\lambda$ & $\\mathrm{thr}_n$ & witness $c$ &
$D(c)$ \\\\
\\midrule
$(1,2,3,4,5)$ & $1/3$ & $1/3$ &
$(\\tfrac1{32},\\tfrac{15}{32},\\tfrac3{32},\\tfrac78)$ &
$97/288{\\approx}.3368$ \\\\
$(1,2,3,4,5,6)$ & $5/14$ & $5/14$ &
$(0,\\tfrac5{16},\\tfrac{15}{16},\\tfrac34,\\tfrac12)$ & $29/80=.3625$ \\\\
$(1,2,3,4,5,7)$ & $1/3$ & $5/14$ &
$(\\tfrac1{16},\\tfrac9{16},0,\\tfrac38,0)$ & $4/11{\\approx}.3636$ \\\\
$(1,2,3,4,5,12)$ & $9/26$ & $5/14$ & non-torsion point &
$5/14=\\mathrm{thr}_6$ \\\\
\\bottomrule
\\end{tabularx}
\\end{table}""",
"""\\caption{Exact witnesses against the covering-radius form
$\\rho(V)\\le\\mathrm{thr}_n$. All values are exact rational lower bounds
for $\\rho$; enclosures from the grid searches are recorded in the
program's data. Witness points are in $c$-coordinates
(\\S\\ref{sec:prelim}). ``$>\\mathrm{thr}$?'' marks whether the witness
\\emph{strictly} exceeds the threshold: only strict rows falsify the
$\\rho$-form; boundary rows (equality) are recorded separately because
they already separate the deepest hole from the center class. The
$\\delta$ entries for non-harmonic instances are computed values, not
harmonic identities.}
\\label{tab:refute}
\\begin{tabularx}{\\textwidth}{@{} l l l >{\\raggedright\\arraybackslash}X l c @{}}
\\toprule
$V$ & $\\delta=\\tfrac12-\\lambda$ & $\\mathrm{thr}_n$ & witness $c$ &
$D(c)$ & $>\\mathrm{thr}$? \\\\
\\midrule
$(1,2,3,4,5)$ & $1/3$ & $1/3$ &
$(\\tfrac{11}{256},\\tfrac{163}{256},\\tfrac{15}{256},\\tfrac{79}{256})$ &
$607/1792{\\approx}.3387$ & yes \\\\
$(1,2,3,4,5)$ & $1/3$ & $1/3$ &
$(\\tfrac{465}{4096},\\tfrac{605}{1024},\\tfrac{1019}{4096},
\\tfrac{315}{4096})$ &
$4183/12288{\\approx}.3404$ & yes \\\\
$(1,2,3,4,5)$ & $1/3$ & $1/3$ &
$(\\tfrac1{32},\\tfrac{15}{32},\\tfrac3{32},\\tfrac78)$ &
$97/288{\\approx}.3368$ & yes \\\\
$(1,2,3,4,5,6)$ & $5/14$ & $5/14$ &
$(0,\\tfrac5{16},\\tfrac{15}{16},\\tfrac34,\\tfrac12)$ & $29/80=.3625$ &
yes \\\\
$(1,2,3,4,5,7)$ & $1/3$ & $5/14$ &
$(\\tfrac1{16},\\tfrac9{16},0,\\tfrac38,0)$ & $4/11{\\approx}.3636$ & yes
\\\\
$(1,2,3,4,5,12)$ & $9/26$ & $5/14$ &
$(\\tfrac1{16},\\tfrac{11}{16},\\tfrac18,\\tfrac7{16},\\tfrac78)$ &
$5/14=\\mathrm{thr}_6$ & no (boundary) \\\\
\\bottomrule
\\end{tabularx}
\\end{table}""")

patch('sec5.tex',
"""Three further exact refutation families fill out the record: the $n=5$
family $(1,2,3,4,m)$ fails deepest-hole for every non-multiple of $5$
up to $m=30$ (twenty exact witnesses, the sharpest being $6/19>5/16$ at
$m=15$); the $n=6$ family fails for every tested $m\\in\\{6,\\dots,12,18,
24\\}$, including all multiples of $6$ through $4n$ (extinction of the
modular deepest-hole law at $n=6$; the narrowest margin is $7/5200$ at
$m=24$, witness value $71/208$); and at $n=6$, $m=9$ the witness sits at
$\\rho\\ge\\mathrm{thr}_6$ exactly. A useful way to read the table: the
$\\rho$-form was already falsified three times over before this session,
and the session's contribution below is on the \\emph{positive} side of the
ladder---certifying where the deepest hole genuinely is the center.""",
"""Three further exact refutation families fill out the record---but it
matters what each refutes. The $n=5$ family $(1,2,3,4,m)$ refutes the
\\emph{modular deepest-hole law}, not the $\\rho$-form: the law's
equality clause (deepest hole at the center when $5\\mid m$) fails
exactly at the multiples $m=5$ ($97/288>1/3$) and $m=15$
($6/19>5/16=\\delta$), while the deepest-hole-at-center property fails
at all twenty non-multiples (exact witnesses, the largest margins
$1/3>3/10$ at $m=6,7,8$). The $n=6$ family extinguishes the modular law
at every tested $m\\in\\{6,\\dots,11,18,24\\}$ (all multiples of $6$
through $4n$ included; the narrowest modular margin is $7/5200$ at
$m=24$, witness value $71/208$ against $\\delta=17/50$); against the
$\\rho$-form itself these witnesses fall short of the threshold
($6/19<1/3$; $71/208<5/14$) and are recorded as deepest-hole evidence,
not covering-radius refutations. At $n=6$, $m=9$ the witness sits at
$\\rho\\ge\\mathrm{thr}_6$ exactly---a boundary row, like the $n=5$
$m=6,7,8$ rows. A useful way to read the table: the $\\rho$-form was
falsified three times over before this session (the three strict rows:
$97/288$, $29/80$, $4/11$), and the session's contribution below is on
the \\emph{positive} side of the ladder---certifying where the deepest
hole genuinely is the center.""")

patch('sec5.tex',
"""but at $d=4$ a global exact
certificate was beyond the program's earlier arrangement machinery.
This session's branch-and-bound method (\\S\\ref{sec:bb}) upgrades three
of the four. At $m=25$ and $m=30$ the trees empty in under five minutes
each ($10{,}414$ and $10{,}646$ certified leaves; $\\rho=4/13$ and
$\\rho=19/62$ exactly), and at $m=20$ the tree empties in sixteen
minutes ($37{,}018$ leaves; $\\rho=13/42$).""",
"""but at $d=4$ a global exact
certificate was beyond the program's earlier arrangement machinery.
(Only $m=10$ is a rung-2 instance; $m=20,25,30$ are the higher
multiples of $5$.) This session's branch-and-bound method
(\\S\\ref{sec:bb}) upgrades three of the four. At $m=25$ and $m=30$ the
trees empty in under five minutes each ($10{,}414$ and $10{,}646$
certified leaves; $\\rho=4/13$ and $\\rho=19/62$ exactly), and at $m=20$
the tree empties in sixteen minutes ($37{,}018$ leaves;
$\\rho=13/42$).""")

patch('sec5.tex',
"""$571/1792=7/22+9/19712$---a quantized minimax loss, not an excess. The
status of $m=10$ is therefore: \\emph{consistent, with the two-sided
evidence substantially strengthened, and exact equality open at $d=4$}.""",
"""$571/1792=7/22+9/19712$---a fixed quantization offset, consistent with
pure minimax loss rather than genuine excess. We stress what this does
and does not establish: the surviving boxes are exactly the ones the
method could not resolve, and a residual volume below $6\\times10^{-7}$
does not exclude a deeper point hiding inside it; the reading is a
diagnostic, not a determination. The status of $m=10$ is therefore:
\\emph{consistent, with the two-sided evidence substantially
strengthened, and exact equality open at $d=4$}.""")

patch('sec5.tex',
"""built on the runner formulation \\eqref{eq:runnerform} rather than on the
arrangement enumeration of the program's records. For a dyadic
$c$-box $C$ define
\\[
U(C)\\;=\\;\\min_{\\sigma\\in[0,1)}\\;\\max_j\\;\\sup_{c\\in C}\\wnc{c_j-v_j\\sigma},
\\]
computed exactly: each inner supremum is piecewise linear in $\\sigma$
with a rational kink set, so the minimum of their maximum is attained at
a rational kink or crossing. By the minimax inequality,
$\\max_C D\\le U(C)$, so a box with exact $U(C)\\le r$ is \\emph{certified}
to contain no point of depth above $r$; boxes that fail are split
dyadically.""",
"""built on the runner formulation \\eqref{eq:runnerform} rather than on the
arrangement enumeration of the program's records. For a dyadic
$c$-box $C$ define
\\[
U(C)\\;=\\;\\min_{\\sigma\\in[0,1)}\\;
\\max\\Bigl(\\wnc{\\sigma},\\;\\max_{2\\le j\\le n}\\;
\\sup_{c\\in C}\\wnc{c_j-v_j\\sigma}\\Bigr),
\\]
the box-supremum analogue of \\eqref{eq:runnerform} itself (the
$\\wnc\\sigma$ term is the depth of the first runner, whose coordinate is
fixed at $0$; the maximum over $j$ runs over $2,\\dots,n$). It is
computed exactly: each inner supremum is piecewise linear in $\\sigma$
with a rational kink set, so the minimum of their maximum is attained at
a rational kink or crossing. By the minimax inequality,
$\\max_C D\\le U(C)$, so a box with exact $U(C)\\le r$ is \\emph{certified}
to contain no point of depth above $r$; boxes that fail are split
dyadically.""")

patch('sec5.tex',
"""it localizes the false. The same control distinguishes the two
$n=5$ worlds sharply: for the harmonic instance the surviving excess
stays above the threshold by a fixed margin, while for the consistent
rung-2 instances the surviving excess is pure minimax loss and decays
with the box size.""",
"""it localizes the false. The same control distinguishes the two
$n=5$ worlds sharply: for the harmonic instance the surviving excess
stays above the threshold by a fixed margin, while for the consistent
$n=5$ multiples ($m=10$ and the now-certified $m=20,25,30$) the
surviving excess is pure minimax loss and decays with the box size.""")

patch('sec5.tex',
"""The archive's content, stated as a single paragraph: the
covering-radius strengthening of the Lonely Runner statement is false
from $n=5$, with three exact witnesses at extremal and non-extremal
instances, an improved harmonic lower bound $\\rho\\ge607/1792$, and an
extinction record at $n=6$;""",
"""The archive's content, stated as a single paragraph: the
covering-radius strengthening of the Lonely Runner statement is false
from $n=5$, with three exact strict witnesses at extremal and
non-extremal instances (plus the audit's improved harmonic record
$\\rho\\ge4183/12288$, superseding $607/1792$), boundary rows at
$n=5$ ($m=6,7,8$) and $n=6$ ($m=9,12$) at exactly the threshold, and a
modular-law extinction record at $n=6$;""")

# ---------------- sec6: four confirmations, minimax characterization
patch('sec6.tex',
"""and, from this session, an independently certified deepest-hole archive
(equality $\\rho=\\delta$ at $n=4$, $m=16$ and at $n=5$,
$m=20,25,30$, with the improved harmonic lower bound
$\\rho\\ge607/1792$ and the $m=10$ residual characterized as pure
minimax loss).""",
"""and, from this session, an independently certified deepest-hole archive
(equality $\\rho=\\delta$ at $n=4$, $m=16$ and at $n=5$,
$m=20,25,30$, with the harmonic lower bound improved twice---to
$\\rho\\ge607/1792$ by the branch-and-bound follow-up and to
$\\rho\\ge4183/12288$ by this paper's audit---and the $m=10$ residual
consistent with pure minimax loss but not excluded to contain an
unresolved deeper point).""")

patch('sec6.tex',
"""only into closures. The four instances are four independent
confirmations of the same three-part statement.""",
"""only into closures. The four instances are four confirmations of the
same three-part statement---three on the pair-sum cascade, one
independent on the zonotope.""")

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
