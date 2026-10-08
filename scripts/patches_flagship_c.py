# patches_flagship_c.py — remaining flagship v2 patches (sec5, sec6, secA, secB)
P = []

P.append(('sec5.tex',
"""and at $n=6$, $m=9$ the witness sits at
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
hole genuinely is the center."""))

# sec5: careful — the three-families sentence begins BEFORE this fragment;
# we must also neutralize the original opening. Handle via second patch
# that removes the now-duplicated opening clause.
P.append(('sec5.tex',
"""Three further exact refutation families fill out the record: the $n=5$
family $(1,2,3,4,m)$ fails deepest-hole for every non-multiple of $5$
up to $m=30$ (twenty exact witnesses, the sharpest being $6/19>5/16$ at
$m=15$); the $n=6$ family fails for every tested $m\\in\\{6,\\dots,12,18,
24\\}$, including all multiples of $6$ through $4n$ (extinction of the
modular deepest-hole law at $n=6$; the narrowest margin is $7/5200$ at
$m=24$, witness value $71/208$); Three further exact refutation families""",
"""Three further exact refutation families"""))

P.append(('sec5.tex',
"""built on the runner formulation \\eqref{eq:runnerform} rather than on the
arrangement enumeration of the program's records. For a dyadic
$c$-box $C$ define
\\[
U(C)\\;=\\;\\min_{\\sigma\\in[0,1)}\\;\\max_j\\;\\sup_{c\\in C}\\wnc{c_j-v_j\\sigma},
\\]
computed exactly: each inner supremum is piecewise linear in $\\sigma$
with a rational kink set, so the minimum of their maximum is attained at
a rational kink or crossing.""",
"""built on the runner formulation \\eqref{eq:runnerform} rather than on
the arrangement enumeration of the program's records. For a dyadic
$c$-box $C$ define
\\[
U(C)\\;=\\;\\min_{\\sigma\\in[0,1)}\\;
\\max\\Bigl(\\wnc{\\sigma},\\;\\max_{2\\le j\\le n}\\;
\\sup_{c\\in C}\\wnc{c_j-v_j\\sigma}\\Bigr),
\\]
the box-supremum analogue of \\eqref{eq:runnerform} itself (the
$\\wnc\\sigma$ term is the depth of the first runner, whose coordinate is
fixed at $0$; the inner maximum runs over $j=2,\\dots,n$). It is
computed exactly: each inner supremum is piecewise linear in $\\sigma$
with a rational kink set, so the minimum of their maximum is attained at
a rational kink or crossing."""))

P.append(('sec6.tex',
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
consistent with pure minimax loss, though a deeper point within the
unresolved volume is not excluded)."""))

# ---------------- secA: Sturm upgrade, decay precision, columns, conjecture
P.append(('secA.tex',
"""Table~\\ref{tab:hstar} lists, for each family and $n$: the
$h^*$-polynomial, real-rootedness (verified at $45$ decimal places on
the exact integer coefficients, with a companion-matrix check agreeing),
log-concavity and Newton's inequalities (exact rational arithmetic), and
the minimum ratio $h_k^2/(h_{k-1}h_{k+1})$, which is the
\\emph{strictness} of log-concavity.""",
"""Table~\\ref{tab:hstar} lists, for each family and $n$: the
$h^*$-polynomial, real-rootedness (verified \\emph{exactly} by
Sturm-sequence sign-change counts on the square-free part of the integer
polynomial---the earlier $45$-decimal-place numerical verification and
the companion-matrix check agreeing---upgraded in this revision from the
numerical form), log-concavity and Newton's inequalities (exact rational
arithmetic), and the minimum ratio $h_k^2/(h_{k-1}h_{k+1})$, which is
the \\emph{strictness} $S$ of log-concavity, together with the Newton
floor $\\bigl((k{+}1)(d{-}k{+}1)\\bigr)/\\bigl(k(d{-}k)\\bigr)$ at the
minimizing $k$ and the excess ratio $S/\\text{floor}$."""))

P.append(('secA.tex',
"""\\caption{$h^*$-polynomials of the LR-zonotopes, $n=3..10$, both core
families. ``RR'' = real-rooted (exact-integer verification at 45 dps);
``LC/Nw'' = log-concave and Newton (exact rational arithmetic);
$S$ = $\\min_k h_k^2/(h_{k-1}h_{k+1})$; $nS$ = $n\\cdot S$. The $n\\le8$
rows reproduce the program's earlier records digit for digit; the
$n=9,10$ rows are new (this session).}
\\label{tab:hstar}
\\begin{tabularx}{\\textwidth}{@{} l c >{\\raggedright\\arraybackslash}X c c r r @{}}
\\toprule
family & $n$ & $h^*(t)$ & RR & LC/Nw & $S$ & $nS$ \\\\""",
"""\\caption{$h^*$-polynomials of the LR-zonotopes, $n=3..10$, both core
families. ``RR'' = real-rooted (exact Sturm verification, this
revision; previously $45$ dps numerical); ``LC/Nw'' = log-concave and
Newton (exact rational arithmetic); $S$ =
$\\min_k h_k^2/(h_{k-1}h_{k+1})$; $nS$ = $n\\cdot S$;
``floor'' = the Newton bound
$\\bigl((k{+}1)(d{-}k{+}1)\\bigr)/\\bigl(k(d{-}k)\\bigr)$ at the
minimizing $k$ (any real-rooted family has $S\\ge\\text{floor}$ there);
$S/f$ = $S/\\text{floor}$. The $n\\le8$ rows reproduce the program's
earlier records digit for digit; the $n=9,10$ rows are new (this
session); the floor and ratio columns are new in this revision.}
\\label{tab:hstar}
\\begin{tabularx}{\\textwidth}{@{} l c >{\\raggedright\\arraybackslash}X c c r r r r @{}}
\\toprule
family & $n$ & $h^*(z)$ & RR & LC/Nw & $S$ & $nS$ & floor & $S/f$ \\\\"""))

# add floor and S/f values to each row
_rows = [
 ("harm. & 3 & $1+7t+4t^2$ & y & y & $12.250$ & $36.8$ \\\\", "harm. & 3 & $1+7z+4z^2$ & y & y & $12.250$ & $36.8$ & $4$ & $3.06$ \\\\"),
 ("rung2 & 3 & $1+11t+6t^2$ & y & y & $20.167$ & $60.5$ \\\\", "rung2 & 3 & $1+11z+6z^2$ & y & y & $20.167$ & $60.5$ & $4$ & $5.04$ \\\\"),
 ("harm. & 4 & $1+18t+35t^2+6t^3$ & y & y & $9.257$ & $37.0$ \\\\", "harm. & 4 & $1+18z+35z^2+6z^3$ & y & y & $9.257$ & $37.0$ & $3$ & $3.09$ \\\\"),
 ("rung2 & 4 & $1+22t+51t^2+10t^3$ & y & y & $9.490$ & $38.0$ \\\\", "rung2 & 4 & $1+22z+51z^2+10z^3$ & y & y & $9.490$ & $38.0$ & $3$ & $3.16$ \\\\"),
 ("harm. & 5 & $1+37t+179t^2+133t^3+10t^4$ & y & y & $6.511$ & $32.6$ \\\\", "harm. & 5 & $1+37z+179z^2+133z^3+10z^4$ & y & y & $6.511$ & $32.6$ & $9/4$ & $2.89$ \\\\"),
 ("rung2 & 5 & $1+45t+239t^2+181t^3+14t^4$ & y & y & $7.013$ & $35.1$ \\\\", "rung2 & 5 & $1+45z+239z^2+181z^3+14z^4$ & y & y & $7.013$ & $35.1$ & $9/4$ & $3.12$ \\\\"),
 ("harm. & 6 & $1+78t+744t^2+1286t^3+399t^4+12t^5$ & y & y & $5.518$\n& $33.1$ \\\\", "harm. & 6 & $1+78z+744z^2+1286z^3+399z^4+12z^5$ & y & y & $5.518$\n& $33.1$ & $2$ & $2.76$ \\\\"),
 ("rung2 & 6 & $1+86t+920t^2+1682t^3+535t^4+16t^5$ & y & y & $5.748$\n& $34.5$ \\\\", "rung2 & 6 & $1+86z+920z^2+1682z^3+535z^4+16z^5$ & y & y & $5.748$\n& $34.5$ & $2$ & $2.87$ \\\\"),
 ("harm. & 7 & $1+147t+2522t^2+8948t^3+7323t^4$\n$+1201t^5+18t^6$ & y & y & $4.335$ & $30.3$ \\\\", "harm. & 7 & $1+147z+2522z^2+8948z^3+7323z^4$\n$+1201z^5+18z^6$ & y & y & $4.335$ & $30.3$ & $16/9$ & $2.44$ \\\\"),
 ("rung2 & 7 & $1+161t+3024t^2+11144t^3+9295t^4$\n$+1551t^5+24t^6$ & y & y & $4.418$ & $30.9$ \\\\", "rung2 & 7 & $1+161z+3024z^2+11144z^3+9295z^4$\n$+1551z^5+24z^6$ & y & y & $4.418$ & $30.9$ & $16/9$ & $2.48$ \\\\"),
 ("harm. & 8 & $1+290t+8317t^2+51458t^3+82947t^4$\n$+35270t^5+3135t^6+22t^7$ & y & y & $3.791$ & $30.3$ \\\\", "harm. & 8 & $1+290z+8317z^2+51458z^3+82947z^4$\n$+35270z^5+3135z^6+22z^7$ & y & y & $3.791$ & $30.3$ & $5/3$ & $2.27$ \\\\"),
 ("rung2 & 8 & $1+298t+9277t^2+60986t^3+102275t^4$\n$+44798t^5+4095t^6+30t^7$ & y & y & $3.829$ & $30.6$ \\\\", "rung2 & 8 & $1+298z+9277z^2+60986z^3+102275t^4$\n$+44798z^5+4095z^6+30z^7$ & y & y & $3.829$ & $30.6$ & $5/3$ & $2.30$ \\\\"),
 ("harm. & 9 & $1+559t+25519t^2+258901t^3+734059t^4$\n$+631993t^5+155305t^6+8035t^7+28t^8$ & y & y & $3.293$ & $29.6$ \\\\", "harm. & 9 & $1+559z+25519z^2+258901z^3+734059z^4$\n$+631993z^5+155305z^6+8035z^7+28z^8$ & y & y & $3.293$ & $29.6$ & $25/16$ & $2.11$ \\\\"),
 ("rung2 & 9 & $1+585t+28885t^2+305035t^3+880149t^4$\n$+764143t^5+188659t^6+9789t^7+34t^8$ & y & y & $3.323$ & $29.9$ \\\\", "rung2 & 9 & $1+585z+28885z^2+305035z^3+880149z^4$\n$+764143z^5+188659z^6+9789z^7+34z^8$ & y & y & $3.323$ & $29.9$ & $25/16$ & $2.13$ \\\\"),
 ("harm. & 10 & $1+1098t+77550t^2+1213626t^3+5529648t^4$\n$+8309694t^5+4180482t^6+626910t^7+19359t^8+32t^9$ & y & y & $2.987$\n& $29.9$ \\\\", "harm. & 10 & $1+1098z+77550z^2+1213626z^3+5529648z^4$\n$+8309694z^5+4180482z^6+626910z^7+19359z^8+32z^9$ & y & y & $2.987$\n& $29.9$ & $3/2$ & $1.99$ \\\\"),
 ("rung2 & 10 & $1+1114t+83790t^2+1377794t^3+6457600t^4$\n$+9866694t^5+5017826t^6+758710t^7+23631t^8+40t^9$ & y & y & $3.004$\n& $30.0$ \\\\", "rung2 & 10 & $1+1114z+83790z^2+1377794z^3+6457600z^4$\n$+9866694z^5+5017826z^6+758710z^7+23631z^8+40z^9$ & y & y & $3.004$\n& $30.0$ & $3/2$ & $2.00$ \\\\"),
]
for old, new in _rows:
    P.append(('secA.tex', old, new))

P.append(('secA.tex',
"""Two quantitative facts were already in the program's records at $n\\le8$
and persist through the new rows. First, the strictness $S$ decays
monotonically in $n$ in both families, with $n\\cdot S$ drifting from
$\\approx37$ at $n=3,4$ to $\\approx29.6$--$30.0$ at $n=9,10$: the decay
is $\\Theta(1/n)$ with a slowly drifting constant, so there is no growing
quantitative margin in the coefficient structure to spend. Second, the
strictness carries no correlation with the covering margins: across the
$n=3$ sweep (twelve instances, $(1,2,m)$, $m=3..14$) the Spearman
coefficient of strictness against the margin is $+0.32$; across the
$n=4$ sweep (eight instances, $(1,2,3,m)$, $m=5..12$) it is $-0.31$.
The sign flips between adjacent rungs; the harmonic family has margin
identically $0$ at every $n$ while its strictness varies by a factor of
four over the table. Both correlations were recomputed in the provenance
run from the persisted instance data and reproduce exactly.""",
"""Two quantitative facts were already in the program's records at $n\\le8$
and persist through the new rows; both are stated here with the precision
the data supports. First, the strictness $S$ decays monotonically in $n$
in both families, with $n\\cdot S$ drifting from $\\approx37$ at
$n=3,4$ (harmonic; the rung-2 family starts at $60.5$) to
$\\approx29.6$--$30.0$ at $n=9,10$: within the measured range the decay
is well described as $\\approx30/n$. That description cannot persist
indefinitely---log-concavity forces $S\\ge1$, and Newton's inequalities
force $S$ above the floor $\\bigl((k{+}1)(d{-}k{+}1)\\bigr)/\\bigl(k(d{-}k)
\\bigr)\\approx1+4/d$ at the middle coefficients, so $nS\\approx30$ would
cross the floor near $n\\approx26$---and much of the measured decline is
the generic Newton-floor decline that any real-rooted family exhibits.
The family-specific content is the excess $S/\\text{floor}$, which declines
from $3.06$ (harmonic) and $5.04$ (rung-2) at $n=3$ to $1.99$ and
$2.00$ at $n=10$ (Table~\\ref{tab:hstar}, Figure~\\ref{fig:nS}); it, too,
declines, so there is no growing quantitative margin in the coefficient
structure to spend. Second, the strictness carries no statistically
distinguishable correlation with the covering margins at the achievable
sample sizes: across the $n=3$ sweep (twelve instances, $(1,2,m)$,
$m=3..14$) the Spearman coefficient of strictness against the margin is
$+0.32$ ($p=0.32$); across the $n=4$ sweep (eight instances,
$(1,2,3,m)$, $m=5..12$) it is $-0.31$ ($p=0.45$). The sign flips
between adjacent rungs and neither coefficient is significant; these
arrays are too small to establish correlation or its absence. The
harmonic family has margin identically $0$ at every $n$ while its
strictness varies by a factor of four over the table. Both coefficients
were recomputed in the provenance run from the persisted instance data
and reproduce exactly."""))

P.append(('secA.tex',
"""\\subsection{The residual question}

What survives is a question of independent arithmetic interest, and we
state it decoupled from the Lonely Runner problem entirely. The
Ehrhart data of the LR-zonotopes is governed by the gcd structure
\\eqref{eq:ehr} of the projected basis---the arithmetic matroid of the
configuration \\cite{dAMoci}. In every instance tested (two core
families through $n=10$, sweeps through $m\\le14$ at $n=3$ and $m\\le12$
at $n=4$), the arithmetic $h^*$-polynomial is real-rooted. The natural
conjecture-shaped question is whether this holds for the projected
bases of the form $(1,\\dots,n-1,m)$ or beyond, and the natural home for
it is the real-rootedness and Lorentzian literature for arithmetic
matroid Tutte specializations \\cite{BrandenHuh,dAMoci,StanleyEC}.""",
"""\\subsection{The residual question}

What survives is a question of independent arithmetic interest, and we
state it decoupled from the Lonely Runner problem entirely. The
Ehrhart data of the LR-zonotopes is governed by the gcd structure
\\eqref{eq:ehr} of the projected basis---the arithmetic matroid of the
configuration \\cite{dAMoci}. In every instance tested (two core
families through $n=10$, sweeps through $m\\le14$ at $n=3$ and $m\\le12$
at $n=4$), the arithmetic $h^*$-polynomial is real-rooted.

\\begin{conjecture}\\label{conj:hstarRR}
The $h^*$-polynomial of the arithmetic matroid of any projected basis
$\\pi([0,1]^n)\\cap\\Lambda$ of the form considered here---in particular,
for every speed vector $v=(v_1,\\dots,v_n)$ with distinct entries and
$\\gcd(v)=1$---is real-rooted.
\\end{conjecture}

The conjecture is stated for the projected-basis zonotopes; the sweeps
suggest it may hold for $(1,\\dots,n-1,m)$ and beyond, and the natural
home for it is the real-rootedness and Lorentzian literature for
arithmetic matroid Tutte specializations
\\cite{BrandenHuh,dAMoci,StanleyEC}."""))
