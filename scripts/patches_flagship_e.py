r"""Bib fixes + main.tex update for flagship v2."""
P = []

# ---- bib.tex: fix Lam-Leung (split into two real papers), remove phantom
# de Bruijn-Schoenberg, add Cusick, Mirsky-Newman caveat
P.append(('bib.tex',
r"""\bibitem{LamLeung}
T.~Y.~Lam and K.~H.~Leung,
\emph{On the cyclotomic polynomial $\Phi_{pq}(X)$},
Amer.\ Math.\ Monthly \textbf{107} (2000), no.~7, 605--612.

\bibitem{MirskyNewman}
L.~Mirsky and D.~J.~Newman,
\emph{A problem in the theory of numbers},
Michigan Math.\ J.\ \textbf{9} (1962), 201--206; see also
H.~Davenport and R.~Rado, \emph{Covering systems of negative
congruences},
J.\ London Math.\ Soc.\ \textbf{38} (1963), 509--516.

\bibitem{deBruijnSchoenberg}
N.~G.~de~Bruijn and T.~J.~Schoenberg,
\emph{On the two vanishing identities of Euler},
Indag.\ Math.\ \textbf{23} (1961), 542--547.

\bibitem{BrandenHuh}""",
r"""\bibitem{LamLeung}
T.~Y.~Lam and K.~H.~Leung,
\emph{On the cyclotomic polynomial $\Phi_{pq}(X)$},
Amer.\ Math.\ Monthly \textbf{103} (1996), 562--564; and
\emph{On vanishing sums of roots of unity},
J.\ Algebra \textbf{224} (2000), 91--109.

\bibitem{MirskyNewman}
L.~Mirsky and D.~J.~Newman,
unpublished (the distinct-moduli theorem for exact covering systems);
see the account in
H.~Davenport and R.~Rado,
\emph{Covering systems of negative congruences},
J.\ London Math.\ Soc.\ \textbf{38} (1963), 509--516.

\bibitem{BrandenHuh}"""))

P.append(('bib.tex',
r"""\bibitem{Wills}
J.~M. Wills,
\emph{Zwei S\"atze \"uber inhomogene diophantische Approximation von
Irrationalzahlen},
Monatsh.\ Math.\ 71 (1967).

\bibitem{Rosenfeld8}""",
r"""\bibitem{Wills}
J.~M. Wills,
\emph{Zwei S\"atze \"uber inhomogene diophantische Approximation von
Irrationalzahlen},
Monatsh.\ Math.\ 71 (1967).

\bibitem{Cusick1973}
T.~W. Cusick,
\emph{View-obstruction problems},
Aequationes Math.\ \textbf{9} (1973), 165--170.

\bibitem{Rosenfeld8}"""))

# ---- main.tex: title (v2), pdftitle, include secC
P.append(('main.tex',
r"""\hypersetup{
  pdftitle={Every Natural Strengthening Is False: Type Mismatch in
            Lossless Reformulations of the Lonely Runner Conjecture},""",
r"""\hypersetup{
  pdftitle={Every Natural Strengthening Examined Is False: Type Mismatch
            in Lossless Reformulations of the Lonely Runner Conjecture
            (v2, audit revision)},"""))

P.append(('main.tex',
r"""\input{sec6}
\input{secA}
\input{secB}

\input{bib}""",
r"""\input{sec6}
\input{secA}
\input{secB}
\input{secC}

\input{bib}"""))
