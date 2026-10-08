r"""Companion v2 part 4: T=10 extrapolation wording, honest-tone reduction,
closed-by-theorem narrowing, sec9 placement table, bib fixes, main.tex."""
P = []

# ---- secsc: T=10 "predicted" -> "extrapolated"
P.append(('secsc.tex',
r"""and unmodeled. The $T=10$ saturation (predicted $\ge 0.9$ at the
$40\%$ budget) is \emph{not measured}.""",
r"""and unmodeled. The $T=10$ saturation (extrapolated $\ge 0.9$ at the
$40\%$ budget from the rung-monotone pattern---two primes, unmodeled
$\epsilon$-modulation, no error control) is \emph{not measured}."""))

# ---- honest-tone reduction: keep load-bearing, neutralize the rest
P.append(('sec1.tex',
r"""\item \textbf{The honest measurement} (\S\ref{sec:measure}). On all""",
r"""\item \textbf{The measurement} (\S\ref{sec:measure}). On all"""))

P.append(('sec5.tex',
r"""\subsection{Reading the decay honestly}""",
r"""\subsection{Interpreting the decay curve}"""))

P.append(('sec5.tex',
r"""report the decline in both directions of honesty---what it does not""",
r"""report the decline in both directions---what it does not"""))

P.append(('secscale.tex',
r"""data; after five iterations of the pattern the honest prior is that""",
r"""data; after five iterations of the pattern the defensible prior is that"""))

P.append(('secscale.tex',
r"""classes is itself a lemma---the honest shape of the remaining distance.""",
r"""classes is itself a lemma---the precise shape of the remaining distance."""))

# ---- closed-by-theorem narrowing at the two broad sites
P.append(('sec1.tex',
r"""pre-committed success criterion and closed by theorem""",
r"""pre-committed success criterion and closed by theorem (for automata
built on the cascade's bijective drift)"""))

P.append(('sec1.tex',
r"""  closed by theorem, not by measurement: the constrained transfer""",
r"""  closed by theorem for the lossless state space, not by measurement:
  the constrained transfer"""))

# ---- sec9: placement counts as a compact display
P.append(('sec9.tex',
r"""where $q=(p-1)/2$ is the antipode and $t=\floor{(p+1)/3}$ (equivalently
$t\equiv\pm3^{-1}\bmod p$, $3t\in\{p-1,p+1\}$). The placement counts are
independent of $p$ and depend only on the class:
$48$ for $(1,1,1)$, $24$ for each other family at $\epsilon\in\{3,7\}$
(total $288$); $528$, $232$, $216$ at $\epsilon\in\{1,5\}$ (total
$2{,}784$). No pairwise $\pm$-distinct covering occurs at any $k\ge2$""",
r"""where $q=(p-1)/2$ is the antipode and $t=\floor{(p+1)/3}$ (equivalently
$t\equiv\pm3^{-1}\bmod p$, $3t\in\{p-1,p+1\}$). The placement counts are
independent of $p$ and depend only on the class and the $\epsilon$-class
of $p$:
\begin{center}
\small
\begin{tabular}{@{}lccc@{}}
\toprule
$\epsilon$ class & families & placements per family & total \\
\midrule
$\epsilon\in\{3,7\}$ & $(1,1,1)$ & $48$ & \multirow{2}{*}{$288$} \\
                     & each other & $24$ & \\
$\epsilon\in\{1,5\}$ & three classes & $528$, $232$, $216$ & $2{,}784$ \\
\bottomrule
\end{tabular}
\end{center}
\noindent No pairwise $\pm$-distinct covering occurs at any $k\ge2$"""))

# ---- bib: same citation repairs as the flagship
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
