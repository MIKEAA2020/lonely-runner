P = []
# secA: insert figure after the decay paragraph
P.append(('secA.tex',
r"""Both coefficients
were recomputed in the provenance run from the persisted instance data
and reproduce exactly.""",
r"""Both coefficients
were recomputed in the provenance run from the persisted instance data
and reproduce exactly.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.92\textwidth]{fig_nS.png}
\caption{The strictness records of Table~\ref{tab:hstar}. Left: the
plateau of $n\cdot S$ (the harmonic family starts at $36.8$, rung-2 at
$60.5$; both drift to $\approx30$). Right: the excess $S/\mathrm{floor}$
over the Newton floor at the minimizing coefficient---the
family-specific content of the decline---which itself decreases from
$3.06$/$5.04$ at $n=3$ to $1.99$/$2.00$ at $n=10$.}
\label{fig:nS}
\end{figure}"""))

# secB: four scripts, control row, repro commands
P.append(('secB.tex',
r"""Three scripts were written and run for this paper; their outputs are
persisted alongside the inherited records and copied into the paper's
artifact bundle.""",
r"""Four scripts were written and run for this paper (the fourth, the
survivor-enclosure follow-up, was added after the branch-and-bound
session); their outputs are persisted alongside the inherited records
and copied into the paper's artifact bundle."""))

P.append(('secB.tex',
r"""$(1,2,3,4,5)$ & $1/3$ & \emph{refused} (control) & --- &
$45{,}355$ $(421)$ \\""",
r"""$(1,2,3,4,5)$ & $1/3$ & \emph{refused} (control) & $17{,}520$
$(421)$ & $45{,}355$ \\"""))

P.append(('secB.tex',
r"""budgets. The artifact bundle accompanying the paper contains the three
scripts, their outputs, the inherited zonotope records
(\texttt{out\_lrc\_zono\_final.json},
\texttt{out\_lrc\_zono\_rescan.json}, and the library
\texttt{lrc\_zono\_lib.py} they import), and the ledger.""",
r"""budgets. One-line reproduction from the artifact bundle's root:
\texttt{python3 hstar\_n9\_n10.py};
\texttt{python3 certify\_rung2\_n5.py};
\texttt{python3 paper\_provenance.py}; and
\texttt{python3 certify\_followup.py}---each writes its
\texttt{out\_*.json} next to itself, deterministically. The artifact
bundle accompanying the paper contains the four scripts, their outputs,
the inherited zonotope records
(\texttt{out\_lrc\_zono\_final.json},
\texttt{out\_lrc\_zono\_rescan.json}, and the library
\texttt{lrc\_zono\_lib.py} they import), and the ledger."""))

# secB: soften the audit-standard opener per gemini2
P.append(('secB.tex',
r"""The audit standard adopted for this paper is stated in one sentence:
\emph{every number cited in the text traces to a script that was
actually run and whose output was persisted.} It was applied
retroactively to all inherited numbers and prospectively to all new
ones. This appendix is the map.""",
r"""The reproducibility standard of this paper is one sentence:
\emph{every number cited in the text traces to a script that was
actually run and whose output was persisted.} It was applied
retroactively to all inherited numbers and prospectively to all new
ones. This appendix is the map."""))
