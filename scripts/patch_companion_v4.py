#!/usr/bin/env python3
"""Companion v4 patch: the harmonic exact value + the n=6 boundary rows
resolved, in the zonotope-strand records and the revision note."""
import sys

BASE = '/home/z/my-project/download/paper_sources_v4/'


def patch(fname, pairs):
    p = BASE + fname
    src = open(p).read()
    n = 0
    for old, new in pairs:
        if new in src:
            continue
        if old not in src:
            print('MISS in %s: %r...' % (fname, old[:70]))
            sys.exit(1)
        src = src.replace(old, new)
        n += 1
    open(p, 'w').write(src)
    print('patched %s (%d edits)' % (fname, n))


secsc = [
    ("""\\emph{The third revision (v3) changes no theorem, lemma, or hypothesis
of this monograph}: it updates the zonotope-strand record of
\\S\\ref{sec:scale} for the companion paper's closure of the last open
certification row (the $n=5$, $m=10$ deepest-hole value, now exactly
$7/22$ by a Lipschitz box certificate), harmonizes the harmonic
lower-bound record with the companion's audited $4183/12288$ witness,
and records the version marker.""",
     """\\emph{The third revision (v3) changes no theorem, lemma, or hypothesis
of this monograph}: it updates the zonotope-strand record of
\\S\\ref{sec:scale} for the companion paper's closure of the last open
certification row (the $n=5$, $m=10$ deepest-hole value, now exactly
$7/22$ by a Lipschitz box certificate), harmonizes the harmonic
lower-bound record with the companion's audited $4183/12288$ witness,
and records the version marker. \\emph{The fourth revision (v4) again
changes no theorem, lemma, or hypothesis}: it updates the
zonotope-strand record for the companion's third certification
session---the harmonic instance's deepest hole closed \\emph{exactly}
($\\rho(1,2,3,4,5)=16/47$, by a balance-family witness whose three tied
pair-dips obey a weighted identity, certified by branch and bound with
two new box certificates), and the two $n=6$ boundary rows resolved as
strict refutations ($\\rho\\ge17/47>5/14$ at $m=12$,
$\\rho\\ge455/1269>5/14$ at $m=9$)---and records the version marker."""),
]
patch('secsc.tex', secsc)

secscale = [
    ("""Its geometric
reading is false: the deepest-hole form fails from $n=5$ with exact
witnesses (the harmonic lower bound improved to $\\rho\\ge 4183/12288$),""",
     """Its geometric
reading is false: the deepest-hole form fails from $n=5$ with exact
witnesses (the harmonic instance closed \\emph{exactly} in the fourth
revision, $\\rho(1,2,3,4,5)=16/47$; the two $n=6$ boundary rows
resolved as strict refutations, $\\rho\\ge17/47$ at $m=12$ and
$\\rho\\ge455/1269$ at $m=9$),"""),
]
patch('secscale.tex', secscale)

sec6 = [
    ("""  onward with exact rational witnesses (at the harmonic $n=5$ instance
  $m=5$, $\\rho\\ge 4183/12288>\\tfrac13$, improved twice---from
  $97/288$ to $607/1792$ by the first certification session's
  follow-up and to $4183/12288$ by the companion paper's external
  audit), the modular law that would""",
     """  onward with exact rational witnesses (at the harmonic $n=5$ instance
  $m=5$, $\\rho=16/47>\\tfrac13$ \\emph{exactly}, closed in the companion
  paper's fourth revision by a balance-family witness and two new box
  certificates---the search record
  $97/288\\to607/1792\\to4183/12288$ retained as provenance; the two
  $n=6$ boundary rows resolved as strict refutations, $17/47$ and
  $455/1269$), the modular law that would"""),
]
patch('sec6.tex', sec6)

cover = [
    ("""        repair, quantified hypotheses, scope levels fixed. Revision
        (v3): zonotope&ndash;strand record updated&mdash;the last open""",
     """        repair, quantified hypotheses, scope levels fixed. Revision
        (v4): zonotope&ndash;strand record updated&mdash;the harmonic
        instance's deepest hole closed exactly (&rho;&thinsp;=&thinsp;16/47),
        the two n&thinsp;=&thinsp;6 boundary rows resolved as strict
        refutations. Revision
        (v3): the last open"""),
]
patch('cover.html', cover)

print('COMPANION V4 PATCH COMPLETE')
