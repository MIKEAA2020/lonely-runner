#!/usr/bin/env python3
"""Flagship v4 patch: the harmonic exact maximum rho(1,2,3,4,5) = 16/47,
the n=6 boundary rows resolved as strict refutations (m=12: 17/47,
m=9: 455/1269), the two new certificates (sliding-sigma gap + triple
identity), and the balance-family record.  Surgical edits with asserts."""
import re
import sys

BASE = '/home/z/my-project/download/paper2_sources_v4/'


def patch(fname, pairs):
    p = BASE + fname
    src = open(p).read()
    for old, new in pairs:
        if new in src:
            continue        # already applied (idempotent re-run)
        if old not in src:
            print('MISS in %s: %r...' % (fname, old[:70]))
            sys.exit(1)
        if src.count(old) != 1:
            print('NON-UNIQUE in %s: %r...' % (fname, old[:70]))
            sys.exit(1)
        src = src.replace(old, new)
    open(p, 'w').write(src)
    print('patched %s (%d edits)' % (fname, len(pairs)))


# ------------------------------------------------------------------ sec1
sec1 = [
    # 1. abstract: the harmonic sentence -> exact value + n=6 rows
    ("""with exact witnesses $\\rho \\ge 607/1792 > 1/3$ at $v=(1,2,3,4,5)$
(improved from $97/288$ during this paper's audit, which further produced
$\\rho \\ge 4183/12288$ at an independently found rational point, adopted
below), $\\rho \\ge 29/80 > 5/14$ at $v=(1,2,3,4,5,6)$, and $\\rho \\ge
4/11 > 5/14$ at the non-extremal instance $v=(1,2,3,4,5,7)$, while its""",
     """with the extremal harmonic instance's covering radius now certified
exactly---$\\rho(1,2,3,4,5)=16/47>1/3$, the audit's $4183/12288$ witness
superseded by an exact balance family whose three tied dips obey the
weighted identity $32\\,r_{(3,5)}+27\\,r_{(4,5)}+35\\,r_{(3,4)}=32$---alongside
exact witnesses $\\rho \\ge 29/80 > 5/14$ at $v=(1,2,3,4,5,6)$, $\\rho \\ge
4/11 > 5/14$ at the non-extremal instance $v=(1,2,3,4,5,7)$, and, resolved
in this revision, the two former $n=6$ boundary rows as further strict
refutations ($\\rho\\ge17/47>5/14$ at $v=(1,2,3,4,5,12)$ and
$\\rho\\ge455/1269>5/14$ at $v=(1,2,3,4,5,9)$), while its"""),
    # 2. abstract: the ladder sentence -> add the third session
    ("""with a new exact Lipschitz box certificate. Every number cited in this abstract is re-asserted by a persisted,""",
     """with a new exact Lipschitz box certificate; a third session (this
revision) then closed the harmonic instance itself,
$\\rho(1,2,3,4,5)=16/47$ exactly, using two further certificates---a
vertex-checked sliding-$\\sigma$ interval test and a triple-identity
certificate that converts the argmax band's weighted dip identity into a
box proof. Every number cited in this abstract is re-asserted by a persisted,"""),
    # 3. intro: the covering-radius paragraph
    ("""It is false from
$n=5$: at the extremal harmonic instance $v=(1,2,3,4,5)$ there is an
exact rational point at depth $97/288>1/3$, at $v=(1,2,3,4,5,6)$ one at
$29/80>5/14$, and at the non-extremal $v=(1,2,3,4,5,7)$ one at
$4/11>5/14$.""",
     """It is false from
$n=5$: at the extremal harmonic instance $v=(1,2,3,4,5)$ the covering
radius is now known \emph{exactly}, $\\rho=16/47>1/3$
(\\S\\ref{sec:bb}), at $v=(1,2,3,4,5,6)$ there is an exact rational point
at depth $29/80>5/14$, at the non-extremal $v=(1,2,3,4,5,7)$ one at
$4/11>5/14$, and---resolved in this revision---at
$v=(1,2,3,4,5,12)$ and $v=(1,2,3,4,5,9)$ points at $17/47$ and
$455/1269$, both $>5/14$."""),
]

# 4. the v4 what-changed paragraph after the v3 one
sec1.append((
    """All computations are exact (rational arithmetic) unless explicitly marked
as a float enclosure; the one exception class is flagged where it occurs
(numerical root finding is supplemented by exact Sturm-sequence
verification, and the statistical statements carry their sample sizes and
$p$-values).""",
    """\\paragraph{What changed in this revision (v4).}
The third certification session closes the archive's oldest open row:
the deepest hole of the harmonic instance itself. The second attempt's
negative control had left $38{,}117$ survivor boxes around the harmonic
deep holes; their centers, re-used as search seeds, led to a
\\emph{balance family} of witnesses
$c(k)=(2k{+}160,\\,3k{-}13,\\,4k{+}91,\\,5k)/235$, $k=20\\dots35$, at each
point of which the three pair-dips of the runner triple $\\{3,4,5\\}$ tie
at $16/47$---a tie \emph{forced} by the weighted identity
$32\\,r_{(3,5)}+27\\,r_{(4,5)}+35\\,r_{(3,4)}=32$, whose $c$-terms cancel
algebraically. A branch and bound at $r=16/47$ initially stranded on the
band: dyadic boxes cannot align with the family's non-dyadic slice edges,
and the single-$\\sigma$ minimax certificate cannot follow the dip
$\\sigma$ as it slides with $c$ along the band. Two new certificates close
it. The \\emph{sliding-$\\sigma$ gap certificate} observes that for a
fixed integer lift vector $m$, a point $c$ admits a valid $\\sigma$ with
all runner terms $\\le r$ iff an interval of affine forms in $c$
intersects $[0,r]\\cup[1-r,1]$; the interval gap is concave in $c$, so
$2^d$ vertex checks certify the whole box---with $\\sigma$ free to slide.
The \\emph{triple-identity certificate} exhibits three pair-dips of a
runner triple $\\{x,y,z\\}$ with weights $w_{(a,b)}=v_c(v_a{+}v_b)$ whose
weighted values sum, over the box, to a constant $\\le(\\sum w)\\,r$
(the $c$-coefficients cancel exactly); the minimum dip value then bounds
$D$ from above. With lattice-aligned splitting at the dip lattice
($q=470$), the tree empties: $782{,}255$ box evaluations, zero
survivors, $5{,}935$~s on one core; $4{,}841$ sampled prunes
re-confirmed exactly. The same session resolved the two $n=6$ boundary
rows as strict refutations: at $v=(1,2,3,4,5,12)$ the balance family of
the triple $\\{2,5,12\\}$ gives $\\rho\\ge17/47>5/14$ (identity
$84\\,r_{(2,5)}+70\\,r_{(2,12)}+34\\,r_{(5,12)}=68$; witness
$(0,\\tfrac{27}{47},0,\\tfrac{25}{94},\\tfrac{22}{47})$), and at
$v=(1,2,3,4,5,9)$ a four-dip chain balance gives
$\\rho\\ge455/1269>5/14$ (witness
$(\\tfrac{13}{423},\\tfrac{22}{47},\\tfrac{58}{423},\\tfrac{695}{846},
\\tfrac{401}{846})$). The $n=5$ boundary rows $m\\in\\{6,7,8\\}$ remain
open. Nothing else in the paper's claims changes.

All computations are exact (rational arithmetic) unless explicitly marked
as a float enclosure; the one exception class is flagged where it occurs
(numerical root finding is supplemented by exact Sturm-sequence
verification, and the statistical statements carry their sample sizes and
$p$-values)."""))

patch('sec1.tex', sec1)

# ------------------------------------------------------------------ sec3
sec3 = [
    ("""$4183/12288$ at
$c=(\\tfrac{465}{4096},\\tfrac{605}{1024},\\tfrac{1019}{4096},
\\tfrac{315}{4096})$): the deepest hole is not the center, and the
covering-radius statement fails at the very instance where the
$\\delta$-form is tight.""",
     """$4183/12288$ at
$c=(\\tfrac{465}{4096},\\tfrac{605}{1024},\\tfrac{1019}{4096},
\\tfrac{315}{4096})$, and the third certification session closed the row
exactly: $\\rho(1,2,3,4,5)=16/47$ at the balance-family witness
$c=(\\tfrac{46}{47},\\tfrac{92}{235},\\tfrac{231}{235},\\tfrac{35}{47})$,
\\S\\ref{sec:bb}): the deepest hole is not the center, and the
covering-radius statement fails at the very instance where the
$\\delta$-form is tight---with the failure now measured exactly."""),
    ("""LR-zonotope & $\\rho\\le\\mathrm{thr}_n$ (covering) & false from $n=5$;
$h^*$ strictness decays, no measurable signal & $607/1792$,
$4183/12288$, $29/80$, $4/11$ (exact) \\\\""",
     """LR-zonotope & $\\rho\\le\\mathrm{thr}_n$ (covering) & false from $n=5$;
$h^*$ strictness decays, no measurable signal & $16/47$ (exact,
certified), $29/80$, $4/11$, $17/47$, $455/1269$ (exact) \\\\"""),
]
patch('sec3.tex', sec3)

# ------------------------------------------------------------------ sec6
sec6 = [
    ("""tree---with the harmonic lower bound improved twice---to
$\\rho\\ge607/1792$ by the branch-and-bound follow-up and to
$\\rho\\ge4183/12288$ by this paper's audit). And""",
     """tree---with the harmonic instance then closed \emph{exactly} in the
third session: $\\rho(1,2,3,4,5)=16/47$, its balance family and weighted
dip identity recorded in \\S\\ref{sec:bb}, and the two former $n=6$
boundary rows resolved as further strict refutations,
$\\rho\\ge17/47>5/14$ at $m=12$ and $\\rho\\ge455/1269>5/14$ at $m=9$).
And"""),
]
patch('sec6.tex', sec6)

print('PART 1 DONE')

# ------------------------------------------------------------------ sec5
sec5 = [
    ("""$(1,2,3,4,5)$ & $1/3$ & $1/3$ &
$(\\tfrac{465}{4096},\\tfrac{605}{1024},\\tfrac{1019}{4096},
\\tfrac{315}{4096})$ &
$4183/12288{\\approx}.3404$ & yes \\\\""",
     """$(1,2,3,4,5)$ & $1/3$ & $1/3$ &
$(\\tfrac{465}{4096},\\tfrac{605}{1024},\\tfrac{1019}{4096},
\\tfrac{315}{4096})$ &
$4183/12288{\\approx}.3404$ & yes \\\\
$(1,2,3,4,5)$ & $1/3$ & $1/3$ &
$(\\tfrac{46}{47},\\tfrac{92}{235},\\tfrac{231}{235},\\tfrac{35}{47})$ &
$16/47{\\approx}.3404=\\rho$ exactly & yes (certified,
\\emph{this revision}) \\\\"""),
    ("""$(1,2,3,4,5,12)$ & $9/26$ & $5/14$ &
$(\\tfrac1{16},\\tfrac{11}{16},\\tfrac18,\\tfrac7{16},\\tfrac78)$ &
$5/14=\\mathrm{thr}_6$ & no (boundary) \\\\""",
     """$(1,2,3,4,5,12)$ & $9/26$ & $5/14$ &
$(0,\\tfrac{27}{47},0,\\tfrac{25}{94},\\tfrac{22}{47})$ &
$17/47{\\approx}.3617$ & yes (boundary resolved,
\\emph{this revision}) \\\\
$(1,2,3,4,5,9)$ & $1/3$ & $5/14$ &
$(\\tfrac{13}{423},\\tfrac{22}{47},\\tfrac{58}{423},\\tfrac{695}{846},
\\tfrac{401}{846})$ &
$455/1269{\\approx}.3586$ & yes (boundary resolved,
\\emph{this revision}) \\\\"""),
    ("""The rung-2 instance at $m=12$ sits exactly at the threshold:
its row is a \\emph{boundary} row (zero slack, not a strict refutation),
recorded because it already separates the deepest hole from the center
($\\delta=9/26<5/14$). Three boundary rows at $n=5$ ($m=6,7,8$, witness
value exactly $1/3=\\mathrm{thr}_5$ against $\\delta=3/10$) and one at
$n=6$ ($m=9$) are recorded in the same spirit in the program data. The
witnesses are non-torsion in every case: the deepest holes are not the
center classes.""",
     """The rung-2 instance at $m=12$ was recorded in earlier revisions as a
\\emph{boundary} row---its witness sat exactly at the threshold, zero
slack, already separating the deepest hole from the center
($\\delta=9/26<5/14$). This revision resolves both $n=6$ boundary rows
as strict refutations: the balance family of the runner triple
$\\{2,5,12\\}$ gives $\\rho\\ge17/47>5/14$ at $m=12$, and a four-dip
chain witness gives $\\rho\\ge455/1269>5/14$ at $m=9$
(\\S\\ref{sec:bb}). Three boundary rows at $n=5$ ($m=6,7,8$, witness
value exactly $1/3=\\mathrm{thr}_5$ against $\\delta=3/10$) remain open.
The witnesses are non-torsion in every case: the deepest holes are not
the center classes."""),
    ("""At $n=6$, $m=9$ the witness sits at
$\\rho\\ge\\mathrm{thr}_6$ exactly---a boundary row, like the $n=5$
$m=6,7,8$ rows. A useful way to read the table: the $\\rho$-form was
falsified three times over before this session (the three strict rows:
$97/288$, $29/80$, $4/11$), and the session's contribution below is on
the \\emph{positive} side of the ladder---certifying where the deepest
hole genuinely is the center.""",
     """At $n=6$, $m=9$ the recorded witness sat at
$\\mathrm{thr}_6$ exactly---a boundary row, like the $n=5$ $m=6,7,8$
rows---until this revision resolved it, and $m=12$ with it, as strict
refutations. A useful way to read the table: the $\\rho$-form was
falsified three times over before the certification sessions (the three
original strict rows: $97/288$, $29/80$, $4/11$), deepened twice by the
audit ($607/1792$, $4183/12288$), and is now measured \\emph{exactly} at
the harmonic instance ($\\rho=16/47$) with two further strict rows at
$n=6$ ($17/47$, $455/1269$); the sessions' contribution below is on the
\\emph{positive} side of the ladder---certifying where the deepest hole
genuinely is the center---and on the exact measurement of the one
instance where it provably is not."""),
    ("""$r=\\delta=\\tfrac13$---where the truth is $\\rho\\ge97/288>\\tfrac13$---the""",
     """$r=\\delta=\\tfrac13$---where the truth is now known exactly,
$\\rho=16/47>\\tfrac13$ (the third session, below)---the"""),
    ("""The archive's content, stated as a single paragraph: the
covering-radius strengthening of the Lonely Runner statement is false
from $n=5$, with three exact strict witnesses at extremal and
non-extremal instances (plus the audit's improved harmonic record
$\\rho\\ge4183/12288$, superseding $607/1792$), boundary rows at
$n=5$ ($m=6,7,8$) and $n=6$ ($m=9,12$) at exactly the threshold, and a
modular-law extinction record at $n=6$; the equivalent
center-depth form survives untouched (it is $\\lrc$); the deepest-hole
modular law is exact at $n=3$, certified at $n=4$ (both sides, all
tested $m$, including $m=16$ by the new method), false at $n=5$ (the
$k=1,3$ multiples fail exactly; $k=2,4,5,6$ are certified, $m=10$
closing in this revision), and extinct at $n=6$; and every number in the paragraph is""",
     """The archive's content, stated as a single paragraph: the
covering-radius strengthening of the Lonely Runner statement is false
from $n=5$, measured \\emph{exactly} at the extremal instance
($\\rho(1,2,3,4,5)=16/47$) and falsified at four further instances with
exact witnesses ($29/80$ and $4/11$; $17/47$ and $455/1269$ at the two
former boundary rows, resolved in this revision; the superseded search
records $97/288\\to607/1792\\to4183/12288$ retained as provenance),
boundary rows at $n=5$ ($m=6,7,8$) still open, and a
modular-law extinction record at $n=6$; the equivalent
center-depth form survives untouched (it is $\\lrc$); the deepest-hole
modular law is exact at $n=3$, certified at $n=4$ (both sides, all
tested $m$, including $m=16$ by the new method), false at $n=5$ (the
$k=1,3$ multiples fail exactly; $k=2,4,5,6$ are certified, $m=10$
closed in the second session, the harmonic instance itself closed
exactly in the third), and extinct at $n=6$; and every number in the paragraph is"""),
]
patch('sec5.tex', sec5)
print('FLAGSHIP V4 PART 2 DONE')

# ------------------------------------------------------------------ sec5: third-session block
p = BASE + 'sec5.tex'
src = open(p).read()
anchor = """attempt is recorded in \\S\\ref{sec:witness} above; its screen is sound
in the direction that matters (no point with $D>r$ can pass under the
threshold), so the search is a genuine refutation instrument and not
merely reassurance.

\\subsection{Summary of the archive}"""
if '\\paragraph{The third session' not in src:
    assert anchor in src and src.count(anchor) == 1
    newblock = """attempt is recorded in \\S\\ref{sec:witness} above; its screen is sound
in the direction that matters (no point with $D>r$ can pass under the
threshold), so the search is a genuine refutation instrument and not
merely reassurance.

\\paragraph{The third session: the harmonic instance, closed exactly.}
The negative control described above was, all along, the deepest open
object in the archive: the harmonic instance's own covering radius had
only ever been bounded below
($97/288\\to607/1792\\to4183/12288$, each improvement found by search).
The third session closed it from both sides. On the search side, the
$30{,}000$ stored survivor centers of the control run served as seeds;
the summit they converge to is not a point but a \\emph{balance family}
$c(k)=(2k{+}160,\\,3k{-}13,\\,4k{+}91,\\,5k)/235$, $k=20\\dots35$, on
which the three pair-dips of the runner triple $\\{3,4,5\\}$---the
crossings of the speed pairs $(3,5)$, $(4,5)$, $(3,4)$---tie at the
common value $16/47$. The tie is forced: with weights
$w_{(a,b)}=v_c(v_a{+}v_b)$ the dip values obey the weighted identity
\\[
32\\,r_{(3,5)}+27\\,r_{(4,5)}+35\\,r_{(3,4)}\\;=\\;32,
\\]
whose $c$-terms cancel algebraically, so
$\\min r\\le 32/94=16/47$ throughout the family's branch region; at the
family the minimum is attained, $D(c(k))=16/47$ exactly. The witness of
Table~\\ref{tab:refute} is $c(35)=(\\tfrac{46}{47},\\tfrac{92}{235},
\\tfrac{231}{235},\\tfrac{35}{47})$.

On the certificate side, the run at $r=16/47$ exposed why the earlier
cascade could not close a band-type argmax, and two new certificates
were added. The \\emph{sliding-$\\sigma$ gap certificate}: for a fixed
integer lift vector $m$, a point $c$ admits a $\\sigma$ with every
runner term $\\le r$ iff the interval
$[\\max_j (c_j-m_j-r)/v_j,\\ \\min_j (c_j-m_j+r)/v_j]$ intersects
$[0,r]\\cup[1-r,1]$; each endpoint is affine in $c$ and the gap is
concave, so its minimum over a box sits at a vertex---$2^d$ vertex
checks certify the whole box, with $\\sigma$ allowed to \\emph{slide}
with $c$ along the band (the single-$\\sigma$ minimax bound cannot do
this). The \\emph{triple-identity certificate}: exhibit three pair-dips
of a runner triple with weights $w_{(a,b)}=v_c(v_a{+}v_b)$ whose
weighted values sum, over the box, to a constant $S\\le(\\sum w)\\,r$
(the $c$-coefficients cancel exactly, verified symbolically per
candidate); at every $c$ the minimum dip is then $\\le S/\\sum w\\le
r$, and the non-pair terms at each dip are checked $\\le r$ by exact
affine interval arithmetic. Splitting is lattice-aligned at the dip
lattice $q=470$: dyadic midpoints alone can never align with the
family's rational slice edges, and boxes straddling them unpruned are
exactly what stranded the first attempts. With these, the tree empties:
$782{,}255$ box evaluations ($3{,}721$ cheap-grid, $45{,}598$
grid-Lipschitz, $61{,}191$ gap, $197{,}330$ triple-identity,
$83{,}288$ exact-minimax prunes; $391{,}127$ splits; maximum depth
$81$), \\emph{zero} survivors, $5{,}935$~s on one core; $4{,}841$
sampled prunes re-confirmed exactly. Combined with the witness,
\\[
\\rho(1,2,3,4,5)\\;=\\;\\tfrac{16}{47}\\quad\\text{exactly.}
\\]
The soundness profile matches the earlier sessions: adversarial boxes
around the known deeper holes refuse the new certificates at their
thresholds, and the $n=4$ and $m=10$ certificates are untouched.

\\paragraph{The two $n=6$ boundary rows, resolved.}
The same balance machinery, applied at $n=6$, resolves the two boundary
rows of Table~\\ref{tab:refute} as strict refutations. At
$v=(1,2,3,4,5,12)$ the runner triple $\\{2,5,12\\}$ carries the identity
$84\\,r_{(2,5)}+70\\,r_{(2,12)}+34\\,r_{(5,12)}=68$ (the same
cancellation; $\\sum w=188$, tie $68/188=17/47$), with the family
$c(k)=(2k,\\,5k{+}25,\\,12k{+}44)/94$ in the $(c_2,c_5,c_{12})$
coordinates; the witness $(0,\\tfrac{27}{47},0,\\tfrac{25}{94},
\\tfrac{22}{47})$ has its three dips at $\\sigma=17/94$, $77/94$,
$87/94$, each with maximum exactly $34/94=17/47>5/14$
(Appendix~\\ref{sec:self}). At $v=(1,2,3,4,5,9)$ the summit is a
\\emph{four-dip chain}---the pairs $(2,4)$, $(2,5)$, $(4,9)$, $(5,9)$
tie near $0.3586$---and the recorded witness
$(\\tfrac{13}{423},\\tfrac{22}{47},\\tfrac{58}{423},\\tfrac{695}{846},
\\tfrac{401}{846})$ certifies $\\rho\\ge455/1269>5/14$
(Appendix~\\ref{sec:self}). Both searches screened $2^{20}$ random
points, the full dyadic grids, climbs, and Nelder--Mead polishes on the
sound $[0,r]\\cup[1-r,1)$ screen before exact re-evaluation; neither
row is an equality. The $n=5$ boundary rows $m\\in\\{6,7,8\\}$ remain
open; a certification attempt at $m=6$ is recorded in the program data
as incomplete (its tree, at the revised budget, had not emptied when
the session closed).

\\subsection{Summary of the archive}"""
    src = src.replace(anchor, newblock)
    open(p, 'w').write(src)
    print('patched sec5.tex (third-session block)')

# ------------------------------------------------------------------ secC
secC = [
    ("""its discovery is itself a datum for \\S\\ref{sec:bb}'s discussion
of search thoroughness.""",
     """its discovery is itself a datum for \\S\\ref{sec:bb}'s discussion
of search thoroughness.

\\emph{The exact harmonic value $16/47$.} Take
$c=(\\tfrac{46}{47},\\tfrac{92}{235},\\tfrac{231}{235},\\tfrac{35}{47})$
(the $k=35$ point of the balance family; in $235$-ths,
$c=(\\tfrac{230}{235},\\tfrac{92}{235},\\tfrac{231}{235},
\\tfrac{175}{235})$). At $\\sigma=\\tfrac4{235}$:
$\\wnc\\sigma=\\tfrac4{235}$;
$\\wnc{c_2-2\\sigma}=\\tfrac{13}{235}$;
$\\wnc{c_3-3\\sigma}=\\tfrac{80}{235}$;
$\\wnc{c_4-4\\sigma}=\\tfrac{20}{235}$;
$\\wnc{c_5-5\\sigma}=\\tfrac{80}{235}$.
At $\\sigma=\\tfrac{19}{235}$ the terms are
$\\tfrac{19}{235},\\tfrac{43}{235},\\tfrac{35}{235},
\\tfrac{80}{235},\\tfrac{80}{235}$; at $\\sigma=\\tfrac{214}{235}$ they
are $\\tfrac{21}{235},\\tfrac{37}{235},\\tfrac{80}{235},
\\tfrac{80}{235},\\tfrac{45}{235}$. Each of the three maxima is exactly
$\\tfrac{80}{235}=\\tfrac{16}{47}$, with opposite-slope binding pairs
$(+3,-5)$, $(+4,-5)$, $(+3,-4)$: hence $D(c)=16/47$, and the weighted
identity $32\\,r_{(3,5)}+27\\,r_{(4,5)}+35\\,r_{(3,4)}=32$ shows no
point of the branch region exceeds it; the branch-and-bound of
\\S\\ref{sec:bb} certifies the matching upper bound.

\\emph{The $17/47$ witness at $m=12$.} Take $v=(1,2,3,4,5,12)$ and
$c=(0,\\tfrac{27}{47},0,\\tfrac{25}{94},\\tfrac{22}{47})$. At
$\\sigma=\\tfrac{17}{94}$ the six terms (over $94$) are
$17,34,3,26,34,28$; at $\\sigma=\\tfrac{77}{94}$ they are
$17,34,11,26,16,34$; at $\\sigma=\\tfrac{87}{94}$ they are
$7,14,19,28,34,34$. Each maximum is exactly
$\\tfrac{34}{94}=\\tfrac{17}{47}>\\tfrac5{14}$, with binding pairs
$(+2,-5)$, $(+2,-12)$, $(+5,-12)$; the identity
$84\\,r_{(2,5)}+70\\,r_{(2,12)}+34\\,r_{(5,12)}=68$ forces the tie value
$68/188=17/47$ on the family.

\\emph{The $455/1269$ witness at $m=9$.} Take $v=(1,2,3,4,5,9)$ and
$c=(\\tfrac{13}{423},\\tfrac{22}{47},\\tfrac{58}{423},
\\tfrac{695}{846},\\tfrac{401}{846})$. At $\\sigma=\\tfrac{247}{1269}$
the six terms are $\\tfrac{247}{1269},\\tfrac{455}{1269},
\\tfrac{49}{423},\\tfrac{455}{1269},\\tfrac{385}{2538},
\\tfrac{5}{18}$: two terms tie at the maximum
$\\tfrac{455}{1269}>\\tfrac5{14}$ with opposite effective slopes, so
$D(c)=455/1269$."""),
]
patch('secC.tex', secC)

# ------------------------------------------------------------------ main + cover
main = [
    ("""  pdftitle={Every Natural Strengthening Examined Is False: Type Mismatch
            in Lossless Reformulations of the Lonely Runner Conjecture
            (v3, the m=10 row certified)},""",
     """  pdftitle={Every Natural Strengthening Examined Is False: Type Mismatch
            in Lossless Reformulations of the Lonely Runner Conjecture
            (v4, the harmonic instance certified: rho=16/47)},"""),
]
patch('main.tex', main)

cover = [
    ("""(4183/12288&thinsp;&gt;&thinsp;1/3 at the harmonic instance,
        from n&thinsp;=&thinsp;5 on). Revision (v3): the m&thinsp;=&thinsp;10
        row certified&mdash;&rho;(1,2,3,4,10)&thinsp;=&thinsp;7/22 exactly, by
        an audit&ndash;driven second attempt pairing a sound search screen
        with an exact Lipschitz box certificate. With the
        program&rsquo;s bounded positive theorem and a certified
        deepest&#8209;hole archive.</div>""",
     """(16/47&thinsp;&gt;&thinsp;1/3 at the harmonic instance&mdash;measured
        exactly, from n&thinsp;=&thinsp;5 on). Revision (v4): the harmonic
        instance certified&mdash;&rho;(1,2,3,4,5)&thinsp;=&thinsp;16/47
        exactly, by a balance&#8209;family witness whose three tied dips
        obey a weighted identity, and a branch&#8209;and&#8209;bounds tree emptied
        by two new certificates (a sliding&#8209;&sigma; interval test and
        a triple&#8209;identity certificate); the two n&thinsp;=&thinsp;6
        boundary rows resolved as further strict refutations
        (17/47 and 455/1269&thinsp;&gt;&thinsp;5/14). With the
        program&rsquo;s bounded positive theorem and a certified
        deepest&#8209;hole archive.</div>"""),
]
patch('cover.html', cover)

print('FLAGSHIP V4 PATCH COMPLETE')
