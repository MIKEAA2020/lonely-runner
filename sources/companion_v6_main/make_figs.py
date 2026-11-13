#!/usr/bin/env python3
"""Figures for the pair-sum LRC paper.
fig_decay.png    -- proved-coverage decay curves (both batteries).
fig_rigidzone.png-- the rigid-zone map with verified ranges.
English labels (paper language), restrained academic palette.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
fm.fontManager.addfont('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

OUT = '/home/z/my-project/scripts/paper'

# ------------------------------------------------------------------ fig 1
# Proved coverage of the n=5 corpora (primitive 5-sets of [1,V]).
V = [16, 24, 32, 40, 48, 56]
same_ports = [96.94, 85.89, 79.52, 74.93, 72.10, 71.13]
same_d3 = [97.26, 87.29, 82.13, 78.03, 75.25, 74.28]
aug_ports = [None, 92.69, 88.10, 82.13, 77.79, 75.57]
aug_d3 = [None, 93.44, 89.64, 84.29, 80.21, 78.16]

fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
ax.plot(V, same_ports, 'o-', color='#7A9AB5', lw=1.6, ms=5,
        label='same battery, ports (P1+P2+P4 at $\\{7,13,17\\}$)')
ax.plot(V, same_d3, 'o--', color='#7A9AB5', lw=1.3, ms=4, mfc='white',
        label='same battery, +D3')
Va = [v for v, x in zip(V, aug_ports) if x is not None]
aug_p = [x for x in aug_ports if x is not None]
aug_d = [x for x in aug_d3 if x is not None]
ax.plot(Va, aug_p, 's-', color='#B5651D', lw=1.6, ms=5,
        label='augmented battery, ports (P4 at $\\{7,13,17,19,37\\}$)')
ax.plot(Va, aug_d, 's--', color='#B5651D', lw=1.3, ms=4, mfc='white',
        label='augmented battery, +D3')
ax.axhline(80, color='#999999', lw=0.9, ls=':')
ax.annotate('80% continuation bar', xy=(16.4, 80.35), fontsize=8.5,
            color='#777777')
ax.axhline(60, color='#999999', lw=0.9, ls=':')
ax.annotate('60% pivot bar', xy=(16.4, 60.35), fontsize=8.5, color='#777777')
ax.annotate('V=48: 80.21% (marginal pass)', xy=(48, 80.21), xytext=(38.2, 85.6),
            fontsize=8.5, color='#B5651D',
            arrowprops=dict(arrowstyle='-', color='#B5651D', lw=0.8))
ax.annotate('V=56: 78.16%', xy=(56, 78.16), xytext=(45.5, 68.0),
            fontsize=8.5, color='#B5651D',
            arrowprops=dict(arrowstyle='-', color='#B5651D', lw=0.8))
ax.set_xlabel('speed bound $V$ (corpus = primitive 5-subsets of $[1,V]$)')
ax.set_ylabel('sets closed by proved lemmas (%)')
ax.set_title('Proved-coverage decline, $n=5$ (corpus sizes 4{,}311 to'
             ' 3{,}712{,}576)', fontsize=11)
ax.set_xlim(14, 58)
ax.set_ylim(55, 100)
ax.set_xticks(V)
ax.grid(True, lw=0.4, alpha=0.35)
ax.legend(fontsize=8.2, loc='lower left', framealpha=0.9)
fig.savefig(OUT + '/fig_decay.png', dpi=300)
plt.close(fig)

# ------------------------------------------------------------------ fig 2
# Rigid-zone map: unit-covering-capable moduli at n=5 (k=5 speeds, T=6)
# and n=4 (k=4, T=5), with verified ranges.
fig, ax = plt.subplots(figsize=(7.2, 2.6), constrained_layout=True)
cap5 = [7, 13, 17, 19, 37]
cap4 = [7, 11, 13]
ax.axvspan(0, 150, color='#3D5A80', alpha=0.08)
ax.axvspan(150, 160, color='#BBBBBB', alpha=0.15)
for N in cap5:
    ax.plot([N], [1], marker='o', ms=9, color='#B5651D')
    ax.annotate(str(N), xy=(N, 1), xytext=(0, 9),
                textcoords='offset points', ha='center', fontsize=9,
                color='#B5651D', fontweight='bold')
for N in cap4:
    ax.plot([N], [0], marker='s', ms=8, color='#3D5A80')
    ax.annotate(str(N), xy=(N, 0), xytext=(0, -16),
                textcoords='offset points', ha='center', fontsize=9,
                color='#3D5A80', fontweight='bold')
ax.annotate('micro-cases $7\\cdot\\{7,11,13\\}=\\{49,77,91\\}$:'
            ' verified safe, obstruction open',
            xy=(70, 1.02), fontsize=8, color='#666666', ha='center')
ax.annotate('composite scan $\\leq 111$; prime scan $\\leq 150$'
            ' (open beyond)',
            xy=(75, -0.14), fontsize=8, color='#888888', ha='center')
ax.set_yticks([0, 1])
ax.set_yticklabels(['$n=4$  ($T=5$, 3 bad sets)', '$n=5$  ($T=6$, 4 bad sets)'],
                   fontsize=9)
ax.set_xlim(0, 160)
ax.set_ylim(-0.35, 1.35)
ax.set_xlabel('modulus $N$')
ax.set_title('Unit-covering-capable moduli (the rigid zone) in the verified'
             ' range', fontsize=11)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig(OUT + '/fig_rigidzone.png', dpi=300)
plt.close(fig)
print('figures written to', OUT)
