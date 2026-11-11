import math, sys
sys.path.insert(0, '.')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from family_Gl import G_
X = (0, 0)
def pos(l):
    P = {'N': (0, 0)}
    for k in range(l + 1):
        r = 1.0 + 1.05 * k
        for i in range(4):
            th = math.radians(90 * i - 45 * k + 45 + 45 * l)
            P[(k, i)] = (r * math.cos(th), r * math.sin(th))
    R = 1.0 + 1.05 * l
    P['S'] = (0, -(R + 1.3))
    P['a'] = (R + 3.0, 1.6); P['b'] = (R + 3.0, 0.0); P['c'] = (R + 3.0, -1.6)
    return P
def draw(ax, l):
    G = G_(l); P = pos(l)
    apcol = {'a': '#e08a2c', 'b': '#3a9a3a', 'c': '#888888'}
    for u, v in G.edges():
        if u in apcol or v in apcol:
            ap = u if u in apcol else v
            ax.plot(*zip(P[u], P[v]), color=apcol[ap], lw=0.5, alpha=0.35, zorder=1)
    for u, v in G.edges():
        if u in apcol or v in apcol: continue
        if 'S' in (u, v):
            w = v if u == 'S' else u
            rad = 0.0 if P[w][1] < 0 else (0.75 if P[w][0] > 0 else -0.75)
            ax.add_patch(FancyArrowPatch(P['S'], P[w], arrowstyle='-', connectionstyle='arc3,rad=%.2f' % rad,
                                         color='#333333', lw=1.0, zorder=2))
            continue
        blue = isinstance(u, tuple) and isinstance(v, tuple) and u[0] == v[0]
        ax.plot(*zip(P[u], P[v]), color='#1f5fd6' if blue else '#222222', lw=2.6 if blue else 1.0, zorder=3)
    for v in G:
        col = apcol.get(v, '#1f5fd6' if isinstance(v, tuple) else 'black')
        if v == X: col = '#d62728'
        ax.scatter(*P[v], s=110 if v in apcol else 55, color=col, zorder=5, edgecolors='white', linewidths=0.8)
        lab = {'N': '$N$', 'S': '$S$', 'a': r'$\alpha$', 'b': r'$\beta$', 'c': r'$\gamma$'}.get(v, '%d%d' % v if isinstance(v, tuple) else v)
        if v == X: lab = '$x$=00'
        dx, dy = (0.18, 0.18)
        ax.text(P[v][0] + dx, P[v][1] + dy, lab, fontsize=15 if v in apcol else 10.5, zorder=6,
                fontweight='bold' if v in apcol or v in ('N', 'S') else None)
    ax.set_title('$G_{%d}$: %d vertices, %d edges' % (l, G.number_of_nodes(), G.number_of_edges()), fontsize=14)
    ax.set_aspect('equal'); ax.axis('off')
fig, axs = plt.subplots(1, 2, figsize=(10.5, 5.2), gridspec_kw={'width_ratios': [1, 1.25]})
draw(axs[0], 1); draw(axs[1], 2)

plt.tight_layout(rect=(0, 0, 1, 1)); plt.savefig('fig_Gl.pdf', bbox_inches='tight'); plt.savefig('fig_Gl.png', dpi=130)
