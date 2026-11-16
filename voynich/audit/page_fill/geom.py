"""Attach pixel geometry (from voynichese word boxes) to ZL3b paragraph lines.
Own implementation: segment boxes into lines by x-reset, then monotone DP
alignment of box-lines to ZL P-type loci on string similarity."""
import difflib
import numpy as np
from parse import load, paragraphs, boxes


def seglines(b):
    lines = []
    cur = [b[0]]
    for prev, bx in zip(b, b[1:]):
        # new line when the word starts left of the previous word's middle
        # or jumps down by more than a line without moving right much
        if bx[1] < prev[1] + 0.5 * prev[3]:
            lines.append(cur)
            cur = [bx]
        else:
            cur.append(bx)
    lines.append(cur)
    return lines


def sim(a, b):
    return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()


def align(boxlines, zl):
    """Needleman-Wunsch: boxlines (list of word lists) vs zl loci (list of word lists).
    Returns list of (i_box or None) for each zl locus."""
    A = ['.'.join(w[0] for w in l) for l in boxlines]
    B = ['.'.join(L['words']) for L in zl]
    n, m = len(A), len(B)
    S = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            S[i, j] = sim(A[i], B[j])
    gap = -0.15
    D = np.zeros((n + 1, m + 1))
    D[:, 0] = gap * np.arange(n + 1)
    D[0, :] = gap * np.arange(m + 1)
    T = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            c = (D[i - 1, j - 1] + S[i - 1, j - 1] - 0.35, D[i - 1, j] + gap, D[i, j - 1] + gap)
            T[i, j] = int(np.argmax(c))
            D[i, j] = c[T[i, j]]
    i, j = n, m
    out = [None] * m
    sc = [0.0] * m
    while i > 0 and j > 0:
        if T[i, j] == 0:
            out[j - 1] = i - 1
            sc[j - 1] = S[i - 1, j - 1]
            i -= 1
            j -= 1
        elif T[i, j] == 1:
            i -= 1
        else:
            j -= 1
    return out, sc


def line_geom(l):
    """Geometry of a box-line. y is taken at the left end (lines slope)."""
    x0 = min(b[1] for b in l)
    x1 = max(b[1] + b[3] for b in l)
    k = l[:4]
    yc = float(np.median([b[2] + b[4] / 2 for b in k]))
    ytop = float(np.median([b[2] for b in k]))
    ybot = float(np.median([b[2] + b[4] for b in k]))
    return dict(x0=x0, x1=x1, yc=yc, ytop=ytop, ybot=ybot,
                ymax=max(b[2] + b[4] for b in l), ymin=min(b[2] for b in l), nw=len(l))


def page_geometry(pg, page):
    """Return list of P-loci dicts with geometry (or None), in ZL order."""
    b = boxes(pg)
    if b is None:
        return None
    bl = seglines(b)
    zl = [L for L in page['loci'] if L['typ'][0] == 'P']
    idx, sc = align(bl, zl)
    out = []
    for L, i, s in zip(zl, idx, sc):
        g = line_geom(bl[i]) if i is not None else None
        out.append(dict(L=L, g=g, score=s, bi=i))
    return out, bl


if __name__ == '__main__':
    P = load()
    for pg, d in P.items():
        if d['meta'].get('I') not in ('S', 'B', 'T'):
            continue
        r = page_geometry(pg, d)
        if r is None:
            continue
        out, bl = r
        miss = [o['L']['num'] for o in out if o['g'] is None]
        low = [(o['L']['num'], round(o['score'], 2)) for o in out if o['g'] is not None and o['score'] < 0.6]
        ys = [o['g']['yc'] for o in out if o['g']]
        dy = np.diff(ys)
        print(pg, d['meta'].get('I'), 'P', len(out), 'boxlines', len(bl), 'miss', miss, 'low', low,
              'nonmono', int((dy <= 0).sum()), 'pitch med %.1f' % np.median(dy) if len(dy) else '')
