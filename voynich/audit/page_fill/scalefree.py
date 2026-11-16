"""Scale-free text-block geometry: height of text block in units of its own width
(robust to per-page image scaling), and in units of mean character width."""
import numpy as np
from parse import load, boxes, paragraphs
from geom import seglines, page_geometry
from fill import R, GS, GT20, G50

P = load()
def blockdims(pg):
    out, bl = page_geometry(pg, P[pg])
    ls = [bl[o['bi']] for o in out if o['bi'] is not None]
    x0 = np.array([min(q[1] for q in l) for l in ls]); x1 = np.array([max(q[1] + q[3] for q in l) for l in ls])
    L = np.median(x0); Rr = np.percentile(x1, 90)
    return float(L), float(Rr)

for pg in G50:
    L, Rr = blockdims(pg)
    r = R[pg]
    r['W'] = Rr - L; r['xL'] = L; r['xR'] = Rr
    r['A'] = (r['y1'] - r['y0']) / r['W']           # text height (first-to-last line) / text width
    r['Hc'] = (r['y1'] - r['y0']) / r['cw']         # text height in character widths
    r['Wc'] = r['W'] / r['cw']

if __name__ == '__main__':
    print('page   n  pitch   W     A=H/W  lines-equivalent(A*W/pitch)  H/cw   W/cw  pitch/cw  pitch/W*100')
    for pg in GT20:
        r = R[pg]
        print(f"{pg:6s} {r['n']:3d} {r['pitch']:5.1f} {r['W']:6.1f} {r['A']:6.3f}   {r['Hc']:6.1f} {r['Wc']:5.1f}  {r['mp']/r['cw']:5.2f}  {100*r['mp']/r['W']:5.2f}")
    A = np.array([R[pg]['A'] for pg in GT20]); print('A sorted:', np.round(np.sort(A), 3))
    W = np.array([R[pg]['W'] for pg in GT20]); print('W: mean %.1f sd %.1f (%.1f%%)' % (W.mean(), W.std(), 100 * W.std() / W.mean()))
    Wc = np.array([R[pg]['Wc'] for pg in GT20]); print('W/cw: mean %.1f sd %.1f (%.1f%%)' % (Wc.mean(), Wc.std(), 100 * Wc.std() / Wc.mean()))
    n = np.array([R[pg]['n'] for pg in GT20]); mp = np.array([R[pg]['mp'] for pg in GT20]); cw = np.array([R[pg]['cw'] for pg in GT20])
    print('corr(n lines, pitch/W) = %.2f ; corr(n, pitch/cw) = %.2f ; corr(n, cw/W) = %.2f ; corr(n, A) = %.2f' % (
        np.corrcoef(n, mp / W)[0, 1], np.corrcoef(n, mp / cw)[0, 1], np.corrcoef(n, cw / W)[0, 1], np.corrcoef(n, A)[0, 1]))
    print('CV of n = %.3f ; CV of pitch/W = %.3f ; CV of A = %.3f' % (n.std() / n.mean(), (mp / W).std() / (mp / W).mean(), A.std() / A.mean()))
