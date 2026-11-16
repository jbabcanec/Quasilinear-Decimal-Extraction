"""Page-fill audit: own re-measurement of claim (b)."""
import json, sys
import numpy as np
from table import build

R = build()
PAGES = list(R)
for pg, r in R.items():
    ys = np.array(r['ys'])
    r['y0'] = float(ys.min()); r['y1'] = float(ys.max())
    r['mp'] = (r['y1'] - r['y0']) / (r['n'] - 1)       # mean pitch (extent / (n-1))

def group(sel):
    return [pg for pg in PAGES if sel(R[pg])]

G50 = group(lambda r: True)
GS = group(lambda r: r['meta']['I'] == 'S')
GT20 = group(lambda r: r['meta']['Q'] == 'T')

def observed(pages, foot, tol=0.5, pitch='pitch'):
    """foot: dict page -> y level.  gap in line units."""
    g = np.array([(foot[pg] - R[pg]['y1']) / R[pg][pitch] for pg in pages])
    return g

def pack_null(pages, foot, rng, mode='shuffle', units='px', tol=0.5, nsim=2000, pool=None, pitch='pitch'):
    """Greedy packing of whole paragraphs.
    units='px'   : a paragraph keeps its home-page pixel height (lines x home pitch)
    units='line' : a paragraph takes (lines x pitch of the page it lands on)
    mode='shuffle': random permutation of the group's paragraphs (wraps around if exhausted)
    mode='draw'   : i.i.d. draws from the group's paragraph-length distribution
    mode='shift'  : true manuscript order, random circular start offset
    Returns array (nsim, npages) of gaps in line units of the target page."""
    para = []   # (lines, home pitch)
    for pg in (pool or pages):
        for L in R[pg]['plen']:
            para.append((L, R[pg][pitch]))
    para = np.array(para, float)
    caps = np.array([foot[pg] - R[pg]['y0'] + R[pg][pitch] for pg in pages])   # px room incl. slot of last line
    pit = np.array([R[pg][pitch] for pg in pages])
    out = np.zeros((nsim, len(pages)))
    npar = len(para)
    for s in range(nsim):
        if mode == 'shuffle':
            seq = para[rng.permutation(npar)]; k = 0
        elif mode == 'shift':
            k = rng.integers(npar); seq = para
        for i in range(len(pages)):
            used = 0.0
            first = True
            while True:
                if mode == 'draw':
                    L, hp = para[rng.integers(npar)]
                else:
                    L, hp = seq[k % npar]
                h = L * (hp if units == 'px' else pit[i])
                if used + h <= caps[i] + 1e-9 or first and False:
                    used += h; first = False
                    if mode != 'draw': k += 1
                else:
                    break       # paragraph does not fit: page closed (in 'draw' mode it is discarded = fresh draw next page; see note)
                if used >= caps[i]: break
            out[s, i] = (caps[i] - used) / pit[i]
    return out

def report(name, pages, foot, tol=0.5, nsim=2000, seed=1, modes=(('shuffle','px'),('shuffle','line'),('shift','line'),('draw','line')), pool=None, quiet=False):
    rng = np.random.default_rng(seed)
    g = observed(pages, foot)
    obs = (g < tol).mean()
    res = dict(name=name, n=len(pages), obs=obs, k=int((g < tol).sum()))
    s = f'{name:44s} n={len(pages):2d} obs {int((g<tol).sum()):2d} = {obs:.0%} |'
    for mode, units in modes:
        sim = pack_null(pages, foot, rng, mode=mode, units=units, nsim=nsim, pool=pool)
        fr = (sim < tol).mean(axis=1)
        p = ((fr >= obs - 1e-12).sum() + 1) / (nsim + 1)
        res[mode + '_' + units] = (fr.mean(), fr.std(), p)
        s += f' {mode}/{units}: {fr.mean():.0%}+-{fr.std():.0%} p={p:.4f} |'
    if not quiet: print(s)
    return res

def foot_quantile(pages, q=90, key=None):
    """foot = q-th percentile of last-line y within each key-group (key: page -> group label)."""
    foot = {}
    if key is None:
        F = np.percentile([R[pg]['y1'] for pg in pages], q)
        return {pg: F for pg in pages}
    groups = {}
    for pg in pages: groups.setdefault(key(pg), []).append(pg)
    for k, pgs in groups.items():
        F = np.percentile([R[pg]['y1'] for pg in pgs], q)
        for pg in pgs: foot[pg] = F
    return foot

if __name__ == '__main__':
    print('== M1: claimed definition (foot = P90 of last-line y in group, tol 0.5 line) ==')
    for nm, G in (('50 text pages (S,B,T)', G50), ('25 S pages', GS), ('23 quire-T pages', GT20)):
        report(nm, G, foot_quantile(G, 90))
    print('gaps (lines) S pages:', {pg: round(float(x), 2) for pg, x in zip(GS, observed(GS, foot_quantile(GS, 90)))})
    print('gaps (lines) other pages:', {pg: round(float(x), 2) for pg, x in zip(G50, observed(G50, foot_quantile(G50, 90))) if pg not in GS})
