"""Reproduce the claim's own statistic with my pipeline (line y = mean box centre over the whole line,
foot = P90 of last-line bottom), then give it (i) the claim-style null (foot treated as known) and
(ii) the like-for-like null (same statistic recomputed on greedy-packed pages)."""
import numpy as np
from parse import load, paragraphs
from geom import page_geometry
P = load()
def pagedata(pg):
    out, bl = page_geometry(pg, P[pg])
    ys = []
    for o in out:
        if o['bi'] is None: ys.append(None); continue
        l = bl[o['bi']]; ys.append(float(np.mean([q[2] + q[4] / 2 for q in l])))
    for i, y in enumerate(ys):
        if y is None: ys[i] = ys[i - 1] if i else ys[i + 1]
    ys = np.array(ys)
    plen = [len(p['lines']) for p in paragraphs(P[pg])]
    pitch = float(np.median(np.diff(np.sort(ys))))
    return dict(pg=pg, I=P[pg]['meta']['I'], Q=P[pg]['meta']['Q'], ys=ys, n=len(ys), plen=plen, pitch=pitch, top=ys.min() - pitch / 2, bottom=ys.max() + pitch / 2)
D = [pagedata(pg) for pg, d in P.items() if d['meta'].get('I') in ('S', 'B', 'T') and any(L['typ'][0] == 'P' for L in d['loci'])]

def stat(bottom, pitch, q=90, tol=0.5):
    ref = np.quantile(bottom, q / 100); return np.mean((ref - bottom) / pitch < tol), ref

def analyse(name, pages, nsim=4000, seed=0, slack=0.3, tol=0.5):
    rng = np.random.default_rng(seed)
    bottom = np.array([p['bottom'] for p in pages]); pitch = np.array([p['pitch'] for p in pages]); top = np.array([p['top'] for p in pages])
    obs, ref = stat(bottom, pitch, tol=tol)
    plen = np.array([L for p in pages for L in p['plen']], float)
    cap = (ref - top) / pitch                      # line slots to the observed foot
    a = np.zeros(nsim); b = np.zeros(nsim)
    for t in range(nsim):
        nb = np.zeros(len(pages))
        for i in range(len(pages)):
            used = 0
            while True:
                L = plen[rng.integers(len(plen))]
                if used + L <= cap[i] + slack: used += L
                else: break
            nb[i] = used
        gap = cap - nb
        a[t] = np.mean(gap < tol)                              # claim-style: foot known
        b[t], _ = stat(top + nb * pitch, pitch, tol=tol)       # like-for-like: foot re-estimated
    print(f'{name:28s} n={len(pages)} obs {obs:.0%} ({int(round(obs*len(pages)))}) | claim-style null {a.mean():.0%}+-{a.std():.0%} p={(np.sum(a>=obs)+1)/(nsim+1):.4f} | like-for-like null {b.mean():.0%}+-{b.std():.0%} p={(np.sum(b>=obs)+1)/(nsim+1):.4f}')

for slack in (0.3, 0.0):
    print('packing slack', slack, 'lines')
    analyse('S (stars/recipes)', [p for p in D if p['I'] == 'S'], slack=slack)
    analyse('S+B+T', D, slack=slack)
    analyse('B+T', [p for p in D if p['I'] != 'S'], slack=slack)
    analyse('quire T only', [p for p in D if p['Q'] == 'T'], slack=slack)
