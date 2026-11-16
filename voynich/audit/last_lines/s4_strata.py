"""Observed vs null (page-pool renewal) by section, hand, layout, paragraph length."""
import numpy as np, sys
from collections import defaultdict, Counter
import zl

rng = np.random.default_rng(7)
rows = zl.read()
paras, an = zl.paragraphs(rows)
pagepool = defaultdict(list)
P = []
for p in paras:
    L = p['lines']
    wc = [len(zl.words(l['raw'])) for l in L]
    P.append(dict(wc=wc, page=L[0]['page'], I=L[0]['I'], H=L[0]['H'], Lg=L[0]['L'],
                  cut=any('<->' in l['raw'] for l in L),
                  cutlast=('<->' in L[-1]['raw']),
                  special=L[-1]['typ'] in ('Pt', 'Pc', 'Pr'),
                  lastpara_on_page=False, n=len(L)))
    pagepool[L[0]['page']] += wc[:-1]
# mark last paragraph on page
lastidx = {}
for i, p in enumerate(P):
    lastidx[p['page']] = i
for i in lastidx.values():
    P[i]['lastpara_on_page'] = True


def null_last(p):
    N = sum(p['wc']); src = pagepool[p['page']]
    rem = N
    while True:
        c = src[rng.integers(len(src))]
        if rem <= c:
            return rem
        rem -= c


R = 400
sims = np.array([[null_last(p) for p in P] for _ in range(R)])  # R x nP


def report(name, mask):
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return
    obs1 = sum(P[i]['wc'][-1] == 1 for i in idx); obs2 = sum(P[i]['wc'][-1] <= 2 for i in idx)
    s1 = (sims[:, idx] == 1).sum(1); s2 = (sims[:, idx] <= 2).sum(1)
    mw = np.mean([np.median(P[i]['wc'][:-1]) for i in idx])
    print(f'{name:38s} n={len(idx):4d} med w/line={mw:4.1f} | 1w obs={obs1:3d} ({100*obs1/len(idx):4.1f}%) null={s1.mean():5.1f} [{np.percentile(s1,2.5):.0f},{np.percentile(s1,97.5):.0f}] '
          f'| <=2w obs={obs2:3d} ({100*obs2/len(idx):4.1f}%) null={s2.mean():5.1f} [{np.percentile(s2,2.5):.0f},{np.percentile(s2,97.5):.0f}]')


n = np.array([p['n'] for p in P])
m4 = n >= 4
report('ALL 4+', m4)
report('ALL 3+', n >= 3)
report('2-3 lines only', (n <= 3))
report('2 lines only', (n == 2))
report('3 lines only', (n == 3))
for I in sorted(set(p['I'] for p in P)):
    report(f'section I={I} (4+)', m4 & np.array([p['I'] == I for p in P]))
for H in sorted(set(p['H'] for p in P)):
    report(f'hand H={H} (4+)', m4 & np.array([p['H'] == H for p in P]))
for Lg in sorted(set(str(p['Lg']) for p in P)):
    report(f'language L={Lg} (4+)', m4 & np.array([str(p['Lg']) == Lg for p in P]))
report('para cut by drawing (4+)', m4 & np.array([p['cut'] for p in P]))
report('para not cut (4+)', m4 & np.array([not p['cut'] for p in P]))
report('last para on page (4+)', m4 & np.array([p['lastpara_on_page'] for p in P]))
report('not last para on page (4+)', m4 & np.array([not p['lastpara_on_page'] for p in P]))
for lo, hi in [(4, 5), (6, 8), (9, 99)]:
    report(f'{lo}-{hi} lines', (n >= lo) & (n <= hi))
# word-count bands of typical line
med = np.array([np.median(p['wc'][:-1]) for p in P])
for lo, hi in [(0, 6), (6, 9), (9, 11), (11, 99)]:
    report(f'typical line {lo}-{hi} words (4+)', m4 & (med >= lo) & (med < hi))
