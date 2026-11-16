"""Real manuscripts (Gallicorpora 15th c.): what do 'lines that stop early' hold,
and what do genuine prose paragraph ends (short line after >=3 full lines) hold?"""
import numpy as np
from collections import defaultdict, Counter
import alto

rows = []  # per line records
for ms, fn in alto.pages():
    for b in alto.read_page(fn):
        if not (b['label'] and b['label'].startswith('MainZone')):
            continue
        L = [l for l in b['lines'] if l['ltype'] == 'DefaultLine' and l['text'].strip()]
        if len(L) < 6:
            continue
        L = sorted(L, key=lambda l: l['y'])
        left = np.percentile([l['x0'] for l in L], 25)
        right = np.percentile([l['x1'] for l in L], 80)
        W = right - left
        r = np.array([(l['x1'] - left) / W for l in L]); s = np.array([(l['x0'] - left) / W for l in L])
        fullshare = np.mean(r >= 0.9)
        for i, l in enumerate(L):
            rows.append(dict(ms=ms, page=fn.split('/')[-1], block=b['id'], i=i, n=len(L), r=r[i], s=s[i],
                             w=len(l['text'].split()), text=l['text'], fullshare=fullshare,
                             prev_full=int(i >= 3 and all(r[i - 3:i] >= 0.9) and all(s[i - 3:i] < 0.2)),
                             next_exists=i + 1 < len(L)))

import pandas as pd
df = pd.DataFrame(rows)
print('lines', len(df), 'MSS', df.ms.nunique(), 'pages', df.page.nunique())
print(df.groupby('ms').agg(lines=('w', 'size'), wpl=('w', 'mean'), fullshare=('fullshare', 'mean'),
                           early=('r', lambda x: np.mean(x < 0.85))).round(2))

for thr in (0.9, 0.85, 0.75):
    e = df[(df.r < thr)]
    print(f'\nNAIVE: all lines with right extent < {thr}: n={len(e)} single-word={np.mean(e.w==1)*100:.1f}%  <=2={np.mean(e.w<=2)*100:.1f}%')
    e2 = e[e.s < 0.15]
    print(f'   and starting at left margin: n={len(e2)} single={np.mean(e2.w==1)*100:.1f}% <=2={np.mean(e2.w<=2)*100:.1f}%')
    print('   by MS:', e2.groupby('ms').agg(n=('w', 'size'), single=('w', lambda x: round(100 * np.mean(x == 1), 1))).to_dict('index'))

print('\nSTRICT prose paragraph ends: block >=75% full lines, line r<thr, starts at left, preceded by 3 full lines')
prose = df[df.fullshare >= 0.75]
full = prose[(prose.r >= 0.9) & (prose.s < 0.15)]
print('prose lines', len(prose), 'full-line words mean', full.w.mean().round(2), 'median', full.w.median())
capdist = full.w.values
for thr in (0.9, 0.85, 0.75):
    pe = prose[(prose.r < thr) & (prose.s < 0.15) & (prose.prev_full == 1)]
    k = int((pe.w == 1).sum()); k2 = int((pe.w <= 2).sum()); n = len(pe)
    # null: remainder distribution P(R=k) ~ P(C>=k)/E[C]; condition on r<thr approximated by rem/cap<thr
    rng = np.random.default_rng(3)
    sims1 = []; sims2 = []
    for _ in range(2000):
        # each paragraph end: draw capacity c (size-biased not needed for remainder: use renewal: pick total uniformly)
        c = rng.choice(capdist, size=4 * n + 50)
        u = rng.integers(1, c + 1)  # remainder uniform on 1..c given capacity of last line (length-biased ignored)
        keep = (u / c) < thr
        u = u[keep][:n]
        sims1.append(np.sum(u == 1)); sims2.append(np.sum(u <= 2))
    print(f' thr={thr}: n={n} single={k} ({100*k/max(n,1):.1f}%) <=2={k2} ({100*k2/max(n,1):.1f}%) | null single={np.mean(sims1):.1f} [{np.percentile(sims1,2.5):.0f},{np.percentile(sims1,97.5):.0f}] '
          f'<=2={np.mean(sims2):.1f} [{np.percentile(sims2,2.5):.0f},{np.percentile(sims2,97.5):.0f}]  median r={pe.r.median():.2f} MSS={pe.ms.nunique()}')
    if thr == 0.85:
        h, _ = np.histogram(pe.r, bins=[-1, .2, .4, .6, .8, 1.0]); print('   r hist', list(h))
        print('   last-line word counts', sorted(Counter(pe.w).items()))
        print('   by MS', pe.groupby('ms').agg(n=('w', 'size'), single=('w', lambda x: int(np.sum(x == 1))), le2=('w', lambda x: int(np.sum(x <= 2)))).to_dict('index'))
        print('   examples of single-word ends:', list(pe[pe.w == 1].text.head(12)))
        pe.to_csv('real_para_ends.csv', index=False)
