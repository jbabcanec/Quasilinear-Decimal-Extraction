import numpy as np, pickle
from collections import defaultdict, Counter
import zl, boxes

rows = zl.read()
bypage = defaultdict(list)
for r in rows:
    bypage[r['page']].append(r)
paras, an = zl.paragraphs(rows)
amap = {}
for pg, rr in bypage.items():
    out, rate = boxes.align_page(pg, rr)
    amap[pg] = out

recs = []
for p in paras:
    L = p['lines']
    pg = L[0]['page']
    if amap[pg] is None:
        continue
    lines = []
    ok = True
    for l in L:
        bx = amap[pg].get(l['num'], [])
        zc = len(zl.words(l['raw']))
        if len(bx) == 0:
            ok = False
        lines.append(dict(num=l['num'], typ=l['typ'], pos=l['pos'], zc=zc, bc=len(bx),
                          x0=min([b[1] for b in bx]) if bx else np.nan,
                          x1=max([b[1] + b[3] for b in bx]) if bx else np.nan,
                          y=np.median([b[2] + b[4] / 2 for b in bx]) if bx else np.nan,
                          ink=sum(b[3] for b in bx), boxes=bx, raw=l['raw']))
    recs.append(dict(page=pg, I=L[0]['I'], H=L[0]['H'], L=L[0]['L'], lines=lines, ok=ok))
pickle.dump(recs, open('box_paras.pkl', 'wb'))
print('paragraphs with boxes', len(recs), 'all loci have boxes', sum(r['ok'] for r in recs))
mis = sum(1 for r in recs for l in r['lines'] if l['zc'] != l['bc'])
tot = sum(len(r['lines']) for r in recs)
print('loci where box count != ZL word count', mis, '/', tot)


def measure(recs, minl, mode='A'):
    n = 0; one = 0; two = 0; ratios = []; singles = []; onez = 0
    for r in recs:
        L = r['lines']
        if mode == 'B':
            L = [l for l in L if l['typ'] not in ('Pt', 'Pc', 'Pr')]
        if len(L) < minl or not all(l['bc'] > 0 for l in L):
            continue
        left = np.median([l['x0'] for l in L[:-1]])
        right = np.median([l['x1'] for l in L[:-1]])
        last = L[-1]
        # extent of last line measured from the paragraph's left margin
        ratio = (last['x1'] - left) / (right - left)
        n += 1
        if last['bc'] <= 1:
            one += 1; singles.append((r['page'], last['num'], last['typ'], last['raw'][:30]))
        if last['zc'] <= 1:
            onez += 1
        if last['bc'] <= 2:
            two += 1
        ratios.append(ratio)
    ratios = np.array(ratios)
    h, _ = np.histogram(ratios, bins=[-5, .2, .4, .6, .8, 1.0, 1.2, 10])
    print(f'boxes mode={mode} min{minl}: n={n} single-box last={one} ({100*one/n:.1f}%) [ZL single={onez}] <=2={two} ({100*two/n:.1f}%) '
          f'extent ratio med={np.median(ratios):.2f} mean={np.mean(ratios):.2f} hist={list(map(int,h))}  share<0.2={np.mean(ratios<0.2):.3f} <0.3={np.mean(ratios<0.3):.3f}')
    return singles, ratios


for mode in 'AB':
    for minl in (4, 3):
        s, ratios = measure(recs, minl, mode)
        if minl == 4:
            print('   ', s)
            print('    deciles of extent ratio', np.round(np.percentile(ratios, [5, 10, 25, 50, 75, 90, 95]), 2))
