"""Geometry attack: use Takahashi word boxes to (a) re-derive lines without trusting any
transliteration's line order, (b) look for interleaving / alternating indentation / split halves,
(c) repeat the test on a third transliteration (Takahashi words)."""
import json
import os
import numpy as np
import pandas as pd
from rapidfuzz import fuzz
from build import build
from stats import encode, pairs_within, permtest, fmt
from run3_attacks import runm

B = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/voynich-units/morphometry_voynichese/voynichese_boxes/'
df, pages, order, recs, paras, orphans = build()


def box_lines(fol):
    d = json.load(open(B + fol + '.js'))
    words = [w[0] for w in d[0]]
    lines, cur = [], []
    for wi, x, y, w, h in d[1]:
        if cur and x < cur[-1][1]:       # x went backwards -> new line
            lines.append(cur)
            cur = []
        cur.append((words[wi], x, y, w, h))
    if cur:
        lines.append(cur)
    return lines


# ZL loci per page (all types) to find pure-paragraph pages
zl_by_page = {}
for r in recs:
    zl_by_page.setdefault(r['folio'], []).append(r)

rows = []
status = {}
for fol in order:
    if not os.path.exists(B + fol + '.js'):
        status[fol] = 'nobox'
        continue
    zl = zl_by_page.get(fol, [])
    if not zl:
        status[fol] = 'nozl'
        continue
    bl = box_lines(fol)
    allP = all(r['typ'].startswith('P') for r in zl)
    if not allP:
        status[fol] = 'mixed'
        continue
    if len(bl) != len(zl):
        status[fol] = 'count_mismatch(%d box vs %d zl)' % (len(bl), len(zl))
        continue
    # verify alignment through first-word similarity
    sims = []
    for r, l in zip(zl, bl):
        zt = r['text'].replace('\x01', '.').replace('\x02', '.').replace(',', '.').replace('.', '')
        bt = ''.join(w[0] for w in l)
        sims.append(fuzz.ratio(zt, bt))
    if np.median(sims) < 80 or np.min(sims) < 40:
        status[fol] = 'lowsim(med %.0f min %.0f)' % (np.median(sims), np.min(sims))
        continue
    status[fol] = 'ok'
    for r, l, s in zip(zl, bl, sims):
        rows.append(dict(folio=fol, num=r['num'], x0=l[0][1], y0=l[0][2] + l[0][4] / 2.0,
                         ytop=l[0][2], h0=l[0][4], xend=l[-1][1] + l[-1][3], tk_w1=l[0][0],
                         nbox=len(l), sim=s,
                         maxgap=max([l[i + 1][1] - (l[i][1] + l[i][3]) for i in range(len(l) - 1)] or [0])))

import collections
c = collections.Counter(s.split('(')[0] for s in status.values())
print('page alignment status:', dict(c))
geo = pd.DataFrame(rows)
m = df.merge(geo, on=['folio', 'num'], how='left')
print('P lines with geometry: %d of %d' % (m.x0.notna().sum(), len(m)))
m['tk_g'] = m.tk_w1.map(lambda w: w[0] if isinstance(w, str) and w and w[0] != '?' else None)
agree = (m.tk_g == m.g_eva)[m.x0.notna()].mean()
print('Takahashi first glyph == ZL first glyph on aligned lines: %.1f%%' % (100 * agree))

# per-page geometry normalisation
m['dy'] = m.groupby('pid').y0.diff()
m['dx'] = m.groupby('pid').x0.diff()
pg_med_dy = m.groupby('folio').dy.transform('median')
m['dyn'] = m.dy / pg_med_dy
has = m.x0.notna()
body = m[(m.li >= 2)].reset_index(drop=True)
hasb = body.x0.notna().values

print('\n== ZL test restricted to geometrically aligned pages ==')
runm(body, mask=hasb, label='G0 ZL glyphs, aligned pages')
print('== third transliteration (Takahashi first word of box-line), same pages ==')
runm(body, col='tk_g', mask=hasb, label='G1 Takahashi glyphs, aligned pages')

print('\n== geometry of consecutive lines inside paragraphs (lines li>=2 and their predecessor) ==')
d = body[hasb & body.dy.notna().values & (body.li >= 3).values]
print('pairs with geometry:', len(d))
print('normalised line step dy/median(page): quantiles', np.round(d.dyn.quantile([0, .01, .05, .25, .5, .75, .95, .99, 1]).values, 2))
print('pairs with dy<=0 (next line not lower):', int((d.dy <= 0).sum()))
print('pairs with step <0.5 or >1.6 of page median:', int(((d.dyn < 0.5) | (d.dyn > 1.6)).sum()))
print('dx (indent change, px) quantiles', np.round(d.dx.quantile([0, .01, .05, .25, .5, .75, .95, .99, 1]).values, 1))
print('pairs with |dx|>30px:', int((d.dx.abs() > 30).sum()), ' >60px:', int((d.dx.abs() > 60).sum()))

# alternating indentation? lag-1 autocorrelation of dx (successive indent changes) & of x0 residual
b2 = body.copy()
b2['x0res'] = b2.x0 - b2.groupby('pid').x0.transform('median')
for lag in (1, 2):
    a, bb = pairs_within(b2.pid.values, lag, mask=b2.x0.notna().values)
    r = np.corrcoef(b2.x0res.values[a], b2.x0res.values[bb])[0, 1]
    print('indent residual autocorrelation lag %d: r=%+.3f (n=%d)' % (lag, r, len(a)))
# line-length alternation
b2['len_res'] = b2.nchar - b2.groupby('pid').nchar.transform('mean')
for lag in (1, 2, 3):
    a, bb = pairs_within(b2.pid.values, lag, mask=(~b2['last']).values)
    r = np.corrcoef(b2.len_res.values[a], b2.len_res.values[bb])[0, 1]
    print('line length (chars, para-centred, last line excl.) autocorr lag %d: r=%+.3f (n=%d)' % (lag, r, len(a)))

print('\n== ZL test on pairs with clean geometry: step within 0.6-1.5 of page median and |indent change|<=15px ==')
g, _ = encode(list(body.g_eva))
_, grp = np.unique(body.pid.values, return_inverse=True)
a, bb = pairs_within(body.pid.values, 1, mask=(g >= 0))
dyn = body.dyn.values[bb]
dx = body.dx.values[bb]
ok = np.isfinite(dyn) & (dyn > 0.6) & (dyn < 1.5) & (np.abs(dx) <= 15)
print('clean pairs      :', fmt(permtest(g, grp, a[ok], bb[ok], nperm=3000, seed=1)))
bad = np.isfinite(dyn) & ~ok
print('irregular pairs  :', fmt(permtest(g, grp, a[bad], bb[bad], nperm=3000, seed=1)))
nog = ~np.isfinite(dyn)
print('no-geometry pairs:', fmt(permtest(g, grp, a[nog], bb[nog], nperm=3000, seed=1)))

# does a big internal gap (drawing) matter? left-margin lines only i.e. lines whose first box is at the page's left text margin
pg_x = m.groupby('folio').x0.transform(lambda s: s.quantile(0.1))
body['x0rel'] = (m.x0 - pg_x)[m.li >= 2].values
x0rel = body.x0rel.values
okm = np.isfinite(dyn) & (np.abs(x0rel[a]) <= 20) & (np.abs(x0rel[bb]) <= 20)
print('both lines start within 20px of page left margin:', fmt(permtest(g, grp, a[okm], bb[okm], nperm=3000, seed=1)))

# list status of non-ok pages
bad_pages = {k: v for k, v in status.items() if v != 'ok'}
print('\nnon-aligned pages (%d):' % len(bad_pages), ', '.join('%s:%s' % kv for kv in list(bad_pages.items())[:60]))
m.to_pickle('geo_merged.pkl')
