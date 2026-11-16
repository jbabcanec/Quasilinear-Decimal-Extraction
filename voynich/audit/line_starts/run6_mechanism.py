"""Mechanism probes: which glyphs drive the avoidance, transition matrix, mid-line controls,
previous-line-ending channel, second-glyph dependence, v101 strata."""
import numpy as np
import pandas as pd
from build import build
from stats import encode, pairs_within, permtest, fmt
from run3_attacks import runm
from zlparse import words_of, first_sign

df, pages, order, recs, paras, orphans = build()
body = df[df.li >= 2].reset_index(drop=True)

print('#### 1. transition table (adjacent non-first lines), observed / expected under within-paragraph permutation')
col = 'g_bench'
g, dic = encode(list(body[col]))
inv = {v: k for k, v in dic.items()}
a, b = pairs_within(body.pid.values, 1, mask=(g >= 0))
G = len(dic)
obs = np.zeros((G, G))
np.add.at(obs, (g[a], g[b]), 1)
exp = np.zeros((G, G))
for pid, idx in body.groupby('pid').indices.items():
    gg = g[idx]
    gg = gg[gg >= 0]
    m = len(gg)
    if m < 2:
        continue
    cnt = np.bincount(gg, minlength=G).astype(float)
    # number of adjacent eligible pairs in this paragraph
    full = g[idx]
    npairs = int(((full[1:] >= 0) & (full[:-1] >= 0)).sum())
    e = np.outer(cnt, cnt)
    e[np.diag_indices(G)] = cnt * (cnt - 1)
    exp += npairs * e / (m * (m - 1))
tot = np.bincount(g[g >= 0], minlength=G)
top = [i for i in np.argsort(-tot)[:9]]
print('glyph      n   diag_obs diag_exp ratio')
for i in top:
    print('  %-4s %5d   %6.0f  %7.1f  %.2f' % (inv[i], tot[i], obs[i, i], exp[i, i], obs[i, i] / exp[i, i] if exp[i, i] else np.nan))
print('obs/exp matrix (row = upper line, col = lower line), top glyphs:')
print('      ' + ' '.join('%5s' % inv[j] for j in top))
for i in top:
    print('  %-4s' % inv[i] + ' '.join('%5.2f' % (obs[i, j] / exp[i, j]) if exp[i, j] > 5 else '    .' for j in top))
print('sum diag obs %.0f exp %.1f ; off-diag obs %.0f exp %.1f' % (np.trace(obs), np.trace(exp), obs.sum() - np.trace(obs), exp.sum() - np.trace(exp)))

print('\n#### 2. coarser classes: does avoidance extend to visually similar signs?')
cls = {'t': 'gallows', 'k': 'gallows', 'p': 'gallows', 'f': 'gallows', 'cth': 'benchgallows', 'ckh': 'benchgallows',
       'cph': 'benchgallows', 'cfh': 'benchgallows', 'ch': 'bench', 'sh': 'bench'}
body['g_cls'] = body.g_bench.map(lambda x: cls.get(x, x) if x is not None else None)
runm(body, col='g_cls', lags=(1,), label='class (gallows merged, benches merged)')
# pairs with different first sign: is the class still avoided?
gc, _ = encode(list(body.g_cls))
_, grp = np.unique(body.pid.values, return_inverse=True)
# gallows-gallows adjacency specifically
isg = body.g_cls.values == 'gallows'
print('  gallows lines: %d; adjacent gallows-gallows pairs: %d; same exact gallows: %d' % (
    isg.sum(), int((isg[a] & isg[b]).sum()), int((isg[a] & isg[b] & (g[a] == g[b])).sum())))
gg2, _ = encode([('G' if x else 'n') for x in isg])
print('  gallows-vs-not as binary feature:', fmt(permtest(gg2, grp, a, b, nperm=2000, seed=3)))

print('\n#### 3. is the 2nd-glyph avoidance independent of the first glyph?')
g2, _ = encode(list(body.g2))
g1, _ = encode(list(body.g_eva))
ok = (g1 >= 0) & (g2 >= 0)
a2, b2 = pairs_within(body.pid.values, 1, mask=ok)
same1 = g1[a2] == g1[b2]
print('  pairs %d ; P(2nd glyph same | 1st same)=%.3f (n=%d) ; P(2nd same | 1st differ)=%.3f (n=%d)' % (
    len(a2), (g2[a2] == g2[b2])[same1].mean(), same1.sum(), (g2[a2] == g2[b2])[~same1].mean(), (~same1).sum()))
# word minus first glyph ("stem")
body['stem'] = body.w1.map(lambda w: w[1:] if isinstance(w, str) and len(w) > 1 else None)
runm(body, col='stem', lags=(1,), label='first word minus its first glyph, identical')
body['stem2'] = body.w1.map(lambda w: w[1:3] if isinstance(w, str) and len(w) > 2 else None)
runm(body, col='stem2', lags=(1,), label='glyphs 2-3 of first word')
body['g3'] = body.w1.map(lambda w: w[2] if isinstance(w, str) and len(w) > 2 else None)
runm(body, col='g3', lags=(1,), label='3rd glyph of first word')
body['lastg1'] = body.w1.map(lambda w: w[-1] if isinstance(w, str) and len(w) > 0 else None)
runm(body, col='lastg1', lags=(1,), label='last glyph of first word')

print('\n#### 4. mid-line horizontal control: do neighbouring words in running text avoid the same first glyph?')
rows = []
for pid, p in enumerate(paras):
    for li, r in enumerate(p):
        w = words_of(r['text'].lstrip('.,\x01\x02'))
        for wi, x in enumerate(w):
            rows.append((pid, li, wi, len(w), first_sign(x, 'eva')))
W = pd.DataFrame(rows, columns=['pid', 'li', 'wi', 'nw', 'g'])
W['g'] = W.g.where(W.g != '?', None)
W['line'] = W.pid * 1000 + W.li
gw, _ = encode(list(W.g))
_, grpw = np.unique(W.pid.values, return_inverse=True)
for lag in (1, 2):
    aw, bw = pairs_within(W.line.values, lag, mask=(gw >= 0))
    # exclude pairs involving the line-first word
    sel = W.wi.values[aw] >= 1
    print('  words at distance %d in same line (line-first word excluded), null = shuffle words in paragraph: %s' % (
        lag, fmt(permtest(gw, grpw, aw[sel], bw[sel], nperm=300, seed=lag))))
    sel0 = W.wi.values[aw] == 0
    print('  line-first word vs word %d of same line: %s' % (lag + 1, fmt(permtest(gw, grpw, aw[sel0], bw[sel0], nperm=300, seed=lag))))
# line wrap: last word of line i -> first word of line i+1
aw, bw = pairs_within(W.pid.values, 1, mask=(gw >= 0))
selwrap = (W.line.values[aw] != W.line.values[bw])
print('  last word of line -> first word of next line: %s' % fmt(permtest(gw, grpw, aw[selwrap], bw[selwrap], nperm=300, seed=5)))

print('\n#### 5. vertical test for every word position (word k of line vs word k of next line)')
for k in range(0, 6):
    key = []
    for pid, p in enumerate(paras):
        for li, r in enumerate(p):
            if li == 0:
                continue
            w = words_of(r['text'].lstrip('.,\x01\x02'))
            x = first_sign(w[k], 'eva') if len(w) > k else None
            key.append(None if x == '?' else x)
    body['wk'] = key
    r = runm(body, col='wk', lags=(1,), quiet=True)[1]
    print('  word %d first glyph: %s' % (k + 1, fmt(r)))
# first word after a drawing interruption (2nd segment) - a "pseudo line start"
r = runm(body, col='seg2g', lags=(1, 2), quiet=True)
for k in r:
    print('  first glyph after first <-> (pseudo line start) lag%d: %s' % (k, fmt(r[k])))

print('\n#### 6. channel through the previous line\'s ending? expected match rate given last glyph of upper line')
g1, d1 = encode(list(body.g_eva))
lg, _ = encode(list(body.lastg))
a, b = pairs_within(body.pid.values, 1, mask=(g1 >= 0))
ok = lg[a] >= 0
a, b = a[ok], b[ok]
G = g1.max() + 1
# distribution of lower-line first glyph given upper-line last glyph (leave-one-out)
tab = np.zeros((lg.max() + 1, G))
np.add.at(tab, (lg[a], g1[b]), 1)
p_match = (tab[lg[a], g1[a]] - (g1[a] == g1[b])) / (tab[lg[a]].sum(1) - 1)
print('  pairs %d obs %.1f%% ; expected if lower first glyph depends only on upper line\'s last glyph: %.1f%%' % (
    len(a), 100 * (g1[a] == g1[b]).mean(), 100 * p_match.mean()))
# same but conditioning on the upper line's last WORD first glyph and section
print('\n#### 7. v101 by ZL page metadata')
from v101parse import v101_table
d = v101_table()
meta = {f: pages[f] for f in pages}
d['H'] = d.folio.map(lambda f: meta.get(f, {}).get('H', '?'))
d['L'] = d.folio.map(lambda f: meta.get(f, {}).get('L', '-'))
d['I'] = d.folio.map(lambda f: meta.get(f, {}).get('I', '?'))
d['Q'] = d.folio.map(lambda f: meta.get(f, {}).get('Q', '?'))
vb = d[d.li >= 2].reset_index(drop=True)
print('  v101 pages without ZL metadata:', sorted(set(vb.folio[vb.H == '?']))[:20])
for colname in ['H', 'L', 'I']:
    for v in sorted(vb[colname].astype(str).unique()):
        mk = (vb[colname].astype(str) == v).values
        if mk.sum() < 30:
            continue
        r = runm(vb, col='g', mask=mk, lags=(1,), quiet=True)
        if 1 in r:
            print('  v101 %s=%s: %s' % (colname, v, fmt(r[1])))
nq = 0
nlo = 0
for v in sorted(vb.Q.astype(str).unique()):
    mk = (vb.Q.astype(str) == v).values
    r = runm(vb, col='g', mask=mk, lags=(1,), quiet=True)
    if 1 in r and r[1]['N'] >= 10:
        nq += 1
        nlo += r[1]['obs'] < r[1]['exp']
        print('  v101 quire %s: %s' % (v, fmt(r[1])))
print('  v101 quires with obs<exp: %d of %d' % (nlo, nq))
