"""Exact-line-index position null, second-glyph conditional test, horizontal comparison,
method calibration on wrapped control texts."""
import numpy as np
import pandas as pd
from build import build
from stats import encode, pairs_within, permtest, fmt, expected
from run3_attacks import runm
from zlparse import words_of, first_sign

rng = np.random.default_rng(11)
df, pages, order, recs, paras, orphans = build()
body = df[df.li >= 2].reset_index(drop=True)

print('#### 1. position-matched nulls using the exact line index')
body['HI'] = body.H.astype(str) + '|' + body.L.astype(str) + '|' + body.I.astype(str)
body['HIli'] = body.HI + '|' + body.li.clip(upper=10).astype(str)
runm(body, group='HIli', label='shuffle within hand/lang/section x exact line index')
body['HIliE'] = body.HI + '|' + body.li.clip(upper=6).astype(str) + '|' + (body.npar - body.li).clip(upper=2).astype(str)
runm(body, group='HIliE', label='... x line index(<=6) x distance-from-end(<=2)')
# per pair-position, with the position-matched null
g, _ = encode(list(body.g_eva))
_, grp = np.unique(body.HIli.values, return_inverse=True)
a, b = pairs_within(body.pid.values, 1, mask=(g >= 0))
for lo, hi in [(2, 2), (3, 3), (4, 5), (6, 99)]:
    sel = (body.li.values[a] >= lo) & (body.li.values[a] <= hi)
    print('   upper line %d-%d, position-matched null: %s' % (lo, hi, fmt(permtest(g, grp, a[sel], b[sel], nperm=1500, seed=2))))

print('\n#### 2. beyond the first glyph: is the 2nd glyph avoided when the 1st glyphs differ?')
g1, _ = encode(list(body.g_eva))
g2, _ = encode(list(body.g2))
ok = (g1 >= 0) & (g2 >= 0)
a, b = pairs_within(body.pid.values, 1, mask=ok)
obs_n = int(((g1[a] != g1[b]) & (g2[a] == g2[b])).sum())
obs_d = int((g1[a] != g1[b]).sum())
idx = np.nonzero(ok)[0]
pidv = body.pid.values[idx]
null = []
for k in range(3000):
    o = np.lexsort((rng.random(len(idx)), pidv))
    p1 = g1.copy()
    p2 = g2.copy()
    p1[idx] = g1[idx][o]
    p2[idx] = g2[idx][o]
    d = p1[a] != p1[b]
    null.append(((d & (p2[a] == p2[b])).sum() / d.sum()))
null = np.array(null)
print('   P(2nd same | 1st differ): obs %.3f (%d/%d) ; shuffle null %.3f sd %.4f ; z=%+.2f' % (
    obs_n / obs_d, obs_n, obs_d, null.mean(), null.std(ddof=1), (obs_n / obs_d - null.mean()) / null.std(ddof=1)))
# conditional the other way: among pairs with same 1st glyph, 2nd glyph
obs_n2 = int(((g1[a] == g1[b]) & (g2[a] == g2[b])).sum())
obs_d2 = int((g1[a] == g1[b]).sum())
print('   P(2nd same | 1st same): obs %.3f (%d/%d)' % (obs_n2 / obs_d2, obs_n2, obs_d2))

print('\n#### 3. horizontal comparison: transition ratios between adjacent mid-line words (first signs)')
rows = []
for pid, p in enumerate(paras):
    for li, r in enumerate(p):
        w = words_of(r['text'].lstrip('.,\x01\x02'))
        for wi, x in enumerate(w):
            rows.append((pid, li, wi, first_sign(x, 'bench'), r['meta'].get('H', '-')))
W = pd.DataFrame(rows, columns=['pid', 'li', 'wi', 'g', 'H'])
W['g'] = W.g.where(W.g != '?', None)
W['line'] = W.pid * 1000 + W.li
gw, dw = encode(list(W.g))
invw = {v: k for k, v in dw.items()}
aw, bw = pairs_within(W.line.values, 1, mask=(gw >= 0))
sel = W.wi.values[aw] >= 1          # neither word is line-first
aw, bw = aw[sel], bw[sel]
Gw = len(dw)
obs = np.zeros((Gw, Gw))
np.add.at(obs, (gw[aw], gw[bw]), 1)
# independence expectation from marginals of mid-line words within each paragraph (pooled approx: product of marginals)
exp = np.zeros((Gw, Gw))
pa = pd.Series(1.0, index=pd.MultiIndex.from_arrays([W.pid.values[aw], gw[aw]])).groupby(level=[0, 1]).sum()
pb = pd.Series(1.0, index=pd.MultiIndex.from_arrays([W.pid.values[bw], gw[bw]])).groupby(level=[0, 1]).sum()
for pid in np.unique(W.pid.values[aw]):
    ca = np.zeros(Gw)
    cb = np.zeros(Gw)
    x = pa.loc[pid]
    ca[x.index.values] = x.values
    y = pb.loc[pid]
    cb[y.index.values] = y.values
    exp += np.outer(ca, cb) / ca.sum()
tops = ['d', 'y', 'o', 'q', 's', 't', 'sh', 'ch']
print('   horizontal (word j -> word j+1, mid-line), obs/exp with paragraph-level marginals:')
print('         ' + ' '.join('%5s' % t for t in tops))
for i in tops:
    print('     %-4s' % i + ' '.join('%5.2f' % (obs[dw[i], dw[j]] / exp[dw[i], dw[j]]) for j in tops))

print('\n#### 4. method calibration on line-wrapped control texts (pseudo-paragraphs = blank-line blocks or 6-line chunks)')
C = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/controls/'
for name in ['naibbe/naibbe_caesar_ciphertext_lines.txt', 'naibbe/greshko_nathist_output_ciphertext.txt',
             'selfcitation/selfcite_seed1.txt', 'selfcitation/selfcite_seed2.txt', 'selfcitation/selfcite_seed3.txt']:
    lines = open(C + name, encoding='utf-8', errors='replace').read().split('\n')
    blocks, cur = [], []
    for l in lines:
        if l.strip() == '':
            if cur:
                blocks.append(cur)
            cur = []
        else:
            cur.append(l.strip())
    if cur:
        blocks.append(cur)
    nb_blank = len(blocks)
    if nb_blank < 20:      # no blank-line structure: chunk
        alll = [l for bl in blocks for l in bl]
        blocks = [alll[i:i + 6] for i in range(0, len(alll), 6)]
    rows = []
    for pid, bl in enumerate(blocks):
        for li, l in enumerate(bl):
            w = l.split()
            rows.append(dict(pid=pid, li=li + 1, g=first_sign(w[0], 'eva') if w else None,
                             w2=first_sign(w[1], 'eva') if len(w) > 1 else None))
    d = pd.DataFrame(rows)
    r = runm(d[d.li >= 2], col='g', lags=(1, 2), quiet=True, nperm=1500)
    print('   %-52s blocks=%d' % (name.split('/')[-1], len(blocks)))
    for k in r:
        print('       first glyph lag%d: %s' % (k, fmt(r[k])))
