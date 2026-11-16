"""Step 1: base replication on ZL3b. Lags 1..6, within-paragraph shuffle null, first line excluded."""
import numpy as np
from build import build
from stats import encode, pairs_within, permtest, fmt

df, *_ = build()
print('P lines in marked paragraphs:', len(df), 'paragraphs:', df.pid.nunique())

def run(d, col='g_eva', lags=(1, 2, 3, 4, 5, 6), group='pid', unit='pid', nperm=4000, label=''):
    d = d.reset_index(drop=True)
    g, _ = encode(list(d[col]))
    unitv = d[unit].values
    _, grp = np.unique(d[group].values, return_inverse=True)
    out = {}
    for lag in lags:
        a, b = pairs_within(unitv, lag, mask=(g >= 0))
        r = permtest(g, grp, a, b, nperm=nperm, seed=lag)
        out[lag] = r
        print('%-34s lag%d %s' % (label, lag, fmt(r)))
    return out

if __name__ == '__main__':
    body = df[df.li >= 2]
    print('non-first lines:', len(body), 'unknown first glyph:', body.g_eva.isna().sum())
    print('\n== BASE: first line excluded, EVA first char, within-paragraph shuffle ==')
    run(body, 'g_eva', label='ZL eva, excl line1')
    print('\n== ch/sh/cth/ckh/cph/cfh as single signs ==')
    run(body, 'g_bench', label='ZL bench, excl line1')
    print('\n== also qo as single sign ==')
    run(body, 'g_bq', lags=(1, 2), label='ZL bench+qo, excl line1')
    print('\n== first line included ==')
    run(df, 'g_eva', lags=(1, 2), label='ZL eva, incl line1')
    print('\n== negative controls: second word first glyph, third word, last glyph, last word first glyph, 2nd glyph of first word ==')
    run(body, 'w2g', lags=(1, 2), label='2nd word first glyph')
    run(body, 'w2g_nc', lags=(1, 2), label='2nd word (comma not a break)')
    run(body, 'w3g', lags=(1, 2), label='3rd word first glyph')
    run(body, 'lastg', lags=(1, 2), label='line-last glyph')
    run(body, 'lastw_first', lags=(1, 2), label='last word first glyph')
    run(body, 'g2', lags=(1, 2), label='2nd glyph of first word')
    run(body, 'w1', lags=(1, 2), label='whole first word')
