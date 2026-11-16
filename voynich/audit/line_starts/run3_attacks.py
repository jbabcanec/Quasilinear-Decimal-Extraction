"""Step 2: attacks on the lag-1 avoidance claim (ZL3b)."""
import sys
import numpy as np
import pandas as pd
from build import build
from stats import encode, pairs_within, permtest, fmt

df, pages, order, recs, paras, orphans = build()
body = df[df.li >= 2].reset_index(drop=True)
NPERM = 3000


def runm(d, col='g_eva', mask=None, lags=(1, 2), group='pid', unit='pid', label='', nperm=NPERM, quiet=False):
    """Rows of d stay in place (adjacency preserved); rows with mask False are ineligible
    (not in pairs, not permuted)."""
    d = d.reset_index(drop=True)
    keys = list(d[col])
    if mask is not None:
        mask = np.asarray(mask)
        keys = [k if m else None for k, m in zip(keys, mask)]
    g, _ = encode(keys)
    unitv = d[unit].values
    _, grp = np.unique(d[group].values.astype(str), return_inverse=True)
    out = {}
    for lag in lags:
        a, b = pairs_within(unitv, lag, mask=(g >= 0))
        if len(a) == 0:
            continue
        r = permtest(g, grp, a, b, nperm=nperm, seed=lag)
        out[lag] = r
        if not quiet:
            print('%-44s lag%d %s' % (label, lag, fmt(r)))
    return out


if __name__ == '__main__':
    b = body
    print('#### A. null variants (pairs = adjacent non-first lines in same paragraph)')
    runm(b, label='A1 shuffle within paragraph')
    runm(b, group='folio', label='A2 shuffle within page')
    b['QS'] = b.Q.astype(str)
    runm(b, group='QS', label='A3 shuffle within quire')
    b['one'] = 'x'
    runm(b, group='one', label='A4 shuffle whole MS')
    b['HL'] = b.H.astype(str) + b.L.astype(str) + b.I.astype(str)
    runm(b, group='HL', label='A5 shuffle within hand x lang x section')
    # position-matched: same hand/lang/section and same position class (line 2 / interior / last)
    pos = np.where(b.li == 2, 'second', np.where(b['last'], 'last', 'mid'))
    b['HLpos'] = b.HL + '|' + pos
    runm(b, group='HLpos', label='A6 shuffle within hand/lang/sect x position')
    b['pagepos'] = b.folio + '|' + pos
    runm(b, group='pagepos', label='A7 shuffle within page x position class')

    print('\n#### A8. paragraph definitions')
    # whole page as one unit, all P lines (incl. first lines, orphans ignored)
    d = df.copy()
    runm(d, unit='folio', group='folio', label='A8a page = unit, all P lines')
    # page as unit but paragraph-first lines made ineligible
    runm(d, unit='folio', group='folio', mask=(d.li >= 2), label='A8b page unit, first lines ineligible')
    # cross-paragraph pairs only: last line of para p and 2nd line of p+1 ... not meaningful; skip
    # paragraph by locus prefix: new paragraph at '@' or '*' prefix
    allP = pd.DataFrame([dict(folio=r['folio'], pre=r['pre'], typ=r['typ'], text=r['text'], pstart=r['pstart'])
                         for r in recs if r['typ'].startswith('P')])
    from zlparse import first_sign, words_of
    allP['g'] = [first_sign((words_of(t.lstrip('.,\x01\x02')) or [''])[0], 'eva') for t in allP.text]
    allP['g'] = allP.g.map(lambda x: None if x in ('?', None) else x)
    newp = (allP.pre.isin(['@', '*'])) | (allP.folio != allP.folio.shift())
    allP['pid2'] = newp.cumsum()
    allP['li2'] = allP.groupby('pid2').cumcount() + 1
    print('   prefix-defined paragraphs:', allP.pid2.nunique())
    runm(allP, col='g', unit='pid2', group='pid2', mask=(allP.li2 >= 2), label='A8c prefix-defined paragraphs (@,*)')

    print('\n#### B. subsets (ineligible lines keep their place; shuffle among eligible lines in paragraph)')
    runm(b, mask=(b.ndraw == 0), label='B1 no drawing interruption (both lines)')
    runm(b, mask=(b.ndraw > 0), label='B2 only drawing-interrupted lines')
    # pages with no interruption at all
    pg_draw = df.groupby('folio').ndraw.sum()
    clean_pages = set(pg_draw[pg_draw == 0].index)
    runm(b, mask=b.folio.isin(clean_pages), label='B3 pages without any <-> at all')
    runm(b, mask=(b.li >= 3), label='B4 exclude paragraph 2nd line')
    runm(b, mask=(~b['last']), label='B5 exclude paragraph last line')
    runm(b, mask=((b.li >= 3) & (~b['last'])), label='B6 interior lines only (no 2nd, no last)')
    runm(b, mask=(b.typ == 'P0'), label='B7 P0 loci only')
    runm(b, mask=(b.nw >= 6), label='B8 both lines >= 6 words (comma=break)')
    runm(b, mask=(b.nw_nocomma >= 6), label='B9 both lines >= 6 words (comma not break)')
    runm(b, mask=(b.nw_nocomma >= 8), label='B9b both lines >= 8 words (comma not break)')
    runm(b, mask=(b.nw < 6), label='B10 both lines < 6 words')
    runm(b, mask=(b.npar >= 9), label='B11 paragraphs with >= 9 lines')
    runm(b, mask=(b.npar <= 5), label='B12 paragraphs with <= 5 lines')
    harsh = (b.li >= 3) & (~b['last']) & (b.typ == 'P0') & (b.ndraw == 0) & (b.nw_nocomma >= 6)
    runm(b, mask=harsh, label='B13 HARSH: interior,P0,no <->,>=6w, EVA')
    runm(b, col='g_bench', mask=harsh, label='B13 HARSH, bench signs')
    runm(b, col='g_bench', mask=harsh, group='folio', label='B13 HARSH, bench, page-shuffle null')
    harsh2 = harsh & b.folio.isin(clean_pages)
    runm(b, col='g_bench', mask=harsh2, label='B14 HARSH + page has no <-> at all')

    print('\n#### C. strata')
    for col, name in [('H', 'hand'), ('L', 'Currier lang'), ('I', 'section'), ('C', 'Currier hand($C)')]:
        for v in sorted(b[col].astype(str).unique()):
            m = (b[col].astype(str) == v)
            if m.sum() < 30:
                continue
            runm(b, mask=m, label='C %s=%s' % (name, v))
    print('\n#### C2. per quire (lag 1 only)')
    rows = []
    for v in sorted(b.Q.astype(str).unique()):
        m = (b.Q.astype(str) == v)
        r = runm(b, mask=m, lags=(1,), quiet=True)
        if 1 in r:
            r = r[1]
            rows.append((v, r['N'], r['obs'], r['exp'], r['rate'], r['exp_rate'], r['z']))
            print('  quire %-2s N=%4d obs=%3d exp=%6.1f  rate=%5.1f%% vs %5.1f%%  z=%+5.2f' % (v, r['N'], r['obs'], r['exp'], 100 * r['rate'], 100 * r['exp_rate'], r['z']))
    print('  quires with obs<exp: %d of %d' % (sum(1 for x in rows if x[2] < x[3]), len(rows)))
    print('\n#### C3. per page sign test (lag 1)')
    lo = hi = eq = 0
    for f, sub in b.groupby('folio', sort=False):
        r = runm(sub, lags=(1,), quiet=True, nperm=50)
        if 1 in r and r[1]['N'] >= 5:
            if r[1]['obs'] < r[1]['exp'] - 1e-9:
                lo += 1
            elif r[1]['obs'] > r[1]['exp'] + 1e-9:
                hi += 1
            else:
                eq += 1
    print('  pages (>=5 pairs) with obs<exp: %d, obs>exp: %d, equal: %d' % (lo, hi, eq))

    print('\n#### D. by position of the pair inside the paragraph (pair = lines li, li+1)')
    for lo_, hi_ in [(2, 2), (3, 3), (4, 4), (5, 6), (7, 9), (10, 99)]:
        # restrict the *upper* line index; implement by counting only those pairs
        d = b.reset_index(drop=True)
        g, _ = encode(list(d.g_eva))
        _, grp = np.unique(d.pid.values, return_inverse=True)
        a, bb = pairs_within(d.pid.values, 1, mask=(g >= 0))
        sel = (d.li.values[a] >= lo_) & (d.li.values[a] <= hi_)
        r = permtest(g, grp, a[sel], bb[sel], nperm=NPERM, seed=1)
        print('  upper line index %2d-%2d: %s' % (lo_, hi_, fmt(r)))
    d = b.reset_index(drop=True)
    g, _ = encode(list(d.g_eva))
    _, grp = np.unique(d.pid.values, return_inverse=True)
    a, bb = pairs_within(d.pid.values, 1, mask=(g >= 0))
    sel = d['last'].values[bb]
    print('  pair ending in paragraph-last line:   %s' % fmt(permtest(g, grp, a[sel], bb[sel], nperm=NPERM, seed=1)))
    print('  pair not involving last line:         %s' % fmt(permtest(g, grp, a[~sel], bb[~sel], nperm=NPERM, seed=1)))

    print('\n#### H. alternative readings: take 2nd alternative in [a:b]')
    df2, *_ = build(alt=1)
    b2 = df2[df2.li >= 2].reset_index(drop=True)
    runm(b2, label='H alt reading 2')
