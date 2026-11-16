import numpy as np, sys
from collections import Counter
import zl

rows = zl.read()


def stats(paras, minlines, comma, label, typical='median_nonlast', verbose=False, lines_of=None):
    n = 0; one = 0; two = 0; ratios = []; lastw = []
    singles = []
    for p in paras:
        L = lines_of(p) if lines_of else p['lines']
        if len(L) < minlines:
            continue
        wc = [len(zl.words(l['raw'], comma)) for l in L]
        if typical == 'median_nonlast':
            typ = np.median(wc[:-1])
        elif typical == 'median_interior':
            typ = np.median(wc[1:-1]) if len(wc) > 2 else np.median(wc[:-1])
        elif typical == 'mean_nonlast':
            typ = np.mean(wc[:-1])
        n += 1
        lastw.append(wc[-1])
        if wc[-1] <= 1:
            one += 1; singles.append((L[-1]['page'], L[-1]['num'], L[-1]['typ'], L[-1]['raw']))
        if wc[-1] <= 2:
            two += 1
        ratios.append(wc[-1] / typ if typ > 0 else np.nan)
    ratios = np.array(ratios)
    h, _ = np.histogram(ratios, bins=[0, .2, .4, .6, .8, 1.0, 1.2, 10])
    print(f'{label:55s} n={n:4d} 1w={one:3d} ({100*one/n:4.1f}%) <=2w={two:3d} ({100*two/n:4.1f}%) '
          f'ratio med={np.nanmedian(ratios):.2f} mean={np.nanmean(ratios):.2f} hist[0-.2,.2-.4,.4-.6,.6-.8,.8-1,1-1.2,1.2+]={list(h)}')
    if verbose:
        for s in singles:
            print('     ', s)
    return n, one, two, ratios, lastw


for ptypes, name in [(('P0', 'P1', 'Pb', 'Pt', 'Pc', 'Pr'), 'A: all P* lines (Pt/Pc/Pr count as last line)'),
                     (('P0', 'P1', 'Pb'), 'B: drop Pt/Pc/Pr lines')]:
    paras, an = zl.paragraphs(rows, ptypes)
    print('==', name, 'paragraphs', len(paras))
    for minl in (4, 3, 2):
        for comma in (False, True):
            stats(paras, minl, comma, f'{name[:1]} min{minl} comma_breaks={comma}', verbose=(minl == 2 and not comma))
    stats(paras, 4, False, 'A/B min4 typical=interior median', typical='median_interior')
    stats(paras, 4, False, 'A/B min4 typical=mean nonlast', typical='mean_nonlast')

# C: merge Pt/Pc/Pr into preceding line
paras, an = zl.paragraphs(rows)


def merged(p):
    L = []
    for l in p['lines']:
        if l['typ'] in ('Pt', 'Pc', 'Pr') and L:
            prev = dict(L[-1]); prev['raw'] = prev['raw'] + '.' + l['raw']; L[-1] = prev
        else:
            L.append(l)
    return L


print('== C: merge Pt/Pc/Pr into preceding line')
for minl in (4, 3):
    stats(paras, minl, False, f'C min{minl}', lines_of=merged)

# distribution of last-line word counts
n, one, two, ratios, lastw = stats(paras, 4, False, 'A min4 again')
print('last-line word-count distribution (A, min4):', sorted(Counter(lastw).items()))
allw = []
for p in paras:
    if len(p['lines']) >= 4:
        allw += [len(zl.words(l['raw'])) for l in p['lines'][:-1]]
print('non-last line word counts: mean', np.mean(allw), 'median', np.median(allw), sorted(Counter(allw).items()))
