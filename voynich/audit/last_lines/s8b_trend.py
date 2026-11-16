"""Word-count of lines by position from paragraph end; control for the penultimate effect."""
import numpy as np
from collections import defaultdict
import zl
rows = zl.read()
paras, an = zl.paragraphs(rows)
def chars(l): 
    w = zl.words(l['raw']); return np.mean([len(x) for x in w]) if w else np.nan
for name, sel in [('ALL', lambda p: True), ('stars S', lambda p: p['lines'][0]['I'] == 'S'), ('herbal H', lambda p: p['lines'][0]['I'] == 'H'),
                  ('bio B', lambda p: p['lines'][0]['I'] == 'B'), ('no drawing cuts', lambda p: not any('<->' in l['raw'] for l in p['lines']))]:
    D = defaultdict(list); C = defaultdict(list); bylast = defaultdict(lambda: defaultdict(list))
    for p in paras:
        L = p['lines']
        if len(L) < 6 or not sel(p):
            continue
        wc = [len(zl.words(l['raw'])) for l in L]
        ref = np.mean(wc[1:-3])  # interior lines excluding first and the last three
        for k in (1, 2, 3):
            D[k].append(wc[-k] - ref); C[k].append(chars(L[-k]))
        D['first'].append(wc[0] - ref)
        C['int'].append(np.nanmean([chars(l) for l in L[1:-3]]))
        b = '1-2' if wc[-1] <= 2 else '3-4' if wc[-1] <= 4 else '5-6' if wc[-1] <= 6 else '7+'
        bylast[b]['pen'].append(wc[-2] - ref); bylast[b]['ante'].append(wc[-3] - ref)
    n = len(D[1])
    if n == 0: continue
    f = lambda v: f'{np.mean(v):+.2f}±{np.std(v)/np.sqrt(len(v)):.2f}'
    print(f'{name:16s} n={n:3d} (6+ lines) words minus interior mean:  first {f(D["first"])}  antepenult {f(D[3])}  penult {f(D[2])}  last {f(D[1])} | chars/word int {np.mean(C["int"]):.2f} ante {np.nanmean(C[3]):.2f} pen {np.nanmean(C[2]):.2f} last {np.nanmean(C[1]):.2f}')
    for b in ('1-2', '3-4', '5-6', '7+'):
        if bylast[b]['pen']:
            print(f'      last={b:4s} n={len(bylast[b]["pen"]):3d}  penult {f(bylast[b]["pen"])}  antepenult {f(bylast[b]["ante"])}')
