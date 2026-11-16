"""Are there hidden unit boundaries (paragraph starts) inside transcribed paragraphs?
Paragraph-initial words overwhelmingly begin with a gallows, esp. p/f. Look for p/f-initial words by
position in NON-first lines, and at the lines beside 'unbroken' stars."""
import re, numpy as np
from collections import Counter, defaultdict
import zl
rows = zl.read(); paras, an = zl.paragraphs(rows)
first_init = Counter(); 
pos_pf = defaultdict(lambda: [0, 0])   # position -> [p/f-initial count, total]
pos_pf_first = defaultdict(lambda: [0, 0])
for p in paras:
    for j, l in enumerate(p['lines']):
        w = zl.words(l['raw'])
        if not w: continue
        if j == 0:
            first_init[w[0][0]] += 1
        for k, x in enumerate(w):
            d = pos_pf_first if j == 0 else pos_pf
            kk = min(k, 9)
            d[kk][1] += 1; d[kk][0] += x[0] in 'pf'
print('paragraph-initial word first letter:', first_init.most_common(8), ' share p/f = %.2f, any gallows(ptkf) = %.2f' % (
    (first_init['p'] + first_init['f']) / sum(first_init.values()), sum(first_init[c] for c in 'ptkf') / sum(first_init.values())))
print('share of p/f-initial words by word position (0-based), NON-first lines:', {k: f'{100*v[0]/v[1]:.2f}% ({v[0]}/{v[1]})' for k, v in sorted(pos_pf.items())})
print('same, FIRST lines:', {k: f'{100*v[0]/v[1]:.1f}%' for k, v in sorted(pos_pf_first.items())})

# unbroken stars in S
lines = open(zl.ZL, encoding='latin-1').read().split('\n')
I = None; prev_star = False; prevloc = None; res = []
for ln in lines:
    m = zl.HDR.match(ln)
    if m:
        I = dict(re.findall(r'\$(\w)=(\S)', m.group(2))).get('I'); prev_star = False; prevloc = None; continue
    if I != 'S': continue
    if ln.startswith('#'):
        if re.search(r'\bstar\b', ln, re.I) and 'points' in ln: prev_star = True
        continue
    m = zl.LOC.match(ln)
    if m:
        if prev_star and '<%>' not in ln:
            w = zl.words(m.group(5)); pw = zl.words(prevloc.group(5)) if prevloc else []
            res.append((m.group(1), m.group(2), w[0][0], [k for k, x in enumerate(w) if x[0] in 'pf'][:3], len(pw), len(w)))
        prev_star = False; prevloc = m
print('\nunbroken stars (star beside a line that is not a transcribed paragraph start):', len(res))
print('  first letter of that line:', Counter(r[2] for r in res).most_common())
print('  p/f-initial word positions in that line (first 3):', Counter(tuple(r[3]) for r in res).most_common(8))
print('  words on previous line (median):', np.median([r[4] for r in res]))
