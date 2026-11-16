"""Stars section: paragraphs are anchored on marginal stars. Check star/paragraph correspondence and last lines."""
import re, numpy as np
from collections import Counter
import zl
lines = open(zl.ZL, encoding='latin-1').read().split('\n')
page = None; I = None
star_follow = Counter(); nstar = 0; npar = 0; par_without_star = 0
prev_star = False
for ln in lines:
    m = zl.HDR.match(ln)
    if m:
        page = m.group(1); I = dict(re.findall(r'\$(\w)=(\S)', m.group(2))).get('I'); prev_star = False; continue
    if I != 'S':
        continue
    if ln.startswith('#'):
        if re.search(r'\bstar\b', ln, re.I) and 'points' in ln:
            nstar += 1; prev_star = True
        continue
    m = zl.LOC.match(ln)
    if m:
        st = '<%>' in ln
        if prev_star:
            star_follow['star then <%>' if st else 'star then NON-start line'] += 1
        elif st:
            par_without_star += 1
        npar += st
        prev_star = False
print('S section: star comments', nstar, 'paragraph starts', npar, dict(star_follow), 'paragraph starts without preceding star comment', par_without_star)
rows = zl.read(); paras, an = zl.paragraphs(rows)
S = [p for p in paras if p['lines'][0]['I'] == 'S']
wc = [[len(zl.words(l['raw'])) for l in p['lines']] for p in S]
print('S paragraphs', len(S), 'line-count dist', sorted(Counter(len(w) for w in wc).items()))
last = [w[-1] for w in wc]; rat = [w[-1] / np.median(w[:-1]) for w in wc]
print('last-line words dist', sorted(Counter(last).items()))
print('ratio hist [0-.2,.2-.4,.4-.6,.6-.8,.8-1,1-1.2,1.2+]', list(np.histogram(rat, bins=[0, .2, .4, .6, .8, 1.0, 1.2, 10])[0]), 'median', np.median(rat))
print('mean words on non-last lines', np.mean([x for w in wc for x in w[:-1]]))
