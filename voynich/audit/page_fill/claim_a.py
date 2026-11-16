"""Claim (a): does the last text line of each page end a paragraph?  Plus a check
that this is not merely transcriber convention: page-final lines should look like
paragraph-final lines (short), page-first lines like paragraph-initial lines (gallows)."""
import numpy as np, collections
from parse import load, paragraphs
P = load()
n = ok = 0; exc = []
for pg, d in P.items():
    pl = [L for L in d['loci'] if L['typ'][0] == 'P']
    if not pl: continue
    n += 1
    if pl[-1]['pend']: ok += 1
    else: exc.append((pg, pl[-1]['num'], pl[-1]['typ'], [L['typ'] for L in pl[-4:]]))
print('pages with P-type loci:', n, ' last P line carries <$>:', ok)
print('exceptions:', exc)
print('pages w/o boxes among those:', [pg for pg in P if any(L['typ'][0]=='P' for L in P[pg]['loci']) and pg in ('fRos','f101r','f101v','f116v')])
# also first P line carries <%>
nf = sum(1 for pg, d in P.items() if [L for L in d['loci'] if L['typ'][0]=='P'] and [L for L in d['loci'] if L['typ'][0]=='P'][0]['pstart'])
print('first P line carries <%>:', nf)

# character counts of lines by role (P0 lines only, multi-line paragraphs)
def nch(L): return sum(len(w) for w in L['words']) + len(L['words']) - 1
roles = collections.defaultdict(list)
gal = collections.defaultdict(list)
for pg, d in P.items():
    ps = paragraphs(d)
    ps = [p for p in ps if all(L['typ'] == 'P0' for L in p['lines'])]
    for k, p in enumerate(ps):
        if len(p['lines']) < 3: continue
        body = np.median([nch(L) for L in p['lines'][1:-1]])
        last = nch(p['lines'][-1]) / body
        key = 'page-final para' if k == len(ps) - 1 else 'other para'
        roles[key].append(last)
        sec = d['meta'].get('I')
        roles[key + ' ' + ('SBT' if sec in 'SBT' else 'other')].append(last)
    for k, p in enumerate(ps):
        w = p['lines'][0]['words']
        if not w: continue
        g = w[0][0] in 'ptkf'
        gal['page-first para' if k == 0 else 'later para'].append(g)
        for L in p['lines'][1:]:
            if L['words']: gal['non-initial line'].append(L['words'][0][0] in 'ptkf')
for k, v in sorted(roles.items()):
    v = np.array(v)
    print(f'{k:28s} n={len(v):4d}  last-line length / body-line length: median {np.median(v):.2f}  mean {v.mean():.2f}  frac>=0.9: {(v>=0.9).mean():.2f}')
for k, v in gal.items():
    print(f'{k:18s} n={len(v):5d} first glyph is gallows (p,t,k,f): {np.mean(v):.2f}')
