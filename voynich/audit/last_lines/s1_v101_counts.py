import re, numpy as np
from collections import Counter
V = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/v101/voyn_101.txt'
LOC = re.compile(r'^<([^.>]+)\.([^>]+)>(.*)$')


def load(numeric_only=True):
    paras = []; cur = []; curpage = None
    odd = Counter()
    for raw in open(V, encoding='latin-1'):
        raw = raw.strip()
        if not raw:
            continue
        m = LOC.match(raw)
        if not m:
            odd['nomatch'] += 1; continue
        page, num, text = m.groups()
        isnum = num.isdigit()
        if numeric_only and not isnum:
            if cur:
                paras.append(cur + ['UNCLOSED']); cur = []
            continue
        if page != curpage and cur:
            paras.append(cur + ['UNCLOSED']); cur = []
        curpage = page
        text = text.strip()
        term = text[-1] if text else ''
        body = text.rstrip('-=').strip()
        cur.append((page, num, body))
        if term == '=':
            paras.append(cur); cur = []
        elif term != '-':
            odd['noterm'] += 1
    return paras, odd


def wc(body, comma=False):
    t = body.replace(',', '.' if comma else '')
    return len([w for w in t.split('.') if w.strip()])


paras, odd = load()
print('paragraph-like units', len(paras), odd, 'unclosed', sum(1 for p in paras if p[-1] == 'UNCLOSED'))
paras = [p for p in paras if p[-1] != 'UNCLOSED']
print('line-count dist', sorted(Counter(len(p) for p in paras).items()))
for minl in (4, 3, 2):
    for comma in (False, True):
        n = one = two = 0; ratios = []; lastw = []; singles = []
        for p in paras:
            if len(p) < minl:
                continue
            w = [wc(b, comma) for _, _, b in p]
            n += 1; one += w[-1] <= 1; two += w[-1] <= 2
            if w[-1] <= 1:
                singles.append(p[-1])
            ratios.append(w[-1] / np.median(w[:-1])); lastw.append(w[-1])
        h, _ = np.histogram(ratios, bins=[0, .2, .4, .6, .8, 1.0, 1.2, 10])
        print(f'v101 min{minl} comma={comma}: n={n} 1w={one} ({100*one/n:.1f}%) <=2w={two} ({100*two/n:.1f}%) '
              f'ratio med={np.median(ratios):.2f} mean={np.mean(ratios):.2f} hist={list(map(int,h))}')
        if not comma:
            print('   last wc dist', sorted(Counter(lastw).items()))
            print('   singles', [(a, b) for a, b, c in singles])
