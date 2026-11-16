"""Hunt for short paragraph endings filed as something else."""
import numpy as np, re
from collections import defaultdict, Counter
import zl, boxes

rows = zl.read()
# 1. what follows a <$> line on the same page?
foll = Counter(); cand = []
for i, r in enumerate(rows[:-1]):
    if r['typ'][0] == 'P' and r['pend']:
        nx = rows[i + 1]
        if nx['page'] != r['page']:
            foll['(page end)'] += 1; continue
        foll[(nx['pos'], nx['typ'], nx['pstart'])] += 1
        if not (nx['typ'][0] == 'P' and nx['pstart']):
            cand.append((r, nx))
print('what follows a paragraph-end locus:'); [print('  ', k, v) for k, v in foll.most_common()]
print('\nnon-paragraph-start loci directly after a <$>: ', len(cand))
for r, nx in cand:
    print(f"   {r['page']}.{r['num']} ({len(zl.words(r['raw']))}w) -> {nx['num']},{nx['pos']}{nx['typ']} ({len(zl.words(nx['raw']))}w): {nx['raw'][:50]}")

# 2. P-type loci that are not inside <%>..<$>
paras, an = zl.paragraphs(rows)
print('\nP-loci outside paragraph markers:', an)

# 3. P-type loci by type and position code, word counts
print('\nPt/Pc/Pr lines: word counts', Counter((r['typ'], len(zl.words(r['raw']))) for r in rows if r['typ'] in ('Pt', 'Pc', 'Pr')))
print('Pt/Pc/Pr without <$>:', [(r['page'], r['num'], r['typ']) for r in rows if r['typ'] in ('Pt', 'Pc', 'Pr') and not r['pend']])
print('paragraph-final loci by position code:', Counter(r['pos'] for r in rows if r['typ'][0] == 'P' and r['pend']))

# 4. all single-word P-type loci anywhere (not just last lines)
sing = [(r['page'], r['num'], r['pos'] + r['typ'], r['pend'], r['pstart'], zl.clean(r['raw'])) for r in rows if r['typ'][0] == 'P' and len(zl.words(r['raw'])) == 1]
print('\nall single-word P-type loci:', len(sing))
for s in sing:
    print('  ', s)

# 5. x-reset inside a locus = locus spans two visual lines (possible merged run-over)
bypage = defaultdict(list)
for r in rows:
    bypage[r['page']].append(r)
print('\nparagraph loci whose boxes contain an x-reset (word lower-left of previous word inside same locus):')
nres = 0; tot = 0; reslast = 0; totlast = 0
for pg, rr in bypage.items():
    out, rate = boxes.align_page(pg, rr)
    if out is None:
        continue
    for r in rr:
        if r['typ'][0] != 'P':
            continue
        bx = out.get(r['num'], [])
        if len(bx) < 2:
            continue
        tot += 1; totlast += r['pend']
        hmed = np.median([b[4] for b in bx])
        resets = [k for k in range(1, len(bx)) if bx[k][1] < bx[k - 1][1] - 30 and bx[k][2] > bx[k - 1][2] + 0.6 * hmed]
        if resets:
            nres += 1; reslast += r['pend']
            tail = len(bx) - resets[-1]
            if r['pend'] or tail <= 2:
                print(f"   {pg}.{r['num']} {r['typ']} pend={r['pend']} boxes={len(bx)} zlwords={len(zl.words(r['raw']))} tail_after_reset={tail} : {[b[0] for b in bx[resets[-1]:]][:4]}")
print('loci with internal reset', nres, '/', tot, ' among paragraph-final loci', reslast, '/', totlast)
