"""Is the pixel 'foot' a property of the page or of the image frame?
Evidence that every image is a content crop scaled to fit ~640 x ~900."""
import numpy as np, collections
from math import comb
from parse import load, boxes
P = load()
rows = []
for pg, d in P.items():
    b = boxes(pg)
    if not b: continue
    xs = [q[1] + q[3] for q in b]; ys = [q[2] + q[4] for q in b]
    rows.append((pg, d['meta'].get('I'), d['meta'].get('Q'), min(q[1] for q in b), max(xs), min(q[2] for q in b), max(ys)))
ymin = np.array([r[5] for r in rows]); xmax = np.array([r[4] for r in rows]); ymax = np.array([r[6] for r in rows])
print('all %d pages: box_ymin<=6 on %d ; hist of ymin for <=12:' % (len(rows), (ymin <= 6).sum()), collections.Counter(ymin[ymin <= 12].tolist()))
print('box_xmax in [633,642]: %d ; box_ymax in [894,906]: %d ; xmax>642: %s' % (((xmax >= 633) & (xmax <= 642)).sum(), ((ymax >= 894) & (ymax <= 906)).sum(), [r[0] for r in rows if r[4] > 642]))
print('ymax>906:', [(r[0], r[6]) for r in rows if r[6] > 906])
print('hist ymax>=870:', sorted(collections.Counter(ymax[ymax >= 870].tolist()).items()))
print('hist xmax>=600:', sorted(collections.Counter(xmax[(xmax >= 600) & (xmax <= 642)].tolist()).items()))
T = [r for r in rows if r[2] == 'T']
a = sum(1 for r in T if r[4] >= 633 and r[6] >= 894); b_ = sum(1 for r in T if r[4] < 633 and r[6] >= 894)
c = sum(1 for r in T if r[4] >= 633 and r[6] < 894); d_ = sum(1 for r in T if r[4] < 633 and r[6] < 894)
print('quire T (text-only stars pages): x-limited&y-limited=%d  y-limited only=%d  x-limited only=%d  neither=%d' % (a, b_, c, d_))
n = a + b_ + c + d_; r1 = a + b_; c1 = a + c
p = sum(comb(c1, k) * comb(n - c1, r1 - k) for k in range(0, a + 1)) / comb(n, r1)
print('Fisher exact (one-sided, as few or fewer both-limited): p = %.5f' % p)
for r in T: print('  ', r, 'aspect if content-fit: %.3f' % (r[6] / r[4]))
# recto/verso
