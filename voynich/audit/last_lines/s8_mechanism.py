"""Mechanism tests: is the scribe re-balancing (penultimate line shortened / last line stretched)?"""
import numpy as np, pickle
from collections import defaultdict, Counter
import zl

recs = pickle.load(open('box_paras.pkl', 'rb'))


def gaps(bx):
    g = [bx[k + 1][1] - (bx[k][1] + bx[k][3]) for k in range(len(bx) - 1)]
    return [x for x in g if -20 < x < 60]   # ignore drawing gaps / resets


rows = []
for r in recs:
    L = r['lines']
    if len(L) < 5 or not all(l['bc'] > 0 for l in L):
        continue
    inter = L[1:-2]           # interior lines excluding first, penultimate and last
    if len(inter) < 2:
        continue
    left = np.median([l['x0'] for l in L[:-1]])
    right = np.median([l['x1'] for l in inter])
    W = right - left
    if W <= 0:
        continue
    pen = L[-2]; last = L[-1]
    g_int = [x for l in inter for x in gaps(l['boxes'])]
    g_last = gaps(last['boxes']); g_pen = gaps(pen['boxes'])
    ww_int = [b[3] for l in inter for b in l['boxes']]
    rows.append(dict(I=r['I'], n=len(L),
                     pen_r=(pen['x1'] - left) / W, pen_w=pen['zc'], int_w=np.median([l['zc'] for l in inter]),
                     int_r=np.median([(l['x1'] - left) / W for l in inter]),
                     int_r_sd=np.std([(l['x1'] - left) / W for l in inter]),
                     last_r=(last['x1'] - left) / W, last_w=last['zc'],
                     last_s=(last['x0'] - left) / W,
                     gap_int=np.mean(g_int) if g_int else np.nan, gap_last=np.mean(g_last) if g_last else np.nan,
                     gap_pen=np.mean(g_pen) if g_pen else np.nan,
                     ww_int=np.mean(ww_int), ww_last=np.mean([b[3] for b in last['boxes']]),
                     ww_pen=np.mean([b[3] for b in pen['boxes']]),
                     cl_int=np.mean([len(b[0]) for l in inter for b in l['boxes']]),
                     cl_last=np.mean([len(b[0]) for b in last['boxes']])))
import pandas as pd
df = pd.DataFrame(rows)
print('paragraphs (5+ lines, boxes on every line):', len(df))
print('penultimate line: right extent rel. interior median: median %.3f mean %.3f ; interior lines by construction 1.0' % (df.pen_r.median(), df.pen_r.mean()))
print('   share of penultimate lines stopping >10%% short of interior margin: %.3f ; interior-line sd of r (typical) %.3f' % ((df.pen_r < 0.9).mean(), df.int_r_sd.median()))
print('penultimate words - interior median words: mean %.2f (se %.2f)' % ((df.pen_w - df.int_w).mean(), (df.pen_w - df.int_w).std() / np.sqrt(len(df))))
for lo, hi in [(1, 2), (3, 4), (5, 6), (7, 99)]:
    d = df[(df.last_w >= lo) & (df.last_w <= hi)]
    print(f'   last line {lo}-{hi} words: n={len(d)}  pen_r median={d.pen_r.median():.3f}  pen_w-int_w mean={(d.pen_w-d.int_w).mean():+.2f}  pen gap/int gap={np.nanmean(d.gap_pen/d.gap_int):.2f}')
print('word gap: last/interior mean ratio %.3f (median %.3f); penult/interior %.3f' % (np.nanmean(df.gap_last / df.gap_int), np.nanmedian(df.gap_last / df.gap_int), np.nanmean(df.gap_pen / df.gap_int)))
print('word width: last/interior %.3f ; penult/interior %.3f ; chars per word last %.2f vs interior %.2f' % ((df.ww_last / df.ww_int).mean(), (df.ww_pen / df.ww_int).mean(), df.cl_last.mean(), df.cl_int.mean()))
print('last line start offset (indent) median %.3f, share starting >10%% right of margin %.3f' % (df.last_s.median(), (df.last_s > 0.1).mean()))
# gap by last-line fill
for lo, hi in [(0, .4), (.4, .6), (.6, .8), (.8, 1.3)]:
    d = df[(df.last_r >= lo) & (df.last_r < hi)]
    print(f'   last_r {lo}-{hi}: n={len(d)} gap_last/gap_int={np.nanmean(d.gap_last/d.gap_int):.2f} ww_last/ww_int={(d.ww_last/d.ww_int).mean():.2f}')
