"""Apply the same geometric criterion to Voynich word boxes as to the real MSS:
a line that stops early (right extent < thr of the 3 preceding lines' margin), starts at the left
margin, preceded by 3 full lines.  No use of ZL paragraph markers for the selection."""
import numpy as np, pickle
from collections import defaultdict, Counter
import zl, boxes

rows = zl.read()
bypage = defaultdict(list)
for r in rows:
    bypage[r['page']].append(r)

recs = []
for pg, rr in bypage.items():
    out, rate = boxes.align_page(pg, rr)
    if out is None:
        continue
    seq = []
    for r in rr:
        if r['typ'][0] != 'P':
            seq.append(None); continue          # break the run at non-paragraph loci
        bx = out.get(r['num'], [])
        if not bx:
            seq.append(None); continue
        seq.append(dict(page=pg, num=r['num'], typ=r['typ'], pend=r['pend'], pstart=r['pstart'], I=r['I'],
                        x0=min(b[1] for b in bx), x1=max(b[1] + b[3] for b in bx), nb=len(bx),
                        zc=len(zl.words(r['raw'])), cut='<->' in r['raw'], raw=r['raw']))
    for i in range(3, len(seq)):
        if seq[i] is None or any(s is None for s in seq[i - 3:i]):
            continue
        prev = seq[i - 3:i]
        left = np.median([p['x0'] for p in prev]); right = np.median([p['x1'] for p in prev]); W = right - left
        if W <= 0:
            continue
        pr = [(p['x1'] - left) / W for p in prev]; ps = [(p['x0'] - left) / W for p in prev]
        prev_full = all(x >= 0.9 for x in pr) and all(abs(x) < 0.15 for x in ps) and not any(p['pend'] for p in prev[:-0] if p is not prev[-1])
        # previous line must not itself be a paragraph end
        prev_full = prev_full and not prev[-1]['pend'] and not prev[-2]['pend']
        cur = seq[i]
        recs.append(dict(page=pg, num=cur['num'], typ=cur['typ'], I=cur['I'], r=(cur['x1'] - left) / W, s=(cur['x0'] - left) / W,
                         nb=cur['nb'], zc=cur['zc'], pend=cur['pend'], prev_full=prev_full, cut=cur['cut'], raw=cur['raw'][:40]))
import pandas as pd
df = pd.DataFrame(recs)
print('candidate lines with 3 preceding P-lines:', len(df), 'with prev 3 full:', int(df.prev_full.sum()))
d = df[df.prev_full]
for thr in (0.9, 0.85, 0.75):
    e = d[(d.r < thr) & (d.s.abs() < 0.15)]
    print(f'thr={thr}: early-stopping lines n={len(e)}  single-box={int((e.nb==1).sum())} ({100*np.mean(e.nb==1):.1f}%)  <=2={int((e.nb<=2).sum())} ({100*np.mean(e.nb<=2):.1f}%) '
          f' ZL-single={int((e.zc==1).sum())}  marked <$> by ZL={100*np.mean(e.pend):.0f}%  median r={e.r.median():.2f}')
    if thr == 0.85:
        h, _ = np.histogram(e.r, bins=[-1, .2, .4, .6, .8, 1.0]); print('   r hist', list(h))
        print('   box-count dist', sorted(Counter(e.nb).items()))
        print('   not marked as paragraph end:', len(e[~e.pend]), ' of which single-box', int((e[~e.pend].nb == 1).sum()), ' <=2', int((e[~e.pend].nb <= 2).sum()))
        print(e[~e.pend][['page', 'num', 'typ', 'r', 'nb', 'zc', 'cut', 'raw']].head(40).to_string())
        print('   single-box rows:'); print(e[e.nb == 1][['page', 'num', 'typ', 'r', 'nb', 'zc', 'pend', 'raw']].to_string())
# relaxed: lines not starting at left margin (right-justified / indented endings)
e = d[(d.s > 0.3)]
print('lines starting far right of margin (s>0.3) after 3 full lines:', len(e)); print(e[['page', 'num', 'typ', 's', 'r', 'nb', 'pend', 'raw']].head(40).to_string())
