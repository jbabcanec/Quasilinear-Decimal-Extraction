"""Geometric re-wrap null from word boxes (own implementation).
Greedy wrap of each paragraph's actual word widths into lines of width W*f, where W is the
paragraph's median non-last line extent.  Also a phase null: f=1 with a random first-line indent."""
import numpy as np, pickle
recs = pickle.load(open('box_paras.pkl', 'rb'))
rng = np.random.default_rng(99)

P = []
for r in recs:
    L = r['lines']
    if len(L) < 4 or not all(l['bc'] > 0 for l in L):
        continue
    if any('<->' in l['raw'] for l in L):
        continue          # drawing gaps make extent != text width
    ww = [b[3] for l in L for b in l['boxes']]
    g = []
    for l in L:
        bx = l['boxes']
        g += [bx[k + 1][1] - (bx[k][1] + bx[k][3]) for k in range(len(bx) - 1)]
    g = [x for x in g if -15 < x < 80]
    gap = max(0.0, float(np.median(g))) if g else 5.0
    W = float(np.median([l['x1'] - l['x0'] for l in L[:-1]]))
    ext = [l['x1'] - l['x0'] for l in L[:-1]]
    P.append(dict(ww=ww, gap=gap, W=W, last=L[-1]['bc'], n=len(L), ext=ext, I=r['I'],
                  last_r=(L[-1]['x1'] - L[-1]['x0']) / W))
print('paragraphs (4+ lines, no drawing cut, boxes everywhere):', len(P))
obs1 = sum(p['last'] == 1 for p in P); obs2 = sum(p['last'] <= 2 for p in P)
print(f'observed: single-box last {obs1} ({100*obs1/len(P):.1f}%), <=2 {obs2} ({100*obs2/len(P):.1f}%)')


def wrap(ww, gap, W, indent=0.0, caps=None):
    """returns (n_lines, words on last line, extent of last line)"""
    x = indent; n = 0; lines = 1; k = 0
    Wc = W if caps is None else caps[0]
    for w in ww:
        need = w if n == 0 else gap + w
        if n > 0 and x + need > Wc:
            lines += 1; x = w; n = 1
            if caps is not None:
                k += 1; Wc = caps[k % len(caps)]
        else:
            x += need; n += 1
    return lines, n, x


for f in (0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3):
    one = two = m = 0; rr = []
    for p in P:
        nl, n, x = wrap(p['ww'], p['gap'], p['W'] * f)
        if nl < 2:
            continue
        m += 1; one += n == 1; two += n <= 2; rr.append(x / (p['W'] * f))
    h, _ = np.histogram(rr, bins=[0, .2, .4, .6, .8, 1.01])
    print(f'rewrap f={f:.1f}: n={m} single={one} ({100*one/m:.1f}%) <=2={two} ({100*two/m:.1f}%)  last-extent hist share={np.round(h/h.sum(),2)}')

# phase null: same width, random indent; also resampled actual line extents as capacities
for label, usecaps in (('phase null f=1, random indent', False), ('phase null, capacities = paragraph\'s own actual line extents (shuffled), random indent', True)):
    ones = []; twos = []; hh = np.zeros(5)
    for rep in range(1000):
        one = two = 0; rr = []
        for p in P:
            caps = None
            if usecaps:
                caps = list(rng.permutation(p['ext']))
            W0 = p['W'] if caps is None else caps[0]
            nl, n, x = wrap(p['ww'], p['gap'], p['W'], indent=rng.uniform(0, W0), caps=caps)
            one += n == 1; two += n <= 2; rr.append(x / p['W'])
        ones.append(one); twos.append(two); hh += np.histogram(rr, bins=[0, .2, .4, .6, .8, 5])[0]
    print(f'{label}: single exp={np.mean(ones):.1f} [{np.percentile(ones,2.5):.0f},{np.percentile(ones,97.5):.0f}] ({100*np.mean(ones)/len(P):.1f}%)  '
          f'<=2 exp={np.mean(twos):.1f} [{np.percentile(twos,2.5):.0f},{np.percentile(twos,97.5):.0f}]  min={np.min(ones)}  last-extent hist share={np.round(hh/hh.sum(),2)}')
h, _ = np.histogram([p['last_r'] for p in P], bins=[0, .2, .4, .6, .8, 5]); print('observed last-extent hist share', np.round(h / h.sum(), 2))
