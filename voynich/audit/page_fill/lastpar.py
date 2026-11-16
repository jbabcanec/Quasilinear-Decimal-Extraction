import numpy as np
from fill import R, GT20, GS
from audit import run, rule_q
# paragraph lengths by position on page, quire T
first = [R[pg]['plen'][0] for pg in GT20]; last = [R[pg]['plen'][-1] for pg in GT20]
mid = [L for pg in GT20 for L in R[pg]['plen'][1:-1]]
allp = [L for pg in GT20 for L in R[pg]['plen']]
print('quire T paragraphs: n=%d mean %.2f lines; page-first mean %.2f; page-last mean %.2f; middle mean %.2f' % (len(allp), np.mean(allp), np.mean(first), np.mean(last), np.mean(mid)))
print('length distribution:', np.bincount(allp))
print('page-last lengths:', last); print('page-first lengths:', first)
# permutation test for last-paragraph mean within page (shuffle order within page)
rng = np.random.default_rng(0)
obs = np.mean(last); sim = [np.mean([rng.choice(R[pg]['plen']) for pg in GT20]) for _ in range(20000)]
print('page-last mean vs random paragraph of same page: obs %.2f null %.2f+-%.2f p(two-sided)=%.3f' % (obs, np.mean(sim), np.std(sim), np.mean(np.abs(np.array(sim) - np.mean(sim)) >= abs(obs - np.mean(sim)))))
# greedy-stopping check: would the next page's first paragraph have fitted in the blank left (pixel frame at face value, P90 foot)?
out, g, mask = run('x', GT20, rule_q(90), modes=(), verbose=False)
fit = 0; tot = 0
for i, pg in enumerate(GT20[:-1]):
    nxt = GT20[i + 1]; L = R[nxt]['plen'][0]
    ok = L <= g[i] + 0.3
    tot += 1; fit += ok
    if g[i] > 1.0: print(f'  {pg}: blank {g[i]:.1f} lines; first paragraph of {nxt} is {L} lines -> {"WOULD HAVE FITTED" if ok else "did not fit"}')
print('pages where the next first paragraph would have fitted in the blank: %d of %d' % (fit, tot))
# P(gap < x lines) under greedy stopping, quire T lengths
plen = np.array(allp, float); gaps = []
for t in range(40000):
    cap = rng.uniform(40, 55); used = 0
    while True:
        L = plen[rng.integers(len(plen))]
        if used + L <= cap: used += L
        else: break
    gaps.append(cap - used)
gaps = np.array(gaps)
print('greedy stopping, gap distribution (lines): P(<0.5)=%.2f P(<1)=%.2f P(<1.5)=%.2f P(<2)=%.2f P(<2.5)=%.2f median %.2f' % tuple([(gaps < x).mean() for x in (0.5, 1, 1.5, 2, 2.5)] + [np.median(gaps)]))
