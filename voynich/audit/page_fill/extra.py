import numpy as np, re
from fill import R, G50, GS, GT20
import scalefree
from audit import run, rule_q, rule_max, leaf, TEXTONLY, frac_full

print('#### D. foot quantile x tolerance scan, pixel frame at face value, like-for-like null (draw mode)')
for nm, G in (('25 S', GS), ('50 SBT', G50)):
    for q in (75, 80, 85, 90, 95, 100):
        for tol in (0.5, 1.0):
            run(f'{nm}: P{q}, tol {tol}', G, rule_q(q), tol=tol, modes=('draw',), nsim=3000)

print('#### E. which pages are "full" (P90 of 25 S pages, tol 0.5) and frame status')
out, g, mask = run('25 S: P90 of group', GS, rule_q(90), modes=(), verbose=False)
for pg, x in zip(GS, g):
    r = R[pg]
    lim = ('H' if r['box_ymax'] >= 894 else '-') + ('W' if r['box_xmax'] >= 633 else '-')
    print(f"  {pg:6s} gap {x:5.2f} lines {'FULL' if x < 0.5 else '    '} frame-limit {lim}  box_ymax {r['box_ymax']} box_xmax {r['box_xmax']} n={r['n']} bif={r['meta']['B']}")
full = np.array(g) < 0.5; Hl = np.array([R[pg]['box_ymax'] >= 894 for pg in GS])
print('  full & H-limited %d, full & not %d, not full & H-limited %d, neither %d' % ((full & Hl).sum(), (full & ~Hl).sum(), (~full & Hl).sum(), (~full & ~Hl).sum()))

print('#### F. leave-one-leaf-out and leave-one-bifolio-out (25 S, P90, tol .5, like-for-like draw null)')
ps = []
for lf in sorted(set(leaf(pg) for pg in GS)):
    G = [pg for pg in GS if leaf(pg) != lf]
    o, _, _ = run('x', G, rule_q(90), modes=('draw',), nsim=1500, verbose=False)
    ps.append((lf, o['k'], round(o['obs'], 2), round(o['draw'][0], 2), round(o['draw'][2], 3)))
print('  leave-one-leaf-out (leaf, full, obs, null, p):', ps)
ps = []
for b in sorted(set(R[pg]['meta']['Q'] + R[pg]['meta']['B'] for pg in GS)):
    G = [pg for pg in GS if R[pg]['meta']['Q'] + R[pg]['meta']['B'] != b]
    o, _, _ = run('x', G, rule_q(90), modes=('draw',), nsim=1500, verbose=False)
    ps.append((b, o['k'], round(o['obs'], 2), round(o['draw'][0], 2), round(o['draw'][2], 3)))
print('  leave-one-bifolio-out:', ps)
# leaf-level clustering of fullness
fl = {}
for pg, x in zip(GS, g): fl.setdefault(leaf(pg), []).append(x < 0.5)
both = [v for v in fl.values() if len(v) == 2]
print('  leaves with both sides: %d ; both full %d ; one full %d ; none %d' % (len(both), sum(all(v) for v in both), sum(any(v) and not all(v) for v in both), sum(not any(v) for v in both)))

print('#### G. elastic line spacing (quire T)')
n = np.array([R[pg]['n'] for pg in GT20]); mp = np.array([R[pg]['mp'] for pg in GT20])
cw = np.array([R[pg]['cw'] for pg in GT20]); W = np.array([R[pg]['W'] for pg in GT20])
Hl = np.array([R[pg]['box_ymax'] >= 894 for pg in GT20])
print('  lines/page: min %d max %d ; among the 11 frame-bottom ("full") pages: %s' % (n.min(), n.max(), sorted(n[Hl].tolist())))
print('  pitch/charwidth: all %.2f..%.2f ; corr with n = %.2f | pitch/textwidth corr with n = %.2f' % ((mp/cw).min(), (mp/cw).max(), np.corrcoef(n, mp/cw)[0,1], np.corrcoef(n, mp/W)[0,1]))
def halves(pg):
    ys = np.array(R[pg]['ys']); k = len(ys)
    i = np.arange(k); h = k // 2
    top = np.polyfit(i[:h], ys[:h], 1)[0]; bot = np.polyfit(i[h:], ys[h:], 1)[0]
    last = np.polyfit(i[-8:], ys[-8:], 1)[0]; rest = np.polyfit(i[:-8], ys[:-8], 1)[0]
    return bot / top, last / rest
hv = np.array([halves(pg) for pg in GT20])
print('  bottom-half/top-half pitch: frame-bottom pages mean %.3f (sd %.3f) ; other pages mean %.3f (sd %.3f)' % (hv[Hl,0].mean(), hv[Hl,0].std(), hv[~Hl,0].mean(), hv[~Hl,0].std()))
print('  last-8-lines/rest pitch:    frame-bottom pages mean %.3f (sd %.3f) ; other pages mean %.3f (sd %.3f)' % (hv[Hl,1].mean(), hv[Hl,1].std(), hv[~Hl,1].mean(), hv[~Hl,1].std()))
# does a page's pitch predict what is needed to end at the foot?  capacity in lines at the page's own pitch
print('  quire T: px extent first-to-last line: mean %.0f sd %.0f (CV %.1f%%) ; lines CV %.1f%% ; if pitch were common, extent CV would equal lines CV' % ((n-1).dot(mp)/len(n)*0+np.mean((n-1)*mp), np.std((n-1)*mp), 100*np.std((n-1)*mp)/np.mean((n-1)*mp), 100*n.std()/n.mean()))
