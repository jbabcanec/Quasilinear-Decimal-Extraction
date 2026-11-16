"""Like-for-like audit of claim (b).
For every definition of the page foot we compute the SAME statistic on the real
pages and on pages simulated by 'whole paragraphs are added until the next one
does not fit' (foot re-estimated from each simulated set, exactly as for the data)."""
import numpy as np, sys
from fill import R, G50, GS, GT20
import scalefree  # adds W, A, cw rulers to R

TEXTONLY = [pg for pg in G50 if R[pg]['meta']['I'] in ('S', 'T') and pg not in ('f66r',)]  # ZL comments: 'text only'; f66r has a drawing in the bottom corner

def leaf(pg):
    import re
    return re.match(r'(f\d+)', pg).group(1)

# ---------- foot rules: (pages, y1 array) -> (foot array, evaluate-mask) ----------
def rule_q(q, key=None):
    def f(pages, y1):
        foot = np.zeros(len(pages)); mask = np.ones(len(pages), bool)
        keys = np.array([key(pg) if key else 'all' for pg in pages])
        for k in set(keys):
            m = keys == k
            foot[m] = np.percentile(y1[m], q)
            if m.sum() < 4: mask[m] = False
        return foot, mask
    return f

def rule_max(key, drop_ref=True):
    def f(pages, y1):
        foot = np.zeros(len(pages)); mask = np.ones(len(pages), bool)
        keys = np.array([key(pg) for pg in pages])
        for k in set(keys):
            m = np.where(keys == k)[0]
            foot[m] = y1[m].max()
            if len(m) < 2: mask[m] = False
            elif drop_ref: mask[m[np.argmax(y1[m])]] = False     # the page that defines the foot is not evaluated
        return foot, mask
    return f

def frac_full(pages, y1, pitch, rule, tol):
    foot, mask = rule(pages, y1)
    g = (foot - y1) / pitch
    return (g[mask] < tol).mean(), g, mask

def simulate(pages, y0, pitch, cap_lines, plen_seq, rng, mode):
    """greedy packing in line units; returns simulated number of lines per page."""
    n = np.zeros(len(pages)); npar = len(plen_seq)
    if mode == 'shuffle': seq = plen_seq[rng.permutation(npar)]; k = 0
    elif mode == 'shift': seq = plen_seq; k = int(rng.integers(npar))
    for i in range(len(pages)):
        used = 0
        while True:
            L = plen_seq[rng.integers(npar)] if mode == 'draw' else seq[k % npar]
            if used + L <= cap_lines[i] + 1e-9:
                used += L
                if mode != 'draw': k += 1
            else:
                if used == 0:             # paragraph longer than the whole page: let it fill the page (rare)
                    used = np.floor(cap_lines[i]);
                    if mode != 'draw': k += 1
                break
        n[i] = used
    return n

def run(name, pages, rule, tol=0.5, ruler='px', nsim=4000, seed=7, modes=('shuffle', 'shift', 'draw'), pool=None, verbose=True, slack=0.0):
    rng = np.random.default_rng(seed)
    if ruler == 'px':
        y0 = np.array([R[pg]['y0'] for pg in pages]); y1 = np.array([R[pg]['y1'] for pg in pages]); pitch = np.array([R[pg]['mp'] for pg in pages])
    elif ruler == 'lines':          # ZL3b line counts only, no pixels
        y0 = np.ones(len(pages)); y1 = np.array([float(R[pg]['n']) for pg in pages]); pitch = np.ones(len(pages))
    elif ruler == 'W':              # heights in units of the page's own text width (scale-free)
        y0 = np.zeros(len(pages)); y1 = np.array([R[pg]['A'] for pg in pages]); pitch = np.array([R[pg]['mp'] / R[pg]['W'] for pg in pages])
    elif ruler == 'cw':             # heights in units of the page's mean character width (scale-free)
        y0 = np.zeros(len(pages)); y1 = np.array([R[pg]['Hc'] for pg in pages]); pitch = np.array([R[pg]['mp'] / R[pg]['cw'] for pg in pages])
    obs, g, mask = frac_full(pages, y1, pitch, rule, tol)
    foot, _ = rule(pages, y1)
    cap = (foot - y0) / pitch + 1 + slack          # line slots down to the estimated foot
    nobs = np.array([R[pg]['n'] for pg in pages])
    cap = np.maximum(cap, 1)
    plen = np.array([L for pg in (pool or pages) for L in R[pg]['plen']], float)
    out = dict(name=name, n_eval=int(mask.sum()), k=int((g[mask] < tol).sum()), obs=obs)
    s = f'{name:58s} eval {int(mask.sum()):2d} full {out["k"]:2d} = {obs:4.0%} |'
    for mode in modes:
        fr = np.zeros(nsim)
        for t in range(nsim):
            ns = simulate(pages, y0, pitch, cap, plen, rng, mode)
            y1s = y0 + (ns - 1) * pitch
            fr[t], _, _ = frac_full(pages, y1s, pitch, rule, tol)
        p = ((fr >= obs - 1e-12).sum() + 1) / (nsim + 1)
        out[mode] = (fr.mean(), fr.std(), p)
        s += f' {mode}: {fr.mean():4.0%}+-{fr.std():3.0%} p={p:.3f} |'
    if verbose: print(s)
    return out, g, mask

if __name__ == '__main__':
    quire = lambda pg: R[pg]['meta']['Q']
    bif = lambda pg: R[pg]['meta']['Q'] + R[pg]['meta']['B']
    sets = (('50 S+B+T', G50), ('25 S', GS), ('23 quire T', GT20), ('30 text-only', TEXTONLY))
    print('#### A. pixel frame taken at face value; foot rule varied; like-for-like null (foot re-estimated in every simulation)')
    for nm, G in sets:
        run(f'{nm}: P90 of group', G, rule_q(90))
    for nm, G in sets:
        run(f'{nm}: P90 of group, tol 1.0 line', G, rule_q(90), tol=1.0)
    for nm, G in sets[:2] + sets[3:]:
        run(f'{nm}: P90 per quire (quires >=4 pages)', G, rule_q(90, quire))
    for nm, G in sets:
        run(f'{nm}: max per quire, ref page dropped', G, rule_max(quire))
    for nm, G in sets:
        run(f'{nm}: max per bifolio, ref page dropped', G, rule_max(bif))
    for nm, G in sets:
        run(f'{nm}: max of the two pages of a leaf, ref dropped', G, rule_max(leaf))
    print('#### B. ZL3b line counts only (no pixels)')
    for nm, G in sets:
        run(f'{nm}: lines, P90 per quire', G, rule_q(90, quire), ruler='lines')
    for nm, G in sets:
        run(f'{nm}: lines, max per bifolio, ref dropped', G, rule_max(bif), ruler='lines')
    for nm, G in sets:
        run(f'{nm}: lines, max per leaf, ref dropped', G, rule_max(leaf), ruler='lines')
    print('#### C. scale-free rulers (text height in units of own text width W / own char width cw)')
    for nm, G in sets[1:]:
        run(f'{nm}: H/W, P90 of group', G, rule_q(90), ruler='W')
        run(f'{nm}: H/W, max per leaf, ref dropped', G, rule_max(leaf), ruler='W')
        run(f'{nm}: H/cw, P90 of group', G, rule_q(90), ruler='cw')
        run(f'{nm}: H/cw, max per leaf, ref dropped', G, rule_max(leaf), ruler='cw')
