"""Permutation machinery for the line-start audit (written from scratch)."""
import numpy as np


def encode(keys):
    """Map hashable keys to ints; None -> -1."""
    d = {}
    out = np.empty(len(keys), dtype=np.int64)
    for i, k in enumerate(keys):
        if k is None or (isinstance(k, float) and k != k):
            out[i] = -1
        else:
            out[i] = d.setdefault(k, len(d))
    return out, d


def pairs_within(unit, lag, mask=None):
    """Index pairs (i, i+lag) where both lines lie in the same unit (consecutive rows)."""
    n = len(unit)
    a = np.arange(0, n - lag)
    b = a + lag
    ok = unit[a] == unit[b]
    if mask is not None:
        ok &= mask[a] & mask[b]
    return a[ok], b[ok]


def expected(g, group, a, b):
    """Exact permutation expectation of #matches over the given pairs when glyphs are
    permuted among eligible (g>=0) lines within each group."""
    elig = g >= 0
    G = int(g.max()) + 1
    grp_ids, gi = np.unique(group, return_inverse=True)
    cnt = np.zeros((len(grp_ids), G))
    np.add.at(cnt, (gi[elig], g[elig]), 1)
    n = cnt.sum(1)
    with np.errstate(divide='ignore', invalid='ignore'):
        same = (cnt * (cnt - 1)).sum(1) / (n * (n - 1))
    same = np.nan_to_num(same)
    ga, gb = gi[a], gi[b]
    e = np.where(ga == gb, same[ga], 0.0)
    diff = ga != gb
    if diff.any():
        e[diff] = (cnt[ga[diff]] * cnt[gb[diff]]).sum(1) / (n[ga[diff]] * n[gb[diff]])
    return e.sum()


def permtest(g, group, a, b, nperm=4000, seed=0, return_null=False):
    """g: int glyph ids (-1 = ineligible, never permuted; pairs touching them must be absent).
    group: permutation stratum for each line. a, b: pair index arrays."""
    rng = np.random.default_rng(seed)
    assert (g[a] >= 0).all() and (g[b] >= 0).all()
    obs = int((g[a] == g[b]).sum())
    N = len(a)
    idx = np.nonzero(g >= 0)[0]
    grp = group[idx]
    # sort idx by group so a lexsort on (rand, grp) permutes within groups
    o0 = np.argsort(grp, kind='stable')
    idx = idx[o0]
    grp = grp[o0]
    gv = g[idx]
    null = np.empty(nperm)
    gp = g.copy()
    for k in range(nperm):
        r = rng.random(len(idx))
        o = np.lexsort((r, grp))
        gp[idx] = gv[o]
        null[k] = (gp[a] == gp[b]).sum()
    exp_exact = expected(g, group, a, b)
    mu, sd = null.mean(), null.std(ddof=1)
    z = (obs - exp_exact) / sd if sd > 0 else float('nan')
    p_lo = (1 + (null <= obs).sum()) / (nperm + 1)
    p_hi = (1 + (null >= obs).sum()) / (nperm + 1)
    res = dict(N=N, obs=obs, rate=obs / N if N else float('nan'), exp=exp_exact,
               exp_rate=exp_exact / N if N else float('nan'), mc_mean=mu, sd=sd, z=z,
               p_lo=p_lo, p_hi=p_hi)
    if return_null:
        res['null'] = null
    return res


def fmt(r):
    return ('N=%5d obs=%5.1f%% exp=%5.1f%% ratio=%.2f z=%+6.2f p_lo=%.4f p_hi=%.4f' %
            (r['N'], 100 * r['rate'], 100 * r['exp_rate'],
             r['rate'] / r['exp_rate'] if r['exp_rate'] else float('nan'), r['z'], r['p_lo'], r['p_hi']))
