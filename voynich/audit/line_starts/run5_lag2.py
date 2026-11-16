"""Is the lag-2 excess an independent 'alternation' signal, or just the arithmetic consequence of
lag-1 avoidance under a within-paragraph permutation null (the paragraph's multiset is fixed, so
matches removed from lag 1 must reappear at other lags)?

Conditional null: permutations of each paragraph's (non-first) lines that have exactly the observed
number of lag-1 matches in that paragraph.  Compare observed lag-k (k>=2) matches against it."""
import sys
import numpy as np
from build import build
from stats import encode

rng = np.random.default_rng(12345)
KEEP = 400


def cond_null(seqs, lags=(2, 3, 4, 5, 6), keep=KEEP, cond='count'):
    obs = {k: 0 for k in lags}
    N = {k: 0 for k in lags}
    mu = {k: 0.0 for k in lags}
    var = {k: 0.0 for k in lags}
    mu_un = {k: 0.0 for k in lags}
    nfail = 0
    for s in seqs:
        s = np.asarray(s)
        m = len(s)
        if m < 3:
            continue
        o1 = int((s[1:] == s[:-1]).sum())
        got = []
        tries = 0
        while sum(len(x) for x in got) < keep and tries < 400:
            tries += 1
            R = rng.random((4000, m)).argsort(1)
            P = s[R]
            c1 = (P[:, 1:] == P[:, :-1]).sum(1)
            got.append(P[c1 == o1])
        P = np.concatenate(got)[:keep]
        if len(P) < 30:
            nfail += 1
            continue
        # unconditional expectation
        _, cnt = np.unique(s, return_counts=True)
        psame = (cnt * (cnt - 1)).sum() / (m * (m - 1))
        for k in lags:
            if m <= k:
                continue
            x = (P[:, k:] == P[:, :-k]).sum(1)
            obs[k] += int((s[k:] == s[:-k]).sum())
            N[k] += m - k
            mu[k] += x.mean()
            var[k] += x.var(ddof=1)
            mu_un[k] += (m - k) * psame
    return obs, N, mu, var, mu_un, nfail


def report(seqs, label):
    obs, N, mu, var, mu_un, nfail = cond_null(seqs)
    print(label, '(paragraphs skipped for sampling failure: %d)' % nfail)
    for k in obs:
        if N[k] == 0:
            continue
        print('  lag%d N=%5d obs=%5.1f%%  uncond.exp=%5.1f%% (ratio %.2f)  cond.exp|lag1=%5.1f%% (ratio %.2f, z=%+.2f)' % (
            k, N[k], 100 * obs[k] / N[k], 100 * mu_un[k] / N[k], obs[k] / mu_un[k],
            100 * mu[k] / N[k], obs[k] / mu[k], (obs[k] - mu[k]) / np.sqrt(var[k])))


if __name__ == '__main__':
    df, *_ = build()
    body = df[df.li >= 2]
    for col in ['g_eva', 'g_bench']:
        g, _ = encode(list(body[col]))
        g = g.copy()
        # unknown first glyphs: give each a unique id so they never match
        unk = np.nonzero(g < 0)[0]
        g[unk] = 10000 + np.arange(len(unk))
        seqs = [g[idx] for idx in body.reset_index(drop=True).groupby('pid').indices.values()]
        report(seqs, 'ZL %s all paragraphs' % col)
    seqs_long = [s for s in seqs if len(s) >= 8]
    report(seqs_long, 'ZL g_bench paragraphs with >=8 non-first lines')
    from v101parse import v101_table
    d = v101_table()
    b = d[d.li >= 2].reset_index(drop=True)
    g, _ = encode(list(b.g))
    unk = np.nonzero(g < 0)[0]
    g[unk] = 10000 + np.arange(len(unk))
    seqs = [g[idx] for idx in b.groupby('pid').indices.values()]
    report(seqs, 'v101 all paragraphs')
