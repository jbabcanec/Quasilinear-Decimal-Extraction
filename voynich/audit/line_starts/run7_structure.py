"""Further structure tests: off-diagonal transition structure beyond same-glyph avoidance,
even/odd lag alternation statistic, bootstrap CI, previous-line-ending channel."""
import numpy as np
import pandas as pd
from build import build
from stats import encode, pairs_within, permtest, fmt

rng = np.random.default_rng(7)
df, pages, order, recs, paras, orphans = build()
body = df[df.li >= 2].reset_index(drop=True)
g, dic = encode(list(body.g_bench))
inv = {v: k for k, v in dic.items()}
gu = g.copy()
unk = np.nonzero(gu < 0)[0]
gu[unk] = 1000 + np.arange(len(unk))
idxs = list(body.groupby('pid').indices.values())
seqs = [gu[i] for i in idxs]
G = len(dic)
KEEP = 300


def table(pairs_a, pairs_b):
    t = np.zeros((G, G))
    ok = (pairs_a < G) & (pairs_b < G)
    np.add.at(t, (pairs_a[ok], pairs_b[ok]), 1)
    return t


# conditional permutations (same lag-1 match count per paragraph) and unconditional ones
cond = []   # per paragraph: array KEEP x m
unc = []
for s in seqs:
    m = len(s)
    if m < 2:
        cond.append(np.tile(s, (KEEP, 1)))
        unc.append(np.tile(s, (KEEP, 1)))
        continue
    o1 = int((s[1:] == s[:-1]).sum())
    got = []
    n = 0
    first = None
    while n < KEEP:
        P = s[rng.random((4000, m)).argsort(1)]
        if first is None:
            first = P[:KEEP]
        c1 = (P[:, 1:] == P[:, :-1]).sum(1)
        sel = P[c1 == o1]
        got.append(sel)
        n += len(sel)
    cond.append(np.concatenate(got)[:KEEP])
    unc.append(first)

obs_t = sum(table(s[:-1], s[1:]) for s in seqs if len(s) > 1)


def perm_tables(perms):
    out = np.zeros((KEEP, G, G))
    for P in perms:
        if P.shape[1] < 2:
            continue
        a = P[:, :-1]
        b = P[:, 1:]
        ok = (a < G) & (b < G)
        k = np.repeat(np.arange(KEEP)[:, None], a.shape[1], 1)
        np.add.at(out, (k[ok], a[ok], b[ok]), 1)
    return out


Tc = perm_tables(cond)
Tu = perm_tables(unc)
Ec = Tc.mean(0)
Eu = Tu.mean(0)
tot = np.bincount(g[g >= 0], minlength=G)
top = [i for i in np.argsort(-tot)[:8]]
print('#### cell counts: observed / expected (unconditional within-paragraph shuffle) [z from shuffle sd]')
sdu = Tu.std(0, ddof=1)
print('         ' + ' '.join('%16s' % inv[j] for j in top))
for i in top:
    print('  %-4s ' % inv[i] + ' '.join('%4d/%5.1f[%+4.1f]' % (obs_t[i, j], Eu[i, j], (obs_t[i, j] - Eu[i, j]) / sdu[i, j] if sdu[i, j] > 0 else 0) for j in top))

# off-diagonal chi-square among cells with expectation >= 5, conditional on the lag-1 match counts
off = ~np.eye(G, dtype=bool) & (Ec >= 5)


def chi(t, E, mask):
    return (((t - E) ** 2 / np.where(E > 0, E, 1))[mask]).sum()


c_obs = chi(obs_t, Ec, off)
c_null = np.array([chi(Tc[k], Ec, off) for k in range(KEEP)])
print('\n#### off-diagonal structure beyond same-sign avoidance (conditional on per-paragraph lag-1 match count)')
print('  cells used: %d ; chi2 obs = %.1f ; null mean %.1f sd %.1f max %.1f ; z=%+.1f ; p=%.4f' % (
    off.sum(), c_obs, c_null.mean(), c_null.std(ddof=1), c_null.max(), (c_obs - c_null.mean()) / c_null.std(ddof=1),
    (1 + (c_null >= c_obs).sum()) / (KEEP + 1)))
# asymmetry: is the off-diagonal table asymmetric (directional, i.e. order matters top->down)?
asym_obs = (np.abs(obs_t - obs_t.T)[off & off.T]).sum() / 2
asym_null = np.array([(np.abs(Tc[k] - Tc[k].T)[off & off.T]).sum() / 2 for k in range(KEEP)])
print('  asymmetry sum|n_ij - n_ji|/2: obs %.0f ; null mean %.1f sd %.1f ; z=%+.1f' % (
    asym_obs, asym_null.mean(), asym_null.std(ddof=1), (asym_obs - asym_null.mean()) / asym_null.std(ddof=1)))
# symmetric part
sym_obs_t = obs_t + obs_t.T
Es = Ec + Ec.T
iu = np.triu(off & off.T, 1)
cs_obs = chi(sym_obs_t, Es, iu)
cs_null = np.array([chi(Tc[k] + Tc[k].T, Es, iu) for k in range(KEEP)])
print('  symmetrised off-diagonal chi2: obs %.1f ; null mean %.1f sd %.1f ; z=%+.1f' % (
    cs_obs, cs_null.mean(), cs_null.std(ddof=1), (cs_obs - cs_null.mean()) / cs_null.std(ddof=1)))
# most deviant off-diagonal pairs (symmetrised)
dev = []
for i in range(G):
    for j in range(i + 1, G):
        if iu[i, j]:
            sd = (Tc[:, i, j] + Tc[:, j, i]).std(ddof=1)
            dev.append(((sym_obs_t[i, j] - Es[i, j]) / sd, inv[i], inv[j], sym_obs_t[i, j], Es[i, j]))
dev.sort()
print('  most depleted unordered pairs:', ['%s-%s %d/%.1f z%+.1f' % (d[1], d[2], d[3], d[4], d[0]) for d in dev[:6]])
print('  most enriched unordered pairs:', ['%s-%s %d/%.1f z%+.1f' % (d[1], d[2], d[3], d[4], d[0]) for d in dev[-6:]])

print('\n#### even-vs-odd lags (2,4,6 minus 3,5) under the conditional null: residual alternation?')


def lagm(P, k):
    if P.shape[1] <= k:
        return np.zeros(P.shape[0])
    return (P[:, k:] == P[:, :-k]).sum(1)


T_obs = 0.0
mu = 0.0
var = 0.0
for s, P in zip(seqs, cond):
    S = s[None, :]
    f = lambda X: lagm(X, 2) + lagm(X, 4) + lagm(X, 6) - lagm(X, 3) - lagm(X, 5)
    T_obs += f(S)[0]
    x = f(P)
    mu += x.mean()
    var += x.var(ddof=1)
print('  T_obs=%.0f null mean=%.1f sd=%.1f z=%+.2f' % (T_obs, mu, np.sqrt(var), (T_obs - mu) / np.sqrt(var)))
# ABA triples: among triples with a!=b and b!=c, how often a==c ; conditional null
T_obs = 0
n_obs = 0
mu = 0.0
var = 0.0
for s, P in zip(seqs, cond):
    if len(s) < 3:
        continue
    f = lambda X: ((X[:, 2:] == X[:, :-2]) & (X[:, 1:-1] != X[:, :-2]) & (X[:, 2:] != X[:, 1:-1])).sum(1)
    T_obs += f(s[None, :])[0]
    x = f(P)
    mu += x.mean()
    var += x.var(ddof=1)
print('  ABA triples: obs=%d cond.null mean=%.1f sd=%.1f z=%+.2f' % (T_obs, mu, np.sqrt(var), (T_obs - mu) / np.sqrt(var)))

print('\n#### page-bootstrap CI for obs/exp ratio (lag 1)')
g1, _ = encode(list(body.g_eva))
_, grp = np.unique(body.pid.values, return_inverse=True)
a, b = pairs_within(body.pid.values, 1, mask=(g1 >= 0))
# per-paragraph obs and exp
from stats import expected
pid_of_pair = body.pid.values[a]
fol = body.folio.values
obs_p = pd.Series((g1[a] == g1[b]).astype(float)).groupby(fol[a]).sum()
exp_p = {}
for f, idx in body.groupby('folio').indices.items():
    sub = body.iloc[idx]
    gg, _ = encode(list(sub.g_eva))
    aa, bb = pairs_within(sub.pid.values, 1, mask=(gg >= 0))
    _, gr = np.unique(sub.pid.values, return_inverse=True)
    exp_p[f] = expected(gg, gr, aa, bb) if len(aa) else 0.0
exp_p = pd.Series(exp_p)
n_p = pd.Series(np.ones(len(a))).groupby(fol[a]).sum()
fl = list(exp_p.index)
o = obs_p.reindex(fl).fillna(0).values
e = exp_p.values
n = n_p.reindex(fl).fillna(0).values
bs = []
for _ in range(4000):
    k = rng.integers(0, len(fl), len(fl))
    bs.append((o[k].sum() / e[k].sum(), o[k].sum() / n[k].sum(), e[k].sum() / n[k].sum()))
bs = np.array(bs)
print('  ratio obs/exp = %.3f, 95%% CI [%.3f, %.3f] ; obs rate CI [%.1f%%, %.1f%%] ; exp rate CI [%.1f%%, %.1f%%]' % (
    o.sum() / e.sum(), *np.percentile(bs[:, 0], [2.5, 97.5]), *(100 * np.percentile(bs[:, 1], [2.5, 97.5])), *(100 * np.percentile(bs[:, 2], [2.5, 97.5]))))

print('\n#### channel through the upper line\'s ending')
lg, _ = encode(list(body.lastg))
ok = lg[a] >= 0
a2, b2 = a[ok], b[ok]
Gn = g1.max() + 1
tab = np.zeros((lg.max() + 1, Gn))
np.add.at(tab, (lg[a2], g1[b2]), 1)
rowsum = tab[lg[a2]].sum(1)
p_match = tab[lg[a2], g1[a2]] / rowsum
print('  pairs %d obs %.1f%% ; expected if lower first glyph depended only on the upper line\'s last glyph (pooled MS): %.1f%% ; pooled-MS independence: %.1f%%' % (
    len(a2), 100 * (g1[a2] == g1[b2]).mean(), 100 * p_match.mean(),
    100 * (np.bincount(g1[a2], minlength=Gn) / len(a2) * np.bincount(g1[b2], minlength=Gn) / len(b2)).sum()))

print('\n#### per-hand diagonal ratios for the main signs (EVA first char)')
inv1 = {v: k for k, v in encode(list(body.g_eva))[1].items()}
for h in ['1', '2', '3']:
    mk = (body.H.astype(str) == h).values
    line = []
    for ch in ['d', 'y', 'o', 'q', 's', 't', 'c']:
        keys = [(k if (k == ch) else ('other' if k is not None else None)) for k in body.g_eva]
        keys = [k if m else None for k, m in zip(keys, mk)]
        gg, dd = encode(keys)
        if ch not in dd:
            continue
        aa, bb = pairs_within(body.pid.values, 1, mask=(gg >= 0))
        # count only ch-ch matches
        code = dd[ch]
        obs = int(((gg[aa] == code) & (gg[bb] == code)).sum())
        # expected ch-ch
        ex = 0.0
        for idx in body.groupby('pid').indices.values():
            x = gg[idx]
            m = (x >= 0).sum()
            if m < 2:
                continue
            c = (x == code).sum()
            npairs = int(((x[1:] >= 0) & (x[:-1] >= 0)).sum())
            ex += npairs * c * (c - 1) / (m * (m - 1))
        line.append('%s %d/%.1f=%.2f' % (ch, obs, ex, obs / ex if ex else float('nan')))
    print('  hand %s: ' % h + ' | '.join(line))
