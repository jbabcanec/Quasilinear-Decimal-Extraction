"""Pair table: every (target word, earlier word on the same page within D lines above).

Works on any token table that has the columns
    folio, line (reading-order line id within folio), k, kr, nline, xc, x0, x1, yc, w, wmed, word
so the real manuscript and every synthetic/control layout go through identical code.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist


def type_distance(words):
    """Levenshtein distance matrix between word types, capped at 3 (uint8)."""
    types = sorted(set(words))
    idx = {t: i for i, t in enumerate(types)}
    D = cdist(types, types, scorer=Levenshtein.distance, score_cutoff=2, dtype=np.uint8)
    # rapidfuzz writes score_cutoff+1 for "beyond cutoff"
    return types, idx, D


def build_pairs(df: pd.DataFrame, dmax: int = 12) -> pd.DataFrame:
    """One row per (a, b) with b on a line 1..dmax above a, same folio."""
    chunks = []
    for folio, g in df.groupby("folio", sort=False):
        g = g.sort_values(["line", "k"])
        words = g["word"].to_numpy()
        types, idx, D = type_distance(words)
        t = np.array([idx[w] for w in words])
        n = len(g)
        if n < 4:
            continue
        # per-target base rate: P(similar to a random other paragraph token on this page)
        tc = np.bincount(t, minlength=len(types)).astype(float)
        base = {}
        for name, thr in (("ex", 0), ("e1", 1), ("e2", 2)):
            M = (D <= thr).astype(float)
            cnt = M @ tc                     # per type: number of tokens within thr
            base[name] = (cnt[t] - 1.0) / (n - 1.0)   # exclude the token itself

        line = g["line"].to_numpy()
        ulines, linepos = np.unique(line, return_inverse=True)   # rank of line on page
        first_tok = np.searchsorted(linepos, np.arange(len(ulines)))
        last_tok = np.append(first_tok[1:], n)
        yline = np.array([g["yc"].to_numpy()[first_tok[i]:last_tok[i]].mean()
                          for i in range(len(ulines))])
        pitch = np.median(np.diff(yline)) if len(yline) > 2 else np.nan

        xc = g["xc"].to_numpy(); x0 = g["x0"].to_numpy(); x1 = g["x1"].to_numpy()
        k = g["k"].to_numpy(); kr = g["kr"].to_numpy(); nl = g["nline"].to_numpy()
        wmed = float(g["wmed"].iloc[0])
        par = g["par"].to_numpy() if "par" in g else np.zeros(n, int)
        gi = g.index.to_numpy()
        seq = np.arange(n)

        A, B, Dd = [], [], []
        for la in range(len(ulines)):
            a_idx = np.arange(first_tok[la], last_tok[la])
            for d in range(1, dmax + 1):
                lb = la - d
                if lb < 0:
                    break
                b_idx = np.arange(first_tok[lb], last_tok[lb])
                aa, bb = np.meshgrid(a_idx, b_idx, indexing="ij")
                A.append(aa.ravel()); B.append(bb.ravel())
                Dd.append(np.full(aa.size, d))
        if not A:
            continue
        A = np.concatenate(A); B = np.concatenate(B); Dd = np.concatenate(Dd)
        dist = D[t[A], t[B]]
        chunks.append(pd.DataFrame(dict(
            folio=folio, a=gi[A], b=gi[B], d=Dd,
            dx=(xc[A] - xc[B]) / wmed,
            # horizontal overlap of the two boxes as a fraction of the narrower one
            ov=np.clip(np.minimum(x1[A], x1[B]) - np.maximum(x0[A], x0[B]), 0, None)
               / np.minimum(x1[A] - x0[A], x1[B] - x0[B]),
            dk=k[A] - k[B], dkr=kr[A] - kr[B], lag=seq[A] - seq[B],
            dy=(yline[linepos[A]] - yline[linepos[B]]) / pitch,
            ka=k[A], kra=kr[A], kb=k[B], krb=kr[B], na=nl[A], nb=nl[B],
            same_par=par[A] == par[B],
            dist=dist,
            p_ex=base["ex"][A], p_e1=base["e1"][A], p_e2=base["e2"][A],
            # source-side propensity and page mean, for a pair-level expectation
            q_e1=base["e1"][B], pbar_e1=float(base["e1"].mean()),
            la=len(words[0]) * 0 + np.char.str_len(words[A].astype(str)),
        )))
    P = pd.concat(chunks, ignore_index=True)
    P["ex"] = P["dist"] == 0
    P["e1"] = P["dist"] <= 1
    P["near"] = P["dist"] == 1          # similar but not identical
    P["e2"] = P["dist"] <= 2
    P["interior"] = (P.ka > 0) & (P.kra > 0) & (P.kb > 0) & (P.krb > 0)
    return P


def lift(P: pd.DataFrame, sim: str = "e1", base: str | None = None) -> float:
    """Observed / expected similarity, expectation = per-target page base rate."""
    base = base or {"ex": "p_ex", "e1": "p_e1", "e2": "p_e2", "near": None}[sim]
    if base is None:
        exp = (P["p_e1"] - P["p_ex"]).sum()
    else:
        exp = P[base].sum()
    return float(P[sim].sum() / exp) if exp > 0 else np.nan


def rulers(P1: pd.DataFrame, rng=None) -> pd.DataFrame:
    """For each target a (pairs with the line directly above only), pick the source that
    each 'ruler' nominates:
        phys : the word whose centre is horizontally closest (the word physically above)
        ord  : the word with the same index counted from the line start
        ordr : the word with the same index counted from the line end
        seq  : the last word of the line above (closest in reading order)
        rand : a random word of the line above
    Returns one row per target with the similarity outcome under each ruler.
    """
    rng = rng or np.random.default_rng(0)
    P1 = P1.copy()
    P1["adx"] = P1["dx"].abs()
    P1["u"] = rng.random(len(P1))
    grp = P1.groupby("a", sort=False)
    out = {}
    pick = {
        "phys": P1.loc[grp["adx"].idxmin()],
        "seq": P1.loc[grp["lag"].idxmin()],
        "rand": P1.loc[grp["u"].idxmin()],
        "ord": P1[P1.dk == 0],
        "ordr": P1[P1.dkr == 0],
    }
    res = None
    for name, sub in pick.items():
        s = sub.set_index("a")[["b", "dist", "dx", "dk", "dkr", "lag"]].add_prefix(name + "_")
        res = s if res is None else res.join(s, how="outer")
    meta = grp[["folio", "ka", "kra", "na", "nb", "p_ex", "p_e1", "p_e2", "same_par", "dy"]].first()
    return res.join(meta)


def ruler_lifts(P1: pd.DataFrame, interior: bool = True, rng=None) -> dict:
    """Observed / expected similarity of the word each ruler nominates.

    A ruler that tends to nominate long or rare words would look bad for reasons that have
    nothing to do with position, so the expectation is pair-level: e_ab = p_a * p_b / pbar,
    with p the page-level propensity of each word to resemble another word of its page.
    Also returns the paired comparison phys-vs-ord on targets where the two rulers nominate
    different words (counts of discordant outcomes).
    """
    rng = rng or np.random.default_rng(0)
    P1 = P1.copy()
    P1["adx"] = P1["dx"].abs()
    P1["u"] = rng.random(len(P1))
    P1["e_ab"] = P1.p_e1 * P1.q_e1 / P1.pbar_e1
    if interior:
        P1 = P1[(P1.ka > 0) & (P1.kra > 0)]
    grp = P1.groupby("a", sort=False)
    pick = {"phys": P1.loc[grp["adx"].idxmin()], "ord": P1[P1.dk == 0], "ordr": P1[P1.dkr == 0],
            "seq": P1.loc[grp["lag"].idxmin()], "rand": P1.loc[grp["u"].idxmin()]}
    out = {k: dict(n=len(v), rate=float(v.e1.mean()), lift=float(v.e1.sum() / v.e_ab.sum()))
           for k, v in pick.items()}
    a = pick["phys"].set_index("a")[["b", "e1", "e_ab"]]
    b = pick["ord"].set_index("a")[["b", "e1", "e_ab"]]
    j = a.join(b, lsuffix="_p", rsuffix="_o", how="inner")
    j = j[j.b_p != j.b_o]
    out["paired"] = dict(n=len(j), phys_only=int((j.e1_p & ~j.e1_o).sum()),
                         ord_only=int((~j.e1_p & j.e1_o).sum()),
                         phys_lift=float(j.e1_p.sum() / j.e_ab_p.sum()),
                         ord_lift=float(j.e1_o.sum() / j.e_ab_o.sum()))
    return out
