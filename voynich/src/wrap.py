"""The line-break test: was the text poured into the lines, or fitted to them?

A scribe writing out a text that already exists (plain or enciphered, copied or dictated)
breaks the line when the next word will not fit.  That leaves two physical fingerprints:

  1. slack     - the blank left at the line end is anywhere between zero and one word, and
  2. coupling  - it is smaller than the word that starts the next line (that word is the
                 one that did not fit), so lines followed by long words end earlier:
                 slope of slack on next-word width is about +0.5.

A scribe who chooses what to write to fill the line leaves neither: the blank is small and
knows nothing about the next line's first word.

All quantities are in units of the page's median word width.  `slack` needs the right
boundary available to the scribe; it is estimated per line from the neighbouring lines of
the same paragraph (see `margins`), and every simulated control is measured by the same
estimator, so estimator bias cancels in the comparison.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

GAP = 0.125            # median inter-word gap, in word widths (measured; see layout.WidthModel)


def line_table(df: pd.DataFrame, wm=None) -> pd.DataFrame:
    """One row per line.  If a WidthModel is given, also carries the glyph-predicted width
    of the first / last word (what the word 'should' occupy, independent of how it was
    actually penned)."""
    df = df.sort_values(["folio", "line", "k"])
    g = df.groupby(["folio", "line"], sort=False)
    L = g.agg(xl=("x0", "min"), E=("x1", "max"), y=("yc", "mean"), n=("k", "size"),
              par=("par", "first"), wmed=("wmed", "first")).reset_index()
    f = df[df.k == 0].set_index(["folio", "line"])
    l = df[df.kr == 0].set_index(["folio", "line"])
    L = L.join(f[["w", "word"]].rename(columns={"w": "w_first", "word": "word_first"}),
               on=["folio", "line"])
    L = L.join(l[["w", "word", "x0"]].rename(columns={"w": "w_last", "word": "word_last",
                                                       "x0": "x0_last"}), on=["folio", "line"])
    for c in ("w_first", "w_last"):
        L[c] = L[c] / L.wmed
    if wm is not None:
        pred = lambda w: wm.coef[-1] + sum(wm.coef[wm.idx[ch]] for ch in w if ch in wm.idx)
        L["p_first"] = L.word_first.map(pred)
        L["p_last"] = L.word_last.map(pred)
    grp = L.groupby("folio", sort=False)
    # the next row must be the same paragraph AND the physically next line
    same = (grp.par.shift(-1) == L.par) & (grp.line.shift(-1) == L.line + 1)
    L["has_next"] = same
    L["w_next"] = grp.w_first.shift(-1).where(same)
    if wm is not None:
        L["p_next"] = grp.p_first.shift(-1).where(same)
    L["par_lines"] = L.groupby(["folio", "par"]).line.transform("size")
    return L


def margins(L: pd.DataFrame, half: int = 3) -> pd.DataFrame:
    """Right boundary available to each line, from its neighbours in the same paragraph.
        M1 : furthest right end among the +-`half` neighbouring lines (own line excluded)
        M2 : second furthest (robust to one overrunning neighbour)
    Paragraph-final lines are not used as neighbours (they end early by design)."""
    L = L.copy()
    M1 = np.full(len(L), np.nan); M2 = np.full(len(L), np.nan)
    for _, g in L.groupby(["folio", "par"], sort=False):
        E = g.E.to_numpy().astype(float)
        usable = np.ones(len(E), bool); usable[-1] = len(E) == 1
        pos = g.index.to_numpy()
        for i in range(len(E)):
            idx = [j for j in range(max(0, i - half), min(len(E), i + half + 1))
                   if j != i and usable[j]]
            if not idx:
                continue
            v = np.sort(E[idx])[::-1]
            M1[L.index.get_loc(pos[i])] = v[0]
            M2[L.index.get_loc(pos[i])] = v[min(1, len(v) - 1)]
    L["M1"], L["M2"] = M1, M2
    L["slack1"] = (L.M1 - L.E) / L.wmed
    L["slack2"] = (L.M2 - L.E) / L.wmed
    return L


def breaks(L: pd.DataFrame, min_par_lines: int = 3) -> pd.DataFrame:
    """Line breaks inside paragraphs with a usable margin estimate."""
    return L[L.has_next & (L.par_lines >= min_par_lines) & L.M1.notna() & L.w_next.notna()].copy()


def _demean(B, col):
    return B[col] - B.groupby("folio")[col].transform("mean")


def wrap_stats(B: pd.DataFrame, slack: str = "slack1", wnext: str = "w_next") -> dict:
    """fit   : share of breaks where the next line's first word (+ a gap) would have fitted
       slope : regression slope of slack on next-word width, folio fixed effects
       med   : median slack;  iqr : interquartile range of slack"""
    x, y = _demean(B, wnext), _demean(B, slack)
    return dict(n=len(B), fit=float((B[wnext] + GAP <= B[slack]).mean()),
                slope=float((x * y).sum() / (x * x).sum()),
                med=float(B[slack].median()),
                iqr=float(B[slack].quantile(.75) - B[slack].quantile(.25)),
                p90=float(B[slack].quantile(.9)))


def boot_ci(B, fn, key, n=1000, seed=0):
    """Percentile bootstrap over folios for one statistic of fn(B)."""
    rng = np.random.default_rng(seed)
    groups = [g for _, g in B.groupby("folio", sort=False)]
    out = []
    for _ in range(n):
        pick = rng.integers(0, len(groups), len(groups))
        parts = []
        for j, i in enumerate(pick):
            gcopy = groups[i].copy(); gcopy["folio"] = j
            parts.append(gcopy)
        out.append(fn(pd.concat(parts))[key])
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def shuffle_next(B, rng, wnext="w_next"):
    """Null for 'the next line's first word is unrelated to this line's ending':
    permute next-word widths among the breaks of the same folio."""
    Bp = B.copy()
    Bp[wnext] = B.groupby("folio")[wnext].transform(lambda s: rng.permutation(s.to_numpy()))
    return Bp


# ---------------------------------------------------------------- simulated scribes
def simulate(pt: pd.DataFrame, Lreal: pd.DataFrame, wm, rng, mode="pour", tol=0.0,
             elastic=0.0, margin="M1", tries=12, pick="first"):
    """Re-lay the REAL word stream on the REAL pages under a stated line-ending rule.

    Each line's left edge and available right boundary are the real ones (boundary =
    max(real line end, neighbour-based margin)).  Words keep their real box widths.

      mode='pour'  the next word goes on this line if it fits (overrun up to `tol`
                   allowed), else it starts the next line.  With `elastic` > 0 the scribe
                   may also squeeze the last word by up to that fraction to make it fit.
      mode='fit'   as 'pour', but at the line end the scribe is free to pick what to write:
                   among the next `tries` words of the stream it takes the first one that
                   fits the remaining space (pick='first') or the one that fills it best
                   (pick='widest'); the others are put back in order.  The word stream is
                   the same multiset, only locally re-ordered by the page.
    Returns a token table (same columns as the manuscript's) for the standard measurements.
    """
    q = pt.sort_values(["folio", "line", "k"])
    rows = []
    for folio, g in q.groupby("folio", sort=False):
        words = g.word.tolist()
        widths = (g.w / g.wmed).tolist()
        s = float(g.wmed.iloc[0])
        Lf = Lreal[Lreal.folio == folio].sort_values("line")
        i = 0
        for r in Lf.itertuples():
            R = np.nanmax([r.E, getattr(r, margin)]) if not np.isnan(getattr(r, margin)) else r.E
            x = r.xl
            placed = 0
            while i < len(words):
                room = (R - x) / s
                w = widths[i]
                if mode == "fit" and placed and room < 2.2:
                    # choose, among upcoming words, the best fill of the remaining room
                    cand = range(i, min(len(words), i + tries))
                    ok = [j for j in cand if widths[j] <= room + tol]
                    if ok:
                        j = max(ok, key=lambda j: widths[j]) if pick == "widest" else ok[0]
                        words.insert(i, words.pop(j)); widths.insert(i, widths.pop(j))
                        w = widths[i]
                        # if another word could still follow, do not force line end yet
                    else:
                        break
                wdraw = w
                if placed and w > room + tol:
                    if elastic and w * (1 - elastic) <= room + tol:
                        wdraw = max(room, w * (1 - elastic))
                    else:
                        break
                rows.append((folio, r.line, words[i], x, x + wdraw * s, r.y, r.par))
                x += wdraw * s + wm.gap(rng) * s
                placed += 1
                i += 1
                if x > R:
                    break
    S = pd.DataFrame(rows, columns=["folio", "line", "word", "x0", "x1", "yc", "par"])
    S["w"] = S.x1 - S.x0
    S["xc"] = (S.x0 + S.x1) / 2
    S["k"] = S.groupby(["folio", "line"]).cumcount()
    S["nline"] = S.groupby(["folio", "line"])["k"].transform("size")
    S["kr"] = S.nline - 1 - S.k
    S["wmed"] = S.folio.map(q.groupby("folio").wmed.first())
    return S
