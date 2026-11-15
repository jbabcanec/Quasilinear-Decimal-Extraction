"""Synthetic layouts on the REAL page geometry.

Every control in this study is a token table with real line extents, real line pitch and a
word-width model fitted to the real boxes, so controls and manuscript pass through exactly
the same measurement code.

    pour(stream)            greedy line-wrap of a pre-existing word stream ("poured" text)
    keep_lines(lines)       a generator's own lines, one per physical line
    generate_onpage(...)    text composed on the page by copying earlier words, with a
                            chosen source rule:
                              'phys'    the word physically overhead (Gaussian in pixels)
                              'ord'     the word with the same count from the line start
                              'diffuse' any word of a nearby line above
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist


# ----------------------------------------------------------------------------- geometry
def page_templates(pt: pd.DataFrame):
    """folio -> dict(wmed, lines=[dict(xl, xr, y, par, par_last, n)])"""
    T = {}
    for folio, g in pt.groupby("folio", sort=False):
        lines = []
        for _, h in g.groupby("line", sort=True):
            lines.append(dict(xl=float(h.x0.min()), xr=float(h.x1.max()), y=float(h.yc.mean()),
                              par=int(h.par.iloc[0]), par_last=bool(h.par_last.iloc[0]),
                              n=len(h)))
        T[folio] = dict(wmed=float(g.wmed.iloc[0]), lines=lines)
    return T


class WidthModel:
    """Box width ~ sum of per-glyph widths (least squares on the real boxes) + noise."""

    def __init__(self, pt: pd.DataFrame):
        chars = sorted(set("".join(pt.word)))
        self.idx = {c: i for i, c in enumerate(chars)}
        X = np.zeros((len(pt), len(chars) + 1))
        for r, w in enumerate(pt.word):
            for c in w:
                X[r, self.idx[c]] += 1
        X[:, -1] = 1
        y = (pt.w / pt.wmed).to_numpy()          # in page-median word widths
        self.coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ self.coef
        self.sigma = float(res.std())
        self.r2 = float(1 - res.var() / y.var())
        self.mean_char = float(np.mean(self.coef[:-1]))
        # gaps between consecutive words of a line, drawing interruptions excluded
        q = pt.sort_values(["folio", "line", "k"])
        gap = (q.x0.shift(-1) - q.x1)[(q.line.shift(-1) == q.line) & (q.folio.shift(-1) == q.folio)]
        gap = (gap / q.wmed.loc[gap.index]).to_numpy()
        self.gaps = gap[(gap > -0.5) & (gap < 1.0)]

    def width(self, word: str, rng, known=True) -> float:
        if known:
            v = self.coef[-1] + sum(self.coef[self.idx[c]] for c in word if c in self.idx)
        else:   # foreign alphabet: every letter gets the mean glyph width
            v = self.coef[-1] + self.mean_char * len(word)
        return float(max(0.15, v + rng.normal(0, self.sigma)))

    def gap(self, rng) -> float:
        return float(rng.choice(self.gaps))


def _finish(rows, T) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=["folio", "line", "word", "x0", "x1", "yc", "par", "src"])
    df["w"] = df.x1 - df.x0
    df["xc"] = (df.x0 + df.x1) / 2
    df["k"] = df.groupby(["folio", "line"]).cumcount()
    df["nline"] = df.groupby(["folio", "line"])["k"].transform("size")
    df["kr"] = df.nline - 1 - df.k
    df["wmed"] = df.folio.map({f: t["wmed"] for f, t in T.items()})
    df["par_last"] = False
    return df


# ------------------------------------------------------------------------------ pouring
def pour(stream, T, wm: WidthModel, rng, known=True, widths=None, tol=0.15):
    """Greedy wrap: put the next word on the line if it fits before the line's real right
    end (+tol word widths), otherwise start the next line.  `stream` is an iterable of
    words; `widths` optionally supplies their widths (in word-width units)."""
    it = iter(stream)
    wit = iter(widths) if widths is not None else None
    rows = []
    pending = None
    try:
        for folio, t in T.items():
            s = t["wmed"]
            for li, L in enumerate(t["lines"]):
                x = L["xl"]
                placed = 0
                while True:
                    if pending is None:
                        w = next(it)
                        pending = (w, (next(wit) if wit is not None else wm.width(w, rng, known)) * s)
                    word, ww = pending
                    if placed and x + ww > L["xr"] + tol * s:
                        break
                    rows.append((folio, li, word, x, x + ww, L["y"], L["par"], -1))
                    x += ww + wm.gap(rng) * s
                    placed += 1
                    pending = None
                    if x > L["xr"]:
                        break
    except StopIteration:
        pass
    return _finish(rows, T)


def keep_lines(lines, T, wm: WidthModel, rng, known=True):
    """One generator line per physical line, left-aligned at the real left edge."""
    it = iter(lines)
    rows = []
    try:
        for folio, t in T.items():
            s = t["wmed"]
            for li, L in enumerate(t["lines"]):
                x = L["xl"]
                for word in next(it):
                    ww = wm.width(word, rng, known) * s
                    rows.append((folio, li, word, x, x + ww, L["y"], L["par"], -1))
                    x += ww + wm.gap(rng) * s
    except StopIteration:
        pass
    return _finish(rows, T)


# --------------------------------------------------------------------- on-page generators
class Vocab:
    """Real Voynich word types with frequencies and their edit-distance-1 neighbourhoods."""

    def __init__(self, words):
        c = Counter(words)
        self.types = np.array(sorted(c))
        self.freq = np.array([c[t] for t in self.types], float)
        self.p = self.freq / self.freq.sum()
        self.idx = {t: i for i, t in enumerate(self.types)}
        D = cdist(list(self.types), list(self.types), scorer=Levenshtein.distance,
                  score_cutoff=1, dtype=np.uint8)
        self.nbr = [np.flatnonzero(row == 1) for row in D]

    def fresh(self, rng):
        return int(rng.choice(len(self.types), p=self.p))

    def mutate(self, i, rng, p_exact):
        nb = self.nbr[i]
        if len(nb) == 0 or rng.random() < p_exact:
            return i
        w = self.freq[nb]
        return int(rng.choice(nb, p=w / w.sum()))


def generate_onpage(T, V: Vocab, wm: WidthModel, rng, rule="phys", lam=0.2, p_exact=0.3,
                    sigma=0.6, line_w=(0.6, 0.3, 0.1)):
    """Compose text on the page.  With probability `lam` a word is a copy (exact with
    prob. p_exact, else an edit-distance-1 variant) of an earlier word chosen by `rule`;
    otherwise it is a fresh draw from the manuscript's word frequencies."""
    rows = []
    mean_w = 1.0
    for folio, t in T.items():
        s = t["wmed"]
        page = []                                  # per line: list of (xc, type index)
        for li, L in enumerate(t["lines"]):
            x = L["xl"]
            cur = []
            while True:
                room = (L["xr"] - x) / s
                if cur and room < 0.55 * mean_w:
                    break
                for attempt in range(6):
                    src = -1
                    ti = None
                    if page and rng.random() < lam:
                        nl = min(len(line_w), len(page))
                        pw = np.array(line_w[:nl]); pw = pw / pw.sum()
                        srcline = page[-1 - int(rng.choice(nl, p=pw))]
                        if srcline:
                            xs = np.array([c[0] for c in srcline])
                            if rule == "phys":
                                pen = x + 0.5 * mean_w * s
                                wgt = np.exp(-0.5 * ((xs - pen) / (sigma * s)) ** 2) + 1e-12
                                j = int(rng.choice(len(xs), p=wgt / wgt.sum()))
                            elif rule == "ord":
                                j = len(cur) + int(rng.choice([-1, 0, 1], p=[0.2, 0.6, 0.2]))
                                j = min(max(j, 0), len(xs) - 1)
                            else:
                                j = int(rng.integers(len(xs)))
                            ti = V.mutate(srcline[j][1], rng, p_exact)
                            src = 1
                    if ti is None:
                        ti = V.fresh(rng)
                    word = str(V.types[ti])
                    ww = wm.width(word, rng)
                    if not cur or ww <= room + 0.15 or attempt == 5:
                        break
                ww *= s
                rows.append((folio, li, word, x, x + ww, L["y"], L["par"], src))
                cur.append((x + ww / 2, ti))
                x += ww + wm.gap(rng) * s
                if x > L["xr"]:
                    break
            page.append(cur)
    return _finish(rows, T)
