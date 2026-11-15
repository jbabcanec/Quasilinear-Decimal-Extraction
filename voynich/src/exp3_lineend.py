"""Experiment 3 - what kind of thing is the last word of a line?
Writes results/exp3_lineend.json.

Follow-ups to the line-break test (exp2):
  filler      is the last word padding added after the next word failed to fit?
  measure     does the last word's length follow the room that was left for it?
  precision   does it fill that room better than a random word that fits would?
  boundary    does the glyph-to-glyph link between neighbouring words cross the line break?
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from layout import WidthModel
from vdata import ROOT, paragraph_text, tokens
from wrap import breaks, line_table, margins, simulate

ROOM_EDGES = [0.3, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9]


def prep(S, wm):
    L = margins(line_table(S, wm))
    L["room_last"] = (L.M1 - L.x0_last) / L.wmed       # room when the last word was begun
    return L, breaks(L)


def fe_slope(B, y, x):
    xx = B[x] - B.groupby("folio")[x].transform("mean")
    yy = B[y] - B.groupby("folio")[y].transform("mean")
    return float((xx * yy).sum() / (xx * xx).sum())


def margin_free_slope(S, wm):
    """No margin estimate at all: within a paragraph, does a line stop further left when
    the next line starts with a wider word?  Wrapping predicts about -0.5."""
    L = line_table(S, wm)
    L = L[L.has_next & (L.par_lines >= 4)].copy()
    L["e"] = L.E / L.wmed
    key = ["folio", "par"]
    L = L[L.groupby(key).e.transform("size") >= 3]
    y = L.e - L.groupby(key).e.transform("mean")
    x = L.w_next - L.groupby(key).w_next.transform("mean")
    ok = y.abs() < 1.5        # lines stopping far short of the paragraph edge (drawings) dropped
    return float((x[ok] * y[ok]).sum() / (x[ok] * x[ok]).sum())


def mi(pairs):
    c = Counter(pairs); n = len(pairs)
    a = Counter(p[0] for p in pairs); b = Counter(p[1] for p in pairs)
    return sum(v / n * np.log2(v * n / (a[p[0]] * b[p[1]])) for p, v in c.items())


def mi_corrected(pairs, rng, m=3000, reps=8, shuffles=30):
    """Plug-in mutual information minus its shuffle bias, on equal-size subsamples."""
    vals = []
    for _ in range(reps):
        sub = [pairs[i] for i in rng.permutation(len(pairs))[:m]]
        x = [p[0] for p in sub]; y = [p[1] for p in sub]
        null = np.mean([mi(list(zip(x, rng.permutation(y)))) for _ in range(shuffles)])
        vals.append(mi(sub) - null)
    return float(np.mean(vals)), float(np.std(vals))


def main():
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"]).reset_index(drop=True)
    wm = WidthModel(pt)
    rng = np.random.default_rng(0)
    L, B = prep(pt, wm)
    out = {}

    # ---- filler: room before the last word vs the next line's first word
    f = {"real": fe_slope(B, "room_last", "w_next")}
    for mode in ("pour", "fit"):
        v = []
        for seed in range(6):
            S = simulate(pt, L, wm, np.random.default_rng(500 + seed), mode=mode, tol=0.15)
            v.append(fe_slope(prep(S, wm)[1], "room_last", "w_next"))
        f["pourer" if mode == "pour" else "fitter"] = [float(np.mean(v)), float(np.min(v)), float(np.max(v))]
    out["filler_slope"] = f

    # ---- robustness: long lines only (short lines beside drawings removed)
    from wrap import wrap_stats
    rob = {"real": wrap_stats(B[B.n >= 6])["slope"]}
    v = []
    for seed in range(4):
        S = simulate(pt, L, wm, np.random.default_rng(700 + seed), mode="pour", tol=0.15)
        Bs = prep(S, wm)[1]
        v.append(wrap_stats(Bs[Bs.n >= 6])["slope"])
    rob["pourer"] = [float(np.mean(v)), float(np.min(v)), float(np.max(v))]
    out["coupling_lines_of_6_plus"] = rob

    # ---- robustness: no margin estimate
    mf = {"real": margin_free_slope(pt, wm)}
    for mode, name in (("pour", "pourer"), ("fit", "fitter")):
        v = [margin_free_slope(simulate(pt, L, wm, np.random.default_rng(900 + seed), mode=mode, tol=0.15), wm)
             for seed in range(4)]
        mf[name] = [float(np.mean(v)), float(np.min(v)), float(np.max(v))]
    out["margin_free_slope"] = mf

    # ---- measure: last word against the room it had
    B["rb"] = pd.cut(B.room_last, ROOM_EDGES)
    g = B.groupby("rb", observed=True)
    inner = pt[(pt.k > 0) & (pt.kr > 0)]
    out["by_room"] = dict(
        room=[float(v) for v in g.room_last.median()], n=[int(v) for v in g.size()],
        glyphs=[float(v) for v in g.word_last.apply(lambda s: s.str.len().mean())],
        slack=[float(v) for v in g.slack1.median()],
        end_m=[float(v) for v in g.word_last.apply(lambda s: s.str.endswith("m").mean())])
    out["end_m"] = dict(last=float(B.word_last.str.endswith("m").mean()),
                        interior=float(inner.word.str.endswith("m").mean()))

    # ---- precision: the real last word vs a random interior word that fits the same room
    pool = (inner.w / inner.wmed).to_numpy()
    cf = []
    for room in B.room_last:
        ok = pool[pool <= room + 0.15]
        if len(ok):
            cf.append(room - rng.choice(ok))
    q = lambda v: [float(np.median(v)), float(np.percentile(v, 75) - np.percentile(v, 25))]
    out["precision"] = dict(real=q(B.slack1), random_fitting_word=q(np.array(cf)))

    # ---- boundary: last glyph of a word vs first glyph of the next
    w, ln, fo, par, k, kr = (pt[c].to_numpy() for c in ("word", "line", "folio", "par", "k", "kr"))
    idx = np.arange(len(w))
    same = lambda a, b: (fo[a] == fo[b]) & (ln[a] == ln[b])
    nxt = lambda a, b: (fo[a] == fo[b]) & (ln[b] == ln[a] + 1) & (par[a] == par[b])
    i, j, i2, j2 = idx[:-1], idx[1:], idx[:-2], idx[2:]
    sets = {
        "neighbours inside a line": (i, j, same(i, j) & (k[i] > 0) & (kr[j] > 0)),
        "one word apart inside a line": (i2, j2, same(i2, j2) & (k[i2] > 0) & (kr[j2] > 0)),
        "second-last word -> last word": (i, j, same(i, j) & (kr[j] == 0)),
        "last word -> first word of next line": (i, j, nxt(i, j) & (kr[i] == 0) & (k[j] == 0)),
        "second-last word -> first word of next line": (i2, j2, nxt(i2, j2) & (kr[i2] == 1) & (k[j2] == 0)),
    }
    bd = {}
    for name, (a, b, m) in sets.items():
        pairs = [(w[x][-1], w[y][0]) for x, y in zip(a[m], b[m])]
        mu, sd = mi_corrected(pairs, rng)
        bd[name] = dict(n=len(pairs), bits=mu, sd=sd)
    out["boundary_mi"] = bd

    # ---- which words pile up at the line end
    cl, ci = Counter(B.word_last), Counter(inner.word)
    rows = sorted(((wd, c, c / len(B), ci.get(wd, 0) / len(inner)) for wd, c in cl.items() if c >= 12),
                  key=lambda r: -r[2] / max(r[3], 1e-5))[:12]
    out["line_final_words"] = [dict(word=r[0], n=r[1], share_last=r[2], share_inside=r[3]) for r in rows]

    (ROOT / "results" / "exp3_lineend.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "line_final_words"}, indent=1))


if __name__ == "__main__":
    main()
