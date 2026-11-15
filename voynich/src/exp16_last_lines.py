"""Experiment 16 - paragraphs do not end on a short line.
Writes results/exp16_last_lines.json.

1. How long is the last line of a paragraph, as a share of that paragraph's normal line?
   A text that simply runs out leaves every length about equally often, so roughly one
   last line in nine holds a single word.
2. Self-control: re-wrap each paragraph's own words (real box widths) into lines a little
   narrower or wider than the real ones.  If the amount of text was not tuned to the
   layout, the last-line lengths should look the same at every width.
3. Real manuscripts: lines that stop early in six French manuscripts.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "control_real"))
from vdata import ROOT, load_zl, paragraph_text, tokens
from wrap import line_table, margins

GAP = 0.125
EDGES = [0, .2, .4, .6, .8, 1.0, 9]


def rewrap(widths, W):
    """Greedy wrap; returns (fill of last line, words in last line, number of lines)."""
    x, n, lines = 0.0, 0, 1
    for w in widths:
        if n and x + w > W + 0.15:
            lines += 1; x, n = 0.0, 0
        x += w + GAP; n += 1
    return (x - GAP) / W, n, lines


def main():
    out = {}
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    L = line_table(pt)
    sec = pt.groupby("folio").illust.first()
    paras = []
    for (f, p), g in L.groupby(["folio", "par"], sort=False):
        g = g.sort_values("line")
        if len(g) < 4:
            continue
        body = g.iloc[1:-1]
        W = float(np.median((body.E - body.xl) / body.wmed))
        tok = pt[(pt.folio == f) & (pt.par == p)]
        last = g.iloc[-1]
        paras.append(dict(folio=f, W=W, widths=(tok.w / tok.wmed).tolist(),
                          fill=float((last.E - last.xl) / last.wmed / W), nwords=int(last.n),
                          body_words=float(np.median(body.n))))
    P = pd.DataFrame(paras)

    def summary(fill, nwords):
        fill = np.asarray(fill); nwords = np.asarray(nwords)
        h, _ = np.histogram(fill, EDGES)
        return dict(n=len(fill), one_word=float((nwords == 1).mean()), up_to_two=float((nwords <= 2).mean()),
                    under_30=float((fill < 0.3).mean()), mid_40_80=float(((fill >= .4) & (fill < .8)).mean()),
                    fifths=(h / len(fill)).round(4).tolist(), mean=float(fill.mean()))
    out["voynich_as_written"] = summary(P.fill, P.nwords)
    for code, nm in (("S", "stars / recipes"), ("H", "plants"), ("B", "bathing")):
        s = P[P.folio.map(sec) == code]
        out[f"voynich_as_written, {nm}"] = summary(s.fill, s.nwords)

    # self-control: the same paragraphs re-wrapped at other line widths
    rw = {}
    for delta in (-0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3):
        r = [rewrap(p.widths, p.W * (1 + delta)) for p in P.itertuples()]
        rw[f"{delta:+.2f}"] = summary([x[0] for x in r], [x[1] for x in r])
    out["voynich_rewrapped"] = rw

    # ZL3b alone (independent transliteration; short 'title' lines counted as last lines)
    pages, loci = load_zl()
    zp, cur = [], []
    for Lc in loci:
        if not Lc["ltype"].startswith("P") or not Lc["words"]:
            if cur: zp.append(cur); cur = []
            continue
        if (Lc["pstart"] or Lc["pos"] in "@*") and cur:
            zp.append(cur); cur = []
        cur.append(len(Lc["words"]))
        if Lc["pend"]:
            zp.append(cur); cur = []
    zp = [p for p in zp if len(p) >= 4]
    ratio = np.array([p[-1] / np.median(p[1:-1]) for p in zp]); lastn = np.array([p[-1] for p in zp])
    out["voynich_zl3b_word_counts"] = summary(ratio, lastn)

    # real manuscripts: lines that stop early
    from control_real import token_table
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/user/Gallicorpora/HTR-MSS-15e-Siecle")
    parts = []
    for ms in ["btv1b100737746", "btv1b55008562q", "btv1b90069505", "btv1b9009670t", "btv1b90595325", "btv1b90616967"]:
        d = token_table(root / "data" / ms)
        d["folio"] = ms + "/" + d.folio; d["par"] = ms + "/" + d.par.astype(str)
        parts.append(d)
    for name, df in (("french_lines_stopping_early", pd.concat(parts, ignore_index=True)),
                     ("voynich_lines_stopping_early", pt)):
        M = margins(line_table(df)); M = M[M.M1.notna()]
        M["fill"] = (M.E - M.xl) / (M.M1 - M.xl)
        S = M[(M.fill < 0.9) & (M.fill > 0)]
        out[name] = summary(S.fill, S.n)

    json.dump(out, open(ROOT / "results" / "exp16_last_lines.json", "w"), indent=1)
    show = lambda k, r: print("%-44s n=%4d  one word %.3f  1-2 words %.3f  under 30%% %.3f  40-80%% %.3f  fifths %s" % (
        k, r["n"], r["one_word"], r["up_to_two"], r["under_30"], r["mid_40_80"], r["fifths"]))
    for k, r in out.items():
        if k == "voynich_rewrapped":
            for d, rr in r.items():
                show("  re-wrapped, line width %s" % d, rr)
        else:
            show(k, r)


if __name__ == "__main__":
    main()
