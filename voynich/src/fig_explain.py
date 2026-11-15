"""Explanatory figure: one real paragraph as written, and the same words wrapped."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

sys.path.insert(0, str(Path(__file__).parent))
from figures import AQUA, BLUE, FIG, GREY, GRID, INK, INK2, ORANGE, SURFACE, head
from layout import WidthModel
from vdata import paragraph_text, tokens
from wrap import line_table, margins, simulate

FOLIO = "f76r"


def panel(ax, S, Lt, lines, title):
    """Word boxes to scale.  The blank is measured the way the analysis measures it: from
    the line's end to the furthest right end among its neighbouring lines."""
    ys = {ln: i for i, ln in enumerate(lines)}
    blanks = []
    for ln in lines:
        g = S[S.line == ln].sort_values("k")
        y = -ys[ln]
        E = g.x1.max()
        M = float(Lt[Lt.line == ln].M1.iloc[0])
        blanks.append(max(M - E, 0) / float(g.wmed.iloc[0]))
        if M > E:
            ax.add_patch(Rectangle((E, y - 0.36), M - E, 0.72, color=ORANGE, linewidth=0))
        for r in g.itertuples():
            last = r.kr == 0
            ax.add_patch(Rectangle((r.x0 + 1, y - 0.30), r.x1 - r.x0 - 2, 0.60, linewidth=0,
                                   facecolor=BLUE if last else GRID))
            ax.text((r.x0 + r.x1) / 2, y, r.word, ha="center", va="center", fontsize=5.6,
                    color="white" if last else INK2)
    ax.set_ylim(-len(lines) + 0.3, 0.9)
    ax.set_xlim(S.x0.min() - 6, S.x1.max() + 40)
    ax.axis("off")
    ax.set_title(f"{title}   (median blank here: {np.median(blanks):.2f} of a word)", loc="left",
                 fontsize=10.5, fontweight="bold", color=INK, pad=4)


def main():
    pt = paragraph_text(tokens()); wm = WidthModel(pt)
    L = margins(line_table(pt, wm))
    g = pt[pt.folio == FOLIO]
    # the longest paragraph on the page, first 11 lines
    par = g.groupby("par").line.nunique().idxmax()
    lines = sorted(g[g.par == par].line.unique())[2:14]     # twelve lines from the middle
    A = g[g.line.isin(lines)]
    Sim = simulate(pt[pt.folio == FOLIO], L[L.folio == FOLIO], wm, np.random.default_rng(3), mode="pour", tol=0.1)
    Sim["par_last"] = False
    Ls = margins(line_table(Sim, wm))
    B = Sim[Sim.line.isin(lines)]
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 8.0))
    panel(axes[0], A, L[L.folio == FOLIO], lines, "As written")
    panel(axes[1], B, Ls, lines, "Same words, wrapped")
    fig.text(0.012, 0.012, "Blue: last word of each line.   Orange: blank up to where neighbouring lines reach.   "
             f"Page {FOLIO}, word boxes to scale.", fontsize=8.8, color=INK2, ha="left", va="bottom")
    head(fig, "Line ends are more even than wrapping would leave them",
         "One paragraph as written, and the same words wrapped. Whole manuscript: 0.13 against 0.32")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.865, bottom=0.05, hspace=0.14)
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "fig0_explain.png", dpi=200)
    print("wrote fig0_explain.png", len(lines), "lines")


if __name__ == "__main__":
    main()
