"""Figures for the report.  Reads results/*.json (and recomputes two light maps)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from vdata import ROOT

# palette: validated categorical slots 1-3 (light surface), text and grid tokens
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8d8c86"
FIG = ROOT / "figures"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 10.5, "text.color": INK,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.linewidth": 1.0,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 1.0, "grid.linestyle": "-",
    "axes.axisbelow": True, "lines.linewidth": 2.0, "lines.solid_capstyle": "round",
    "legend.frameon": False, "legend.fontsize": 9.5,
})


def head(fig, title, sub, y=0.985):
    fig.text(0.012, y, title, fontsize=13.5, fontweight="bold", va="top", ha="left", color=INK)
    fig.text(0.012, y - 0.052 * 6.0 / fig.get_figheight(), sub, fontsize=10, va="top", ha="left",
             color=INK2)


def tidy(ax, grid="y"):
    ax.grid(False)
    ax.grid(True, axis=grid)
    ax.tick_params(length=0)


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / name, dpi=200)
    plt.close(fig)
    print("wrote", name)


def fig_slack(E):
    rows = [("The manuscript", E["slack_values_real"], BLUE),
            ("Simulated scribe wrapping a finished text", E["sim_examples"]["pourer"]["slack"], ORANGE),
            ("Simulated scribe choosing words that fit", E["sim_examples"]["fitter"]["slack"], AQUA)]
    fig, axes = plt.subplots(3, 1, figsize=(7.6, 6.6), sharex=True)
    edges = np.arange(-0.5, 2.001, 0.1)
    for ax, (name, v, col) in zip(axes, rows):
        v = np.asarray(v, float)
        h, _ = np.histogram(np.clip(v, -0.5, 2.0), edges)
        share = 100 * h / len(v)
        ax.bar(edges[:-1] + 0.05, share, width=0.078, color=col, linewidth=0)
        tidy(ax)
        ax.set_ylim(0, 34)
        ax.set_yticks([0, 10, 20, 30])
        ax.set_ylabel("% of line ends")
        med, q1, q3 = np.median(v), np.percentile(v, 25), np.percentile(v, 75)
        ax.text(0.995, 0.93, name, transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
                fontweight="bold", color=INK)
        ax.text(0.995, 0.74, f"median {med:.2f}   middle half spans {q3 - q1:.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color=INK2)
    axes[-1].set_xlabel("Blank space left at the end of the line, in average word widths "
                        "(values past 2 are stacked at 2)")
    head(fig, "The manuscript's lines are filled far more evenly than wrapping allows",
         "Blank left at the line end, measured against the neighbouring lines' right edge")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.875, bottom=0.085, hspace=0.16)
    save(fig, "fig1_slack.png")


def _binned(slack, wnext, edges):
    slack, wnext = np.asarray(slack, float), np.asarray(wnext, float)
    xs, ys = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (wnext >= lo) & (wnext < hi)
        if m.sum() >= 30:
            xs.append(np.median(wnext[m])); ys.append(np.median(slack[m]))
    return np.array(xs), np.array(ys)


def fig_coupling(E):
    edges = np.array([0.2, 0.6, 0.8, 0.95, 1.1, 1.25, 1.45, 1.7, 2.6])
    series = [("Manuscript", E["slack_values_real"], E["wnext_values_real"], BLUE),
              ("Wrapped text (simulated)", E["sim_examples"]["pourer"]["slack"],
               E["sim_examples"]["pourer"]["wnext"], ORANGE),
              ("Fitted text (simulated)", E["sim_examples"]["fitter"]["slack"],
               E["sim_examples"]["fitter"]["wnext"], AQUA)]
    fig, ax = plt.subplots(figsize=(7.6, 4.9))
    ends = []
    for name, s, w, col in series:
        x, y = _binned(s, w, edges)
        ax.plot(x, y, color=col, marker="o", markersize=7, markeredgecolor=SURFACE,
                markeredgewidth=1.5, label=name)
        ends.append((x[-1], y[-1], name))
    tidy(ax)
    ax.set_xlim(0.3, 2.45); ax.set_ylim(0, 1.0)
    ax.set_xlabel("Width of the word that starts the next line (average word widths)")
    ax.set_ylabel("Typical blank left at the end of this line")
    lab = {"Manuscript": "Manuscript", "Wrapped text (simulated)": "Wrapped text",
           "Fitted text (simulated)": "Fitted text"}
    offs = {"Manuscript": 0.03, "Wrapped text (simulated)": 0.0, "Fitted text (simulated)": -0.03}
    for x, y, name in ends:
        ax.text(x + 0.045, y + offs[name], lab[name], va="center", ha="left", fontsize=9.5, color=INK)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), handlelength=1.6)
    head(fig, "The line end knows nothing about the word that follows it",
         "If a line ended because the next word did not fit, bigger next words would leave bigger blanks")
    fig.subplots_adjust(left=0.09, right=0.865, top=0.855, bottom=0.115)
    save(fig, "fig2_coupling.png")


def fig_forest(E):
    real = E["real"]
    rows = [("All paragraph text", real["slope"], real["slope_ci"], real["n"])]
    order = ["herbal", "pharmaceutical", "balneological", "stars / recipes", "text only",
             "Currier A", "Currier B"] + sorted(k for k in E["subsets"] if k.startswith("scribe"))
    for k in order:
        if k in E["subsets"] and E["subsets"][k]["slope_ci"][1] - E["subsets"][k]["slope_ci"][0] < 0.8:
            v = E["subsets"][k]          # (subsets too small to say anything are left to the table)
            rows.append(((k[0].upper() + k[1:]).replace("Scribe", "Hand"), v["slope"], v["slope_ci"], v["n"]))
    pour, fit = E["sims"]["pourer"]["slope"], E["sims"]["fitter"]["slope"]
    fig, ax = plt.subplots(figsize=(7.6, 0.36 * len(rows) + 2.0))
    ys = np.arange(len(rows))[::-1]
    for (name, s, ci, n), y in zip(rows, ys):
        ax.plot(ci, [y, y], color=BLUE, linewidth=2)
        ax.plot([s], [y], "o", color=BLUE, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.axvspan(pour[2], pour[3], color=ORANGE, alpha=0.18, linewidth=0)
    ax.axvline(pour[0], color=ORANGE, linewidth=2)
    ax.axvspan(fit[2], fit[3], color=AQUA, alpha=0.18, linewidth=0)
    ax.axvline(fit[0], color=AQUA, linewidth=2)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{r[0]}  (n={r[3]:,})" for r in rows], color=INK)
    tidy(ax, "x")
    ax.set_xlim(-0.35, 0.72); ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("Coupling between a line's blank end and the next line's first word")
    top = len(rows) - 0.42
    ax.text(pour[0], top, f"wrapped text {pour[0]:.2f}", ha="center", va="bottom", fontsize=9.5, color=INK)
    ax.text(fit[0] - 0.012, top, f"fitted text {fit[0]:.2f}", ha="right", va="bottom", fontsize=9.5, color=INK)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color=BLUE, marker="o", markersize=7, label="Manuscript, 95% interval"),
                        Line2D([], [], color=ORANGE, label="Wrapped text, 20 simulations"),
                        Line2D([], [], color=AQUA, label="Fitted text, 20 simulations")],
               loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=3, handlelength=1.6, fontsize=9.5)
    head(fig, "Every section and every hand sits with the fitted text",
         "Line-break coupling by part of the manuscript, against the two simulated scribes")
    fig.subplots_adjust(left=0.3, right=0.975, top=0.86, bottom=0.175)
    save(fig, "fig3_sections.png")


def fig_lastword(D):
    x = np.array(D["room_mid"]);
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.2))
    ax = axes[0]
    ax.plot(x, D["glyphs"], color=BLUE, marker="o", markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.5)
    ax.axhline(D["glyphs_interior"], color=GREY, linewidth=1)
    ax.text(1.88, D["glyphs_interior"] - 0.12, f"words inside the line: {D['glyphs_interior']:.1f}",
            fontsize=9, color=INK2, va="top", ha="right")
    ax.set_ylim(0, 7); ax.set_ylabel("Glyphs in the last word of the line")
    ax.set_title("Its length follows the space", loc="left", fontsize=10.5, fontweight="bold", color=INK)
    ax = axes[1]
    ax.plot(x, 100 * np.array(D["end_m"]), color=BLUE, marker="o", markersize=7, markeredgecolor=SURFACE,
            markeredgewidth=1.5)
    ax.axhline(100 * D["end_m_interior"], color=GREY, linewidth=1)
    ax.text(x[-1], 100 * D["end_m_interior"] + 0.6, f"words inside the line: {100 * D['end_m_interior']:.1f}%",
            fontsize=9, color=INK2, va="bottom", ha="right")
    ax.set_ylim(0, 30); ax.set_ylabel("% of last words ending in the glyph m")
    ax.set_title("Its ending does not", loc="left", fontsize=10.5, fontweight="bold", color=INK)
    for ax in axes:
        tidy(ax)
        ax.set_xlim(0.3, 1.9)
        ax.set_xlabel("Room left when the last word was begun\n(average word widths)")
    head(fig, "The last word of a line is cut to measure",
         "Last word of each line, by how much room remained before it was written")
    fig.subplots_adjust(left=0.08, right=0.985, top=0.78, bottom=0.17, wspace=0.26)
    save(fig, "fig4_lastword.png")


def fig_kernel(K):
    c = np.array(K["centres"])
    fig, ax = plt.subplots(figsize=(7.6, 4.9))
    ser = [("Manuscript, 1-2 lines above", K["real_near"], BLUE),
           ("Manuscript, 6-20 lines above", K["real_far"], AQUA),
           ("Simulated 5% copier, 1-2 lines above", K["syn_near"], ORANGE)]
    for name, v, col in ser:
        ax.plot(c, v, color=col, marker="o", markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.5,
                label=name)
    ax.axhline(1, color=GREY, linewidth=1)
    tidy(ax)
    ax.set_ylim(0.8, 2.3); ax.set_xlim(-4.4, 4.4)
    ax.set_xlabel("Horizontal offset of the earlier word (average word widths; 0 = directly overhead)")
    ax.set_ylabel("Similarity, relative to two random words of the page")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), handlelength=1.6)
    i = int(np.argmax(K["syn_near"]))
    ax.annotate("what copying 1 word in 20\nfrom overhead looks like", xy=(c[i] - 0.12, K["syn_near"][i]),
                xytext=(-4.1, 1.95), fontsize=9.5, color=INK, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=INK2, linewidth=1))
    head(fig, "No extra likeness to the word directly overhead",
         "The mild hump at 0 is the same on lines far above, where copying cannot be the cause")
    fig.subplots_adjust(left=0.09, right=0.985, top=0.855, bottom=0.115)
    save(fig, "fig5_kernel.png")


def fig_calibration(X):
    cal, real, bound = X["calibration"], X["real"], X["bound"]
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    xs = [0.0] * len(cal["null"]) + [float(l) for l, runs in cal["phys"].items() for _ in runs]
    ys = [r["C"] for r in cal["null"]] + [r["C"] for runs in cal["phys"].values() for r in runs]
    ax.axhspan(real["C_ci"][0], real["C_ci"][1], color=BLUE, alpha=0.16, linewidth=0)
    ax.axhline(real["C"], color=BLUE, linewidth=2, label="Manuscript, with 95% interval")
    lx = np.array([0, 0.105])
    ax.plot(100 * lx, bound["intercept"] + bound["slope"] * lx, color=GREY, linewidth=1)
    ax.plot(100 * np.array(xs), ys, "o", color=ORANGE, markersize=7, markeredgecolor=SURFACE,
            markeredgewidth=1.5, label="Simulated pages with a known share of overhead copies")
    tidy(ax)
    ax.set_xlim(-0.4, 10.6); ax.set_ylim(0.8, 2.6)
    ax.set_xlabel("Share of words that are copies of the word physically overhead (%)")
    ax.set_ylabel("Overhead-similarity contrast C")
    ax.legend(loc="upper left", handlelength=1.6)
    ub = 100 * bound["lam_upper95"]
    ax.annotate(f"the manuscript allows at most {ub:.1f}%", xy=(ub, real["C_ci"][1]), xytext=(2.9, 1.5),
                fontsize=9.5, color=INK, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=INK2, linewidth=1))
    head(fig, "The test would have caught even 1-2% of words copied from overhead",
         "Contrast C on synthetic pages with known copying, against the manuscript's value")
    fig.subplots_adjust(left=0.09, right=0.985, top=0.845, bottom=0.125)
    save(fig, "fig6_calibration.png")


def extra_data():
    """Two light recomputations the figures need (cached in results/fig_extra.json)."""
    p = ROOT / "results" / "fig_extra.json"
    if p.exists():
        return json.loads(p.read_text())
    from kernel import DX_EDGES, lift_map
    from layout import Vocab, WidthModel, generate_onpage, page_templates
    from pairs import build_pairs
    from vdata import paragraph_text, tokens
    from wrap import breaks, line_table, margins
    pt = paragraph_text(tokens()); wm = WidthModel(pt)
    # last word vs the room it had
    L = margins(line_table(pt, wm)); B = breaks(L)
    B["room"] = (B.M1 - B.x0_last) / B.wmed
    edges = [0.3, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9]
    B["rb"] = pd.cut(B.room, edges)
    g = B.groupby("rb", observed=True)
    inner = pt[(pt.k > 0) & (pt.kr > 0)]
    D = dict(room_mid=[float(v) for v in g.room.median()], n=[int(v) for v in g.size()],
             glyphs=[float(v) for v in g.word_last.apply(lambda s: s.str.len().mean())],
             end_m=[float(v) for v in g.word_last.apply(lambda s: s.str.endswith("m").mean())],
             glyphs_interior=float(inner.word.str.len().mean()),
             end_m_interior=float(inner.word.str.endswith("m").mean()),
             end_m_last=float(B.word_last.str.endswith("m").mean()))
    # kernel maps
    bands = [(1, 2), (6, 20)]
    P = build_pairs(pt, dmax=20)
    m = lift_map(P[P.interior], "e1", "dx", DX_EDGES, bands, nboot=2)
    T = page_templates(pt); V = Vocab(pt.word)
    S = generate_onpage(T, V, wm, np.random.default_rng(42), rule="phys", lam=0.05)
    Ps = build_pairs(S, dmax=20)
    ms = lift_map(Ps[Ps.interior], "e1", "dx", DX_EDGES, bands, nboot=2)
    sel = slice(1, -1)                                   # drop the two wide outer bins
    K = dict(centres=((DX_EDGES[:-1] + DX_EDGES[1:]) / 2)[sel].tolist(),
             real_near=m["lift"][0][sel].tolist(), real_far=m["lift"][1][sel].tolist(),
             syn_near=ms["lift"][0][sel].tolist(), syn_far=ms["lift"][1][sel].tolist())
    out = dict(lastword=D, kernel=K)
    p.write_text(json.dumps(out, indent=1))
    return out


def main():
    E = json.loads((ROOT / "results" / "exp2_wrap.json").read_text())
    X = json.loads((ROOT / "results" / "exp1_kernel.json").read_text())
    extra = extra_data()
    fig_slack(E); fig_coupling(E); fig_forest(E)
    fig_lastword(extra["lastword"]); fig_kernel(extra["kernel"]); fig_calibration(X)


if __name__ == "__main__":
    main()
