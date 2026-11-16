"""Do the traces the copyist strategies would leave show up in the manuscript?

Reads the study's own data layer (../../src/vdata.py, read only).

A. Displaced repeats (R1).  A copyist who dodges a repeated line opening by moving the
   break by one word leaves the dodged word as the last word of the line above or the
   second word of the line below.  measure.displacement_profile on the ZL3b paragraphs,
   next to the same profile from simulated copyists whose R1 matches the Voynich.

B. Elastic writing near the paragraph end (R2).  Per line: written width of the words
   over the width their glyphs predict (layout.WidthModel), and the mean gap, by distance
   from the paragraph's last line.  A copyist who loosens or tightens the last few lines
   to tune the last line leaves those lines looser or denser than the lines before.

Writes results/voynich_traces.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from measure import displacement_profile, r1          # noqa: E402


def zl_paragraphs():
    """As src/exp13_line_starts.zl_paragraphs: lists of lines (word lists) from ZL3b."""
    from vdata import load_zl
    pages, loci = load_zl()
    out, cur, key = [], [], None
    for L in loci:
        if not L["ltype"].startswith("P") or L["ltype"] == "Pt" or not L["words"]:
            if cur: out.append(cur)
            cur, key = [], None
            continue
        if L["pstart"] or L["pos"] in "@*" or key != L["folio"]:
            if cur: out.append(cur)
            cur, key = [], L["folio"]
        cur.append(L["words"])
        if L["pend"]:
            out.append(cur); cur, key = [], None
    if cur: out.append(cur)
    return out


def part_a():
    Z = zl_paragraphs()
    paras = [[[w[0] for w in ln] for ln in lines[1:]] for lines in Z]
    fake = [[dict(glyph=ln[0]) for ln in [["#"]] + p] for p in paras]       # r1() drops the first line
    out = {"voynich_r1": r1(fake)["lag1"], "voynich": displacement_profile(paras, shuffles=500)}
    # the same profile from copyists that reproduce R1 (coarse glyphs, pooled over four texts)
    from copyist import Scribe, lay_text
    from texts import load
    sims = {
        "no goal (elastic whole line 8%)": dict(squeeze=.08, stretch=.08, zone=None),
        "re-flow, look-ahead 1, elastic whole line 8%": dict(squeeze=.08, stretch=.08, zone=None, avoid_repeat=1000),
        "re-flow, look-ahead 1, elastic last-3 20%, blank cap 0.6": dict(squeeze=.2, stretch=.2, avoid_repeat=1000,
                                                                          blank_cap=.6),
        "re-flow, look-ahead 3, elastic last-3 20%": dict(squeeze=.2, stretch=.2, avoid_repeat=1000, lookahead=3),
        "re-flow, whole paragraph, elastic whole line 5%": dict(squeeze=.05, stretch=.05, zone=None,
                                                                avoid_repeat=1000, lookahead=None),
        "spelling variants only (max), look-ahead 2": dict(variants="max", avoid_repeat=1000, lookahead=2),
    }
    for name, kw in sims.items():
        pooled, lay_all = [], []
        for t in ("caesar", "pliny", "genesis", "recipes"):
            lay, _ = lay_text(load(t, max_words=25000), Scribe(coarse_glyphs=True, **kw), seed=1)
            pooled += [[ln["wg"] for ln in par[1:]] for par in lay]; lay_all += lay
        out["sim: " + name] = dict(r1=r1(lay_all)["lag1"], profile=displacement_profile(pooled))
    return out


def part_b():
    from layout import WidthModel
    from vdata import paragraph_text, tokens
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    wm = WidthModel(pt)
    pred = lambda w: wm.coef[-1] + sum(wm.coef[wm.idx[c]] for c in w if c in wm.idx)
    pt = pt.assign(pw=pt.word.map(pred), aw=pt.w / pt.wmed)
    q = pt.assign(gap=(pt.groupby(["folio", "line"]).x0.shift(-1) - pt.x1) / pt.wmed)
    rows = []
    for (f, p), g in q.groupby(["folio", "par"], sort=False):
        lines = [h for _, h in g.groupby("line", sort=True)]
        if len(lines) < 5:
            continue
        for i, h in enumerate(lines):
            if i == 0 or len(h) < 4:
                continue
            gaps = h.gap.dropna(); gaps = gaps[(gaps > -0.5) & (gaps < 1.0)]
            rows.append((f"{f}/{p}", len(lines) - 1 - i, float(np.log(h.aw.sum() / h.pw.sum())),
                         float(gaps.mean()) if len(gaps) else np.nan, len(h)))
    R = np.array([(r[1], r[2], r[3]) for r in rows]); par = np.array([r[0] for r in rows])
    # paragraph fixed effects: each line relative to its paragraph's lines of distance >= 4
    out = {"width_model_r2": wm.r2, "width_model_sigma": wm.sigma, "lines": len(rows)}
    base = {}
    for k in np.unique(par):
        m = (par == k) & (R[:, 0] >= 4)
        if m.sum():
            base[k] = (R[m, 1].mean(), np.nanmean(R[m, 2]))
    for d, name in ((0, "last line"), (1, "1 before last"), (2, "2 before last"), (3, "3 before last")):
        m = (R[:, 0] == d) & np.array([k in base for k in par])
        dens = R[m, 1] - np.array([base[k][0] for k in par[m]])
        gap = R[m, 2] - np.array([base[k][1] for k in par[m]])
        gap = gap[~np.isnan(gap)]
        out[name] = dict(n=int(m.sum()), word_width_vs_earlier=float(dens.mean()),
                         se=float(dens.std() / np.sqrt(len(dens))),
                         loose_over_5pct=float((dens > 0.05).mean()), tight_over_5pct=float((dens < -0.05).mean()),
                         gap_vs_earlier=float(gap.mean()), gap_se=float(gap.std() / np.sqrt(len(gap))))
    early = R[R[:, 0] >= 4, 1]
    out["line_density_sd_all"] = float(R[R[:, 0] >= 1, 1].std())
    return out


def main():
    out = {"A_displaced_repeats": part_a()}
    try:
        out["B_density_near_paragraph_end"] = part_b()
    except Exception as e:                                   # word boxes not available
        out["B_density_near_paragraph_end"] = {"error": repr(e)}
    json.dump(out, open(HERE / "results" / "voynich_traces.json", "w"), indent=1)
    A = out["A_displaced_repeats"]
    print("A. displaced repeats.  Voynich R1: %.3f vs %.3f shuffled (ratio %.2f)" % (
        A["voynich_r1"]["same"], A["voynich_r1"]["shuffled"], A["voynich_r1"]["ratio"]))
    keys = list(A["voynich"].keys())
    print("%-58s" % "" + "".join("%13s" % k.replace("same line, ", "own ").replace("line below, ", "below ")
                                 .replace(" word", "")[:12] for k in keys))
    print("%-58s" % "Voynich (lift; z)" + "".join("%6.2f %+5.1f " % (A["voynich"][k]["lift"], A["voynich"][k]["z"]) for k in keys))
    for name, r in A.items():
        if name.startswith("sim"):
            print("%-58s" % (name[5:] + " [R1 %.2f]" % r["r1"]["ratio"])[:58].ljust(58) +
                  "".join("%6.2f %+5.1f " % (r["profile"][k]["lift"], r["profile"][k]["z"]) for k in keys))
    B = out["B_density_near_paragraph_end"]
    print("\nB. written/predicted word width and gap, relative to the paragraph's earlier lines")
    for k, v in B.items():
        print("  ", k, {a: round(b, 4) for a, b in v.items()} if isinstance(v, dict) else v)


if __name__ == "__main__":
    main()
