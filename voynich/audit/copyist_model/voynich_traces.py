"""Do the traces the copyist strategies would leave show up in the manuscript?

Reads the study's own data layer (../../src/vdata.py, read only).

A. Displaced repeats (R1).  A copyist who dodges a repeated line opening by moving the
   break by one word leaves the dodged word as the last word of the line above or the
   second word of the line below.  measure.displacement_profile on the ZL3b paragraphs,
   next to the same profile from simulated copyists whose R1 matches the Voynich.

B. Elastic writing (R1-R3).  Per line: written width of the words over the width their
   glyphs predict (layout.WidthModel).  The strategies that work by squeezing and
   stretching need whole lines to swing by 3.5-4.7% (sd); how much do Voynich lines swing?

C. Elastic tuning of the last line (R2).  Is a paragraph written looser where its last
   line came out short and denser where it came out long, as the simulated copyists do?

B and C see letter spacing and size only.  A copyist who adjusts by abbreviating changes
the glyphs themselves and leaves nothing for them to see.

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


def _voynich_lines():
    """Word and line tables from the boxes: written and glyph-predicted widths."""
    import pandas as pd
    from layout import WidthModel
    from vdata import paragraph_text, tokens
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    wm = WidthModel(pt)
    pred = lambda w: wm.coef[-1] + sum(wm.coef[wm.idx[c]] for c in w if c in wm.idx)
    pt = pt.assign(pw=pt.word.map(pred).clip(lower=0.1), aw=pt.w / pt.wmed)
    gap = (pt.groupby(["folio", "line"]).x0.shift(-1) - pt.x1) / pt.wmed
    g0 = float(np.nanmedian(gap[(gap > -0.5) & (gap < 1.0)]))
    pt["gap"] = gap
    ok = (gap > -0.3) & (gap < 1.0)                                     # an ordinary word space follows
    pt["lr"] = np.log(pt.aw.clip(lower=0.05) / pt.pw)                   # word alone
    pt["lra"] = np.log((pt.aw + gap.where(ok, 0.0)).clip(lower=0.05) / (pt.pw + np.where(ok, g0, 0.0)))  # word + gap
    rows = []
    for (f, p), g in pt.groupby(["folio", "par"], sort=False):
        lines = [h for _, h in g.groupby("line", sort=True)]
        for i, h in enumerate(lines):
            rows.append(dict(folio=f, par=f"{f}/{p}", line=h.line.iloc[0], i=i, d=len(lines) - 1 - i, n=len(lines),
                             nw=len(h), dens=float(np.log(h.aw.sum() / h.pw.sum())),
                             xl=h.x0.min() / h.wmed.iloc[0], E=h.x1.max() / h.wmed.iloc[0],
                             broken=bool((h.gap > 1.0).any())))
    L = pd.DataFrame(rows)
    return pt.merge(L[["folio", "line", "i", "d", "n", "nw", "broken"]], on=["folio", "line"]), L, wm


def part_b(pt, L, shuffles=100, seed=0):
    """Elastic writing shows as whole lines written wider or narrower than their glyphs
    predict.  Line-to-line sd of log(written / predicted), and how much of it is left once
    the same words' deviations are scattered at random over the lines of the paragraph
    (first, last and inner words kept apart): the coherent, per-line component."""
    rng = np.random.default_rng(seed)
    body = pt[(pt.i >= 1) & (pt.d >= 1) & (pt.nw >= 5) & ~pt.broken].copy()
    body["pos"] = np.where(body.k == 0, "f", np.where(body.kr == 0, "l", "m"))
    out = {}
    for col, name in (("lr", "words only"), ("lra", "words and the gaps after them")):
        for grp, gname in (("par", "within paragraph"), ("folio", "within page")):
            body["r"] = body[col] - body.groupby([grp, "pos"])[col].transform("mean")
            obs = float(body.groupby(["folio", "line"]).r.mean().std())
            null = []
            for _ in range(shuffles):
                body["rs"] = body.groupby([grp, "pos"]).r.transform(lambda v: rng.permutation(v.to_numpy()))
                null.append(body.groupby(["folio", "line"]).rs.mean().std())
            out[f"{name}, {gname}"] = dict(line_sd=obs, scattered_sd=float(np.mean(null)),
                                           coherent_sd=float(np.sqrt(max(0.0, obs ** 2 - np.mean(null) ** 2))))
        body["r"] = body[col] - body.groupby(["par", "pos"])[col].transform("mean")
        z = body[body.kr <= 2].groupby(["folio", "line"]).r.mean()
        null = []
        for _ in range(shuffles):
            body["rs"] = body.groupby(["par", "pos"]).r.transform(lambda v: rng.permutation(v.to_numpy()))
            null.append(body[body.kr <= 2].groupby(["folio", "line"]).rs.mean().std())
        out[f"{name}, last three words of the line, within paragraph"] = dict(
            line_sd=float(z.std()), scattered_sd=float(np.mean(null)),
            coherent_sd=float(np.sqrt(max(0.0, z.std() ** 2 - np.mean(null) ** 2))))
    out["lines"] = int(body.groupby(["folio", "line"]).ngroups)
    return out


def part_c(L):
    """A copyist who tunes the last line by loosening or tightening the lines above it
    leaves looser paragraphs where the last line came out short-to-middling and denser
    ones where it came out long: paragraph density falls as last-line fill rises.
    The same statistic from simulated copyists that reproduce R2."""
    import pandas as pd
    B = L[(L.i >= 1) & (L.d >= 1) & (L.nw >= 4) & ~L.broken].copy()
    B["dp"] = B.dens - B.groupby("folio").dens.transform("mean")
    rows = []
    for par, g in L.groupby("par", sort=False):
        g = g.sort_values("i")
        b = B[B.par == par]
        if len(g) < 4 or len(b) < 2:
            continue
        body = g.iloc[1:-1]; last = g.iloc[-1]
        rows.append((float((last.E - last.xl) / np.median(body.E - body.xl)), float(b.dp.mean())))
    V = np.array(rows)

    def stats(R):
        slope, icpt = np.polyfit(R[:, 0], R[:, 1], 1)
        res = R[:, 1] - (slope * R[:, 0] + icpt)
        se = float(np.sqrt(res.var() / (len(R) * R[:, 0].var())))
        bins = {}
        for lo, hi in ((0, .3), (.3, .45), (.45, .6), (.6, .75), (.75, .9), (.9, 9)):
            m = (R[:, 0] >= lo) & (R[:, 0] < hi)
            bins[f"{lo:.2f}-{hi:.2f}"] = dict(n=int(m.sum()), density=float(R[m, 1].mean()) if m.any() else None)
        return dict(paragraphs=len(R), corr=float(np.corrcoef(R[:, 0], R[:, 1])[0, 1]), slope=float(slope),
                    slope_se=se, density_sd=float(R[:, 1].std()), by_last_line_fill=bins)
    out = {"voynich": stats(V)}
    from copyist import Scribe, lay_text
    from texts import load
    el = lambda x, z=None: dict(squeeze=x, stretch=x, zone=z)
    sims = {"no goal, elastic whole line 8%": dict(el(.08)),
            "R2 goal, whole paragraph, elastic whole line 8%, casting-off error 3%":
                dict(el(.08), widow="shape", lookahead=None, plan_noise=.03),
            "R2 goal, last 4 lines, elastic whole line 10%": dict(el(.10), widow="shape", lookahead=4),
            "R2 goal, last 4 lines, elastic last-3 words 30%": dict(el(.30, 3), widow="shape", lookahead=4)}
    for name, kw in sims.items():
        R = []
        for t in ("caesar", "pliny", "genesis", "recipes"):
            lay, _ = lay_text(load(t, max_words=25000), Scribe(**kw), seed=1)
            for par in lay:
                if len(par) >= 4:
                    R.append((par[-1]["E"] / np.median([l["E"] for l in par[1:-1]]),
                              float(np.mean([np.log(l["E"] / l["N0"]) for l in par[1:-1]]))))
        out["sim: " + name] = stats(np.array(R))
    return out


def main():
    out = {"A_displaced_repeats": part_a()}
    try:
        pt, L, wm = _voynich_lines()
        out["width_model"] = dict(r2=wm.r2, sigma=wm.sigma)
        out["B_line_density_swing"] = part_b(pt, L)
        out["C_paragraph_density_vs_last_line"] = part_c(L)
    except Exception as e:                                   # word boxes not available
        out["B_C_error"] = repr(e)
    json.dump(out, open(HERE / "results" / "voynich_traces.json", "w"), indent=1)
    A = out["A_displaced_repeats"]
    print("A. displaced repeats.  Voynich R1: %.3f vs %.3f shuffled (ratio %.2f)" % (
        A["voynich_r1"]["same"], A["voynich_r1"]["shuffled"], A["voynich_r1"]["ratio"]))
    keys = list(A["voynich"].keys())
    print("   lift over chance (z) of 'begins with the glyph that opens line i', by position")
    print("%-58s" % "" + "".join("%13s" % k.replace("same line, ", "own ").replace("line below, ", "below ")
                                 .replace(" word", "")[:12] for k in keys))
    print("%-58s" % "Voynich" + "".join("%6.2f %+5.1f " % (A["voynich"][k]["lift"], A["voynich"][k]["z"]) for k in keys))
    for name, r in A.items():
        if name.startswith("sim"):
            print("%-58s" % (name[5:] + " [R1 %.2f]" % r["r1"]["ratio"])[:58].ljust(58) +
                  "".join("%6.2f %+5.1f " % (r["profile"][k]["lift"], r["profile"][k]["z"]) for k in keys))
    if "B_line_density_swing" in out:
        print("\nB. line-to-line swing of written/predicted width (Voynich; simulated elastic copyists need 0.035-0.047)")
        for k, v in out["B_line_density_swing"].items():
            if isinstance(v, dict):
                print("   %-62s line sd %.4f  scattered %.4f  coherent %.4f" % (k, v["line_sd"], v["scattered_sd"], v["coherent_sd"]))
        print("\nC. paragraph density against last-line fill")
        for k, v in out["C_paragraph_density_vs_last_line"].items():
            print("   %-72s n=%4d corr %+.3f slope %+.4f (se %.4f) | by fill: %s" % (
                k[:72], v["paragraphs"], v["corr"], v["slope"], v["slope_se"],
                " ".join("%+.3f" % b["density"] if b["density"] is not None else "  -  " for b in v["by_last_line_fill"].values())))


if __name__ == "__main__":
    main()
