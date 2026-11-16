"""R2 in detail: what is the least copyist freedom that takes one-word paragraph endings
from about 12% (a text that simply runs out) to about 2%, and two-words-or-fewer endings
from about 24% to 9-11% - and what would each way of doing it leave on the page?

Ways a copyist can get rid of a stray last word, each run separately:
  pull up     squeeze the line(s) before it so the stray word fits        (squeeze only)
  push down   stretch the line(s) before it so another word moves down    (stretch only)
  stop short  end an earlier line a word or two early, leaving a blank    (no elasticity)
  both        squeeze or stretch, whichever is less work
  abbreviate  the same through optional abbreviation (changes the glyphs, not the spacing)
over the last 2, 3 or 4 lines or the whole paragraph, with the goal 'single' (no one-word
last line) or 'two' (also, less pressingly, no two-word one).

Traces, in the simulation and (same definitions, word boxes) in the manuscript:
  over-full      last line as long as a full line
  full, dense    full last lines (>= 85%) written denser than the paragraph's other lines
  pen short      line before the last stops early (< 90% of the paragraph's normal line):
                 in all paragraphs, and in those whose last line holds at most 0.40 of a
                 normal line's words against those where it holds 0.55 or more
  pen loose      line before the last written looser than the paragraph's other lines,
                 same split
  slope          paragraph density against last-line fill

    python3 r2_focus.py          -> results/r2_focus.json, prints the table
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from copyist import Scribe, lay_text
from measure import r2, r3
from texts import load

HERE = Path(__file__).parent
TEXT_SET = ("caesar", "pliny", "genesis", "recipes")
_cache = {}


def trace_stats(P):
    """P: list of paragraphs (>= 4 lines); each a list of (width, log density, words) per
    line, first line included."""
    over, full_d, part_d, pall = [], [], [], []
    ps = {"small": [], "big": []}; pl = {"small": [], "big": []}; R = []
    for par in P:
        w = np.array([x[0] for x in par]); d = np.array([x[1] for x in par]); nw = par[-1][2]
        body = np.median(w[1:-1]); fill = w[-1] / body
        over.append(fill >= 1.0)
        others = d[1:-2] if len(par) > 4 else d[1:-1]
        ref = np.nanmean(others)
        (full_d if fill >= 0.85 else part_d).append(d[-1] - np.nanmean(d[1:-1]))
        rel = nw / np.median([x[2] for x in par[1:-1]])     # last line's words over a normal line's
        key = "small" if rel <= 0.40 else "big" if rel >= 0.55 else None
        if key:
            ps[key].append(w[-2] < 0.9 * body); pl[key].append(d[-2] - ref)
        pall.append(w[-2] < 0.9 * body)
        R.append((fill, np.nanmean(d[1:-1])))
    R = np.array(R); R = R[~np.isnan(R[:, 1])]
    slope = float(np.polyfit(R[:, 0], R[:, 1], 1)[0])
    res = R[:, 1] - np.polyval(np.polyfit(R[:, 0], R[:, 1], 1), R[:, 0])
    m = lambda v: float(np.nanmean(v)) if len(v) else float("nan")
    se = lambda v: float(np.nanstd(v) / np.sqrt(max(len(v), 1)))
    return dict(paragraphs=len(P), over_full=m(over), pen_short_all=m(pall),
                full_last_density=m(full_d), other_last_density=m(part_d),
                full_last_density_se=se(full_d),
                pen_short_small=m(ps["small"]), pen_short_big=m(ps["big"]), n_small=len(ps["small"]), n_big=len(ps["big"]),
                pen_loose_small=m(pl["small"]), pen_loose_big=m(pl["big"]),
                pen_loose_diff_se=float(np.hypot(se(pl["small"]), se(pl["big"]))),
                slope=slope, slope_se=float(np.sqrt(res.var() / (len(R) * R[:, 0].var()))))


def run_one(job):
    label, kw, text = job
    if text not in _cache:
        _cache[text] = load(text, max_words=25000)
    lay, T = lay_text(_cache[text], Scribe(**kw), seed=1)
    a, c = r2(lay), r3(lay, T.unit)
    P = [[(ln["E"], float(np.log(ln["E"] / ln["N0"])), ln["nwords"]) for ln in par] for par in lay if len(par) >= 4]
    out = dict(one_word=a["one_word"], up_to_two=a["up_to_two"], under_30=a["under_30"], blank=c["median"],
               coupling=c["coupling"])
    out.update(trace_stats(P))
    return label, text, out


def strategies():
    S = [("nothing (text simply runs out)", {})]
    hz = [(2, "last 2 lines"), (3, "last 3 lines"), (4, "last 4 lines"), (None, "whole paragraph")]
    for goal in ("single", "two"):
        for L, hn in hz:
            g = dict(widow=goal, lookahead=L)
            for x in (.05, .08, .10, .12, .15):
                S.append((f"{goal} | pull up: squeeze {x:.0%} | {hn}", dict(g, squeeze=x, stretch=0, zone=None)))
                S.append((f"{goal} | push down: stretch {x:.0%} | {hn}", dict(g, squeeze=0, stretch=x, zone=None)))
            for x in (.03, .05, .06, .08, .10):
                S.append((f"{goal} | both: elastic {x:.0%} | {hn}", dict(g, squeeze=x, stretch=x, zone=None)))
            for x in (.15, .20, .30):
                S.append((f"{goal} | both: elastic last 3 words {x:.0%} | {hn}", dict(g, squeeze=x, stretch=x, zone=3)))
            for x in (.5, 1.0):
                S.append((f"{goal} | abbreviate: table x{x} | {hn}", dict(g, abbr=x, zone=None)))
            S.append((f"{goal} | stop short: blank up to 2.5 words | {hn}", dict(g, blank_cap=2.5)))
            S.append((f"{goal} | stop short + elastic last 3 words 15% | {hn}",
                      dict(g, blank_cap=2.5, squeeze=.15, stretch=.15, zone=3)))
    for L, hn in hz[1:]:
        g = dict(widow="two+", lookahead=L)
        for x in (.06, .08, .10, .12):
            S.append((f"two+ | both: elastic {x:.0%} | {hn}", dict(g, squeeze=x, stretch=x, zone=None)))
        S.append((f"two+ | both: elastic last 3 words 30% | {hn}", dict(g, squeeze=.3, stretch=.3, zone=3)))
        S.append((f"two+ | abbreviate: table x1.0 | {hn}", dict(g, abbr=1.0, zone=None)))
    return S


def voynich():
    sys.path.insert(0, str(HERE.parents[1] / "src"))
    from voynich_traces import _voynich_lines
    pt, L, wm = _voynich_lines()
    L["ill"] = L.folio.map(pt.groupby("folio").illust.first())
    ok = (L.nw >= 3) & ~L.broken
    L["dp"] = (L.dens - L[ok & (L.i >= 1) & (L.d >= 1)].groupby("folio").dens.transform("mean").reindex(L.index)
               .fillna(L.folio.map(L[ok].groupby("folio").dens.mean()))).where(ok)       # density relative to the page
    res = {}
    for name, sect in (("all sections", None), ("stars section (no drawings in the text)", ("S",)),
                       ("stars, jars and text-only pages", ("S", "P", "T")), ("plant pages", ("H",))):
        P = []
        for par, g in L.groupby("par", sort=False):
            g = g.sort_values("i")
            if len(g) < 4 or (sect and g.ill.iloc[0] not in sect):
                continue
            P.append([(float(r.E - r.xl), float(r.dp), int(r.nw)) for r in g.itertuples()])
        out = trace_stats(P)
        nw = np.array([p[-1][2] for p in P])
        out.update(one_word=float((nw == 1).mean()), up_to_two=float((nw <= 2).mean()),
                   words_per_line=float(np.median([x[2] for p in P for x in p[1:-1]])))
        res[name] = out
    return res


def main():
    S = strategies()
    jobs = [(lab, kw, t) for lab, kw in S for t in TEXT_SET]
    print(len(jobs), "layouts ...", flush=True)
    with Pool(4) as pool:
        res = pool.map(run_one, jobs, chunksize=4)
    agg = {}
    for lab, text, s in res:
        agg.setdefault(lab, []).append(s)
    out = {"sim": {lab: {k: float(np.nanmean([s[k] for s in v])) for k in v[0]} for lab, v in agg.items()}}
    try:
        out["voynich"] = voynich()
    except Exception as e:
        out["voynich_error"] = repr(e)
    json.dump(out, open(HERE / "results" / "r2_focus.json", "w"), indent=1)
    show(out)


def show(out):
    print("last line: one word / two or fewer (%) | last line over-full (%), full last lines denser than others (log) | "
          "line before last short (%): all paragraphs, short endings, long endings | line before last looser than the "
          "paragraph (log): short endings, long endings | slope of paragraph density on last-line fill")
    hdr = "%-62s %5s %5s | %5s %6s | %5s %6s %6s | %7s %7s | %7s" % (
        "", "1 wd", "<=2", "over", "dense", "short", "(sh.)", "(long)", "loose", "(long)", "slope")
    print(hdr)
    row = lambda lab, s: print("%-62s %5.1f %5.1f | %5.1f %+6.3f | %5.1f %6.1f %6.1f | %+7.3f %+7.3f | %+7.3f" % (
        lab[:62], 100 * s["one_word"], 100 * s["up_to_two"], 100 * s["over_full"], s["full_last_density"] - s["other_last_density"],
        100 * s["pen_short_all"], 100 * s["pen_short_small"], 100 * s["pen_short_big"], s["pen_loose_small"], s["pen_loose_big"],
        s["slope"]))
    for name, v in out.get("voynich", {}).items():
        row("VOYNICH %s (%d par., %.0f words/line)" % (name, v["paragraphs"], v["words_per_line"]), v)
        print("%-62s n short/long endings %d/%d; se: dense %.3f, loose difference %.3f, slope %.3f" % (
            "", v["n_small"], v["n_big"], v["full_last_density_se"], v["pen_loose_diff_se"], v["slope_se"]))
    for lab, s in out["sim"].items():
        row(lab, s)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "show":
        show(json.load(open(HERE / "results" / "r2_focus.json")))
    else:
        main()
