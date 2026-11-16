"""Devil's-advocate test: can a COPYIST of a pre-existing text produce the Voynich's three
layout regularities?  Lays real Latin texts out with copyist.Scribe under graded freedoms
and measures R1-R3 as the manuscript study does.

    python3 run_experiments.py            # all sweeps -> results/sweeps.json, prints tables
    python3 run_experiments.py r1         # one family: r3 | r1 | r2 | all3

Every row is the mean over the texts in TEXT_SET (paragraphs cut at sentence ends to 4-8
lines, 25,000 words each); [lo-hi] is the range across texts for the headline number.
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from copyist import Scribe, lay_text
from measure import summary
from texts import load

HERE = Path(__file__).parent
TEXT_SET = ("caesar", "pliny", "genesis", "recipes")
MAX_WORDS = 25000
_cache = {}


def paragraphs(name, mode="matched"):
    if (name, mode) not in _cache:
        _cache[(name, mode)] = load(name, mode, max_words=MAX_WORDS)
    return _cache[(name, mode)]


def run_one(job):
    fam, label, kw, text, mode = job
    sc = Scribe(**kw)
    lay, T = lay_text(paragraphs(text, mode), sc, seed=1)
    s = summary(lay, T.unit)
    s["W_letters"] = T.W
    return fam, label, text, s


# ------------------------------------------------------------------ freedoms
def el(x, zone=3, sq=None, st=None):
    return dict(squeeze=x if sq is None else sq, stretch=x if st is None else st, zone=zone)


FREEDOMS_R3 = [
    ("none (plain wrapper)", {}),
    ("squeeze last word 30%", el(0, 1, sq=.30, st=0)),
    ("elastic last 3 words 5%", el(.05)), ("elastic last 3 words 10%", el(.10)),
    ("elastic last 3 words 15%", el(.15)), ("elastic last 3 words 20%", el(.20)),
    ("elastic last 3 words 30%", el(.30)),
    ("elastic whole line 2%", el(.02, None)), ("elastic whole line 3%", el(.03, None)),
    ("elastic whole line 5%", el(.05, None)), ("elastic whole line 8%", el(.08, None)),
    ("abbreviation table x0.5, last 3 words", dict(abbr=.5)),
    ("abbreviation table x1.0, last 3 words", dict(abbr=1.0)),
    ("abbreviation table x0.25, whole line", dict(abbr=.25, zone=None)),
    ("abbreviation table x0.5, whole line", dict(abbr=.5, zone=None)),
    ("abbreviation table x1.0, whole line", dict(abbr=1.0, zone=None)),
    ("drop <=1 letter in half the words, whole line", dict(abbr_generic=(1, .5), zone=None)),
    ("drop <=2 letters in half the words, whole line", dict(abbr_generic=(2, .5), zone=None)),
    ("drop <=2 letters in half the words, last 3", dict(abbr_generic=(2, .5))),
    ("split words at syllables", dict(split="syllable")),
    ("split words anywhere", dict(split="any")),
    ("line filler up to 0.5 word", dict(filler=.5)), ("line filler up to 1 word", dict(filler=1.0)),
    ("abbr x0.5 + elastic 10%, last 3 words", dict(abbr=.5, **el(.10))),
    ("abbr x0.5 + elastic 5%, last 3 words", dict(abbr=.5, **el(.05))),
    ("squeeze-only last 3 words 20% + filler 0.5", dict(filler=.5, **el(0, 3, sq=.2, st=0))),
]

FREEDOMS_R1 = [
    ("none", {}),
    ("elastic last 3 words 10%", el(.10)), ("elastic last 3 words 15%", el(.15)),
    ("elastic last 3 words 20%", el(.20)), ("elastic last 3 words 30%", el(.30)),
    ("elastic whole line 2%", el(.02, None)), ("elastic whole line 3%", el(.03, None)),
    ("elastic whole line 5%", el(.05, None)), ("elastic whole line 8%", el(.08, None)),
    ("elastic whole line 10%", el(.10, None)),
    ("abbr table x0.25, whole line", dict(abbr=.25, zone=None)),
    ("abbr table x0.5, whole line", dict(abbr=.5, zone=None)),
    ("abbr table x1.0, whole line", dict(abbr=1.0, zone=None)),
    ("abbr table x1.0, last 3 words", dict(abbr=1.0)),
    ("abbr x0.5 + elastic 10%, last 3 words", dict(abbr=.5, **el(.10))),
]
LOOK = [(1, "1"), (2, "2"), (3, "3"), (None, "par")]

FREEDOMS_R2 = [
    ("none", {}),
    ("elastic whole line 3%", el(.03, None)), ("elastic whole line 5%", el(.05, None)),
    ("elastic whole line 8%", el(.08, None)), ("elastic whole line 10%", el(.10, None)),
    ("elastic whole line 15%", el(.15, None)),
    ("tighten only: squeeze whole line 8%", el(0, None, sq=.08, st=0)),
    ("tighten only: squeeze whole line 15%", el(0, None, sq=.15, st=0)),
    ("loosen only: stretch whole line 8%", el(0, None, sq=0, st=.08)),
    ("loosen only: stretch whole line 15%", el(0, None, sq=0, st=.15)),
    ("elastic last 3 words 20%", el(.20)), ("elastic last 3 words 30%", el(.30)),
    ("abbr table x0.5, whole line", dict(abbr=.5, zone=None)),
    ("abbr table x1.0, whole line", dict(abbr=1.0, zone=None)),
]


def jobs_for(family):
    J = []
    add = lambda fam, label, kw, texts=TEXT_SET, mode="matched": J.extend(
        (fam, label, kw, t, mode) for t in texts)
    if family in ("r3", "all"):
        for lab, kw in FREEDOMS_R3:
            add("R3", lab, kw)
    if family in ("r1", "all"):
        for lab, kw in FREEDOMS_R1:
            add("R1 baseline (no goal)", lab, kw)
            add("R1 baseline (no goal), coarse", lab, dict(kw, coarse_glyphs=True))
            for L, ln in LOOK:
                add(f"R1 re-flow alone, look-ahead {ln}", lab, dict(kw, lookahead=L, avoid_repeat=1000))
                if ln in ("1", "3", "par"):
                    add(f"R1 re-flow alone, coarse glyphs, look-ahead {ln}", lab,
                        dict(kw, lookahead=L, avoid_repeat=1000, coarse_glyphs=True))
        for lev in ("modest", "generous", "max"):
            for L, ln in LOOK[:2]:
                add("R1 spelling variants alone", f"{lev}, look-ahead {ln}",
                    dict(variants=lev, lookahead=L, avoid_repeat=1000))
                add("R1 spelling variants alone, coarse glyphs", f"{lev}, look-ahead {ln}",
                    dict(variants=lev, lookahead=L, avoid_repeat=1000, coarse_glyphs=True))
        for sp in ("syllable", "any"):
            add("R1 baseline (no goal)", f"split {sp}", dict(split=sp))
            for L, ln in (LOOK[0], LOOK[2]):
                add("R1 word-splitting alone", f"split {sp}, look-ahead {ln}",
                    dict(split=sp, lookahead=L, avoid_repeat=1000))
                add("R1 word-splitting alone", f"split {sp}, coarse glyphs, look-ahead {ln}",
                    dict(split=sp, lookahead=L, avoid_repeat=1000, coarse_glyphs=True))
        for lab, kw in (("elastic whole 3% + modest variants", dict(variants="modest", **el(.03, None))),
                        ("elastic last-3 10% + modest variants", dict(variants="modest", **el(.10))),
                        ("elastic last-3 10% + generous variants", dict(variants="generous", **el(.10))),
                        ("elastic whole 5% + lag-2 seeking", dict(seek_lag2=300, **el(.05, None))),
                        ("elastic whole 8% + lag-2 seeking", dict(seek_lag2=300, **el(.08, None)))):
            for L, ln in (LOOK[0], LOOK[2]):
                add("R1 combinations", f"{lab}, look-ahead {ln}", dict(kw, lookahead=L, avoid_repeat=1000))
                add("R1 combinations, coarse glyphs", f"{lab}, look-ahead {ln}",
                    dict(kw, lookahead=L, avoid_repeat=1000, coarse_glyphs=True))
    if family in ("r2", "all"):
        for lab, kw in FREEDOMS_R2:
            add("R2 baseline (no goal)", lab, kw)
            for L, ln in ((2, "2"), (3, "3"), (4, "4"), (None, "par")):
                add(f"R2 shape goal, re-flow last {ln} lines", lab, dict(kw, lookahead=L, widow="shape"))
            for goal in ("single", "short"):
                for L, ln in ((3, "3"), (None, "par")):
                    add(f"R2 goal '{goal}', re-flow last {ln} lines", lab, dict(kw, lookahead=L, widow=goal))
        for lab, kw in (("elastic whole line 5%", el(.05, None)), ("elastic whole line 8%", el(.08, None))):
            for L, ln in ((3, "3"), (None, "par")):
                add(f"R2 shape goal, NATURAL paragraphs, re-flow last {ln} lines", lab,
                    dict(kw, lookahead=L, widow="shape"), texts=("caesar", "pliny", "recipes"), mode="natural")
    if family in ("extra", "all"):
        # word-splitting combined with a little elasticity (the usual pairing in Latin books)
        for lab, kw in (("split syllable + elastic last-3 5%", dict(split="syllable", **el(.05))),
                        ("split syllable + elastic last-3 10%", dict(split="syllable", **el(.10))),
                        ("split syllable + abbr x0.5 last-3", dict(split="syllable", abbr=.5)),
                        ("split syllable + generous variants", dict(split="syllable", variants="generous"))):
            add("R1 extra: baseline (no goal)", lab, kw)
            add("R1 extra: baseline (no goal), coarse", lab, dict(kw, coarse_glyphs=True))
            for L, ln in (LOOK[0], LOOK[2]):
                add("R1 extra: splitting combined", f"{lab}, look-ahead {ln}",
                    dict(kw, lookahead=L, avoid_repeat=1000))
                add("R1 extra: splitting combined, coarse glyphs", f"{lab}, look-ahead {ln}",
                    dict(kw, lookahead=L, avoid_repeat=1000, coarse_glyphs=True))
        # how much does the result depend on the cap on deliberate blanks?
        for cap in (0.0, 0.3, 0.6, 1.0, 1.5):
            for lab, kw in (("elastic last-3 15%", el(.15)), ("elastic whole 5%", el(.05, None)),
                            ("elastic last-3 20%", el(.20))):
                for L, ln in (LOOK[0], LOOK[2]):
                    add(f"R1 extra: blank cap {cap}, coarse glyphs", f"{lab}, look-ahead {ln}",
                        dict(kw, lookahead=L, avoid_repeat=1000, blank_cap=cap, coarse_glyphs=True))
        # paragraphs as they stand in the files (longer than the Voynich's)
        for lab, kw in (("elastic whole 5%", el(.05, None)), ("elastic whole 8%", el(.08, None)),
                        ("elastic last-3 20%", el(.20))):
            for L, ln in (LOOK[0], LOOK[2]):
                add("R1 extra: NATURAL paragraphs, coarse glyphs", f"{lab}, look-ahead {ln}",
                    dict(kw, lookahead=L, avoid_repeat=1000, coarse_glyphs=True),
                    texts=("caesar", "pliny", "recipes"), mode="natural")
    if family in ("all3", "all"):
        for lab, kw in (("elastic whole 8%", el(.08, None)), ("elastic whole 10%", el(.10, None)),
                        ("abbr table x1.0 whole line", dict(abbr=1.0, zone=None)),
                        ("abbr x0.5 + elastic whole 5%", dict(abbr=.5, **el(.05, None)))):
            for L, ln in ((3, "3"), (None, "par")):
                add("ALL THREE at once (coarse glyphs)", f"{lab}, look-ahead {ln}",
                    dict(kw, lookahead=L, avoid_repeat=1000, widow="shape", coarse_glyphs=True))
                add("ALL THREE at once (coarse glyphs, lag-2 seeking)", f"{lab}, look-ahead {ln}",
                    dict(kw, lookahead=L, avoid_repeat=1000, seek_lag2=300, widow="shape", coarse_glyphs=True))
    return J


def aggregate(results):
    by = {}
    for fam, label, text, s in results:
        by.setdefault(fam, {}).setdefault(label, {})[text] = s
    out = {}
    for fam, labs in by.items():
        out[fam] = {}
        for label, per in labs.items():
            keys = next(iter(per.values())).keys()
            mean = {k: float(np.nanmean([per[t][k] for t in per])) for k in keys}
            rng = {k: [float(np.nanmin([per[t][k] for t in per])), float(np.nanmax([per[t][k] for t in per]))]
                   for k in ("r1_ratio", "r2_one_word", "r2_under_30", "r3_median", "r3_coupling")}
            out[fam][label] = dict(mean=mean, range=rng, per_text=per)
    return out


def show(agg):
    for fam, labs in agg.items():
        print("\n== " + fam)
        for label, r in labs.items():
            m, g = r["mean"], r["range"]
            if fam.startswith("R3"):
                print("  %-46s blank median %.2f [%.2f-%.2f] iqr %.2f p90 %.2f | coupling %+.2f [%+.2f,%+.2f] | "
                      "cannot reach margin %2.0f%% | density sd %.3f  splits %2.0f%%" % (
                          label, m["r3_median"], *g["r3_median"], m["r3_iqr"], m["r3_p90"], m["r3_coupling"],
                          *g["r3_coupling"], 100 * m["r3_cannot_reach_margin"], m["t_adj_sd"], 100 * m["t_split_rate"]))
            elif fam.startswith("R1") or fam.startswith("ALL"):
                extra = ""
                if fam.startswith("ALL"):
                    extra = " | R2 one %.3f <30%% %.3f mid %.2f over %.3f" % (
                        m["r2_one_word"], m["r2_under_30"], m["r2_mid_40_80"], m["r2_over_full"])
                print("  %-46s same %.3f vs shuffled %.3f ratio %.2f [%.2f-%.2f] | lag2 %.3f (x%.2f) lag3 %.3f (x%.2f) | "
                      "R3 med %.2f cpl %+.2f | last-word %.3f 2nd-word-below %.3f | dens sd %.3f var %2.0f%% frag %2.0f%%%s" % (
                          label, m["r1_same"], m["r1_shuffled"], m["r1_ratio"], *g["r1_ratio"], m["r1_lag2_same"],
                          m["r1_lag2_ratio"], m["r1_lag3_same"], m["r1_lag3_ratio"], m["r3_median"], m["r3_coupling"],
                          m["t_disp_last"], m["t_disp_second"], m["t_adj_sd"], 100 * m["t_variant_rate"],
                          100 * m["t_fragment_starts"], extra))
            else:
                print("  %-40s one word %.3f [%.3f-%.3f] | <30%% %.3f [%.3f-%.3f] | 40-80%% %.2f | over-full %.3f | "
                      "R3 med %.2f cpl %+.2f | last-3 lines: mean %+.3f loose>5%% %.2f tight>5%% %.2f ; "
                      "earlier: %+.3f %.2f %.2f" % (
                          label, m["r2_one_word"], *g["r2_one_word"], m["r2_under_30"], *g["r2_under_30"],
                          m["r2_mid_40_80"], m["r2_over_full"], m["r3_median"], m["r3_coupling"], m["t_pen_mean"],
                          m["t_pen_loose5"], m["t_pen_tight5"], m["t_early_mean"], m["t_early_loose5"],
                          m["t_early_tight5"]))


def main():
    family = sys.argv[1] if len(sys.argv) > 1 else "all"
    J = jobs_for(family)
    print(f"{len(J)} layouts ...", flush=True)
    with Pool(4) as pool:
        res = pool.map(run_one, J, chunksize=4)
    agg = aggregate(res)
    out = HERE / "results" / ("sweeps.json" if family == "all" else f"sweeps_{family}.json")
    json.dump(agg, open(out, "w"), indent=1)
    show(agg)
    print("\nwritten", out)


if __name__ == "__main__":
    main()
