"""Experiment 11 - do the words follow the pictures, or the scribe?
Writes results/exp11_sections.json.

Pages are grouped by scribe (hand, from ZL3b) and by what is drawn on them.  Two samples of
250 words are drawn from two groups (from different pages) and their vocabulary overlap is
measured: the share of words they have in common, counting repeats (sum over word types of
the smaller of the two frequencies).  Averaged over many draws.

If the words follow the pictures, pages with the same kind of picture should share
vocabulary even when different scribes wrote them.  If the words follow the scribe, one
scribe's pages should share vocabulary whatever is drawn on them.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from rapidfuzz.distance import Levenshtein

sys.path.insert(0, str(Path(__file__).parent))
from vdata import ROOT, paragraph_text, tokens

NAMES = {"H": "plants", "P": "jars and plant parts", "B": "bathing figures", "S": "stars / recipes",
         "T": "text only"}
SAMPLE, DRAWS = 250, 400


def sample(pages, rng, exclude=()):
    """Whole pages in random order until SAMPLE words; returns words and the pages used."""
    order = [p for p in rng.permutation(list(pages)) if p not in exclude]
    words, used = [], []
    for p in order:
        words += pages[p]; used.append(p)
        if len(words) >= SAMPLE:
            return words[:SAMPLE], used
    return None, used


def overlap(a, b):
    ca, cb = Counter(a), Counter(b)
    return sum(min(ca[w], cb[w]) for w in ca) / SAMPLE


def compare(A, B, rng, same):
    vals = []
    for _ in range(DRAWS):
        a, used = sample(A, rng)
        b, _ = sample(B, rng, exclude=used if same else ())
        if a is not None and b is not None:
            vals.append(overlap(a, b))
    return (float(np.mean(vals)), float(np.std(vals) / np.sqrt(len(vals))), len(vals)) if len(vals) > 20 else None


def main():
    rng = np.random.default_rng(0)
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    groups = {}
    for (hand, ill), g in pt.groupby(["hand", "illust"]):
        if hand in "1235" and ill in NAMES and len(g) >= 500:
            groups[f"scribe {hand}, {NAMES[ill]}"] = dict(
                hand=hand, ill=ill, words=len(g), lang=g.lang.mode().iloc[0],
                pages={f: list(h.word) for f, h in g.groupby("folio")})
    keys = list(groups)
    res = {}
    for i, a in enumerate(keys):
        for b in keys[i:]:
            r = compare(groups[a]["pages"], groups[b]["pages"], rng, same=(a == b))
            if r:
                ga, gb = groups[a], groups[b]
                kind = ("same scribe, same pictures" if a == b else
                        "same scribe, different pictures" if ga["hand"] == gb["hand"] else
                        "different scribe, same pictures" if ga["ill"] == gb["ill"] else
                        "different scribe, different pictures")
                res[f"{a}  vs  {b}"] = dict(overlap=r[0], se=r[1], kind=kind,
                                            same_language=ga["lang"] == gb["lang"])
    summary = {}
    for kind in ("same scribe, same pictures", "same scribe, different pictures",
                 "different scribe, same pictures", "different scribe, different pictures"):
        for lang in (True, False):
            v = [r["overlap"] for r in res.values() if r["kind"] == kind and r["same_language"] == lang]
            if v:
                summary[f"{kind}{'' if lang else ' (different Currier language)'}"] = dict(
                    mean=float(np.mean(v)), n=len(v), min=float(min(v)), max=float(max(v)))

    # how each group writes: neighbour repeats and word length
    style = {}
    for k, g in groups.items():
        sub = pt[(pt.hand == g["hand"]) & (pt.illust == g["ill"])]
        pairs = [(a, b) for _, ln in sub.groupby(["folio", "line"]) for a, b in zip(ln.word[:-1], ln.word[1:])]
        style[k] = dict(words=g["words"], language=g["lang"],
                        identical=float(np.mean([a == b for a, b in pairs])),
                        near=float(np.mean([a != b and Levenshtein.distance(a, b, score_cutoff=1) <= 1 for a, b in pairs])),
                        mean_len=float(sub.word.str.len().mean()),
                        top10=float(sub.word.value_counts().head(10).sum() / len(sub)))
    json.dump(dict(pairs=res, summary=summary, style=style),
              open(ROOT / "results" / "exp11_sections.json", "w"), indent=1)
    for k, r in sorted(res.items(), key=lambda kv: -kv[1]["overlap"]):
        print("%.3f ± %.3f  %-38s %s%s" % (r["overlap"], r["se"], r["kind"], k, "" if r["same_language"] else "  [A vs B]"))
    print()
    for k, v in summary.items():
        print("%-70s mean %.3f  (%d comparisons, %.3f to %.3f)" % (k, v["mean"], v["n"], v["min"], v["max"]))
    print()
    for k, s in style.items():
        print("%-34s words %5d lang %s  identical %.1f%%  near %.1f%%  mean length %.2f  top-10 share %.2f" % (
            k, s["words"], s["language"], 100 * s["identical"], 100 * s["near"], s["mean_len"], s["top10"]))


if __name__ == "__main__":
    main()


def key_contrast(boots=300, draws=40):
    """Page-bootstrap interval for the decisive comparison inside Currier language B:
    scribe 2's plant pages against (a) plant pages by scribes 3 and 5, (b) scribe 2's own
    bathing and text-only pages."""
    rng = np.random.default_rng(1)
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    pages = lambda m: {f: list(h.word) for f, h in pt[m].groupby("folio")}
    A = pages((pt.hand == "2") & (pt.illust == "H"))
    B = pages(pt.hand.isin(["3", "5"]) & (pt.illust == "H") & (pt.lang == "B"))
    C = pages((pt.hand == "2") & pt.illust.isin(["B", "T"]))

    def ov(X, Y):
        v = []
        for _ in range(draws):
            a, _u = sample(X, rng); b, _u = sample(Y, rng)
            if a is not None and b is not None:
                v.append(overlap(a, b))
        return np.mean(v)

    def resample(X):
        keys = list(X); pick = rng.integers(0, len(keys), len(keys))
        return {f"{keys[i]}#{j}": X[keys[i]] for j, i in enumerate(pick)}
    same_pic, same_scribe = ov(A, B), ov(A, C)
    d = []
    for _ in range(boots):
        a, b, c = resample(A), resample(B), resample(C)
        d.append(ov(a, b) - ov(a, c))
    out = dict(pages=dict(scribe2_plants=len(A), other_scribes_plants=len(B), scribe2_other=len(C)),
               same_pictures_other_scribes=float(same_pic), same_scribe_other_pictures=float(same_scribe),
               difference=float(same_pic - same_scribe),
               ci=[float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
               share_positive=float(np.mean(np.array(d) > 0)))
    p = ROOT / "results" / "exp11_sections.json"
    r = json.load(open(p)); r["key_contrast"] = out; json.dump(r, open(p, "w"), indent=1)
    print(out)


if __name__ == "__main__":
    key_contrast()
