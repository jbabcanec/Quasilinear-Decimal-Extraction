"""Experiment 15 - is "the same word twice in a row" really about being in a row?
Writes results/exp15_repeats.json.

For each text: how often a word is identical (or within one letter) to the next word,
divided by how often that happens when the words of every line are shuffled inside the
line.  1.0 = neighbours repeat exactly as often as any two words of the line would.
Real language should be far below 1 for identical words (nobody writes "the the").
Also: the same rates for words 2, 3, 4 apart inside a line.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
from rapidfuzz.distance import Levenshtein

sys.path.insert(0, str(Path(__file__).parent))
import exp8_genres as G8
import exp9_indic as G9
from exp6_llm_context import load
from vdata import ROOT


def counts(lines, k=1):
    ident = near = n = 0
    for ln in lines:
        for a, b in zip(ln[:-k], ln[k:]):
            n += 1
            if a == b:
                ident += 1
            elif Levenshtein.distance(a, b, score_cutoff=1) <= 1:
                near += 1
    return ident, near, n


def analyse(lines, rng, shuffles=60):
    i, nr, n = counts(lines)
    sh = np.array([counts([list(rng.permutation(l)) for l in lines])[:2] for _ in range(shuffles)])
    ei, en = sh[:, 0].mean(), sh[:, 1].mean()
    out = dict(pairs=n, identical=i / n, identical_shuffled=ei / n, near=nr / n, near_shuffled=en / n,
               ratio_identical=float(i / ei) if ei > 0 else None,
               ratio_near=float(nr / en) if en > 0 else None,
               z_identical=float((i - ei) / max(sh[:, 0].std(), 1e-9)),
               z_near=float((nr - en) / max(sh[:, 1].std(), 1e-9)))
    out["by_distance"] = [dict(zip(("identical", "near"), (c[0] / c[2], c[1] / c[2])))
                          for c in (counts(lines, k) for k in range(1, 6))]
    return out


def main():
    rng = np.random.default_rng(0)
    out = {}
    corp = {"Voynich": [ln for p in load("voynich") for ln in p],
            "Timm-Schinner generator": [ln for p in load("timm") for ln in p],
            "Naibbe cipher of Latin": [ln for p in load("naibbe") for ln in p],
            "Six real French manuscripts": [ln for p in load("french") for ln in p]}
    for label, (kind, globs) in G8.TEXTS.items():
        w = G8.read_latin(globs)[:60000]
        corp[f"Latin: {label}"] = [w[i:i + 9] for i in range(0, len(w), 9)]
    for label, (kind, globs) in G9.TEXTS.items():
        w = G9.read(globs)[:60000]
        if len(w) > 5000:
            corp[f"Sanskrit: {label}"] = [w[i:i + 9] for i in range(0, len(w), 9)]
    raw = ""
    for g in ("4_rellit/buddh/*dhAraNI*.txt", "4_rellit/buddh/sa_mahAmantrAnusAriNI.txt",
              "4_rellit/buddh/sa_mahAmAyUrIvidyArAjJI.txt"):
        for f in sorted(G9.G.glob(g)):
            raw += "\n" + f.read_text(errors="replace")
    w = [x for s in re.findall(r"tadyath[āa]\b(.{20,2500}?)sv[āa]h[āa]", raw, flags=re.S | re.I)
         for x in G9.words_of(s)]
    corp["Sanskrit: spell strings"] = [w[i:i + 9] for i in range(0, len(w), 9)]
    for name, lines in corp.items():
        r = analyse(lines, rng)
        out[name] = r
        print("%-40s identical %.2f%% vs %.2f%% shuffled (x%s, z %+.1f) | near %.2f%% vs %.2f%% (x%.2f, z %+.1f)" % (
            name, 100 * r["identical"], 100 * r["identical_shuffled"],
            "%.2f" % r["ratio_identical"] if r["ratio_identical"] is not None else "n/a", r["z_identical"],
            100 * r["near"], 100 * r["near_shuffled"], r["ratio_near"] or 0, r["z_near"]), flush=True)
    json.dump(out, open(ROOT / "results" / "exp15_repeats.json", "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
