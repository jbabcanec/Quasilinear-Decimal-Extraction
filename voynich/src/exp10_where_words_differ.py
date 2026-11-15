"""Experiment 10 - when two neighbouring words differ by one letter, where is that letter?
Writes results/exp10_where.json.

Rhyming sound-play changes the start (ili mili tili), counting changes the end (xii xiii),
grammar changes the end (magna magnam).  Where does the Voynich change?
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from rapidfuzz.distance import Levenshtein

sys.path.insert(0, str(Path(__file__).parent))
import exp8_genres as G8
import exp9_indic as G9
from exp6_llm_context import load
from vdata import ROOT


def roman(n):
    out = ""
    for v, s in ((1000, "m"), (500, "d"), (100, "c"), (50, "l"), (10, "x"), (5, "v"), (1, "i")):
        while n >= v:
            out += s; n -= v
    return out


def where(lines):
    c = Counter(); kinds = Counter(); exp = np.zeros(3); n = 0
    for ln in lines:
        for a, b in zip(ln[:-1], ln[1:]):
            if a == b or min(len(a), len(b)) < 3:
                continue
            ops = Levenshtein.editops(a, b)
            if len(ops) != 1:
                continue
            tag, i, j = ops[0].tag, ops[0].src_pos, ops[0].dest_pos
            long = a if len(a) >= len(b) else b
            pos = i if tag != "insert" else j
            L = len(long)
            c["start" if pos == 0 else "end" if pos >= L - 1 else "middle"] += 1
            kinds["swap a letter" if tag == "replace" else "add or drop a letter"] += 1
            exp += np.array([1, L - 2, 1]) / L
            n += 1
    return dict(pairs=n, start=c["start"] / n, middle=c["middle"] / n, end=c["end"] / n,
                chance=dict(zip(("start", "middle", "end"), (exp / n).tolist())),
                swap=kinds["swap a letter"] / n)


def main():
    out = {}
    out["Voynich"] = where([ln for page in load("voynich") for ln in page])
    raw = ""
    for g in ("4_rellit/buddh/*dhAraNI*.txt", "4_rellit/buddh/sa_mahAmantrAnusAriNI.txt",
              "4_rellit/buddh/sa_mahAmAyUrIvidyArAjJI.txt"):
        for f in sorted(G9.G.glob(g)):
            raw += "\n" + f.read_text(errors="replace")
    import re
    w = [x for s in re.findall(r"tadyath[āa]\b(.{20,2500}?)sv[āa]h[āa]", raw, flags=re.S | re.I)
         for x in G9.words_of(s)]
    out["Spell strings (Sanskrit)"] = where([w[i:i + 9] for i in range(0, len(w), 9)])
    lat = []
    for label, (kind, globs) in G8.TEXTS.items():
        words = G8.read_latin(globs)
        lat += [words[i:i + 9] for i in range(0, len(words), 9)]
    out["Eighteen Latin texts"] = where(lat)
    for key, label in (("timm", "Timm-Schinner generator"), ("naibbe", "Naibbe cipher of Latin"),
                       ("french", "Six real French manuscripts")):
        out[label] = where([ln for page in load(key) for ln in page])
    nums = [roman(n) for n in range(1, 4000)]
    out["Counting in Roman numerals (generated)"] = where([nums[i:i + 9] for i in range(0, len(nums), 9)])
    (ROOT / "results" / "exp10_where.json").write_text(json.dumps(out, indent=1))
    print("%-40s %6s  %6s %6s %6s   (by chance: %s)   %s" % ("", "pairs", "start", "middle", "end", "start/middle/end", "letter swapped"))
    for k, r in out.items():
        ch = r["chance"]
        print("%-40s %6d  %5.0f%% %5.0f%% %5.0f%%   (%.0f%% / %.0f%% / %.0f%%)   %5.0f%%" % (
            k, r["pairs"], 100 * r["start"], 100 * r["middle"], 100 * r["end"],
            100 * ch["start"], 100 * ch["middle"], 100 * ch["end"], 100 * r["swap"]))


if __name__ == "__main__":
    main()
