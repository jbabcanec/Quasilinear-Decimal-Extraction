"""Experiment 12 - do the words fall into two families that take turns?
Writes results/exp12_alternation.json.

In writing where each unit is a letter or syllable, units split into two families (vowels
and consonants) that alternate.  Here the units are whole words.  The split is found on
half the lines (spectral bipartition of the neighbour graph, most anti-clustered direction)
and scored on the other half:

    alternation = share of neighbouring pairs that cross the split  /  share expected if
                  the two families were mixed at random

1.0 = no alternation.  Letters of Latin score far above 1.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from exp6_llm_context import load
from vdata import ROOT

K = 400


def split_score(lines, rng, reps=10):
    vals = []
    for _ in range(reps):
        order = rng.permutation(len(lines))
        tr = [lines[i] for i in order[: len(lines) // 2]]
        te = [lines[i] for i in order[len(lines) // 2:]]
        cnt = Counter(w for ln in tr for w in ln)
        types = [w for w, _ in cnt.most_common(K)]
        ix = {w: i for i, w in enumerate(types)}
        A = np.zeros((len(types), len(types)))
        for ln in tr:
            for a, b in zip(ln[:-1], ln[1:]):
                if a in ix and b in ix and a != b:
                    A[ix[a], ix[b]] += 1; A[ix[b], ix[a]] += 1
        d = A.sum(1); m = d.sum()
        if m == 0:
            continue
        Dm = np.where(d > 0, 1 / np.sqrt(np.maximum(d, 1e-9)), 0)
        B = Dm[:, None] * (A - np.outer(d, d) / m) * Dm[None, :]
        vals_, vecs = np.linalg.eigh(B)
        side = vecs[:, 0] > 0                      # most negative eigenvalue: anti-clustering
        cls = {w: bool(side[i]) for w, i in ix.items()}
        cross = same = 0
        seq = []
        for ln in te:
            for a, b in zip(ln[:-1], ln[1:]):
                if a in cls and b in cls and a != b:
                    cross += cls[a] != cls[b]; same += cls[a] == cls[b]
                    seq.append(cls[a]); seq.append(cls[b])
        p = np.mean(seq)
        vals.append((cross / (cross + same)) / (2 * p * (1 - p)))
    return float(np.mean(vals)), float(np.std(vals))


def main():
    rng = np.random.default_rng(0)
    out = {}
    corp = {"Voynich": [ln for p in load("voynich") for ln in p],
            "Six real French manuscripts": [ln for p in load("french") for ln in p],
            "Latin (Pliny), words": [ln for p in load("latin") for ln in p],
            "Naibbe cipher of Latin": [ln for p in load("naibbe") for ln in p],
            "Timm-Schinner generator": [ln for p in load("timm") for ln in p]}
    letters = [c for ln in corp["Latin (Pliny), words"] for w in ln for c in w]
    corp["Latin, letter by letter"] = [letters[i:i + 9] for i in range(0, 34000, 9)]
    # Latin written as syllable-like pairs of letters: a cipher unit of two letters
    corp["Latin, two letters at a time"] = [[letters[j] + letters[j + 1] for j in range(i, i + 18, 2)]
                                           for i in range(0, 68000, 18)]
    sh = [list(rng.permutation(ln)) for ln in corp["Voynich"]]
    corp["Voynich, words shuffled within lines"] = sh
    for name, lines in corp.items():
        mu, sd = split_score(lines, rng)
        out[name] = dict(alternation=mu, sd=sd)
        print("%-40s alternation %.3f ± %.3f" % (name, mu, sd), flush=True)
    json.dump(out, open(ROOT / "results" / "exp12_alternation.json", "w"), indent=1)


if __name__ == "__main__":
    main()
