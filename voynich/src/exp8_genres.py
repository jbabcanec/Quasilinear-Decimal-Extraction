"""Experiment 8 - which kind of real text does the Voynich most resemble?
Writes results/exp8_genres.json.

Each text is cut into windows of 7,500 words (lines of 9 words; the Voynich keeps its own
lines) and described by eight numbers that do not depend on the alphabet:

  repeat       neighbouring words identical
  near         neighbouring words within one letter of each other, not identical
  order        word-order strength (see exp7)
  phrases      repeated word triples, relative to the same lines shuffled
  vocab        distinct words / words
  top10        share of text taken by the ten commonest words
  recur        words already seen in the previous 30 words
  lencv        spread of word length (sd / mean)

The numbers are standardised across the real texts; the distance between the Voynich and a
text is the Euclidean distance between their standardised profiles.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from rapidfuzz.distance import Levenshtein

sys.path.insert(0, str(Path(__file__).parent))
from exp6_llm_context import load
from exp7_absences import chance_onesided
from vdata import EXT, ROOT

LL = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/user/cltk/lat_text_latin_library")
N = 7500
TEXTS = {   # label: (kind, [files])
    "Psalms": ("prayer / song", ["bible/psalms.txt"]),
    "Leviticus": ("ritual law", ["bible/leviticus.txt"]),
    "Numbers": ("census and lists", ["bible/numbers.txt"]),
    "Chronicles I": ("genealogy lists", ["bible/chronicles1.txt"]),
    "Matthew": ("gospel narrative", ["bible/matthew.txt"]),
    "Rule of Benedict": ("monastic rule", ["benedict.txt"]),
    "Apicius": ("cooking recipes", ["apicius/apicius*.txt"]),
    "Cato, On Farming": ("farm handbook, recipes, charms", ["cato/cato.agri.txt"]),
    "Notitia Dignitatum": ("official list", ["notitia1.txt", "notitia2.txt"]),
    "Hyginus": ("star myths", ["hyginus/hyginus*.txt"]),
    "Manilius": ("astrology in verse", ["manilius*.txt"]),
    "Vergil": ("epic verse", ["vergil/aen*.txt"]),
    "Caesar": ("war narrative", ["caesar/gall*.txt"]),
    "Pliny, Natural History": ("encyclopedia", ["pliny.nh*.txt"]),
    "Gaius, Institutes": ("law textbook", ["gaius*.txt"]),
    "Vitruvius": ("building manual", ["vitruvius*.txt"]),
    "Isidore": ("etymologies", ["isidore/*.txt"]),
    "Frankish Annals": ("year-by-year chronicle", ["annalesregnifrancorum.txt"]),
}
FEATS = ["repeat", "near", "order", "phrases", "vocab", "top10", "recur", "lencv"]


def read_latin(globs):
    words = []
    for g in globs:
        for f in sorted(LL.glob(g)):
            words += re.findall(r"[a-z]+", f.read_text(errors="replace").lower())
    return words


def windows(lines, n=N, max_windows=4):
    out, cur, count = [], [], 0
    for ln in lines:
        cur.append(ln); count += len(ln)
        if count >= n:
            out.append(cur); cur, count = [], 0
            if len(out) == max_windows:
                break
    return out


def features(lines, rng):
    stream = [w for ln in lines for w in ln]
    n = len(stream)
    pairs = [(a, b) for ln in lines for a, b in zip(ln[:-1], ln[1:])]
    rep = np.mean([a == b for a, b in pairs])
    near = np.mean([a != b and Levenshtein.distance(a, b, score_cutoff=1) <= 1 for a, b in pairs])
    c = Counter((a, b) for ln in lines for a, b in zip(ln[1:-2], ln[2:-1]) if a != b)
    obs = exp = tot = 0.0
    for (a, b), n_ab in c.items():
        n_ba = c.get((b, a), 0)
        if a > b and n_ba:
            continue
        m = n_ab + n_ba
        if m >= 4:
            obs += abs(n_ab - n_ba); exp += chance_onesided(m) * m; tot += m
    order = ((obs / tot) - (exp / tot)) / (1 - exp / tot) if tot else np.nan

    def trip(ls):
        t = Counter(tuple(ln[i:i + 3]) for ln in ls for i in range(len(ln) - 2))
        return sum(v for v in t.values() if v >= 2) / max(1, sum(t.values()))
    sh = np.mean([trip([list(rng.permutation(ln)) for ln in lines]) for _ in range(10)])
    phrases = (trip(lines) + 1e-4) / (sh + 1e-4)
    cnt = Counter(stream)
    last, recur = {}, 0
    for i, w in enumerate(stream):
        if w in last and i - last[w] <= 30:
            recur += 1
        last[w] = i
    L = np.array([len(w) for w in stream])
    return dict(repeat=float(rep), near=float(near), order=float(order), phrases=float(phrases),
                vocab=len(cnt) / n, top10=sum(v for _, v in cnt.most_common(10)) / n,
                recur=recur / n, lencv=float(L.std() / L.mean()))


def profile(lines, rng):
    ws = windows(lines)
    fs = [features(w, rng) for w in ws]
    return {k: float(np.nanmean([f[k] for f in fs])) for k in FEATS} | {"windows": len(ws)}


def main():
    rng = np.random.default_rng(0)
    real = {}
    for label, (kind, globs) in TEXTS.items():
        words = read_latin(globs)
        if len(words) < N:
            print("skipped (too short):", label, len(words)); continue
        lines = [words[i:i + 9] for i in range(0, len(words), 9)]
        real[label] = profile(lines, rng) | {"kind": kind, "words": len(words)}
    other = {}
    for key, label in (("voynich", "Voynich"), ("french", "Six real French manuscripts"),
                       ("naibbe", "Naibbe cipher of Latin"), ("timm", "Timm-Schinner generator")):
        other[label] = profile([ln for page in load(key) for ln in page], rng)
    M = np.array([[real[t][f] for f in FEATS] for t in real])
    mu, sd = M.mean(0), M.std(0)
    z = lambda p: (np.array([p[f] for f in FEATS]) - mu) / sd
    zv = z(other["Voynich"])
    dist = {t: float(np.linalg.norm(z(real[t]) - zv)) for t in real}
    for t in ("Six real French manuscripts", "Naibbe cipher of Latin", "Timm-Schinner generator"):
        dist[t] = float(np.linalg.norm(z(other[t]) - zv))
    out = dict(features=FEATS, real=real, other=other, voynich_z=dict(zip(FEATS, zv.tolist())),
               distance_to_voynich=dict(sorted(dist.items(), key=lambda kv: kv[1])),
               real_range={f: [float(M[:, i].min()), float(M[:, i].max())] for i, f in enumerate(FEATS)})
    (ROOT / "results" / "exp8_genres.json").write_text(json.dumps(out, indent=1))
    print("%-30s %-30s " % ("text", "kind") + " ".join("%8s" % f for f in FEATS) + "   dist")
    rows = [(t, real[t]["kind"], real[t]) for t in real] + [(t, "", other[t]) for t in other]
    for t, kind, p in sorted(rows, key=lambda r: dist.get(r[0], -1)):
        print("%-30s %-30s " % (t, kind) + " ".join("%8.3f" % p[f] for f in FEATS) + "   %5.2f" % dist.get(t, 0))
    print("Voynich z-scores:", {f: round(v, 1) for f, v in zip(FEATS, zv)})


if __name__ == "__main__":
    main()
