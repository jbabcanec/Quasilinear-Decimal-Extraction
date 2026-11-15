"""Experiment 7 - two things real text has: word order, and repeated phrases.
Writes results/exp7_absences.json.

For each corpus, using only pairs of neighbouring words inside a line (line-first and
line-last words excluded, so line-position effects cannot fake an order):

  order      For two words that occur next to each other in both orders or either, how
             one-sided is the order?  0 = a coin flip (no word order), 1 = always the same
             way round.  Corrected for the one-sidedness that small counts give by chance.
  phrases    How many neighbouring word pairs / triples recur, relative to the same lines
             with their words shuffled?  1 = no more repeated phrases than chance.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from math import comb
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from exp6_llm_context import load
from vdata import ROOT


def chance_onesided(n):
    """E|2X - n| / n for X ~ Binomial(n, 1/2)."""
    return sum(comb(n, x) * abs(2 * x - n) for x in range(n + 1)) / (2 ** n * n)


def order_strength(lines, min_n=6):
    c = Counter()
    for ln in lines:
        for a, b in zip(ln[1:-2], ln[2:-1]):
            if a != b:
                c[(a, b)] += 1
    obs = exp = tot = 0.0
    pairs = 0
    for (a, b), n_ab in c.items():
        if a < b or (b, a) not in c:
            n_ba = c.get((b, a), 0)
            if a > b and n_ba:
                continue
            n = n_ab + n_ba
            if n >= min_n:
                obs += abs(n_ab - n_ba); exp += chance_onesided(n) * n; tot += n; pairs += 1
    raw, ch = obs / tot, exp / tot
    return dict(pairs=pairs, tokens=int(tot), raw=raw, chance=ch, strength=(raw - ch) / (1 - ch))


def repeats(lines, k):
    c = Counter(tuple(ln[i:i + k]) for ln in lines for i in range(1, len(ln) - k))
    n = sum(c.values())
    return sum(v for v in c.values() if v >= 2) / max(n, 1)


def phrase_ratio(lines, k, rng, reps=20):
    obs = repeats(lines, k)
    sh = np.mean([repeats([[ln[0]] + list(rng.permutation(ln[1:-1])) + [ln[-1]] if len(ln) > 3 else ln
                           for ln in lines], k) for _ in range(reps)])
    return dict(observed=obs, shuffled=float(sh), ratio=obs / sh if sh > 0 else None)


def main():
    rng = np.random.default_rng(0)
    out = {}
    names = {"voynich": "Voynich", "french": "Six real French manuscripts", "latin": "Latin (Pliny)",
             "naibbe": "Naibbe cipher of Latin", "timm": "Timm-Schinner generator"}
    for key, label in names.items():
        lines = [ln for page in load(key) for ln in page]
        out[label] = dict(words=sum(len(l) for l in lines), order=order_strength(lines),
                          pairs=phrase_ratio(lines, 2, rng), triples=phrase_ratio(lines, 3, rng))
        r = out[label]
        print("%-28s words %6d | order strength %.2f (raw %.2f, chance %.2f, %d word pairs) | repeated pairs x%.2f  triples x%s" % (
            label, r["words"], r["order"]["strength"], r["order"]["raw"], r["order"]["chance"], r["order"]["pairs"],
            r["pairs"]["ratio"], ("%.2f" % r["triples"]["ratio"]) if r["triples"]["ratio"] else "n/a"), flush=True)
    (ROOT / "results" / "exp7_absences.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()


def voynich_by_gap():
    """Word-order strength for Voynich neighbours written close together vs far apart.
    Closely written 'words' are often pieces of one word, which have a fixed order."""
    from collections import Counter
    from vdata import paragraph_text, tokens
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"]).reset_index(drop=True)
    g = pt.groupby(["folio", "line"])
    pt["nxt"] = g.word.shift(-1)
    pt["gap"] = (g.x0.shift(-1) - pt.x1) / pt.wmed
    q = pt[(pt.k > 0) & (pt.kr > 1) & pt.nxt.notna() & (pt.word != pt.nxt) & (pt.gap < 0.6)]
    med = q.gap.median()
    out = {}
    for name, sub in (("written close together", q[q.gap < med]), ("written far apart", q[q.gap >= med])):
        c = Counter(zip(sub.word, sub.nxt))
        obs = exp = tot = 0.0; pairs = 0
        for (a, b), n_ab in c.items():
            if a < b or (b, a) not in c:
                n_ba = c.get((b, a), 0)
                if a > b and n_ba:
                    continue
                n = n_ab + n_ba
                if n >= 6:
                    obs += abs(n_ab - n_ba); exp += chance_onesided(n) * n; tot += n; pairs += 1
        raw, ch = obs / tot, exp / tot
        out[name] = dict(pairs=pairs, raw=raw, chance=ch, strength=(raw - ch) / (1 - ch))
        print("Voynich neighbours %-24s order strength %.2f (raw %.2f, chance %.2f, %d word pairs)" % (
            name, out[name]["strength"], raw, ch, pairs))
    p = ROOT / "results" / "exp7_absences.json"
    d = json.loads(p.read_text()); d["voynich_by_gap"] = out; p.write_text(json.dumps(d, indent=1))


if __name__ == "__main__":
    voynich_by_gap()
