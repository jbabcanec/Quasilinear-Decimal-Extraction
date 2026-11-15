"""Experiment 13 - consecutive lines avoid starting the same way.
Writes results/exp13_line_starts.json.

For each paragraph (first line excluded, since paragraph openings have their own gallows
rule), compare how often two lines begin with the same glyph against the same lines
shuffled within the paragraph.  Lag 1 = the line directly below, lag 2 = two lines below.

Measured twice for the Voynich, on two independent sources of line structure and text:
the ZL3b transliteration alone, and the word boxes (Takahashi's transliteration).
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from exp6_llm_context import load
from vdata import ROOT, load_zl, paragraph_text, tokens


def zl_paragraphs():
    """Paragraphs from ZL3b only: lists of lines (word lists), with page metadata."""
    pages, loci = load_zl()
    out, cur, key = [], [], None
    for L in loci:
        if not L["ltype"].startswith("P") or L["ltype"] == "Pt" or not L["words"]:
            if cur: out.append((key, cur))
            cur, key = [], None
            continue
        new = L["pstart"] or L["pos"] in "@*" or key is None or key[0] != L["folio"]
        if new:
            if cur: out.append((key, cur))
            meta = pages.get(L["folio"], {})
            cur, key = [], (L["folio"], meta.get("L", "?"), meta.get("H", "?"), meta.get("I", "?"))
        cur.append(L["words"])
        if L["pend"]:
            out.append((key, cur)); cur, key = [], None
    if cur: out.append((key, cur))
    return out


def box_paragraphs():
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    F = pt[pt.k == 0]
    out = []
    for (f, p), g in F.groupby(["folio", "par"], sort=False):
        run, prev = [], None
        for ln, w in zip(g.line, g.word):
            if prev is not None and ln != prev + 1:
                out.append(((f,), run)); run = []
            run.append([w]); prev = ln
        out.append(((f,), run))
    return out


def test(seqs, key=lambda w: w[0], lags=(1, 2, 3), shuffles=1000, seed=0):
    """seqs: list of sequences of line-opening words (paragraph-first line already removed)."""
    rng = np.random.default_rng(seed)
    seqs = [[key(w) for w in s] for s in seqs if len(s) >= 2]
    res = {}
    for lag in lags:
        n = sum(max(0, len(s) - lag) for s in seqs)
        if n < 50:
            continue
        obs = sum(a == b for s in seqs for a, b in zip(s[:-lag], s[lag:])) / n
        null = np.array([sum(a == b for s in seqs for t in [list(rng.permutation(s))]
                             for a, b in zip(t[:-lag], t[lag:])) for _ in range(shuffles)]) / n
        res[f"lag{lag}"] = dict(n=n, same=obs, shuffled=float(null.mean()), sd=float(null.std()),
                                z=float((obs - null.mean()) / null.std()))
    return res


def fmt(r):
    return " | ".join("%s: %.3f vs %.3f (z %+.1f, n=%d)" % (k, v["same"], v["shuffled"], v["z"], v["n"])
                      for k, v in r.items())


def main():
    out = {}
    Z = zl_paragraphs()
    firsts = lambda paras: [[ln[0] for ln in lines[1:]] for _, lines in paras]
    strip = lambda w: w[1:] if len(w) > 3 and w[0] in "yds" else w
    out["Voynich (ZL3b)"] = test(firsts(Z))
    out["Voynich (word boxes, Takahashi)"] = test([[ln[0] for ln in lines[1:]] for _, lines in box_paragraphs()])
    out["Voynich (ZL3b), first two glyphs"] = test(firsts(Z), key=lambda w: w[:2], lags=(1, 2))
    out["Voynich (ZL3b), y/d/s prefix removed"] = test(firsts(Z), key=lambda w: strip(w)[0], lags=(1, 2))
    out["Voynich (ZL3b), second word of the line"] = test(
        [[ln[1] for ln in lines[1:] if len(ln) > 1] for _, lines in Z], lags=(1, 2))
    out["Voynich (ZL3b), last word of the line, its last glyph"] = test(
        [[ln[-1] for ln in lines[1:]] for _, lines in Z], key=lambda w: w[-1], lags=(1, 2))
    for lang in "AB":
        out[f"Voynich, Currier {lang}"] = test(firsts([p for p in Z if p[0][1] == lang]), lags=(1, 2))
    for hand in "1235":
        out[f"Voynich, scribe {hand}"] = test(firsts([p for p in Z if p[0][2] == hand]), lags=(1, 2))
    names = {"H": "plants", "S": "stars / recipes", "B": "bathing figures", "P": "jars"}
    for ill, nm in names.items():
        out[f"Voynich, {nm} pages"] = test(firsts([p for p in Z if p[0][3] == ill]), lags=(1, 2))
    for key, label in (("french", "Six real French manuscripts"), ("timm", "Timm-Schinner generator"),
                       ("latin", "Latin (Pliny), wrapped"), ("naibbe", "Naibbe cipher of Latin, wrapped")):
        out[label] = test([[ln[0] for ln in page] for page in load(key)])
    # which openings follow which (ZL3b)
    c = Counter((a[0], b[0]) for s in firsts(Z) for a, b in zip(s[:-1], s[1:]))
    tot = sum(c.values()); pa = Counter(); pb = Counter()
    for (a, b), n in c.items():
        pa[a] += n; pb[b] += n
    main_g = [g for g, _ in pa.most_common(8)]
    out["transition_lift"] = {a: {b: c[(a, b)] * tot / (pa[a] * pb[b]) for b in main_g} for a in main_g}
    out["opening_share"] = {g: pa[g] / tot for g in main_g}
    json.dump(out, open(ROOT / "results" / "exp13_line_starts.json", "w"), indent=1)
    for k, r in out.items():
        if k not in ("transition_lift", "opening_share"):
            print("%-52s %s" % (k, fmt(r)))
    print("next line (column) given this line (row), relative to chance:")
    print("     " + "     ".join(main_g))
    for a in main_g:
        print("  %s  " % a + "  ".join("%4.2f" % out["transition_lift"][a][b] for b in main_g))


if __name__ == "__main__":
    main()


def extensions():
    """Where the rule applies: paragraph openers, restarts after a drawing, and lines that
    open with an added y/d/s letter versus lines that open with an ordinary word.
    Every null re-orders whole lines (or relabels positions), so nothing is ever compared
    with itself."""
    from collections import defaultdict
    rng = np.random.default_rng(0)
    Z = zl_paragraphs()
    out = {}
    by = defaultdict(list)
    for key, lines in Z:
        by[key[0]].append(lines[0][0])
    out["paragraph openers on a page"] = test([v for v in by.values() if len(v) >= 2], lags=(1,), shuffles=2000)["lag1"]

    seen = Counter(w for _, lines in Z for ln in lines for w in ln[1:])
    pref = lambda w: len(w) > 2 and w[0] in "yds" and seen.get(w[1:], 0) >= 3
    seqs = [[ln[0] for ln in lines[1:]] for _, lines in Z if len(lines) >= 3]

    def rate(order_fn, keep, left, right):
        hit = n = 0
        for sq in seqs:
            t = order_fn(sq)
            for prev, cur in zip(t[:-1], t[1:]):
                if keep(cur):
                    n += 1; hit += left(cur) == right(prev)
        return hit / n, n

    def versus(keep, left, right=lambda w: w[0], shuffles=1000):
        obs, n = rate(lambda sq: sq, keep, left, right)
        null = [rate(lambda sq: list(rng.permutation(sq)), keep, left, right)[0] for _ in range(shuffles)]
        return dict(n=n, same=obs, shuffled=float(np.mean(null)), z=float((obs - np.mean(null)) / np.std(null)))
    out["lines opening with an ordinary word"] = versus(lambda w: not pref(w), lambda w: w[0])
    out["lines opening with an added y/d/s: the added letter"] = versus(pref, lambda w: w[0])
    out["lines opening with an added y/d/s: the word underneath"] = versus(pref, lambda w: w[1])

    # restarts after a drawing, stacked one above the other (word boxes)
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    g = pt.groupby(["folio", "line"])
    pt["gap_before"] = (pt.x0 - g.x1.shift(1)) / pt.wmed
    S = pt[(pt.gap_before > 0.6) & (pt.k > 0)].reset_index(drop=True)
    pairs = []
    for f, h in S.groupby("folio"):
        L = h[["line", "x0", "wmed"]].assign(i=h.index).to_dict("records")
        for a in L:
            for b in L:
                if b["line"] == a["line"] + 1 and abs(b["x0"] - a["x0"]) < 1.5 * a["wmed"]:
                    pairs.append((a["i"], b["i"]))
    ia = np.array([p[0] for p in pairs]); ib = np.array([p[1] for p in pairs])
    lab = S.word.str[0].to_numpy(); fo = S.folio.to_numpy()
    obs = float((lab[ia] == lab[ib]).mean())
    null = []
    for _ in range(2000):
        l2 = lab.copy()
        for f in set(fo[ia]):
            m = np.flatnonzero(fo == f); l2[m] = rng.permutation(lab[m])
        null.append((l2[ia] == l2[ib]).mean())
    out["restarts after a drawing, stacked"] = dict(n=len(pairs), same=obs, shuffled=float(np.mean(null)),
                                                    z=float((obs - np.mean(null)) / np.std(null)))
    p = ROOT / "results" / "exp13_line_starts.json"
    d = json.load(open(p)); d["extensions"] = out; json.dump(d, open(p, "w"), indent=1)
    for k, r in out.items():
        print("%-56s n=%4d  same %.3f vs shuffled %.3f  (z %+.1f)" % (k, r["n"], r["same"], r["shuffled"], r["z"]))


if __name__ == "__main__":
    extensions()
