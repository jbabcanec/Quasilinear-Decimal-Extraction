"""Experiment 5 - does the text run on across line breaks?
Writes results/exp5_continuity.json and figures/fig7_continuity.png.

The 'link' between two consecutive words is the mutual information (bits) between the last
glyph of the first and the first glyph of the second, minus what the same pairs give when
partners are shuffled within the page (so page dialect and sampling bias cancel).

Measured for the Voynich by the physical gap between the two words (tight ... wide, a
drawing in the way, a line break), for six real 15th-century French manuscripts (inside a
line / across a line break), and for two published models of the text.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "control_real"))
from vdata import EXT, ROOT, paragraph_text, tokens

FRENCH = ["btv1b100737746", "btv1b55008562q", "btv1b90069505", "btv1b9009670t",
          "btv1b90595325", "btv1b90616967"]


def _mi(a, b):
    n = len(a)
    c = Counter(zip(a, b)); ca = Counter(a); cb = Counter(b)
    return sum(v / n * np.log2(v * n / (ca[x] * cb[y])) for (x, y), v in c.items())


def link(left, right, page, rng, m=3000, reps=8, shuffles=25):
    """Excess mutual information over within-page shuffles, on subsamples of m pairs."""
    left, right, page = np.asarray(left), np.asarray(right), np.asarray(page)
    m = min(m, len(left))
    vals = []
    for _ in range(reps):
        s = rng.permutation(len(left))[:m]
        a, b, p = left[s], right[s], page[s]
        base = np.argsort(p, kind="stable")
        null = []
        for _ in range(shuffles):
            perm = np.lexsort((rng.random(m), p))
            bs = np.empty_like(b); bs[base] = b[perm]
            null.append(_mi(a, bs))
        vals.append(_mi(a, b) - np.mean(null))
    return dict(bits=float(np.mean(vals)), sd=float(np.std(vals)), n=int(len(left)))


def frame(df):
    """Arrays and the standard pair masks for a token table sorted in reading order."""
    d = df.sort_values(["folio", "line", "k"]).reset_index(drop=True)
    A = {c: d[c].to_numpy() for c in ("word", "line", "folio", "par", "k", "kr")}
    idx = np.arange(len(d)); i, j = idx[:-1], idx[1:]
    same = (A["folio"][i] == A["folio"][j]) & (A["line"][i] == A["line"][j])
    brk = ((A["folio"][i] == A["folio"][j]) & (A["line"][j] == A["line"][i] + 1)
           & (A["par"][i] == A["par"][j]) & (A["kr"][i] == 0) & (A["k"][j] == 0))
    return d, A, i, j, same, brk


def measure(A, i, j, mask, rng, fl=lambda a: a[-1], fr=lambda b: b[0], keep=None, m=3000):
    a, b = A["word"][i[mask]], A["word"][j[mask]]
    pg = A["folio"][i[mask]]
    if keep is not None:
        ok = np.array([keep(x, y) for x, y in zip(a, b)])
        a, b, pg = a[ok], b[ok], pg[ok]
    return link([fl(x) for x in a], [fr(y) for y in b], pg, rng, m=m)


def lines_table(lines, per_page=30):
    rows = []
    for n, words in enumerate(lines):
        for k, w in enumerate(words):
            rows.append((f"p{n // per_page}", n % per_page, n // per_page, w, k, len(words) - 1 - k))
    return pd.DataFrame(rows, columns=["folio", "line", "par", "word", "k", "kr"])


def main():
    rng = np.random.default_rng(0)
    out = {}

    # ------------------------------------------------------------------ Voynich
    pt = paragraph_text(tokens())
    d, A, i, j, same, brk = frame(pt)
    gap = (d.x0.to_numpy()[j] - d.x1.to_numpy()[i]) / d.wmed.to_numpy()[i]
    inner = same & (A["k"][i] > 0) & (A["kr"][j] > 0)
    q = np.quantile(gap[inner & (gap < 0.6)], [0, .2, .4, .6, .8, 1.0])
    v = {}
    names = ["tightest fifth of gaps", "second fifth", "middle fifth", "fourth fifth", "widest fifth of gaps"]
    for nm, lo, hi in zip(names, q[:-1], q[1:]):
        v[nm] = measure(A, i, j, inner & (gap >= lo) & (gap <= hi) & (gap < 0.6), rng) | {"gap": [float(lo), float(hi)]}
    v["a drawing in the way"] = measure(A, i, j, same & (gap > 0.6), rng)
    v["across a line break"] = measure(A, i, j, brk, rng)
    out["voynich_by_gap"] = v

    # is continuity merely hidden by line-start prefixes or line-end forms?
    seen = Counter(d.word[d.k > 0])
    strip = lambda w: w[1:] if len(w) > 2 and w[0] in "yds" and seen.get(w[1:], 0) >= 3 else w
    wide = inner & (gap >= 0.2) & (gap < 0.6)
    starts = A["word"][j[brk]]
    rep = {"share_starting_yds": float(np.mean([w[0] in "yds" for w in starts])),
           "remainder_known_at_line_start": float(np.mean([strip(w) != w for w in starts if w[0] in "yds"])),
           "remainder_known_inside_line": float(np.mean([strip(w) != w for w in d.word[d.k > 0] if w[0] in "yds"]))}
    tests = {
        "as written": dict(),
        "to the second glyph of the next word": dict(fr=lambda b: b[1] if len(b) > 1 else "_"),
        "after removing a y/d/s prefix": dict(fr=lambda b: strip(b)[0]),
        "only next words without y/d/s": dict(keep=lambda a, b: b[0] not in "yds"),
        "...and previous word not ending in m": dict(keep=lambda a, b: b[0] not in "yds" and a[-1] != "m"),
    }
    for nm, kw in tests.items():
        rep[nm] = {"wide gap inside a line": measure(A, i, j, wide, rng, m=1200, **kw),
                   "across a line break": measure(A, i, j, brk, rng, m=1200, **kw)}
    out["voynich_repairs"] = rep

    # ------------------------------------------------------------------ real manuscripts
    from control_real import token_table
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/user/Gallicorpora/HTR-MSS-15e-Siecle")
    fr = {}
    parts = []
    for ms in FRENCH:
        t = token_table(root / "data" / ms)
        t["word"] = t.word.map(lambda s: re.sub(r"[^\w]", "", s.lower()) or "_")
        t["folio"] = ms + "/" + t.folio
        t["par"] = ms + "/" + t.par.astype(str)
        parts.append(t)
        _, Af, fi, fj, fsame, fbrk = frame(t)
        fin = fsame & (Af["k"][fi] > 0) & (Af["kr"][fj] > 0)
        m = int(fbrk.sum())
        fr[ms] = {"inside a line": measure(Af, fi, fj, fin, rng, m=m),
                  "across a line break": measure(Af, fi, fj, fbrk, rng, m=m)}
    _, Af, fi, fj, fsame, fbrk = frame(pd.concat(parts, ignore_index=True))
    fin = fsame & (Af["k"][fi] > 0) & (Af["kr"][fj] > 0)
    fr["pooled"] = {"inside a line": measure(Af, fi, fj, fin, rng),
                    "across a line break": measure(Af, fi, fj, fbrk, rng)}
    out["french"] = fr

    # ------------------------------------------------------------------ models of the text
    mod = {}
    files = {"Timm-Schinner generator": EXT / "controls" / "selfcitation" / "selfcite_seed1.txt",
             "Naibbe cipher of Latin": EXT / "controls" / "naibbe" / "naibbe_caesar_ciphertext_lines.txt"}
    for nm, path in files.items():
        lines = [ln.split() for ln in path.read_text().splitlines() if len(ln.split()) >= 3]
        _, Am, mi_, mj, msame, mbrk = frame(lines_table(lines))
        mod[nm] = {"inside a line": measure(Am, mi_, mj, msame, rng),
                   "across a line break": measure(Am, mi_, mj, mbrk, rng)}
    out["models"] = mod

    (ROOT / "results" / "exp5_continuity.json").write_text(json.dumps(out, indent=1))
    for k, r in v.items():
        print("voynich  %-26s %.3f ± %.3f  (n=%d)" % (k, r["bits"], r["sd"], r["n"]))
    for k, r in rep.items():
        if isinstance(r, dict):
            print("repair   %-40s wide gap %.3f | line break %.3f" % (k, r["wide gap inside a line"]["bits"], r["across a line break"]["bits"]))
    for k, r in fr.items():
        print("french   %-16s inside %.3f ± %.3f | across %.3f ± %.3f" % (k, r["inside a line"]["bits"], r["inside a line"]["sd"], r["across a line break"]["bits"], r["across a line break"]["sd"]))
    for k, r in mod.items():
        print("model    %-26s inside %.3f | across %.3f" % (k, r["inside a line"]["bits"], r["across a line break"]["bits"]))


if __name__ == "__main__":
    main()
