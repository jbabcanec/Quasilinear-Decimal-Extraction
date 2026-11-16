"""Experiment 17 - whole pages: do paragraphs spill over, and how full are the pages?
Writes results/exp17_pages.json.

1. ZL3b: does the last text line of every page end a paragraph?  If text simply ran on, a
   page break would fall on a paragraph end about once in (lines per paragraph).
2. Word boxes, text-heavy pages: blank left at the foot of the page, in lines, against what
   packing the same whole paragraphs in shuffled order would leave.  A scribe who refuses
   to split a paragraph but has no control over its length must leave part of a
   paragraph's height blank at the foot of most pages.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from vdata import ROOT, folio_sort_key, load_zl, paragraph_text, tokens

NAMES = {"H": "plants", "S": "stars / recipes", "B": "bathing", "T": "text only", "P": "jars"}


def main():
    rng = np.random.default_rng(0)
    out = {}
    pages, loci = load_zl()
    info, par_lines, cur = {}, [], 0
    for L in loci:
        if L["ltype"].startswith("P") and L["words"]:
            d = info.setdefault(L["folio"], dict(first=L, last=L)); d["last"] = L
            cur += 1
            if L["pend"]:
                par_lines.append((pages[L["folio"]].get("I"), cur)); cur = 0
    Z = pd.DataFrame([dict(sec=pages[f].get("I"), ends=bool(d["last"]["pend"])) for f, d in info.items()])
    PL = pd.DataFrame(par_lines, columns=["sec", "n"])
    sp = {}
    for code, nm in NAMES.items():
        z = Z[Z.sec == code]; pl = PL[PL.sec == code]
        if len(z):
            sp[nm] = dict(pages=int(len(z)), pages_ending_on_a_paragraph_end=int(z.ends.sum()),
                          expected_share_if_text_ran_on=float(1 / pl.n.mean()))
    out["page_breaks"] = sp
    out["page_breaks_all"] = dict(pages=int(len(Z)), ending_on_a_paragraph_end=int(Z.ends.sum()))

    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    g = pt.groupby(["folio", "line"]).agg(y=("yc", "mean"), par=("par", "first"), sec=("illust", "first")).reset_index()

    def analyse(secs):
        pg = []
        for f, h in g[g.sec.isin(secs)].groupby("folio"):
            h = h.sort_values("y")
            if len(h) < 20:
                continue
            pitch = float(np.median(np.diff(h.y)))
            top, bottom = h.y.iloc[0] - pitch / 2, h.y.iloc[-1] + pitch / 2
            starts = np.append(np.sort(h.groupby("par").y.min().to_numpy()) - pitch / 2, bottom)
            pg.append(dict(key=folio_sort_key(f), top=top, bottom=bottom, pitch=pitch, heights=np.diff(starts)))
        pg.sort(key=lambda p: p["key"])
        ref = np.quantile([p["bottom"] for p in pg], 0.9)
        blank = np.array([(ref - p["bottom"]) / p["pitch"] for p in pg])
        cap = [ref - p["top"] for p in pg]; pitch = float(np.mean([p["pitch"] for p in pg]))
        allh = np.concatenate([p["heights"] for p in pg])

        def pack(order):
            res, i = [], 0
            for c in cap:
                used = 0.0
                while i < len(order) and used + order[i] <= c + 0.3 * pitch:
                    used += order[i]; i += 1
                if i >= len(order):
                    break
                res.append((c - used) / pitch)
            return np.array(res)
        sims = [pack(rng.permutation(allh)) for _ in range(4000)]
        tight = np.array([np.mean(s < 0.5) for s in sims]); med = np.array([np.median(s) for s in sims])
        return dict(pages=len(pg), paragraphs=int(len(allh)), lines_per_paragraph=float(np.mean(allh) / pitch),
                    blank_median=float(np.median(blank)), under_half_a_line=float(np.mean(blank < 0.5)),
                    packed_blank_median=float(med.mean()), packed_under_half_a_line=float(tight.mean()),
                    packed_under_half_a_line_sd=float(tight.std()),
                    p=float(np.mean(tight >= np.mean(blank < 0.5))))
    out["page_fill"] = {"stars / recipes": analyse(["S"]), "bathing and text-only": analyse(["B", "T"]),
                        "all text-heavy pages": analyse(["S", "B", "T"])}
    json.dump(out, open(ROOT / "results" / "exp17_pages.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
