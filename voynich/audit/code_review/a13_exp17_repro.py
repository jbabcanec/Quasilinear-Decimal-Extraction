"""Verbatim copy of exp17's analyse() with the shared rng, no file writes: does the JSON reproduce?"""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np, pandas as pd
from vdata import ROOT, folio_sort_key, load_zl, paragraph_text, tokens
rng = np.random.default_rng(0)
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
    return dict(pages=len(pg), under_half_a_line=float(np.mean(blank < 0.5)), packed_under_half_a_line=float(tight.mean()),
                sd=float(tight.std()), p=float(np.mean(tight >= np.mean(blank < 0.5))), blank_median=float(np.median(blank)), packed_blank_median=float(med.mean()))
res = {"stars / recipes": analyse(["S"]), "bathing and text-only": analyse(["B", "T"]), "all text-heavy pages": analyse(["S", "B", "T"])}
J = json.load(open(ROOT / "results" / "exp17_pages.json"))["page_fill"]
for k, v in res.items():
    print(k, {a: round(b, 4) for a, b in v.items()}, "| JSON p", J[k]["p"], "tight", round(J[k]["packed_under_half_a_line"], 4))
# number of sub-20-line pages in these sections
n = g.groupby("folio").agg(n=("line", "size"), sec=("sec", "first"))
print("S/B/T pages with paragraph text:", (n.sec.isin(list("SBT"))).sum(), " with <20 lines (excluded):", n[n.sec.isin(list("SBT")) & (n.n < 20)].to_dict("index"))
