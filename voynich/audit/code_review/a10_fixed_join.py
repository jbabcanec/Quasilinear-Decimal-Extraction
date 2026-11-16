"""Rebuild the token table with ONE change: a new geometric line also starts when the next
box begins left of the previous box's midpoint (original rule: x0 < prev.x0 - 3, which
cannot see a one-word line followed by a line starting at the same margin).
Then rerun the 'as written' part of exp16 unchanged."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
from collections import defaultdict
import numpy as np, pandas as pd
import vdata
from vdata import load_zl, load_boxes, _align_chars, _add_structure, folio_sort_key, BOX_DIR, paragraph_text
from wrap import line_table, margins

def build(fixed):
    pages, loci = load_zl()
    by_folio = defaultdict(list)
    for L in loci: by_folio[L["folio"]].append(L)
    rows = []
    for path in sorted(BOX_DIR.glob("*.js"), key=lambda p: folio_sort_key(p.stem)):
        folio = path.stem; boxes = load_boxes(folio)
        if not boxes: continue
        meta = pages.get(folio, {}); zl = by_folio.get(folio, [])
        locus, frac = _align_chars([b["word"] for b in boxes], [L["words"] for L in zl])
        gline, g = [], 0
        for i, b in enumerate(boxes):
            if i:
                p = boxes[i - 1]
                if b["x0"] < p["x0"] - 3 or (fixed and b["x0"] < p["x0"] + 0.5 * p["w"]):
                    g += 1
            gline.append(g)
        for i, b in enumerate(boxes):
            j = int(locus[i]); L = zl[j] if j >= 0 else None
            rows.append(dict(folio=folio, seq=i, word=b["word"], x0=b["x0"], x1=b["x0"] + b["w"], xc=b["x0"] + b["w"] / 2,
                y0=b["y0"], yc=b["y0"] + b["h"] / 2, w=b["w"], h=b["h"], gline=gline[i], zl_locus=j, zl_frac=frac[i],
                zl_num=L["num"] if L else "", ltype=L["ltype"] if L else "", zl_pos=L["pos"] if L else "",
                zl_pstart=bool(L and L["pstart"]), zl_pend=bool(L and L["pend"]),
                lang=meta.get("L", "?"), hand=meta.get("H", "?"), illust=meta.get("I", "?"), quire=meta.get("Q", "?")))
    return _add_structure(pd.DataFrame(rows))

EDGES = [0, .2, .4, .6, .8, 1.0, 9]
def exp16_as_written(df, keep_pt=False):
    m = df["is_par"] & df["line_agree"]
    if not keep_pt: m &= df["ltype"] != "Pt"
    pt = df[m].sort_values(["folio", "line", "k"])
    L = line_table(pt)
    paras = []
    for (f, p), g in L.groupby(["folio", "par"], sort=False):
        g = g.sort_values("line")
        if len(g) < 4: continue
        body = g.iloc[1:-1]; W = float(np.median((body.E - body.xl) / body.wmed)); last = g.iloc[-1]
        paras.append(dict(fill=float((last.E - last.xl) / last.wmed / W), nwords=int(last.n)))
    P = pd.DataFrame(paras)
    h, _ = np.histogram(P.fill, EDGES)
    M = margins(line_table(pt)); M = M[M.M1.notna()].copy(); M["fill"] = (M.E - M.xl) / (M.M1 - M.xl)
    S = M[(M.fill < 0.9) & (M.fill > 0)]
    return dict(lines=len(L), n=len(P), one_word=int((P.nwords == 1).sum()), one_word_share=round(float((P.nwords == 1).mean()), 4),
                up_to_two=round(float((P.nwords <= 2).mean()), 4), under_30=round(float((P.fill < .3).mean()), 4),
                mid_40_80=round(float(((P.fill >= .4) & (P.fill < .8)).mean()), 4), fifths=(h / len(P)).round(3).tolist(),
                stop_early_n=len(S), stop_early_one_word=round(float((S.n == 1).mean()), 4))

orig = build(False); fix = build(True)
fix.to_pickle("tokens_fixed.pkl")
print("original rule            :", exp16_as_written(orig))
print("fixed rule               :", exp16_as_written(fix))
print("fixed rule, Pt lines kept:", exp16_as_written(fix, keep_pt=True))
P = fix[fix.is_par]
lines = P.groupby(["folio", "line"]).agg(n=("k", "size"), agree=("line_agree", "first"))
print("fixed: P lines by word count -> share surviving line_agree")
print(lines.assign(nb=lines.n.clip(upper=8)).groupby("nb").agg(lines=("n", "size"), agree=("agree", "mean")).round(3).to_string())
