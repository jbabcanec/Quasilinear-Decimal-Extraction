"""exp16: is the 'last line' used for each paragraph really the paragraph's last ZL line?"""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np, pandas as pd
from collections import Counter, defaultdict
from vdata import tokens, paragraph_text, load_zl
from wrap import line_table
df = tokens()
pages, loci = load_zl()
by_folio = defaultdict(list)
for L in loci: by_folio[L["folio"]].append(L)
pt = paragraph_text(df).sort_values(["folio", "line", "k"])
L = line_table(pt)
# majority zl locus per surviving line
zl_of = pt.groupby(["folio", "line"]).zl_locus_f.agg(lambda s: s.value_counts().index[0])
rows = []
for (f, p), g in L.groupby(["folio", "par"], sort=False):
    g = g.sort_values("line")
    if len(g) < 4: continue
    last = g.iloc[-1]
    j = int(zl_of.loc[(f, last.line)])
    zl = by_folio[f]
    # walk forward in ZL from locus j to the end of its paragraph (first pend at or after j, P-type loci only)
    k = j; trail = []
    while not zl[k]["pend"]:
        k += 1
        if k >= len(zl) or zl[k]["pstart"]: k = None; break
        trail.append(zl[k])
    # all df lines of this paragraph (before the paragraph_text filter)
    allp = df[(df.folio == f) & (df.par == p)].groupby("line").agg(n=("k", "size"), agree=("line_agree", "first"), lt=("ltype", "first"))
    dropped_after = allp[allp.index > last.line]
    rows.append(dict(folio=f, par=p, used_last_n=int(last.n), used_is_pend=bool(zl[j]["pend"]),
                     zl_last_words=(len(zl[k]["words"]) if k is not None else -1),
                     zl_last_type=(zl[k]["ltype"] if k is not None else "?"),
                     zl_last_num=(zl[k]["num"] if k is not None else "?"),
                     n_trailing_zl=len(trail), n_dropped_after=len(dropped_after),
                     dropped_after=";".join(f"{r.lt}:{r.n}:{'A' if r.agree else 'x'}" for r in dropped_after.itertuples())))
R = pd.DataFrame(rows)
print("paragraphs analysed (>=4 surviving lines):", len(R))
print("used last line carries <$> :", R.used_is_pend.sum(), " does not:", (~R.used_is_pend).sum())
bad = R[~R.used_is_pend]
print("true (ZL) last line word count for those where the used last line is NOT the real last:")
print(bad.zl_last_words.value_counts().sort_index().to_string())
print(bad.zl_last_type.value_counts().to_string())
print(bad.head(60).to_string())
R.to_csv("a03_lastline.csv", index=False)
# corrected one-word rate
true_one = int(((R.used_is_pend) & (R.used_last_n == 1)).sum() + ((~R.used_is_pend) & (R.zl_last_words == 1)).sum())
print("exp16 one-word count as coded:", int((R.used_last_n == 1).sum()), "of", len(R))
print("one-word true last lines among the same paragraphs (ZL count where dropped):", true_one, "of", len(R), "=", true_one / len(R))
true_le2 = int(((R.used_is_pend) & (R.used_last_n <= 2)).sum() + ((~R.used_is_pend) & (R.zl_last_words <= 2) & (R.zl_last_words > 0)).sum())
print("<=2 words: coded", int((R.used_last_n <= 2).sum()), " corrected", true_le2, true_le2 / len(R))
