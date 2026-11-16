"""One-word lines: are they dropped by the geometric line rule (x must jump left by >3px)?"""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np, pandas as pd
from collections import Counter, defaultdict
from vdata import tokens, paragraph_text, load_zl
df = tokens()
pages, loci = load_zl()
byf = defaultdict(list)
for L in loci: byf[L["folio"]].append(L)
# ZL: P-type loci by word count (excluding Pt)
folios_with_boxes = set(df.folio)
zl_p = [L for L in loci if L["ltype"].startswith("P") and L["words"] and L["folio"] in folios_with_boxes]
print("ZL P-type loci on boxed folios:", len(zl_p), " Pt:", sum(L["ltype"] == "Pt" for L in zl_p))
c = Counter(min(len(L["words"]), 8) for L in zl_p if L["ltype"] != "Pt")
print("ZL P loci (not Pt) by ZL word count:", sorted(c.items()))
# which ZL loci are represented by a surviving line in paragraph_text?
pt = paragraph_text(df)
surv = set(zip(pt.folio, pt.zl_locus_f))
idx = {id(L): i for f, v in byf.items() for i, L in enumerate(v)}
tab = defaultdict(lambda: [0, 0])
miss1 = []
for L in zl_p:
    if L["ltype"] == "Pt": continue
    n = min(len(L["words"]), 8)
    ok = (L["folio"], idx[id(L)]) in surv
    tab[n][0] += 1; tab[n][1] += ok
    if n == 1 and not ok: miss1.append(L)
print("ZL words -> (loci, surviving in paragraph_text, share)")
for n in sorted(tab): print("  ", n, tab[n], round(tab[n][1] / tab[n][0], 3))
tot = sum(v[0] for v in tab.values()); s = sum(v[1] for v in tab.values())
print("overall survival", s, tot, s / tot)
# mechanism for the dropped one-word loci
reasons = Counter()
for L in miss1:
    f = L["folio"]; j = idx[id(L)]
    g = df[df.folio == f]
    t = g[g.zl_locus_f == j]
    if t.empty:
        reasons["no box aligned to it"] += 1; continue
    gl = t.gline.iloc[0]
    others = g[(g.gline == gl) & (g.zl_locus_f != j)]
    if len(others):
        nxt = others[others.seq > t.seq.max()]
        prv = others[others.seq < t.seq.min()]
        reasons["shares geometric line with NEXT locus (x did not jump left)" if len(nxt) and not len(prv) else
                "shares geometric line with PREVIOUS locus" if len(prv) and not len(nxt) else "both"] += 1
    else:
        reasons["other"] += 1
print("dropped one-word ZL loci:", len(miss1)); print(reasons)
# is it at paragraph end?
print("of the dropped one-word loci, paragraph-final (<$>):", sum(L["pend"] for L in miss1),
      " also paragraph-initial:", sum(L["pend"] and L["pstart"] for L in miss1))
# show a few
for L in [L for L in miss1 if L["pend"] and not L["pstart"]][:8]:
    f = L["folio"]; j = idx[id(L)]
    g = df[df.folio == f]; t = g[g.zl_locus_f == j]
    if t.empty: print(f, L["num"], L["words"], "no box"); continue
    s0 = t.seq.min()
    print(f, L["num"], L["ltype"], L["words"])
    print(g[(g.seq >= s0 - 2) & (g.seq <= s0 + 2)][["seq", "word", "x0", "y0", "gline", "zl_locus_f", "zl_num", "line", "line_agree"]].to_string(index=False))
