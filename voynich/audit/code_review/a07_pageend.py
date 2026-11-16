"""Is '<$>' on the last line of a page evidence (the line looks like a paragraph end) or convention?"""
import sys
sys.dont_write_bytecode = True
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/audit/code_review")
import ivtff
pages, loci = ivtff.read()
P = ivtff.paragraphs(loci, include_pt=True, break_on_other=False)
# page-last P locus
last_on_page = {}
first_on_page = {}
for L in loci:
    if L["ltype"].startswith("P") and L["words"]:
        last_on_page[L["folio"]] = L["fileline"]; first_on_page.setdefault(L["folio"], L["fileline"])
notend = [f for f in last_on_page if not [L for L in loci if L["fileline"] == last_on_page[f]][0]["pend"]]
print("pages with P text:", len(last_on_page), " last P line lacks <$>:", notend)
for f in notend:
    L = [L for L in loci if L["fileline"] == last_on_page[f]][0]; print("   ", f, L["num"], L["pos"], L["ltype"], L["raw"][:80])
rows = []
for p in P:
    if len(p) < 4: continue
    n = [len(L["words"]) for L in p]
    ratio = n[-1] / np.median(n[1:-1])
    sec = pages[p[0]["folio"]].get("I")
    rows.append((sec, p[-1]["fileline"] == last_on_page[p[0]["folio"]], ratio, n[-1], p[0]["fileline"] == first_on_page[p[0]["folio"]], p[0]["words"][0][0]))
import pandas as pd
R = pd.DataFrame(rows, columns=["sec", "page_final", "ratio", "nlast", "page_first", "g0"])
for name, m in (("all sections", R.sec.notna()), ("S,B,T text-heavy", R.sec.isin(list("SBT"))), ("S", R.sec == "S"), ("H plants", R.sec == "H")):
    a, b = R[m & R.page_final], R[m & ~R.page_final]
    print("%-18s paragraphs(>=4 lines) ending at page foot: n=%3d mean last/body %.2f, share>=0.8: %.2f, <0.3: %.2f | ending mid-page: n=%3d mean %.2f, >=0.8: %.2f, <0.3: %.2f"
          % (name, len(a), a.ratio.mean(), (a.ratio >= .8).mean(), (a.ratio < .3).mean(), len(b), b.ratio.mean(), (b.ratio >= .8).mean(), (b.ratio < .3).mean()))
# gallows-initial paragraph starts: page-first vs others (all paragraphs)
rows = []
for p in P:
    f = p[0]["folio"]
    rows.append((pages[f].get("I"), p[0]["fileline"] == first_on_page[f], p[0]["words"][0][0] in "ptkf"))
G = pd.DataFrame(rows, columns=["sec", "page_first", "gallows"])
for name, m in (("all", G.sec.notna()), ("S,B,T", G.sec.isin(list("SBT")))):
    print(name, "first line of page opens with a gallows: %.2f (n=%d); other paragraph openings: %.2f (n=%d)"
          % (G[m & G.page_first].gallows.mean(), (m & G.page_first).sum(), G[m & ~G.page_first].gallows.mean(), (m & ~G.page_first).sum()))
# non-paragraph-initial lines opening with gallows, for reference
allp = [L for p in P for L in p[1:]]
print("lines that do not open a paragraph, opening with a gallows: %.3f (n=%d)" % (np.mean([L["words"][0][0] in "ptkf" for L in allp]), len(allp)))
