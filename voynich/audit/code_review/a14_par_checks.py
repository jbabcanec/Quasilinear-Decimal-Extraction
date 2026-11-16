import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/audit/code_review")
from collections import Counter
import numpy as np
from vdata import load_zl
import ivtff
pages, loci = load_zl()
P = [L for L in loci if L["ltype"].startswith("P")]
print("P loci:", len(P), " @/* without <%>:", sum(L["pos"] in "@*" and not L["pstart"] for L in P),
      " <%> on + or = loci:", sum(L["pstart"] and L["pos"] not in "@*" for L in P))
print("  examples @/* without <%>:", [(L["folio"], L["num"], L["pos"], L["ltype"]) for L in P if L["pos"] in "@*" and not L["pstart"]][:12])
# open/close balance in file order (P loci only)
open_ = False; bad = Counter()
for L in loci:
    if not L["ltype"].startswith("P"):
        if open_: bad["non-P locus inside an open paragraph"] += 1
        continue
    if L["pstart"]:
        if open_: bad["<%> while open"] += 1
        open_ = True
    elif not open_: bad["P line outside <%>..<$>"] += 1
    if L["pend"]: open_ = False
print(bad)
# single-word paragraph-final lines in ZL (all paragraph sizes), incl / excl Pt
pg2, loci2 = ivtff.read()
for inc in (True, False):
    Q = ivtff.paragraphs(loci2, include_pt=inc, break_on_other=False)
    n1 = sum(len(p[-1]["words"]) == 1 and len(p) > 1 for p in Q)
    n1h = sum(len(p[-1]["hard"]) == 1 and len(p) > 1 for p in Q)
    big = [p for p in Q if len(p) >= 4]
    b1 = sum(len(p[-1]["words"]) == 1 for p in big); b1h = sum(len(p[-1]["hard"]) == 1 for p in big)
    print("Pt included" if inc else "Pt excluded", "paragraphs", len(Q), "multi-line ending on one word:", n1, "(commas joined:", n1h, ")",
          "| >=4 lines:", len(big), "one-word endings", b1, "= %.3f" % (b1 / len(big)), "(commas joined:", b1h, "= %.3f)" % (b1h / len(big)))
