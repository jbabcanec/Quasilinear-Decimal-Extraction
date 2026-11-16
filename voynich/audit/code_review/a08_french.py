import sys, re, unicodedata
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/control_real")
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import Counter
import numpy as np, pandas as pd
import control_real as cr
root = Path("/home/user/Gallicorpora/HTR-MSS-15e-Siecle/data")
NS = cr.NS
for ms in sorted(p for p in root.glob("*") if p.is_dir()):
    files = sorted(ms.glob("*.chocomufin.xml"))
    zl, ll, ncol = Counter(), Counter(), 0
    nsbad = 0
    for f in files:
        r = ET.parse(f).getroot()
        if not r.tag.startswith(NS): nsbad += 1
        labels = {t.get("ID"): t.get("LABEL") for t in r.iter(f"{NS}OtherTag")}
        mz = 0
        for b in r.iter(f"{NS}TextBlock"):
            lab = labels.get(b.get("TAGREFS", ""), "?")
            zl[lab] += 1
            mz += lab == "MainZone"
            for tl in b.iter(f"{NS}TextLine"):
                ll[(lab.split(":")[0], labels.get(tl.get("TAGREFS", ""), "?"))] += 1
        ncol += mz
    df = cr.token_table(ms)
    nl = df.groupby(["folio", "line"]).ngroups if len(df) else 0
    print(ms.name, "pages", len(files), "ns_bad", nsbad, "MainZone blocks", ncol, "lines kept", nl)
    print("    zones:", dict(zl)); print("    lines:", dict(ll))
