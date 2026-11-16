import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/control_real")
from pathlib import Path
import numpy as np, pandas as pd
from control_real import token_table
from wrap import line_table, margins
root = Path("/home/user/Gallicorpora/HTR-MSS-15e-Siecle/data")
parts = []
for ms in ["btv1b100737746", "btv1b55008562q", "btv1b90069505", "btv1b9009670t", "btv1b90595325", "btv1b90616967"]:
    d = token_table(root / ms); d["folio"] = ms + "/" + d.folio; d["par"] = ms + "/" + d.par.astype(str); parts.append(d)
df = pd.concat(parts, ignore_index=True)
L = line_table(df)
M = margins(L); M = M[M.M1.notna()].copy()
M["fill"] = (M.E - M.xl) / (M.M1 - M.xl)
S = M[(M.fill < 0.9) & (M.fill > 0)]
print("lines", len(L), "early-stopping", len(S), "one word", (S.n == 1).sum(), (S.n == 1).mean())
txt = df.groupby(["folio", "line"]).word.agg(" ".join)
S = S.join(txt.rename("text"), on=["folio", "line"])
# position of each early-stopping line in its region
L["rank"] = L.groupby(["folio", "par"]).cumcount(); L["size"] = L.groupby(["folio", "par"]).line.transform("size")
S = S.join(L.set_index(["folio", "line"])[["rank", "size"]], on=["folio", "line"])
S["region_last"] = S["rank"] == S["size"] - 1
print("early-stopping lines that are the LAST line of their zone:", S.region_last.sum(), " one-word among those:", ((S.n == 1) & S.region_last).sum())
print(S[S.n == 1][["folio", "line", "rank", "size", "fill", "text"]].to_string())
print(S.groupby(S.folio.str.split("/").str[0]).agg(n=("n", "size"), one=("n", lambda s: (s == 1).sum())).to_string())
