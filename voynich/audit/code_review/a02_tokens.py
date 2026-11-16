import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np, pandas as pd
import vdata
from vdata import tokens, build_tokens, paragraph_text, load_zl
df = tokens()
print(df.shape, list(df.columns))
if "--rebuild" in sys.argv:
    df2 = build_tokens()
    print("rebuild equal to cache:", df2.equals(df), df2.shape)
    df2.to_pickle("/home/user/Quasilinear-Decimal-Extraction/voynich/audit/code_review/tokens_rebuilt.pkl")
print("unaligned tokens:", (df.zl_locus < 0).sum(), (df.zl_locus < 0).mean())
print("zl_frac<0.5:", (df.zl_frac < 0.5).mean(), " <1:", (df.zl_frac < 1).mean())
print("line_agree false share of tokens:", (~df.line_agree).mean())
lines = df.groupby(["folio", "line"]).agg(n=("k", "size"), agree=("line_agree", "first"), ltype=("ltype", "first"),
        first_unal=("zl_locus", lambda s: s.iloc[0] < 0), any_unal=("zl_locus", lambda s: (s < 0).any()),
        is_par=("is_par", "first"), par=("par", "first")).reset_index()
print("lines:", len(lines), "agree:", lines.agree.mean())
print("lines by ltype:\n", lines.groupby("ltype").agg(n=("n", "size"), agree=("agree", "mean")).to_string())
P = lines[lines.is_par]
print("P lines:", len(P), "agree", P.agree.mean(), " first token unaligned:", P.first_unal.sum())
print("P lines by word count: agree share")
print(P.assign(nb=P.n.clip(upper=8)).groupby("nb").agg(lines=("n", "size"), agree=("agree", "mean")).to_string())
pt = paragraph_text(df)
print("pt tokens", len(pt), "lines", pt.groupby(["folio", "line"]).ngroups, "pars", pt.par.nunique())
