import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np, pandas as pd
from layout import WidthModel
from vdata import paragraph_text, tokens
from wrap import breaks, line_table, margins, simulate, wrap_stats, boot_ci
import exp2_wrap as e2
for name, df in (("cache", tokens()), ("fixed join", pd.read_pickle("tokens_fixed.pkl"))):
    pt = paragraph_text(df)
    wm = WidthModel(pt)
    real, S0, L, B = e2.measure(pt, wm)
    print(name, {k: round(v, 4) if isinstance(v, float) else v for k, v in real.items()})
    if name == "cache":
        # sanity of the estimator: independent recomputation of the slope with paragraph-level max as margin
        Lx = L.copy()
        Lx["Mpar"] = Lx.groupby(["folio", "par"]).E.transform("max")
        Bx = Lx[Lx.has_next & (Lx.par_lines >= 3) & Lx.w_next.notna()].copy()
        Bx["s"] = (Bx.Mpar - Bx.E) / Bx.wmed
        x = Bx.w_next - Bx.groupby("folio").w_next.transform("mean"); y = Bx.s - Bx.groupby("folio").s.transform("mean")
        print("  slope with paragraph-max margin:", round(float((x * y).sum() / (x * x).sum()), 4), "n", len(Bx))
        # within-paragraph demeaning instead of folio
        x = B.w_next - B.groupby(["folio", "par"]).w_next.transform("mean"); y = B.slack1 - B.groupby(["folio", "par"]).slack1.transform("mean")
        print("  slope with paragraph fixed effects:", round(float((x * y).sum() / (x * x).sum()), 4))
        S = simulate(pt, L, wm, np.random.default_rng(1000), mode="pour", tol=0.15)
        st, _, Ls, Bs = e2.measure(S, wm)
        print("  wrapper seed0:", {k: round(v, 4) if isinstance(v, float) else v for k, v in st.items()})
        print("  sim lines used", len(Ls), "of", len(L), "; sim words", len(S), "of", len(pt))
