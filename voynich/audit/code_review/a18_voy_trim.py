import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np
from layout import WidthModel
from vdata import paragraph_text, tokens
from wrap import breaks, line_table, margins, simulate, wrap_stats
pt = paragraph_text(tokens()); wm = WidthModel(pt)
L = margins(line_table(pt, wm)); B = breaks(L)
print("Voynich real: all %+.3f (n=%d), blank<1.5 %+.3f (n=%d), |blank|<1 %+.3f, share blank>=1.5: %.3f, share blank<0: %.3f"
      % (wrap_stats(B)["slope"], len(B), wrap_stats(B[B.slack1 < 1.5])["slope"], (B.slack1 < 1.5).sum(),
         wrap_stats(B[B.slack1.abs() < 1])["slope"], (B.slack1 >= 1.5).mean(), (B.slack1 < 0).mean()))
for seed in range(2):
    S = simulate(pt, L, wm, np.random.default_rng(1000 + seed), mode="pour", tol=0.15)
    Bs = breaks(margins(line_table(S, wm)))
    print("  wrapper seed", seed, "all %+.3f, blank<1.5 %+.3f, |blank|<1 %+.3f" % (wrap_stats(Bs)["slope"], wrap_stats(Bs[Bs.slack1 < 1.5])["slope"], wrap_stats(Bs[Bs.slack1.abs() < 1])["slope"]))
# by section, herbal vs not
ill = pt.groupby("folio").illust.first(); B["ill"] = B.folio.map(ill)
for k, b in B.groupby("ill"):
    if len(b) > 100: print("  ", k, len(b), "%+.3f" % wrap_stats(b)["slope"], "trim %+.3f" % wrap_stats(b[b.slack1 < 1.5])["slope"])
