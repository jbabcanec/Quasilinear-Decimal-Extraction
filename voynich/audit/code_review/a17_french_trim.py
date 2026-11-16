"""French control: real vs own wrapper, with the same trim (blank < 1.5 word widths) on both."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/control_real")
from pathlib import Path
import numpy as np
import control_real as cr
from wrap import breaks, line_table, margins, simulate, wrap_stats
root = Path("/home/user/Gallicorpora/HTR-MSS-15e-Siecle/data")
for ms in ["btv1b100737746", "btv1b55008562q", "btv1b90069505", "btv1b9009670t", "btv1b90595325", "btv1b90616967"]:
    df = cr.token_table(root / ms)
    L = margins(line_table(df)); B = breaks(L)
    stub = cr.Gap(cr.SPACE / max(1e-9, float(df.word.map(cr.nchar).median())))
    full, trim = [], []
    for seed in range(4):
        S = simulate(df, L, stub, np.random.default_rng(seed), mode="pour", tol=0.15); S["par_last"] = False
        Bs = breaks(margins(line_table(S)))
        full.append(wrap_stats(Bs)["slope"]); trim.append(wrap_stats(Bs[Bs.slack1 < 1.5])["slope"])
    big = B[B.slack1 >= 1.5]
    print("%s real: all %+.2f (n=%d) | blank<1.5: %+.2f (n=%d; %d lines removed) || wrapper: all %+.2f | blank<1.5: %+.2f"
          % (ms, wrap_stats(B)["slope"], len(B), wrap_stats(B[B.slack1 < 1.5])["slope"], (B.slack1 < 1.5).sum(), len(big), np.mean(full), np.mean(trim)))
