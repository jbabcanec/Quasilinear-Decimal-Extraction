"""exp17 audit: (1) reproduce, (2) stuck paragraphs in pack(), (3) calibrated null where the
reference foot is estimated from the simulated pages exactly as it is from the real ones."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np, pandas as pd
from vdata import folio_sort_key, load_zl, paragraph_text, tokens
df = tokens()
pt = paragraph_text(df).sort_values(["folio", "line", "k"])
g = pt.groupby(["folio", "line"]).agg(y=("yc", "mean"), par=("par", "first"), sec=("illust", "first")).reset_index()
pages, loci = load_zl()

def build(secs):
    pg = []
    for f, h in g[g.sec.isin(secs)].groupby("folio"):
        h = h.sort_values("y")
        if len(h) < 20:
            print("   excluded (<20 lines):", f, len(h)); continue
        pitch = float(np.median(np.diff(h.y)))
        top, bottom = h.y.iloc[0] - pitch / 2, h.y.iloc[-1] + pitch / 2
        starts = np.append(np.sort(h.groupby("par").y.min().to_numpy()) - pitch / 2, bottom)
        pg.append(dict(f=f, key=folio_sort_key(f), top=top, bottom=bottom, pitch=pitch, heights=np.diff(starts), nl=len(h)))
    pg.sort(key=lambda p: p["key"])
    return pg

def analyse(secs, nsim=4000, seed=0):
    rng = np.random.default_rng(seed)
    pg = build(secs)
    ref = np.quantile([p["bottom"] for p in pg], 0.9)
    blank = np.array([(ref - p["bottom"]) / p["pitch"] for p in pg])
    cap = [ref - p["top"] for p in pg]; pitch = float(np.mean([p["pitch"] for p in pg]))
    allh = np.concatenate([p["heights"] for p in pg])
    obs = np.mean(blank < 0.5)
    print(" pages", len(pg), "paragraph heights", len(allh), "ref", round(ref, 1), "obs tight", obs,
          " obs with blank<=0:", np.mean(blank <= 0), " obs in (0,0.5):", np.mean((blank > 0) & (blank < 0.5)),
          " blank<-0.3:", np.mean(blank < -0.3))
    print(" blank (lines) sorted:", np.round(np.sort(blank), 1).tolist())
    print(" bottoms px:", np.round(sorted(p["bottom"] for p in pg)).tolist())
    print(" tallest paragraphs (lines):", np.round(np.sort(allh)[-6:] / pitch, 1).tolist(), " capacities (lines):", np.round(np.sort(cap)[:6] / pitch, 1).tolist(), "...", np.round(np.sort(cap)[-3:] / pitch, 1).tolist())

    def pack(order, cap, tol=0.3):
        res, used_l, i, stuck = [], [], 0, 0
        for c in cap:
            used = 0.0
            while i < len(order) and used + order[i] <= c + tol * pitch:
                used += order[i]; i += 1
            if i >= len(order):
                break
            if used == 0.0: stuck += 1
            res.append((c - used) / pitch); used_l.append(used)
        return np.array(res), np.array(used_l), stuck
    tight, tight_cal, stucks, npg = [], [], [], []
    tops = np.array([p["top"] for p in pg])
    for _ in range(nsim):
        s, used, st = pack(rng.permutation(allh), cap)
        tight.append(np.mean(s < 0.5)); stucks.append(st); npg.append(len(s))
        # calibrated: treat the simulated pages as data, estimate ref the same way
        bot = tops[:len(used)] + used
        r2 = np.quantile(bot, 0.9)
        tight_cal.append(np.mean((r2 - bot) / pitch < 0.5))
    tight, tight_cal = np.array(tight), np.array(tight_cal)
    print(" as coded:   packed tight %.3f sd %.3f  p=%.4f" % (tight.mean(), tight.std(), np.mean(tight >= obs)))
    print(" sims with a page left EMPTY because the next paragraph is taller than the page: %.3f ; mean empty pages per sim %.2f; pages filled per sim %.1f of %d"
          % (np.mean(np.array(stucks) > 0), np.mean(stucks), np.mean(npg), len(pg)))
    print(" calibrated: same 90th-percentile reference re-estimated on each simulated book: tight %.3f sd %.3f  p=%.4f"
          % (tight_cal.mean(), tight_cal.std(), np.mean(tight_cal >= obs)))
    return pg

for secs in (["S"], ["B", "T"], ["S", "B", "T"]):
    print(secs); analyse(secs)
