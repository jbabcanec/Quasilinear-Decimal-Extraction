"""exp17 calibrated null, larger simulation, two symmetric definitions of 'the foot'."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np, pandas as pd
from vdata import folio_sort_key, paragraph_text, tokens
pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
g = pt.groupby(["folio", "line"]).agg(y=("yc", "mean"), par=("par", "first"), sec=("illust", "first")).reset_index()
def run(secs, nsim=20000, seed=7):
    rng = np.random.default_rng(seed)
    pg = []
    for f, h in g[g.sec.isin(secs)].groupby("folio"):
        h = h.sort_values("y")
        if len(h) < 20: continue
        pitch = float(np.median(np.diff(h.y)))
        top, bottom = h.y.iloc[0] - pitch / 2, h.y.iloc[-1] + pitch / 2
        starts = np.append(np.sort(h.groupby("par").y.min().to_numpy()) - pitch / 2, bottom)
        pg.append(dict(key=folio_sort_key(f), top=top, bottom=bottom, pitch=pitch, heights=np.diff(starts)))
    pg.sort(key=lambda p: p["key"])
    tops = np.array([p["top"] for p in pg]); bots = np.array([p["bottom"] for p in pg]); pit = np.array([p["pitch"] for p in pg])
    pitch = float(pit.mean()); allh = np.concatenate([p["heights"] for p in pg])
    ref = np.quantile(bots, 0.9); cap = ref - tops
    stat = {"q90": lambda b, pp: np.mean((np.quantile(b, 0.9) - b) / pp < 0.5),
            "max": lambda b, pp: np.mean((b.max() - b) / pp < 0.5)}
    obs = {k: f(bots, pit) for k, f in stat.items()}
    coded, cal = [], {k: [] for k in stat}
    for _ in range(nsim):
        order = rng.permutation(allh); i = 0; used_l = []
        for c in cap:
            used = 0.0
            while i < len(order) and used + order[i] <= c + 0.3 * pitch:
                used += order[i]; i += 1
            if i >= len(order): break
            used_l.append(used)
        used_l = np.array(used_l); b = tops[:len(used_l)] + used_l
        coded.append(np.mean((cap[:len(used_l)] - used_l) / pitch < 0.5))
        for k, f in stat.items(): cal[k].append(f(b, pitch))
    coded = np.array(coded)
    print(secs, "observed (as coded) %.2f | null as coded %.3f±%.3f p=%.4f" % (obs["q90"], coded.mean(), coded.std(), np.mean(coded >= obs["q90"])))
    for k in stat:
        c = np.array(cal[k]); print("    same estimator on both sides, foot=%s: observed %.2f  null %.3f±%.3f  p=%.4f" % (k, obs[k], c.mean(), c.std(), np.mean(c >= obs[k])))
for secs in (["S"], ["B", "T"], ["S", "B", "T"]): run(secs)
