"""Experiment 2 - the line-break test.  Writes results/exp2_wrap.json.

Question: when a line ends, is that because the next word of a pre-existing text did not
fit (poured), or because the scribe chose what to write so as to fill the line (fitted)?
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from layout import WidthModel
from vdata import ROOT, paragraph_text, tokens
from wrap import GAP, boot_ci, breaks, line_table, margins, shuffle_next, simulate, wrap_stats

NSEED = 20
SCRIBES = {            # simulated scribes: (mode, overrun tolerance, elastic squeeze, pick)
    "pourer": ("pour", 0.15, 0.0, "first"),
    "squeezer": ("pour", 0.15, 0.30, "first"),
    "fitter": ("fit", 0.15, 0.0, "first"),
}


def natural(S, L, wm):
    """Glyph-predicted ('natural') widths: what the words should occupy however penned."""
    pred = lambda w: wm.coef[-1] + sum(wm.coef[wm.idx[ch]] for ch in w if ch in wm.idx)
    S = S.copy()
    S["pw"] = S.word.map(pred)
    S["res"] = S.w / S.wmed - S.pw
    nat = S.groupby(["folio", "line"]).agg(N=("pw", "sum"), n=("pw", "size"))
    L = L.join((nat.N + (nat.n - 1) * GAP).rename("Nat"), on=["folio", "line"])
    L["nslack"] = (L.M1 - L.xl) / L.wmed - L.Nat
    return S, L


def measure(S, wm):
    """All line-end statistics for one token table."""
    L = margins(line_table(S, wm))
    S, L = natural(S, L, wm)
    B = breaks(L)
    st = wrap_stats(B, "slack1", "w_next")
    st["slope_natural"] = wrap_stats(B, "nslack", "p_next")["slope"]
    st["slope_m2"] = wrap_stats(B, "slack2", "w_next")["slope"]
    ww = S.w / S.wmed
    m = S.nline >= 6
    inner = ww[m & (S.k > 1) & (S.kr > 1)].mean()
    st["first_rel"] = float(ww[m & (S.k == 0)].mean() / inner)
    st["last_rel"] = float(ww[m & (S.kr == 0)].mean() / inner)
    return st, S, L, B


def main():
    pt = paragraph_text(tokens())
    wm = WidthModel(pt)
    out = {"n_tokens": len(pt), "gap": GAP, "width_model_r2": wm.r2, "width_model_sigma": wm.sigma}

    real, S0, L, B = measure(pt, wm)
    meta = pt.groupby("folio")[["illust", "lang", "hand"]].first()
    B = B.join(meta, on="folio")
    for key, sl, wn in (("slope", "slack1", "w_next"), ("slope_natural", "nslack", "p_next"),
                        ("med", "slack1", "w_next"), ("iqr", "slack1", "w_next")):
        k = "slope" if key.startswith("slope") else key
        real[key + "_ci"] = boot_ci(B, lambda b, sl=sl, wn=wn: wrap_stats(b, sl, wn), k, n=1000)
    rng = np.random.default_rng(0)
    sh = np.array([[wrap_stats(shuffle_next(B, rng), "slack1")["slope"],
                    wrap_stats(shuffle_next(B, rng, "p_next"), "nslack", "p_next")["slope"]]
                   for _ in range(1000)])
    real["shuffle_slope_sd"] = float(sh[:, 0].std())
    real["shuffle_slope_natural_sd"] = float(sh[:, 1].std())
    real["p_slope_vs_shuffle"] = float((np.abs(sh[:, 0]) >= abs(real["slope"])).mean())
    out["real"] = real
    out["slack_values_real"] = np.round(B.slack1.to_numpy(), 3).tolist()
    out["wnext_values_real"] = np.round(B.w_next.to_numpy(), 3).tolist()

    # subsets of the manuscript
    sub = {}
    names = {"H": "herbal", "S": "stars / recipes", "B": "balneological", "T": "text only",
             "P": "pharmaceutical"}
    for code, nm in names.items():
        b = B[B.illust == code]
        if len(b) >= 100:
            sub[nm] = wrap_stats(b) | {"slope_ci": boot_ci(b, wrap_stats, "slope", n=500)}
    for lang in "AB":
        b = B[B.lang == lang]
        sub["Currier " + lang] = wrap_stats(b) | {"slope_ci": boot_ci(b, wrap_stats, "slope", n=500)}
    for hand in sorted(B.hand.unique()):
        b = B[B.hand == hand]
        if len(b) >= 100:
            sub["scribe " + str(hand)] = wrap_stats(b) | {"slope_ci": boot_ci(b, wrap_stats, "slope", n=500)}
    out["subsets"] = sub

    # simulated scribes on the real pages
    sims = {}
    example = {}
    for name, (mode, tol, el, pick) in SCRIBES.items():
        runs = []
        for seed in range(NSEED):
            S = simulate(pt, L, wm, np.random.default_rng(1000 + seed), mode=mode, tol=tol,
                         elastic=el, pick=pick)
            st, _, _, Bs = measure(S, wm)
            runs.append(st)
            if seed == 0:
                example[name] = (np.round(Bs.slack1.to_numpy(), 3).tolist(),
                                 np.round(Bs.w_next.to_numpy(), 3).tolist())
        sims[name] = {k: [float(np.mean([r[k] for r in runs])), float(np.std([r[k] for r in runs])),
                          float(np.min([r[k] for r in runs])), float(np.max([r[k] for r in runs]))]
                      for k in ("slope", "slope_natural", "slope_m2", "med", "iqr", "p90", "fit",
                                "first_rel", "last_rel")}
        print(name, {k: round(v[0], 3) for k, v in sims[name].items()}, flush=True)
    out["sims"] = sims
    out["sim_examples"] = {k: {"slack": v[0], "wnext": v[1]} for k, v in example.items()}

    # ---- loophole 1: elastic penmanship (is the last word physically squeezed / stretched?)
    el = {}
    q = S0.sort_values(["folio", "line", "k"])
    q["gap_after"] = (q.groupby(["folio", "line"]).x0.shift(-1) - q.x1) / q.wmed
    for kr in range(6):
        m = (q.kr == kr) & (q.nline >= 7)
        gm = m & (q.gap_after < 1)
        el[f"kr{kr}"] = dict(n=int(m.sum()), res_mean=float(q.res[m].mean()), res_sd=float(q.res[m].std()),
                             gap_mean=float(q.gap_after[gm].mean()) if kr else None)
    Lb = L.join(S0[S0.kr == 0].set_index(["folio", "line"])[["res", "pw"]]
                .rename(columns={"res": "res_last", "pw": "pw_last"}), on=["folio", "line"])
    Lb = Lb[Lb.M1.notna() & (Lb.par_lines >= 3) & Lb.has_next]
    Lb["over"] = Lb.pw_last - (Lb.M1 - Lb.x0_last) / Lb.wmed
    x = Lb.over - Lb.groupby("folio").over.transform("mean")
    y = Lb.res_last - Lb.groupby("folio").res_last.transform("mean")
    el["slope_residual_on_overfill"] = float((x * y).sum() / (x * x).sum())
    out["elasticity"] = el

    # ---- loophole 2: words split across the line break
    vocab = Counter(pt.word)
    q = pt.sort_values(["folio", "line", "k"])
    w, ln, fo, par = (q[c].to_numpy() for c in ("word", "line", "folio", "par"))
    within = (ln[1:] == ln[:-1]) & (fo[1:] == fo[:-1])
    across = (ln[1:] == ln[:-1] + 1) & (fo[1:] == fo[:-1]) & (par[1:] == par[:-1])
    joined = np.array([a + b for a, b in zip(w[:-1], w[1:])])
    att = np.array([vocab.get(s, 0) >= 2 for s in joined])
    short = np.array([len(s) <= 8 for s in joined])
    out["splitting"] = dict(
        within_line=float(att[within].mean()), across_break=float(att[across].mean()),
        within_line_short=float(att[within & short].mean()),
        across_break_short=float(att[across & short].mean()),
        n_within=int(within.sum()), n_across=int(across.sum()))

    # ---- glyph counts by position
    gl = pt.word.str.len()
    m = pt.nline >= 6
    out["glyphs_by_position"] = dict(first=float(gl[m & (pt.k == 0)].mean()),
                                     interior=float(gl[m & (pt.k > 1) & (pt.kr > 1)].mean()),
                                     last=float(gl[m & (pt.kr == 0)].mean()))
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "exp2_wrap.json").write_text(json.dumps(out, indent=1))
    show = {k: v for k, v in real.items()}
    print("REAL", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in show.items()})
    for k, v in sub.items():
        print("  %-18s n=%d slope %+.3f %s  med %.3f iqr %.3f" % (
            k, v["n"], v["slope"], np.round(v["slope_ci"], 3), v["med"], v["iqr"]))
    print("elasticity", {k: v for k, v in el.items() if k.startswith("slope")}, el["kr0"], el["kr3"])
    print("splitting", out["splitting"])


if __name__ == "__main__":
    main()
