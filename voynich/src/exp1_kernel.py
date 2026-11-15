"""Experiment 1 - is a word copied from the word physically above it?
Writes results/exp1_kernel.json.

Measures the manuscript, then proves the measurement on synthetic pages (real geometry)
where the answer is known, and uses those to convert the manuscript's number into an upper
bound on the share of words that can be copies of the word overhead.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from kernel import ADX_EDGES, D_BANDS, DX_EDGES, lift_map, line_decay, local_contrast
from layout import Vocab, WidthModel, generate_onpage, keep_lines, page_templates, pour
from pairs import build_pairs, ruler_lifts
from vdata import EXT, ROOT, paragraph_text, tokens

LAMS = (0.01, 0.02, 0.05, 0.10)
N_NULL, N_SEED = 12, 4


def summarise(S, full=False, nboot=1):
    P = build_pairs(S, dmax=20)
    r = {"n_tokens": int(len(S)), "words_per_line": float(S.groupby(["folio", "line"]).size().mean())}
    c = local_contrast(P, "e1", nboot=max(nboot, 2))
    r["C"] = c["C"]
    r["parts"] = {k: c[k] for k in ("near_line_overhead", "near_line_aside",
                                    "far_line_overhead", "far_line_aside")}
    if nboot > 2:
        r["C_ci"] = [c["lo"], c["hi"]]
    r["rulers"] = ruler_lifts(P[P.d == 1])
    ld = line_decay(P[P.interior], nboot=2)["lift"]
    r["decay"] = [float(v) for v in ld]
    if full:
        return r, P
    return r


def main():
    pt = paragraph_text(tokens())
    T = page_templates(pt); wm = WidthModel(pt); V = Vocab(pt.word)
    out = {}

    # ------------------------------------------------------------------ the manuscript
    real, P = summarise(pt, full=True, nboot=2000)
    real["events_overhead"] = local_contrast(P, "e1", nboot=2)["events_overhead"]
    for sim in ("ex", "e2"):
        c = local_contrast(P, sim, nboot=2000)
        real["C_" + sim] = [c["C"], c["lo"], c["hi"]]
    c = local_contrast(P[P.interior], "e1", nboot=2000)
    real["C_interior"] = [c["C"], c["lo"], c["hi"]]
    c = local_contrast(P, "e1", near=(0, 0.5), nboot=2000)
    real["C_tight"] = [c["C"], c["lo"], c["hi"]]
    c = local_contrast(P, "e1", d_near=(1, 1), nboot=2000)
    real["C_line1"] = [c["C"], c["lo"], c["hi"]]
    meta = pt.groupby("folio")[["lang", "illust"]].first()
    Pm = P.join(meta, on="folio")
    for lang in "AB":
        c = local_contrast(Pm[Pm.lang == lang], "e1", nboot=1000)
        real["C_lang" + lang] = [c["C"], c["lo"], c["hi"]]
    m = lift_map(P[P.interior], "e1", "dx", DX_EDGES, nboot=500)
    real["map"] = {k: np.round(m[k], 4).tolist() for k in ("lift", "lo", "hi")}
    real["map_dx_edges"] = DX_EDGES.tolist(); real["map_d_bands"] = D_BANDS
    ld = line_decay(P[P.interior], nboot=500)
    real["decay_ci"] = [np.round(ld["lo"], 4).tolist(), np.round(ld["hi"], 4).tolist()]
    out["real"] = real
    print("REAL  C=%.3f [%.3f, %.3f]  rulers phys %.3f ord %.3f  paired %s" % (
        real["C"], *real["C_ci"], real["rulers"]["phys"]["lift"], real["rulers"]["ord"]["lift"],
        real["rulers"]["paired"]), flush=True)
    del P, Pm

    # --------------------------------------------------- synthetic pages, known answer
    cal = {"null": [], "phys": {}, "ord": {}, "diffuse": {}}
    for seed in range(N_NULL):
        S = generate_onpage(T, V, wm, np.random.default_rng(100 + seed), lam=0.0)
        cal["null"].append(summarise(S))
        print("null", seed, round(cal["null"][-1]["C"], 3), flush=True)
    for rule in ("phys", "ord", "diffuse"):
        for lam in LAMS:
            runs = []
            for seed in range(N_SEED):
                S = generate_onpage(T, V, wm, np.random.default_rng(200 + seed), rule=rule, lam=lam)
                runs.append(summarise(S))
            cal[rule][str(lam)] = runs
            print(rule, lam, [round(r["C"], 3) for r in runs],
                  [(round(r["rulers"]["phys"]["lift"], 2), round(r["rulers"]["ord"]["lift"], 2)) for r in runs],
                  flush=True)
    out["calibration"] = cal

    # upper bound on overhead copying: C(lam) is linear in lam for the 'phys' rule
    null_C = np.array([r["C"] for r in cal["null"]])
    xs = [0.0] * len(null_C) + [lam for lam in LAMS for _ in range(N_SEED)]
    ys = list(null_C) + [r["C"] for lam in LAMS for r in cal["phys"][str(lam)]]
    slope, icpt = np.polyfit(xs, ys, 1)
    out["bound"] = dict(null_mean=float(null_C.mean()), null_sd=float(null_C.std(ddof=1)),
                        slope=float(slope), intercept=float(icpt),
                        lam_point=float(max(0.0, (real["C"] - icpt) / slope)),
                        lam_upper95=float((real["C_ci"][1] - icpt) / slope))
    print("BOUND", out["bound"], flush=True)

    # ------------------------------------------------------------- other texts as controls
    ctl = {}
    q = pt.sort_values(["folio", "line", "k"])
    S = pour(q.word.tolist(), T, wm, np.random.default_rng(3), widths=(q.w / q.wmed).tolist())
    ctl["Voynich word stream, re-wrapped"] = summarise(S, nboot=500)
    for seed in (1, 2, 3):
        lines = [ln.split() for ln in (EXT / "controls" / "selfcitation" / f"selfcite_seed{seed}.txt")
                 .read_text().splitlines() if ln.strip()]
        S = keep_lines(lines, T, wm, np.random.default_rng(seed))
        ctl[f"Timm-Schinner generator, own lines (seed {seed})"] = summarise(S, nboot=500)
    words = (EXT / "controls" / "naibbe" / "greshko_nathist_output_ciphertext.txt").read_text().split()
    S = pour(words, T, wm, np.random.default_rng(4))
    ctl["Naibbe cipher of Latin, wrapped"] = summarise(S, nboot=500)
    lat = []
    for f in sorted((EXT / "voynich-units" / "herbal").glob("pliny*.tess")):
        for ln in f.read_text(errors="replace").splitlines():
            lat.extend(re.findall(r"[a-z]+", re.sub(r"<[^>]*>", " ", ln.lower())))
    S = pour(lat, T, wm, np.random.default_rng(5), known=False)
    ctl["Latin (Pliny), wrapped"] = summarise(S, nboot=500)
    for k, v in ctl.items():
        print("%-48s C=%.3f %s rulers phys %.2f ord %.2f decay d1 %.2f d5 %.2f" % (
            k, v["C"], np.round(v["C_ci"], 3), v["rulers"]["phys"]["lift"], v["rulers"]["ord"]["lift"],
            v["decay"][0], v["decay"][4]), flush=True)
    out["controls"] = ctl

    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "exp1_kernel.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
