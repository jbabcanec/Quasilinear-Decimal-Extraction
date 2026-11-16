"""The three regularities, measured on a simulated layout exactly as the manuscript study
measures them on the Voynich (src/exp13_line_starts.py, exp16_last_lines.py, wrap.py),
plus the traces each copying strategy would leave.

Voynich values for reference:
  R1  same opening glyph as the line above 9.6% vs 16.0% shuffled (ratio 0.60);
      two lines down 20.9% vs 16.2% (ratio 1.29)
  R2  last line a single word 0.4% (boxes) - 2.2% (ZL3b); under 30% of a normal line
      5.0% - 7.2%; 40-80% of a normal line 66%; over-full 4.6%
  R3  median blank 0.13, middle-half spread 0.33, 90th centile 0.81;
      coupling to the next line's first word 0.05 (95% interval -0.03 to 0.14)
"""
from __future__ import annotations

from collections import Counter

import numpy as np

VOYNICH = dict(r1_ratio=0.60, r1_same=0.096, r1_shuffled=0.160, r1_lag2_ratio=1.29,
               r2_one_word=(0.004, 0.022), r2_under30=(0.050, 0.072), r2_mid=0.66, r2_overfull=0.046,
               r3_median=0.133, r3_iqr=0.333, r3_p90=0.81, r3_coupling=0.05, r3_coupling_ci=(-0.03, 0.14))


def r1(layout, lags=(1, 2, 3), shuffles=200, seed=0):
    """Paragraph-first lines set aside; same glyph at lag k against the same lines shuffled
    within the paragraph (exact expectation; sd from permutations)."""
    rng = np.random.default_rng(seed)
    seqs = [[ln["glyph"] for ln in par[1:]] for par in layout]
    seqs = [s for s in seqs if len(s) >= 2]
    out = {}
    for lag in lags:
        n = sum(max(0, len(s) - lag) for s in seqs)
        obs = sum(a == b for s in seqs for a, b in zip(s[:-lag], s[lag:]))
        exp = sum((len(s) - lag) * sum(c * (c - 1) for c in Counter(s).values()) / (len(s) * (len(s) - 1))
                  for s in seqs if len(s) > lag)
        null = []
        for _ in range(shuffles):
            k = 0
            for s in seqs:
                if len(s) > lag:
                    t = [s[i] for i in rng.permutation(len(s))]
                    k += sum(a == b for a, b in zip(t[:-lag], t[lag:]))
            null.append(k)
        sd = float(np.std(null)) or 1.0
        out[f"lag{lag}"] = dict(n=n, same=obs / n, shuffled=exp / n, ratio=obs / exp if exp else np.nan,
                                z=(obs - exp) / sd)
    return out


def r2(layout, min_lines=4):
    fill, nw, tight = [], [], []
    for par in layout:
        if len(par) < min_lines:
            continue
        body = np.median([ln["E"] for ln in par[1:-1]])
        last = par[-1]
        fill.append(last["E"] / body); nw.append(last["nwords"])
        tight.append(last["N0"] > last["W"] * 1.0001)
    fill, nw = np.array(fill), np.array(nw)
    return dict(n=len(fill), one_word=float((nw == 1).mean()), under_30=float((fill < 0.3).mean()),
                mid_40_80=float(((fill >= 0.4) & (fill < 0.8)).mean()), over_full=float((fill >= 1.0).mean()),
                squeezed_last=float(np.mean(tight)), mean_fill=float(fill.mean()))


def r3(layout, unit, half=3, page=4):
    """Blank at the line end, against the right boundary estimated from the neighbouring
    lines of the paragraph (wrap.margins: furthest line end among the three lines either
    side, paragraph-final lines not used), and its coupling to the width of the next
    line's first word (slope within pseudo-pages of `page` paragraphs)."""
    slack, wnext, grp, true_blank, structural = [], [], [], [], []
    for pi, par in enumerate(layout):
        if len(par) < 3:
            continue
        E = np.array([ln["E"] for ln in par])
        for i in range(len(par) - 1):
            idx = [j for j in range(max(0, i - half), min(len(par) - 1, i + half + 1)) if j != i]
            if not idx:
                continue
            slack.append((E[idx].max() - E[i]) / unit)
            wnext.append(par[i + 1]["first_len"] / unit)
            grp.append(pi // page)
            true_blank.append((par[i]["W"] - E[i]) / unit)
            structural.append(max(0.0, par[i]["W"] - par[i]["hi"]) / unit)
    slack, wnext, grp = np.array(slack), np.array(wnext), np.array(grp)
    x, y = wnext.copy(), slack.copy()
    for g in np.unique(grp):
        m = grp == g
        x[m] -= x[m].mean(); y[m] -= y[m].mean()
    structural = np.array(structural)
    return dict(n=len(slack), median=float(np.median(slack)),
                iqr=float(np.percentile(slack, 75) - np.percentile(slack, 25)),
                p90=float(np.percentile(slack, 90)), coupling=float((x * y).sum() / (x * x).sum()),
                true_median=float(np.median(true_blank)),
                cannot_reach_margin=float((structural > 0.25).mean()))


def traces(layout):
    """What the strategies leave behind.
       disp_last    the last word of a line begins with the glyph that begins the line
       disp_second  the second word of a line begins with the glyph that begins the line above
                    (a repeat avoided by re-flow has to land in one of these two places)
       adj_sd       spread of written width / natural width over full lines (elastic writing
                    and abbreviation show up as lines denser or looser than their neighbours)
       pen_*        the same for the three lines before a paragraph's last line
    """
    dl = dln = ds = dsn = 0
    adj, pen, early, splits, variants, frag, nlines, wpl = [], [], [], 0, 0, 0, 0, []
    for par in layout:
        for i, ln in enumerate(par):
            if not ln["final"]:
                a = ln["E"] / ln["N0"] - 1.0
                adj.append(a)
                (pen if i >= len(par) - 4 else early).append(a)
                splits += ln["split_end"]; nlines += 1
                if 0 < i:
                    wpl.append(ln["nwords"])
            if i >= 1:
                variants += ln["glyph"] != ln["default_glyph"]; frag += ln["frag"]
                if ln["last"] and ln["nwords"] >= 3:
                    dln += 1; dl += ln["last"] == ln["glyph"]
                if i >= 2 and ln["second"]:
                    dsn += 1; ds += ln["second"] == par[i - 1]["glyph"]
    adj, pen, early = np.array(adj), np.array(pen), np.array(early)
    n1 = sum(len(p) - 1 for p in layout)
    return dict(disp_last=dl / max(dln, 1), disp_second=ds / max(dsn, 1),
                adj_sd=float(adj.std()), adj_mean=float(adj.mean()),
                adj_over5=float((np.abs(adj) > 0.05).mean()),
                pen_mean=float(pen.mean()) if len(pen) else np.nan,
                early_mean=float(early.mean()) if len(early) else np.nan,
                pen_loose5=float((pen > 0.05).mean()) if len(pen) else np.nan,
                pen_tight5=float((pen < -0.05).mean()) if len(pen) else np.nan,
                early_loose5=float((early > 0.05).mean()) if len(early) else np.nan,
                early_tight5=float((early < -0.05).mean()) if len(early) else np.nan,
                split_rate=splits / max(nlines, 1), variant_rate=variants / max(n1, 1),
                fragment_starts=frag / max(n1, 1), words_per_line=float(np.mean(wpl)) if wpl else np.nan)


def displacement_profile(paras, shuffles=300, seed=0, min_words=5):
    """Where do the missing line-opening repeats go?

    `paras`: paragraphs as lists of lines, each line the list of its words' first glyphs,
    paragraph-first line already removed.  A copyist who dodges a repeat by moving the line
    break by one word leaves the dodged word either as the LAST word of the line above
    (line squeezed) or as the SECOND word of the line below (line above ended early).  So
    relative to the opening glyph g of line i, this reports how often the same glyph opens
      - the last / 2nd-last / 3rd-last word of line i itself, and
      - the 2nd / 3rd / 4th word of line i+1 (and, for reference, its 1st = R1 itself),
    each against the same statistic with the opening glyph taken from another line of the
    paragraph (exact expectation; sd from shuffles).  Lines of fewer than `min_words` words
    are skipped so the positions never coincide."""
    rng = np.random.default_rng(seed)
    paras = [p for p in paras if len(p) >= 2]
    spec = [("same line, last word", 0, -1), ("same line, 2nd-last word", 0, -2),
            ("same line, 3rd-last word", 0, -3), ("line below, 1st word (R1)", 1, 0),
            ("line below, 2nd word", 1, 1), ("line below, 3rd word", 1, 2), ("line below, 4th word", 1, 3)]
    out = {}
    for name, dl, pos in spec:
        obs = n = 0; exp = 0.0; null = np.zeros(shuffles)
        for p in paras:
            first = [ln[0] for ln in p]
            cnt = Counter(first)
            tg, src = [], []                           # target glyph, index of the line whose opening it is tested against
            for i in range(len(p) - dl):
                t = p[i + dl]
                if len(t) < min_words and pos != 0:
                    continue
                tg.append(t[pos]); src.append((i, i + dl))
            if not tg:
                continue
            m = len(p)
            for g, (i, t) in zip(tg, src):
                obs += g == first[i]; n += 1
                # expectation: opening glyph of a random line other than the target's own line
                exp += (cnt[g] - (first[t] == g)) / (m - 1)
            fa = np.array([hash(x) for x in first]); ga = np.array([hash(x) for x in tg])
            ti = np.array([t for _, t in src])
            for s in range(shuffles):
                # opening glyph of a random OTHER line of the paragraph
                j = (ti + rng.integers(1, m, len(ti))) % m
                null[s] += (fa[j] == ga).sum()
        sd = float(null.std()) or 1.0
        out[name] = dict(n=n, rate=obs / max(n, 1), expected=exp / max(n, 1),
                         lift=obs / exp if exp else np.nan, z=(obs - exp) / sd)
    return out


def summary(layout, unit):
    a, b, c, t = r1(layout), r2(layout), r3(layout, unit), traces(layout)
    out = dict(r1_same=a["lag1"]["same"], r1_shuffled=a["lag1"]["shuffled"], r1_ratio=a["lag1"]["ratio"],
               r1_z=a["lag1"]["z"], r1_n=a["lag1"]["n"], r1_lag2_ratio=a["lag2"]["ratio"], r1_lag2_z=a["lag2"]["z"],
               r1_lag2_same=a["lag2"]["same"], r1_lag3_same=a["lag3"]["same"], r1_lag3_ratio=a["lag3"]["ratio"])
    out.update({"r2_" + k: v for k, v in b.items()})
    out.update({"r3_" + k: v for k, v in c.items()})
    out.update({"t_" + k: v for k, v in t.items()})
    return out
