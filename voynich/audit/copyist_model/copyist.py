"""A copyist laying out a pre-existing text, with adjustable scribal freedoms.

THE MODEL
---------
Geometry.  A word's width is its number of letters; words are separated by a gap of
`gap` letters (0.7, i.e. about 0.125 of a word, the Voynich's measured gap).  Every line of
a paragraph has the same width W, set so that a line holds about nine words.  All reported
lengths are in units of the median written word.

The text is fixed: the copyist may not add, drop or reorder words.  What he may do:

 (a) abbreviate.  Each word has a number of letters that abbreviation can save
     (scribal.abbr_saving, or a generic "up to k letters in a share p of words").  `abbr`
     is the share of that saving the scribe is free to use or not; his habitual level sits
     at `abbr_default` of it, so he can both contract further and write out in full.
 (b) split a word across the line break (`split`: at syllable boundaries, or anywhere, with
     at least two letters on each side).  The next line then begins with the remainder.
 (c) write elastically: squeeze by up to `squeeze`, stretch by up to `stretch`.
     (a) and (c) act on the last `zone` words of the line (None = the whole line).
 (d) plan `lookahead` lines ahead (1 = this line's break only, seeing the word that would
     open the next line; None = the whole paragraph).  Planning is exact: dynamic
     programming over break points, re-run after every line (receding horizon).
 (e) end a line with a filler stroke of up to `filler` word widths.
 (f) begin a word with another glyph where spelling allows it (scribal.initial_variants).

So a line holding atoms a..b has a natural width N0 and can be written at any width in
[lo, hi].  It is admissible if lo <= W (it can be made to fit).  A line that cannot be
brought out to the margin (hi < W) leaves a planned blank W - hi.  To keep line ends
"within R3", a break that leaves a planned blank above `blank_cap` is admissible only if no
fuller line fits: the copyist never stops early just to move a word, unless he can fill
the space.

What he wants (costs, lexicographic in effect):
   paragraph ending  `widow`  none | single (no one-word last line) | short (also none
                              under 30% of a line) | shape (also prefer 40-80%)
   line openings     `avoid_repeat` > 0: do not begin a line with the glyph that begins
                              the line above;  `seek_lag2` > 0: prefer the glyph two above
   line ends         a small blank, then little squeezing/stretching, few splits, the
                     default spelling.

Execution.  A line that can reach the margin is ended at W - u, u exponential with mean
`jitter` word widths (clipped to [lo, hi]); this is the ragged right edge of real writing.
`jitter` is the one noise parameter; 0.28 makes a scribe with ample freedom show the
Voynich's median blank (0.13) under the manuscript study's own margin estimator
(measure.r3), so a simulated median above that means lines that could not be fitted.

Opening glyphs.  Latin spreads its word-initial letters more evenly than the Voynich
spreads its line-opening glyphs (chance of a repeat about 7% against 16%), which makes
repeats easier to dodge.  `coarse_glyphs` merges the Latin letters into eight classes
with the Voynich's shares (19, 18, 17, 16, 15, 7, 4, 3%) as a harder test.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from scribal import abbr_saving, initial_variants, split_points

BIG = 1e12


@dataclass
class Scribe:
    gap: float = 0.7
    words_per_line: float = 9.0
    abbr: float = 0.0
    abbr_default: float = 0.5
    abbr_generic: tuple | None = None      # (k letters, share p of words)
    split: str = "none"
    min_frag: int = 2
    squeeze: float = 0.0
    stretch: float = 0.0
    zone: int | None = 3
    lookahead: int | None = 1
    filler: float = 0.0
    variants: str = "none"
    avoid_repeat: float = 0.0
    seek_lag2: float = 0.0
    widow: str = "none"
    blank_cap: float = 0.3
    jitter: float = 0.28
    coarse_glyphs: bool = False            # merge letters into 8 classes with the Voynich's shares
    split_cost: float = 2.0
    variant_cost: float = 1.0

    def label(self):
        d = asdict(self); base = asdict(Scribe())
        return ", ".join(f"{k}={v}" for k, v in d.items() if v != base[k]) or "plain wrapper"


class Text:
    """A text prepared for one scribe: atoms, savings, line width."""

    def __init__(self, paras, sc: Scribe, seed: int = 0):
        rng = np.random.default_rng(seed + 12345)
        self.sc = sc
        self.paras = []
        lens, lens_def = [], []
        for par in paras:
            L, S, WI, WID, G, VAR, TXT = [], [], [], [], [], [], []
            for wid, w in enumerate(par):
                if sc.abbr_generic:
                    k, p = sc.abbr_generic
                    sv = float(min(k, max(0, len(w) - 2))) if rng.random() < p else 0.0
                else:
                    sv = sc.abbr * abbr_saving(w)
                pts = split_points(w, sc.split, sc.min_frag)
                cuts = (0,) + pts + (len(w),)
                for j in range(len(cuts) - 1):
                    pc = w[cuts[j]:cuts[j + 1]]
                    L.append(float(len(pc))); S.append(sv * len(pc) / len(w))
                    WI.append(j == 0); WID.append(wid); G.append(w[cuts[j]])
                    VAR.append(initial_variants(w, sc.variants) if j == 0 else ())
                    TXT.append(w[cuts[j]:])            # what a line starting here begins with
                lens.append(len(w)); lens_def.append(len(w) - sc.abbr_default * sv)
            n = len(L)
            cl = np.concatenate([[0.0], np.cumsum(L)]).tolist()
            cs = np.concatenate([[0.0], np.cumsum(S)]).tolist()
            cw = np.concatenate([[0], np.cumsum(WI)]).tolist()
            wstart = {}
            for k in range(n):
                wstart.setdefault(WID[k], k)
            self.paras.append(dict(n=n, cl=cl, cs=cs, cw=cw, WI=WI, WID=WID, G=G, VAR=VAR, TXT=TXT,
                                   wstart=wstart, words=par))
        if sc.coarse_glyphs:
            from collections import Counter
            freq = Counter(P["G"][k] for P in self.paras for k in range(P["n"]) if P["WI"][k])
            tot = sum(freq.values())
            room = [0.19, 0.18, 0.17, 0.16, 0.155, 0.07, 0.04, 0.035]
            cmap = {}
            for g, c in freq.most_common():
                j = int(np.argmax(room)); cmap[g] = "ABCDEFGH"[j]; room[j] -= c / tot
            self.glyph_map = cmap
            for P in self.paras:
                P["G"] = [cmap.get(g, g) for g in P["G"]]
                P["VAR"] = [tuple(v for v in dict.fromkeys(cmap.get(x, x) for x in var) if v != g)
                            for var, g in zip(P["VAR"], P["G"])]
        self.unit = float(np.median(lens_def))          # one "word width"
        self.mean_word = float(np.mean(lens_def))
        self.W = (sc.words_per_line + 0.5) * (self.mean_word + sc.gap) - sc.gap


def _geom(P, sc: Scribe, W: float, unit: float):
    """Returns line(a, b) -> (N0, lo, hi) for one paragraph."""
    cl, cs, cw, WID, wstart = P["cl"], P["cs"], P["cw"], P["WID"], P["wstart"]
    d, gap, zone, sq, st = sc.abbr_default, sc.gap, sc.zone, sc.squeeze, sc.stretch
    fill = sc.filler * unit

    def line(a, b):
        letters = cl[b] - cl[a]; sv = cs[b] - cs[a]
        gaps = cw[b] - cw[a + 1]
        N0 = letters - d * sv + gap * gaps
        if zone is None:
            zsv, Z0 = sv, N0
        else:
            z = max(a, wstart[max(WID[a], WID[b - 1] - zone + 1)])
            zsv = cs[b] - cs[z]
            Z0 = (cl[b] - cl[z]) - d * zsv + gap * (cw[b] - cw[z + 1])
        return N0, N0 - (1 - d) * zsv - sq * Z0, N0 + d * zsv + st * Z0 + fill
    return line


def lay_paragraph(P, sc: Scribe, W: float, unit: float, rng):
    n, WI, WID, G, VAR = P["n"], P["WI"], P["WID"], P["G"], P["VAR"]
    line = _geom(P, sc, W, unit)
    cand_memo, memo = {}, {}
    w_rep, w_lag2 = sc.avoid_repeat, sc.seek_lag2
    use_glyphs = w_rep > 0 or w_lag2 > 0
    max_shrink = 1.0 - sc.squeeze

    def candidates(a):
        """[(b, cost, final)] admissible ends of a line starting at atom a."""
        if a in cand_memo:
            return cand_memo[a]
        out, b = [], a + 1
        while b <= n:
            N0, lo, hi = line(a, b)
            if lo <= W + 1e-9:
                out.append((b, N0, lo, hi))
            elif (N0 * max_shrink - (P["cs"][b] - P["cs"][a])) > W:
                break
            b += 1
        if not out:                                    # a single over-long atom: take it anyway
            N0, lo, hi = line(a, a + 1); out = [(a + 1, N0, lo, hi)]
        nonfinal = [c for c in out if c[0] < n]
        fullest = max((c[0] for c in nonfinal), default=None)
        final_fits = out[-1][0] == n
        res = []
        for b, N0, lo, hi in out:
            if b == n:
                E = min(N0, W)
                nwords = WID[n - 1] - WID[a] + 1
                fillr = E / W
                c = 100.0 * abs(E - N0) / W
                if sc.widow != "none":
                    if nwords == 1:
                        c += 1e5
                    if sc.widow in ("short", "shape") and fillr < 0.30:
                        c += 3e4
                    if sc.widow == "shape":
                        c += 2e3 * max(0.0, 0.42 - fillr, fillr - 0.80)
                res.append((b, c, True))
            else:
                E = min(W, hi)
                blank = (W - E) / unit
                if blank > sc.blank_cap and not (b == fullest and not final_fits):
                    continue
                c = 100.0 * blank + 100.0 * abs(E - N0) / W + (0.0 if WI[b] else sc.split_cost)
                res.append((b, c, False))
        cand_memo[a] = res
        return res

    def glyph_options(b):
        return ((G[b], 0.0),) + tuple((g, sc.variant_cost) for g in VAR[b])

    def V(a, gc, gp, depth):
        """Least cost of finishing from a line that starts at atom a with glyph gc (line
        above: gp), looking `depth` lines ahead.  Returns (cost, b, next glyph)."""
        if depth == 0:
            return 0.0, None, None
        nd = depth - 1 if depth > 0 else -1            # depth -1 = to the end of the paragraph
        key = (a, gc, gp, depth)
        if key in memo:
            return memo[key]
        best = (BIG, None, None)
        for b, c, final in candidates(a):
            if final:
                if c < best[0]:
                    best = (c, b, None)
                continue
            if use_glyphs:
                for gn, vc in glyph_options(b):
                    t = c + vc + (w_rep if gn == gc else 0.0) - (w_lag2 if gn == gp else 0.0)
                    t += V(b, gn, gc, nd)[0]
                    if t < best[0]:
                        best = (t, b, gn)
            else:
                t = c + V(b, None, None, nd)[0]
                if t < best[0]:
                    best = (t, b, G[b])
        memo[key] = best
        return best

    depth = sc.lookahead if sc.lookahead is not None else -1
    lines, a = [], 0
    gc = G[0] if use_glyphs else None
    gp = None
    glyph_written = G[0]
    while a < n:
        _, b, gn = V(a, gc, gp, depth)
        N0, lo, hi = line(a, b)
        final = b == n
        if final:
            E = min(N0, W)
        else:
            u = rng.exponential(sc.jitter * unit) if sc.jitter > 0 else 0.0
            E = min(max(W - u, lo), min(hi, W))
        # first word (or fragment) of this line, as a transliteration would see it
        k = a + 1
        while k < b and not WI[k]:
            k += 1
        first_len = P["cl"][k] - P["cl"][a]
        second = G[k] if k < b else ""
        k2 = b - 1
        while k2 > a and not WI[k2]:
            k2 -= 1
        wg = [glyph_written] + [G[j] for j in range(a + 1, b) if WI[j]]   # first glyph of each word
        lines.append(dict(a=a, b=b, glyph=glyph_written, default_glyph=G[a], frag=not WI[a], wg=wg,
                          first_len=first_len, second=second, last=G[k2] if k2 > a else "",
                          nwords=WID[b - 1] - WID[a] + 1, N0=N0, lo=lo, hi=hi, E=E, final=final,
                          split_end=(not final) and (not WI[b]), W=W))
        if final:
            break
        glyph_written = gn if gn is not None else G[b]
        gp, gc = (gc, gn) if use_glyphs else (None, None)
        a = b
    return lines


def lay_text(paras, sc: Scribe, seed: int = 0):
    """-> (list of paragraphs, each a list of line dicts; Text)."""
    T = Text(paras, sc, seed)
    rng = np.random.default_rng(seed)
    return [lay_paragraph(P, sc, T.W, T.unit, rng) for P in T.paras], T
