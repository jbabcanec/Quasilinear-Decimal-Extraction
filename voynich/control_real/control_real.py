"""Control: the Voynich line-end measurements on real medieval manuscripts of known text.

Reads line-level ALTO ground truth (baseline + transcription per line), gives every word a
synthetic box whose width is proportional to its character count, and passes the result
through the same code as the Voynich pages (src/wrap.py), including the simulated wrapper
and fitter on each manuscript's own words and geometry.

Usage: python3 control_real/control_real.py <dataset root> [<dataset root> ...]
Writes results/control_real.json
"""
from __future__ import annotations

import json
import re
import sys
import unicodedata
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from wrap import boot_ci, breaks, line_table, margins, simulate, wrap_stats   # noqa: E402

SPACE = 0.6            # a word space counts as this many characters of width
NS = "{http://www.loc.gov/standards/alto/ns-v4#}"


class Gap:
    """Stand-in for layout.WidthModel: a constant inter-word gap, in word widths."""
    def __init__(self, g): self.g = g
    def gap(self, rng): return self.g


def nchar(word: str) -> int:
    return sum(1 for c in word if not unicodedata.combining(c))


def read_page(path: Path):
    """-> list of regions; each region is a list of (xl, xr, y, text) sorted by y."""
    root = ET.parse(path).getroot()
    labels = {t.get("ID"): t.get("LABEL") for t in root.iter(f"{NS}OtherTag")}
    regions = []
    for block in root.iter(f"{NS}TextBlock"):
        if labels.get(block.get("TAGREFS", ""), "") != "MainZone":
            continue
        lines = []
        for tl in block.iter(f"{NS}TextLine"):
            if labels.get(tl.get("TAGREFS", ""), "DefaultLine") != "DefaultLine":
                continue
            pts = [float(v) for v in re.split(r"[ ,]+", (tl.get("BASELINE") or "").strip()) if v]
            text = " ".join(s.get("CONTENT", "") for s in tl.iter(f"{NS}String")).strip()
            if len(pts) < 4 or not text:
                continue
            xs, ys = pts[0::2], pts[1::2]
            lines.append((min(xs), max(xs), float(np.mean(ys)), text))
        lines.sort(key=lambda r: r[2])
        if len(lines) >= 4:
            regions.append(lines)
    return regions


def token_table(ms_dir: Path) -> pd.DataFrame:
    rows, par = [], 0
    for page in sorted(ms_dir.glob("*.xml")):
        ln = 0
        for region in read_page(page):
            par += 1
            for xl, xr, y, text in region:
                words = text.split()
                units = [nchar(w) for w in words]
                total = sum(units) + SPACE * (len(words) - 1)
                if total <= 0 or xr <= xl:
                    continue
                px = (xr - xl) / total
                x = xl
                for w, u in zip(words, units):
                    rows.append((page.stem, ln, par, w, x, x + u * px, y))
                    x += (u + SPACE) * px
                ln += 1
            ln += 5                      # regions are never "next line" to each other
    df = pd.DataFrame(rows, columns=["folio", "line", "par", "word", "x0", "x1", "yc"])
    df["w"] = df.x1 - df.x0
    df["xc"] = (df.x0 + df.x1) / 2
    df["k"] = df.groupby(["folio", "line"]).cumcount()
    df["nline"] = df.groupby(["folio", "line"])["k"].transform("size")
    df["kr"] = df.nline - 1 - df.k
    df["wmed"] = df.groupby("folio").w.transform("median")
    df["par_last"] = False
    return df


def margin_free(L):
    B = L[L.has_next & (L.par_lines >= 4)].copy()
    B["e"] = B.E / B.wmed
    key = ["folio", "par"]
    y = B.e - B.groupby(key).e.transform("mean")
    x = B.w_next - B.groupby(key).w_next.transform("mean")
    ok = y.abs() < 1.5
    return float((x[ok] * y[ok]).sum() / (x[ok] * x[ok]).sum())


def measure(df):
    L = margins(line_table(df))
    B = breaks(L)
    st = wrap_stats(B)
    st["margin_free"] = margin_free(L)
    return st, L, B


def main():
    out = {}
    for root in sys.argv[1:]:
        for ms in sorted(p for p in Path(root).glob("data/*") if p.is_dir()):
            df = token_table(ms)
            if df.empty or df.groupby(["folio", "line"]).ngroups < 150:
                continue
            st, L, B = measure(df)
            if st["n"] < 150:
                continue
            g = df.groupby(["folio", "line"])
            ends = L[L.has_next & (L.par_lines >= 4)]
            res = dict(n_pages=int(df.folio.nunique()), n_lines=int(g.ngroups),
                       words_per_line=float(g.size().mean()),
                       sample=[" ".join(h.word) for _, h in list(g)[10:13]],
                       real=st | {"slope_ci": boot_ci(B, wrap_stats, "slope", n=600)},
                       chars_first=float(df[df.k == 0].word.map(nchar).mean()),
                       chars_last=float(df[df.kr == 0].word.map(nchar).mean()),
                       chars_inside=float(df[(df.k > 0) & (df.kr > 0)].word.map(nchar).mean()))
            stub = Gap(SPACE / max(1e-9, float(df.word.map(nchar).median())))
            for mode, name in (("pour", "wrapper"), ("fit", "fitter")):
                runs = []
                for seed in range(8):
                    S = simulate(df, L, stub, np.random.default_rng(seed), mode=mode, tol=0.15)
                    S["par_last"] = False
                    runs.append(measure(S)[0])
                res[name] = {k: [float(np.mean([r[k] for r in runs])), float(np.min([r[k] for r in runs])),
                                 float(np.max([r[k] for r in runs]))] for k in ("slope", "med", "iqr", "margin_free")}
            out[ms.name] = res
            r = res["real"]
            print("%-18s pages %2d lines %4d w/line %4.1f | n=%4d blank med %.2f iqr %.2f coupling %+.2f [%+.2f,%+.2f] mfree %+.2f"
                  " | wrapper %.2f/%.2f/%+.2f  fitter %.2f/%.2f/%+.2f | %s" % (
                      ms.name, res["n_pages"], res["n_lines"], res["words_per_line"], r["n"], r["med"], r["iqr"],
                      r["slope"], *r["slope_ci"], r["margin_free"],
                      res["wrapper"]["med"][0], res["wrapper"]["iqr"][0], res["wrapper"]["slope"][0],
                      res["fitter"]["med"][0], res["fitter"]["iqr"][0], res["fitter"]["slope"][0],
                      res["sample"][0][:60]), flush=True)
    (ROOT / "results" / "control_real.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
