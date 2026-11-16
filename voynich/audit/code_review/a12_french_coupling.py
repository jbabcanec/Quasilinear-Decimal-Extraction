"""French control: (a) duplicate pages left? (b) coupling when a removed Heading/DropCapital
line is not silently bridged (line numbers skip), (c) last-word length with punctuation stripped."""
import sys, re, hashlib
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/control_real")
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np, pandas as pd
import control_real as cr
from wrap import breaks, line_table, margins, wrap_stats, boot_ci
NS = cr.NS
root = Path("/home/user/Gallicorpora/HTR-MSS-15e-Siecle/data")
MSS = ["btv1b100737746", "btv1b55008562q", "btv1b90069505", "btv1b9009670t", "btv1b90595325", "btv1b90616967"]

def read_page2(path):
    r = ET.parse(path).getroot()
    labels = {t.get("ID"): t.get("LABEL") for t in r.iter(f"{NS}OtherTag")}
    regions = []
    for block in r.iter(f"{NS}TextBlock"):
        if labels.get(block.get("TAGREFS", ""), "") != "MainZone": continue
        lines = []
        for tl in block.iter(f"{NS}TextLine"):
            lab = labels.get(tl.get("TAGREFS", ""), "DefaultLine")
            pts = [float(v) for v in re.split(r"[ ,]+", (tl.get("BASELINE") or "").strip()) if v]
            text = " ".join(s.get("CONTENT", "") for s in tl.iter(f"{NS}String")).strip()
            if len(pts) < 4 or not text: continue
            xs, ys = pts[0::2], pts[1::2]
            lines.append((min(xs), max(xs), float(np.mean(ys)), text, lab == "DefaultLine"))
        lines.sort(key=lambda r: r[2])
        if sum(l[4] for l in lines) >= 4: regions.append(lines)
    return regions

def table2(ms_dir):
    rows, par = [], 0
    for page in sorted(ms_dir.glob("*.chocomufin.xml")):
        ln = 0
        for region in read_page2(page):
            par += 1
            for xl, xr, y, text, keep in region:
                if not keep:
                    ln += 3; continue          # a removed line breaks adjacency
                words = text.split(); units = [cr.nchar(w) for w in words]
                total = sum(units) + cr.SPACE * (len(words) - 1)
                if total <= 0 or xr <= xl: continue
                px = (xr - xl) / total; x = xl
                for w, u in zip(words, units):
                    rows.append((page.stem, ln, par, w, x, x + u * px, y)); x += (u + cr.SPACE) * px
                ln += 1
            ln += 5
    df = pd.DataFrame(rows, columns=["folio", "line", "par", "word", "x0", "x1", "yc"])
    df["w"] = df.x1 - df.x0; df["xc"] = (df.x0 + df.x1) / 2
    df["k"] = df.groupby(["folio", "line"]).cumcount(); df["nline"] = df.groupby(["folio", "line"])["k"].transform("size")
    df["kr"] = df.nline - 1 - df.k; df["wmed"] = df.groupby("folio").w.transform("median"); df["par_last"] = False
    return df

strip = lambda w: re.sub(r"[^\w]", "", w)
for ms in MSS:
    files = sorted((root / ms).glob("*.chocomufin.xml"))
    hashes = [hashlib.md5(" ".join(t[3] for reg in cr.read_page(f) for t in reg).encode()).hexdigest() for f in files]
    d0 = cr.token_table(root / ms); d1 = table2(root / ms)
    out = []
    for d in (d0, d1):
        L = margins(line_table(d)); B = breaks(L)
        st = wrap_stats(B); out.append((st["n"], round(st["slope"], 3), round(st["med"], 2)))
    # trimmed: drop breaks after a line that stops early by more than 1.5 word widths (paragraph ends), for information
    L = margins(line_table(d0)); B = breaks(L); Bt = B[B.slack1 < 1.5]
    d0["nc"] = d0.word.map(lambda w: len(strip(w)))
    dd = d0[d0.nc > 0]
    last = dd[dd.kr == 0].nc.mean(); inside = dd[(dd.k > 0) & (dd.kr > 0)].nc.mean(); first = dd[dd.k == 0].nc.mean()
    print(ms, "pages", len(files), "distinct page texts", len(set(hashes)), "| as coded (n, slope, med):", out[0],
          "| removed lines break adjacency:", out[1], "| slack<1.5 only: n=%d slope %.3f" % (len(Bt), wrap_stats(Bt)["slope"]),
          "| chars w/o punctuation last/inside %.2f/%.2f = %+.0f%%, first %+.0f%%" % (last, inside, 100 * (last / inside - 1), 100 * (first / inside - 1)))
