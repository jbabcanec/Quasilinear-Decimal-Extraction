"""Data layer: word boxes (voynichese.com geometry) joined to the ZL3b transliteration.

Every analysis in this study reads the token table built here:

    one row per word token that has a pixel box, with
      - geometry  : x0, x1, xc, y0, yc, w, h   (pixels on the voynichese.com page image)
      - text      : word (Takahashi EVA, as boxed)
      - structure : folio, line id, index in line (k), index from line end (kr),
                    paragraph id, line-in-paragraph, paragraph first/last line flags
      - metadata  : Currier language, hand, illustration type, quire (from ZL3b page headers)

Lines are defined twice, independently:
  * geometrically  (a new line starts when x jumps back to the left), and
  * textually      (the ZL3b locus each boxed word aligns to, by character alignment).
`line` is the textual one where the alignment is confident, else the geometric one.
"""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "external" / "voynich-units" / "data"
BOX_DIR = EXT / "voynich-units" / "morphometry_voynichese" / "voynichese_boxes"
ZL_PATH = EXT / "voynich-units" / "ZL3b.txt"
CACHE = ROOT / "cache"

_PAGE_RE = re.compile(r"^<(f\d+[rv]\d*)>\s*<!(.*)>")
_LOCUS_RE = re.compile(r"^<(f\d+[rv]\d*)\.(\d+[a-z]?),([@+*=~&/!])(\w\w?)(?:;\w)?>\s*(.*)$")


def folio_sort_key(folio: str):
    m = re.match(r"f(\d+)([rv])(\d*)", folio)
    return (int(m.group(1)), m.group(2), int(m.group(3) or 0))


def clean_zl_text(raw: str):
    """IVTFF line body -> (list of words, paragraph_start, paragraph_end).

    Alternative readings [a:b] take the first; ligature braces are dropped; inline
    comments and special-code glyphs are removed; '<->' (drawing intrusion) and both
    '.' and ',' are word breaks.
    """
    pstart = "<%>" in raw
    pend = "<$>" in raw
    s = raw.replace("<%>", "").replace("<$>", "").replace("<->", ".")
    s = re.sub(r"<![^>]*>", "", s)
    s = re.sub(r"<[^>]*>", ".", s)
    s = re.sub(r"\[([^:\]]*):[^\]]*\]", r"\1", s)
    s = re.sub(r"@\d+;", "?", s)
    s = s.replace("{", "").replace("}", "").replace("'", "")
    words = [w for w in re.split(r"[.,\s]+", s) if w]
    return words, pstart, pend


@lru_cache(maxsize=1)
def load_zl():
    """Return (pages, loci). pages: folio -> header vars. loci: list of dicts in file order."""
    pages, loci = {}, []
    for raw in ZL_PATH.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("#") or not raw.strip():
            continue
        m = _PAGE_RE.match(raw)
        if m and "." not in raw.split(">")[0]:
            pages[m.group(1)] = dict(re.findall(r"\$(\w)=(\w+)", m.group(2)))
            continue
        m = _LOCUS_RE.match(raw)
        if not m:
            continue
        folio, num, pos, ltype, body = m.groups()
        words, pstart, pend = clean_zl_text(body)
        loci.append(dict(folio=folio, num=num, pos=pos, ltype=ltype, words=words,
                         pstart=pstart, pend=pend))
    return pages, loci


def load_boxes(folio: str):
    data = json.loads((BOX_DIR / f"{folio}.js").read_text())
    vocab = [item[0] for item in data[0]]
    return [dict(word=vocab[e[0]], x0=e[1], y0=e[2], w=e[3], h=e[4]) for e in data[1]]


def _align_chars(box_words, zl_words_by_locus):
    """Character-level alignment of the boxed word stream to the ZL locus stream.

    Returns, per boxed word, the index of the ZL locus that most of its characters align
    to (or -1), and the fraction of its characters that aligned.
    """
    a_chars, a_owner = [], []
    for i, w in enumerate(box_words):
        a_chars.extend(w)
        a_owner.extend([i] * len(w))
    b_chars, b_owner = [], []
    for j, words in enumerate(zl_words_by_locus):
        for w in words:
            b_chars.extend(w)
            b_owner.extend([j] * len(w))
    votes = [Counter() for _ in box_words]
    sm = SequenceMatcher(a=a_chars, b=b_chars, autojunk=False)
    for blk in sm.get_matching_blocks():
        for o in range(blk.size):
            votes[a_owner[blk.a + o]][b_owner[blk.b + o]] += 1
    locus = np.full(len(box_words), -1)
    frac = np.zeros(len(box_words))
    for i, c in enumerate(votes):
        if c:
            j, n = c.most_common(1)[0]
            locus[i] = j
            frac[i] = n / max(1, len(box_words[i]))
    return locus, frac


def build_tokens() -> pd.DataFrame:
    pages, loci = load_zl()
    by_folio = defaultdict(list)
    for L in loci:
        by_folio[L["folio"]].append(L)

    rows = []
    for path in sorted(BOX_DIR.glob("*.js"), key=lambda p: folio_sort_key(p.stem)):
        folio = path.stem
        boxes = load_boxes(folio)
        if not boxes:
            continue
        meta = pages.get(folio, {})
        zl = by_folio.get(folio, [])
        locus, frac = _align_chars([b["word"] for b in boxes], [L["words"] for L in zl])

        # geometric lines: x jumps back left
        gline, g = [], 0
        for i, b in enumerate(boxes):
            if i and b["x0"] < boxes[i - 1]["x0"] - 3:
                g += 1
            gline.append(g)

        for i, b in enumerate(boxes):
            j = int(locus[i])
            L = zl[j] if j >= 0 else None
            rows.append(dict(
                folio=folio, seq=i, word=b["word"],
                x0=b["x0"], x1=b["x0"] + b["w"], xc=b["x0"] + b["w"] / 2,
                y0=b["y0"], yc=b["y0"] + b["h"] / 2, w=b["w"], h=b["h"],
                gline=gline[i], zl_locus=j, zl_frac=frac[i],
                zl_num=L["num"] if L else "", ltype=L["ltype"] if L else "",
                zl_pos=L["pos"] if L else "",
                zl_pstart=bool(L and L["pstart"]), zl_pend=bool(L and L["pend"]),
                lang=meta.get("L", "?"), hand=meta.get("H", "?"),
                illust=meta.get("I", "?"), quire=meta.get("Q", "?"),
            ))
    df = pd.DataFrame(rows)
    return _add_structure(df)


def _add_structure(df: pd.DataFrame) -> pd.DataFrame:
    """Final line ids, in-line indices, paragraph structure, per-page scale."""
    out = []
    for folio, g in df.groupby("folio", sort=False):
        g = g.copy()
        # A physical line = maximal run that shares both the geometric line and the ZL locus.
        # Where they disagree we split, so a 'line' is never wider than either definition.
        zl = g["zl_locus"].to_numpy().copy()
        # fill unaligned tokens from neighbours inside the same geometric line
        gl = g["gline"].to_numpy()
        for i in range(len(zl)):
            if zl[i] < 0:
                prev = zl[i - 1] if i and gl[i - 1] == gl[i] else -1
                nxt = -1
                for k in range(i + 1, len(zl)):
                    if gl[k] != gl[i]:
                        break
                    if zl[k] >= 0:
                        nxt = zl[k]
                        break
                zl[i] = prev if prev >= 0 else nxt
        g["zl_locus_f"] = zl
        line = np.zeros(len(g), dtype=int)
        cur = 0
        for i in range(1, len(g)):
            if gl[i] != gl[i - 1] or zl[i] != zl[i - 1]:
                cur += 1
            line[i] = cur
        g["line"] = line
        g["line_agree"] = True
        # flag lines where the geometric and textual segmentations disagree
        for ln, h in g.groupby("line"):
            same_g = (g["gline"] == h["gline"].iloc[0]).sum()
            same_z = (g["zl_locus_f"] == h["zl_locus_f"].iloc[0]).sum()
            if same_g != len(h) or same_z != len(h):
                g.loc[h.index, "line_agree"] = False
        g["k"] = g.groupby("line").cumcount()
        g["nline"] = g.groupby("line")["k"].transform("size")
        g["kr"] = g["nline"] - 1 - g["k"]
        g["wmed"] = float(np.median(g["w"]))
        out.append(g)
    df = pd.concat(out, ignore_index=True)

    # A line has one locus type: the commonest among its aligned words.  (Words that did not
    # align to ZL inherit their line's type, so lines are kept or dropped whole.)
    lt = (df[df["ltype"] != ""].groupby(["folio", "line"])["ltype"]
          .agg(lambda s: s.value_counts().index[0]).rename("ltype_line"))
    df = df.join(lt, on=["folio", "line"])
    df["ltype"] = df["ltype_line"].fillna("")
    df = df.drop(columns="ltype_line")

    # paragraph structure from ZL markers, carried on physical lines
    df["is_par"] = df["ltype"].str.startswith("P")
    first = df.groupby(["folio", "line"]).agg(
        pstart=("zl_pstart", "max"), pend=("zl_pend", "max"),
        ltype=("ltype", "first"), zl_pos=("zl_pos", "first")).reset_index()
    par_id, cur = [], -1
    prev_folio, prev_end = None, True
    for r in first.itertuples():
        new = (r.folio != prev_folio) or prev_end or r.pstart or r.zl_pos in "@*"
        if new:
            cur += 1
        par_id.append(cur)
        prev_folio, prev_end = r.folio, bool(r.pend)
    first["par"] = par_id
    first["line_in_par"] = first.groupby("par").cumcount()
    first["par_nlines"] = first.groupby("par")["line"].transform("size")
    first["par_first"] = first["line_in_par"] == 0
    first["par_last"] = first["line_in_par"] == first["par_nlines"] - 1
    df = df.merge(first[["folio", "line", "par", "line_in_par", "par_nlines",
                         "par_first", "par_last"]], on=["folio", "line"], how="left")
    return df


def tokens(refresh: bool = False) -> pd.DataFrame:
    CACHE.mkdir(exist_ok=True)
    p = CACHE / "tokens.pkl"
    if p.exists() and not refresh:
        return pd.read_pickle(p)
    df = build_tokens()
    df.to_pickle(p)
    return df


def paragraph_text(df: pd.DataFrame, min_frac: float = 0.0) -> pd.DataFrame:
    """Running paragraph text only (ZL locus type P*), lines where both line definitions agree."""
    m = df["is_par"] & df["line_agree"] & (df["ltype"] != "Pt")
    return df[m].copy()


if __name__ == "__main__":
    df = tokens(refresh=True)
    print(df.shape, df["folio"].nunique(), "folios")
    print("aligned to a ZL locus:", (df.zl_locus >= 0).mean().round(4),
          " mean char-align frac:", df.zl_frac.mean().round(4))
    print("line definitions agree:", df.line_agree.mean().round(4))
    print(df["ltype"].value_counts().head(12).to_dict())
    pt = paragraph_text(df)
    print("paragraph tokens:", len(pt), " lines:", pt.groupby(["folio", "line"]).ngroups,
          " paragraphs:", pt["par"].nunique())
    print(pt.groupby("lang").size().to_dict(), pt.groupby("illust").size().to_dict())
    print(pt.groupby(["folio", "line"]).size().describe().round(2).to_dict())
