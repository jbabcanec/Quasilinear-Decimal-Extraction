"""Line-level reading of HTR ground truth (ALTO or PAGE XML) and the three statistics.

The XML is treated as untrusted data: it is parsed with entity resolution and network
access switched off, and nothing from the datasets is ever imported or executed.

Units
  page    one XML file (after de-duplication, see run_all.py)
  column  one main-text block of a page; a block whose lines form two side-by-side
          stacks is split into two columns; lines are sorted top to bottom
  line    (xl, xr, y, text): the two ends and mean height of the baseline, and the text

S1  share of line pairs (lag 1, lag 2) inside a column whose first letter is the same,
    against the same lines in shuffled order inside the column.
S2  lines that stop early (right end below 90% of the way from their own left edge to
    the furthest right end among the +-3 neighbouring lines of the column): the share
    holding exactly one word, and the share shorter than 30%.
S3  share of neighbouring words inside a line that are identical, divided by its
    expectation when the words of each line are put in random order.
"""
from __future__ import annotations

import hashlib
import re
import unicodedata
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from lxml import etree

PARSER = dict(resolve_entities=False, no_network=True, huge_tree=False, recover=True,
              remove_comments=True, remove_pis=True, load_dtd=False, dtd_validation=False)

MAIN_ZONE = re.compile(r"^(MainZone|MainText|Column_\d+|text|paragraph)\b", re.I)
BAD_ZONE = re.compile(r"margin|number|running|title|quire|stamp|seal|graphic|music|damage|dropcap|drop-cap"
                      r"|decoration|artefact|table|heading|header|footer|page-number|catch|signature"
                      r"|picture|inscription|caption|footnote|gloss|formula|figure", re.I)
MAIN_LINE = re.compile(r"^(DefaultLine|default|CustomLine|HalfLine|paragraph)\b", re.I)
BAD_LINE = re.compile(r"heading|rubric|dropcap|interlin|music|number|margin|tironian|catch|signature", re.I)


def local(el) -> str:
    return etree.QName(el).localname if isinstance(el.tag, str) else ""


def nums(s: str) -> list[float]:
    return [float(v) for v in re.findall(r"-?\d+(?:\.\d+)?", s or "")]


def parse(path: Path):
    """-> (kind, image name, zones); a zone is (label, [(line label, xl, xr, y, text, sig)])."""
    try:
        root = etree.parse(str(path), etree.XMLParser(**PARSER)).getroot()
    except Exception:
        return None
    if root is None:
        return None
    kind = local(root)
    if kind == "alto":
        return ("alto",) + _alto(root)
    if kind == "PcGts":
        return ("page",) + _page(root)
    if kind == "TEI":
        got = _tei(root)
        return ("tei",) + got if got[1] else None
    return None


def _geom(base: list[float], poly: list[float]):
    """Line ends from the baseline (preferred) or from the outline."""
    src = base if len(base) >= 4 else poly
    if len(src) < 4:
        return None
    xs, ys = src[0::2], src[1::2]
    if len(base) >= 4:
        return min(xs), max(xs), float(np.mean(ys))
    return min(xs), max(xs), float(max(ys))


def _alto(root):
    labels = {t.get("ID"): (t.get("LABEL") or "") for t in root.iter() if local(t) == "OtherTag"}

    def lab(el):
        refs = (el.get("TAGREFS") or "").split()
        got = [labels[r] for r in refs if r in labels]
        return got[0] if got else ""
    image = ""
    for el in root.iter():
        if local(el) == "fileName":
            image = (el.text or "").strip()
            break
    zones = []
    for block in root.iter():
        if local(block) != "TextBlock":
            continue
        lines = []
        for tl in block.iter():
            if local(tl) != "TextLine":
                continue
            base = nums(tl.get("BASELINE"))
            poly = []
            for p in tl.iter():
                if local(p) == "Polygon":
                    poly = nums(p.get("POINTS"))
                    break
            if len(base) < 4 and len(poly) < 4 and tl.get("HPOS") and tl.get("WIDTH"):
                x, w = float(tl.get("HPOS")), float(tl.get("WIDTH"))
                yb = float(tl.get("VPOS") or 0) + float(tl.get("HEIGHT") or 0)
                poly = [x, yb, x + w, yb]
            g = _geom(base, poly)
            text = " ".join(s.get("CONTENT", "") for s in tl.iter() if local(s) == "String").strip()
            if g is None or not text:
                continue
            lines.append((lab(tl), *g, text, tl.get("BASELINE") or ""))
        zones.append((lab(block), lines))
    return image, zones


def _page(root):
    image = ""
    zones = []
    for pg in root.iter():
        if local(pg) == "Page":
            image = pg.get("imageFilename") or ""
            break

    def struct(el):
        m = re.search(r"structure\s*\{[^}]*type:([^;}]+)", el.get("custom") or "")
        return (m.group(1).strip() if m else "") or (el.get("type") or "")
    for reg in root.iter():
        if local(reg) != "TextRegion":
            continue
        lines = []
        for tl in reg:                                   # direct children only (nested regions separately)
            if local(tl) != "TextLine":
                continue
            base, poly, text = [], [], None
            for ch in tl:
                if local(ch) == "Baseline":
                    base = nums(ch.get("points"))
                elif local(ch) == "Coords":
                    poly = nums(ch.get("points"))
                elif local(ch) == "TextEquiv":
                    for u in ch:
                        if local(u) == "Unicode":
                            text = (u.text or "").strip()
            if text is None:                              # word-level only
                ws = []
                for w in tl:
                    if local(w) == "Word":
                        for u in w.iter():
                            if local(u) == "Unicode" and u.text:
                                ws.append(u.text.strip())
                                break
                text = " ".join(ws)
            g = _geom(base, poly)
            if g is None or not text:
                continue
            lines.append((struct(tl), *g, text, " ".join("%g" % v for v in base)))
        zones.append((struct(reg), lines))
    return image, zones


TEI_SKIP = {"expan", "reg", "note", "corr", "figure", "fw"}     # editorial layers: keep what is on the page


def _tei(root):
    """TEI transcriptions exported from Transkribus (Koenigsfelden): every <lb facs="#id"/> opens
    a line whose outline is the <zone xml:id="id"> of the facsimile; a <p facs="#region"> without
    @rend or @type is a block of main text (dorsal notes and marginal numbers carry one)."""
    xid = "{http://www.w3.org/XML/1998/namespace}id"
    zone = {z.get(xid): nums(z.get("points")) for z in root.iter() if local(z) == "zone" and z.get(xid)}
    zones = []
    for body in root.iter():
        if local(body) != "body":
            continue
        for par in body.iter():
            if local(par) != "p" or not par.get("facs"):
                continue
            label = "MainZone" if not (par.get("rend") or par.get("type")) else "MarginTextZone:" + (
                par.get("rend") or par.get("type"))
            lines: list[list] = []

            def add(t):
                if t and lines:
                    lines[-1][1].append(t)

            def walk(el, skip):
                if local(el) == "lb" and not skip:
                    lines.append([(el.get("facs") or "").lstrip("#"), []])
                if not skip:
                    add(el.text)
                for ch in el:
                    walk(ch, skip or local(ch) in TEI_SKIP)
                    if not skip:
                        add(ch.tail)
            walk(par, False)
            out = []
            for zid, parts in lines:
                text = re.sub(r"\s+", " ", "".join(parts)).strip()
                g = _geom([], zone.get(zid, []))
                if g is None or not text:
                    continue
                xs, ys = zone[zid][0::2], zone[zid][1::2]
                out.append(("DefaultLine", min(xs), max(xs), float(np.mean(ys)), text,
                            " ".join("%g" % v for v in zone[zid])))
            zones.append((label, out))
    return "", zones


def tei_header(path: Path):
    """-> (language, year) of a Koenigsfelden TEI file, or None."""
    try:
        root = etree.parse(str(path), etree.XMLParser(**PARSER)).getroot()
    except Exception:
        return None
    lang, year = "", None
    for el in root.iter():
        k = local(el)
        if k == "textLang" and not lang:
            lang = (el.text or "").strip()
        elif k == "origDate" and year is None:
            m = re.search(r"\d{4}", el.get("when") or el.get("notBefore") or el.get("from") or el.text or "")
            year = int(m.group(0)) if m else None
        elif k == "body":
            break
    return lang, year


# ------------------------------------------------------------------------------ columns
def split_columns(lines):
    """One block -> list of columns.  A block is cut where the left edges jump by more than
    a fifth of its width and the lines to the left end before the lines to the right begin."""
    if len(lines) < 2:
        return [lines]
    lines = sorted(lines, key=lambda r: r[1])
    width = max(r[2] for r in lines) - lines[0][1]
    if width <= 0:
        return [lines]
    for i in range(1, len(lines)):
        if lines[i][1] - lines[i - 1][1] > 0.2 * width:
            left, right = lines[:i], lines[i:]
            if len(left) >= 2 and len(right) >= 2 and \
                    np.median([r[2] for r in left]) <= np.median([r[1] for r in right]) + 0.03 * width:
                return split_columns(left) + split_columns(right)
    return [lines]


def order_column(lines):
    """Sort top to bottom and join pieces of one written line that were segmented apart
    (same height, side by side)."""
    lines = sorted(lines, key=lambda r: r[3])
    if len(lines) < 3:
        return [list(r) for r in lines]
    ys = np.array([r[3] for r in lines])
    step = float(np.median(np.diff(ys))) if len(ys) > 2 else 0.0
    out = [list(lines[0])]
    for r in lines[1:]:
        p = out[-1]
        apart = min(r[2], p[2]) - max(r[1], p[1]) <= 0.1 * min(r[2] - r[1], p[2] - p[1])
        if step > 0 and r[3] - p[3] < 0.35 * step and apart:
            a, b = (p, r) if p[1] <= r[1] else (r, p)
            out[-1] = [p[0], min(p[1], r[1]), max(p[2], r[2]), (p[3] + r[3]) / 2, a[4] + " " + b[4], p[5] + "|" + r[5]]
        else:
            out.append(list(r))
    return out


def columns_of(zones, zone_ok, line_ok, min_lines=4):
    cols = []
    for zlabel, lines in zones:
        if not zone_ok(zlabel, lines):
            continue
        keep = [r for r in lines if line_ok(r[0])]
        for col in split_columns(keep):
            col = order_column(col)
            if len(col) >= min_lines:
                cols.append((zlabel, col))
    return cols


def default_zone_ok(unlabelled: bool = False):
    """Main text only.  `unlabelled`: also accept blocks that carry no label at all
    (datasets without a zone vocabulary) provided they are not small."""
    def ok(label, lines):
        if label and BAD_ZONE.search(label):
            return False
        if label and MAIN_ZONE.search(label):
            return True
        return bool(unlabelled and not label and len(lines) >= 8)
    return ok


def default_line_ok(unlabelled: bool = True):
    def ok(label):
        if label and BAD_LINE.search(label):
            return False
        if label and MAIN_LINE.search(label):
            return True
        return bool(unlabelled and not label)
    return ok


# --------------------------------------------------------------------------------- text
KEEP = {"⁊", "&"}                       # the Tironian et and the ampersand are words, not punctuation


def is_sign(c: str) -> bool:
    """Letters, numerals, private-use (MUFI) characters, and the two signs for 'et'.
    Punctuation, paragraph marks and other symbols are not signs."""
    k = unicodedata.category(c)
    return c in KEEP or k[0] in "LN" or k == "Co"


def first_letter(text: str) -> str:
    """First letter of the line: lower-cased (long s = s), diacritics and abbreviation marks
    dropped, leading punctuation skipped."""
    for c in unicodedata.normalize("NFD", text):
        if is_sign(c):
            return c.casefold()[:1]
    return ""


def words(text: str) -> list[str]:
    """Words of a line: lower-cased, punctuation removed, abbreviation marks kept."""
    out = []
    for tok in unicodedata.normalize("NFC", text).casefold().split():
        w = "".join(c for c in tok if is_sign(c) or unicodedata.category(c)[0] == "M")
        if w:
            out.append(w)
    return out


def geom_signature(cols_or_zones) -> str:
    """Identifies a page by the coordinates of its lines, so that copies of one page that
    differ only in their characters (normalised re-exports, ALTO next to PAGE) are found."""
    sig = sorted(re.sub(r"[^0-9]+", " ", r[5]).strip() for _, lines in cols_or_zones for r in lines)
    return hashlib.sha1("\n".join(sig).encode()).hexdigest()


def text_signature(cols) -> str:
    t = "\n".join("".join(words(r[4])) for _, col in cols for r in col)
    return hashlib.sha1(t.encode()).hexdigest()


# --------------------------------------------------------------------------- statistics
def s1(seqs, lags=(1, 2), shuffles=1000, seed=0):
    """seqs: one sequence of first letters per column (reading order).
    Columns of equal length are stacked so that every shuffle is one array operation."""
    rng = np.random.default_rng(seed)
    code: dict[str, int] = {}
    by_len = defaultdict(list)
    for s in seqs:
        s = [code.setdefault(c, len(code)) for c in s if c]
        if len(s) >= 2:
            by_len[len(s)].append(s)
    stacks = {n: np.array(v, dtype=np.int32) for n, v in by_len.items()}
    out = {}
    for lag in lags:
        n = sum(a.shape[0] * (L - lag) for L, a in stacks.items() if L > lag)
        if n < 30:
            continue
        obs = sum(int((a[:, lag:] == a[:, :-lag]).sum()) for L, a in stacks.items() if L > lag)
        exact = 0.0                      # expectation under shuffling, in closed form (check on the shuffles)
        for L, a in stacks.items():
            if L > lag:
                for row in a:
                    k = np.bincount(row)
                    exact += (L - lag) * float((k * (k - 1)).sum()) / (L * (L - 1))
        null = np.zeros(shuffles)
        for L, a in stacks.items():
            if L <= lag:
                continue
            chunk = max(1, int(2e6 // (a.size + 1)))
            done = 0
            while done < shuffles:
                m = min(chunk, shuffles - done)
                idx = np.argsort(rng.random((m,) + a.shape), axis=-1)
                p = np.take_along_axis(np.broadcast_to(a, (m,) + a.shape), idx, axis=-1)
                null[done:done + m] += (p[..., lag:] == p[..., :-lag]).sum(axis=(1, 2))
                done += m
        sd = float(null.std())
        out[f"lag{lag}"] = dict(n=int(n), same=obs / n, shuffled=float(null.mean()) / n, sd=sd / n, expected=exact / n,
                                z=float((obs - null.mean()) / sd) if sd > 0 else float("nan"))
    return out


def s2(columns, half=3, early=0.9, short=0.3):
    """columns: lists of (xl, xr, n words).  Returns counts over lines that stop early.
    `fifths`: how far the early lines run (below 20%, 20-40, 40-60, 60-80, 80-90% of the way).
    `clear_*`: the same count restricted to lines that stop below 80%, which leaves out the
    ordinary unevenness of a right margin that was not justified."""
    n_lines = n_early = one = under = two = clear = clear_one = 0
    hist = [0, 0, 0, 0, 0]
    for col in columns:
        xr = np.array([r[1] for r in col], float)
        for i, (xl, e, nw) in enumerate(col):
            idx = [j for j in range(max(0, i - half), min(len(col), i + half + 1)) if j != i]
            if not idx:
                continue
            m1 = xr[idx].max()
            if m1 <= xl:
                continue
            n_lines += 1
            fill = (e - xl) / (m1 - xl)
            if 0 < fill < early:
                n_early += 1
                one += nw == 1
                two += nw <= 2
                under += fill < short
                hist[min(4, int(fill / 0.2))] += 1
                if fill < 0.8:
                    clear += 1
                    clear_one += nw == 1
    return dict(lines=n_lines, early=n_early, early_share=n_early / n_lines if n_lines else None,
                one_word=one / n_early if n_early else None, one_word_n=int(one),
                up_to_two=two / n_early if n_early else None,
                under_30=under / n_early if n_early else None, under_30_n=int(under),
                fifths=[round(h / n_early, 4) for h in hist] if n_early else None,
                clear=clear, clear_share=clear / n_lines if n_lines else None,
                clear_one_word=clear_one / clear if clear else None)


def s3(lines, mc_shuffles=0, seed=0):
    """lines: lists of words.  The shuffled rate is exact: for a line of n words in which
    word w occurs c_w times, a random order has sum_w c_w(c_w-1)/n identical neighbours."""
    pairs = ident = 0
    expect = 0.0
    for ln in lines:
        n = len(ln)
        if n < 2:
            continue
        pairs += n - 1
        ident += sum(a == b for a, b in zip(ln[:-1], ln[1:]))
        expect += sum(c * (c - 1) for c in Counter(ln).values()) / n
    out = dict(pairs=pairs, identical=ident / pairs if pairs else None, identical_n=int(ident),
               shuffled=expect / pairs if pairs else None, expected_n=float(expect),
               ratio=float(ident / expect) if expect > 0 else None)
    if mc_shuffles:
        rng = np.random.default_rng(seed)
        tot = 0
        for _ in range(mc_shuffles):
            for ln in lines:
                if len(ln) >= 2:
                    p = rng.permutation(len(ln))
                    tot += sum(ln[i] == ln[j] for i, j in zip(p[:-1], p[1:]))
        out["shuffled_mc"] = tot / mc_shuffles / pairs if pairs else None
    return out


def rhyme(cols_last2, lags=(1, 2), shuffles=200, seed=0):
    """How often lines 1 or 2 apart end on the same two letters, against shuffled order.
    Used only to tell verse from prose."""
    r = s1(cols_last2, lags=lags, shuffles=shuffles, seed=seed)
    return {k: (v["same"], v["shuffled"]) for k, v in r.items()}
