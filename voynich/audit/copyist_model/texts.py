"""Real Latin texts for the copyist model, as paragraphs of words.

Source: the Latin Library plain-text files (cltk/lat_text_latin_library).

Two ways of getting paragraphs:
  natural  - the blank-line blocks of the file (chapters, recipes, sections), kept when they
             are prose of a usable length;
  matched  - the text is cut at sentence ends (verse numbers for the Vulgate, which has no
             punctuation) into paragraphs whose target length is drawn to give 4-8 lines of
             about nine words, the Voynich range.

Words are lower-cased and stripped to letters; numerals, section numbers, headings and the
Latin Library header and footer are dropped.  Punctuation takes no width.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

LATIN = Path("/home/user/cltk/lat_text_latin_library")

TEXTS = {
    "caesar": ["caesar/gall1.txt", "caesar/gall2.txt", "caesar/gall3.txt", "caesar/gall4.txt",
               "caesar/gall5.txt", "caesar/gall6.txt", "caesar/gall7.txt"],
    "pliny": ["pliny.nh2.txt", "pliny.nh3.txt", "pliny.nh4.txt", "pliny.nh5.txt"],
    "genesis": ["bible/genesis.txt"],
    "recipes": ["apicius/apicius1.txt", "apicius/apicius2.txt", "apicius/apicius3.txt",
                "apicius/apicius4.txt", "apicius/apicius5.txt", "cato/cato.agri.txt"],
    "isidore_plants": ["isidore/17.txt"],
}

FOOTER = re.compile(r"the latin library|the classics page|^\s*(caesar|pliny the elder|apicius|"
                    r"cato|isidore|vulgate|christian latin|medieval latin)\s*$", re.I)
WORD = re.compile(r"[a-z]+")
END = re.compile(r"[.?!;:]$")


def _blocks(path: Path):
    """Blank-line separated blocks of a Latin Library file, header and footer removed."""
    raw = path.read_text(encoding="utf-8", errors="ignore").replace("\r", "")
    raw = raw.replace("æ", "ae").replace("œ", "oe").replace("Æ", "Ae").replace("Œ", "Oe")
    blocks, cur = [], []
    for line in raw.split("\n"):
        if line.strip():
            cur.append(line.strip())
        elif cur:
            blocks.append(" ".join(cur)); cur = []
    if cur:
        blocks.append(" ".join(cur))
    out = []
    for i, b in enumerate(blocks):
        if FOOTER.search(b):
            continue
        toks = b.split()
        letters = [t for t in toks if re.search(r"[A-Za-z]", t)]
        if not letters:
            continue                                   # navigation rows of numbers
        if i < 3 and len(letters) < 25:
            continue                                   # file title
        if sum(t.isupper() for t in letters) > 0.6 * len(letters):
            continue                                   # headings in capitals
        out.append(b)
    return out


def _sentences(block: str, verse: bool):
    """-> list of sentences, each a list of clean lower-case words."""
    block = re.sub(r"\[[^\]]{0,12}\]", " ", block)      # [ 12 ] section numbers
    block = re.sub(r"\([^)]{0,40}\)", " ", block)       # (Metam. 5,341) references
    sents, cur = [], []
    for tok in block.split():
        if re.fullmatch(r"[\d.,:;]+", tok):
            if verse and cur:
                sents.append(cur); cur = []
            continue
        if tok.isupper() and len(tok) > 1 and re.fullmatch(r"[IVXLCDM.]+", tok):
            continue                                   # Roman numerals
        w = "".join(WORD.findall(tok.lower()))
        if not w:
            continue
        cur.append(w)
        if not verse and END.search(tok) and len(w) > 2:   # 'C.' 'p.' are not sentence ends
            sents.append(cur); cur = []
    if cur:
        sents.append(cur)
    return sents


def load(name: str, mode: str = "matched", seed: int = 0, lines=(4, 8), words_per_line: float = 9.0,
         max_words: int | None = 40000):
    """-> list of paragraphs, each a list of words."""
    verse = name == "genesis"
    rng = np.random.default_rng(seed)
    paras = []
    stream = []                                        # sentences, for mode='matched'
    for f in TEXTS[name]:
        for b in _blocks(LATIN / f):
            s = _sentences(b, verse)
            if not s:
                continue
            if mode == "natural":
                w = [x for snt in s for x in snt]
                if 3 * words_per_line <= len(w) <= 30 * words_per_line:
                    paras.append(w)
            else:
                stream.extend(s)
    if mode != "natural":
        cur, target = [], None
        for s in stream:
            if target is None:
                target = rng.integers(lines[0], lines[1] + 1) * words_per_line - rng.uniform(0, words_per_line)
            cur.extend(s)
            if len(cur) >= target:
                paras.append(cur); cur, target = [], None
        if len(cur) >= 2 * words_per_line:
            paras.append(cur)
    if max_words:
        out, n = [], 0
        for p in paras:
            if n >= max_words:
                break
            out.append(p); n += len(p)
        paras = out
    return paras


if __name__ == "__main__":
    for nm in TEXTS:
        for mode in ("matched", "natural"):
            P = load(nm, mode)
            n = sum(map(len, P))
            L = np.array([len(p) for p in P]) if P else np.array([0])
            ml = np.mean([len(w) for p in P for w in p]) if P else 0
            print("%-15s %-8s paragraphs %4d  words %6d  words/par median %5.0f (10-90%%: %3.0f-%3.0f)  "
                  "letters/word %.2f" % (nm, mode, len(P), n, np.median(L), np.percentile(L, 10),
                                         np.percentile(L, 90), ml))
            if P:
                print("     ", " ".join(P[min(3, len(P) - 1)][:22]))
