"""How much of the simulator's freedom did real medieval Latin copyists use?

Reads line-by-line graphematic transcriptions (abbreviation signs kept as written) of
Latin manuscripts of the 12th-16th century: HTR-United/CREMMA-Medieval-LAT (ALTO XML,
main text zones).  For each manuscript:

  abbreviated words      share of words carrying at least one abbreviation sign
  signs per word         and, at about 1.5 letters saved per sign, the share of letters saved
  at the line end        the same share for the last word of a line against the others
                         (do scribes abbreviate more to fit the line?)
  line-to-line swing     sd across lines of signs-per-word beyond what scattering the
                         same words at random over the page's lines gives, turned into a
                         line-width swing (1.5 letters per sign): the abbreviation
                         analogue of the simulator's "density sd"
  split words            share of line breaks that fall inside a word (lower bound: the
                         two halves join into a known Latin word and are not both words)
  first glyph            share of lines that open with a sign Latin spelling does not force
                         (tironian et, con-sign, p/q with stroke ...), and R1 on line openings

    python3 calibrate_real_mss.py [dataset root]     -> results/real_mss.json
"""
from __future__ import annotations

import json
import re
import sys
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
NS = "{http://www.loc.gov/standards/alto/ns-v4#}"
ROOT = Path("/home/user/HTR-United/CREMMA-Medieval-LAT")
LATIN = Path("/home/user/cltk/lat_text_latin_library")
SIGNS = set("ꝑꝓꝗꝙꝯꝰ⁊&ꝝꝵħłꝉẜ÷ꝫꝭꝓꝕꝧꝟ")
PUNCT = re.compile(r"[.,;:!?·⸱\-/()\[\]'\"¶]+")


def nsigns(word: str) -> int:
    return sum(1 for c in word if c in SIGNS or unicodedata.category(c) == "Mn" or 0x363 <= ord(c) <= 0x36F
               or 0x1DC0 <= ord(c) <= 0x1DFF)


def plain(word: str) -> str:
    return "".join(c for c in unicodedata.normalize("NFD", word.lower()) if "a" <= c <= "z")


def pages(ms_dir: Path):
    files = sorted(ms_dir.glob("*.chocomufin.xml")) or sorted(ms_dir.glob("*.xml"))
    for f in files:
        root = ET.parse(f).getroot()
        labels = {t.get("ID"): t.get("LABEL") for t in root.iter(f"{NS}OtherTag")}
        for block in root.iter(f"{NS}TextBlock"):
            if labels.get(block.get("TAGREFS", ""), "") != "MainZone":
                continue
            lines = []
            for tl in block.iter(f"{NS}TextLine"):
                if labels.get(tl.get("TAGREFS", ""), "DefaultLine") != "DefaultLine":
                    continue
                txt = " ".join(s.get("CONTENT", "") for s in tl.iter(f"{NS}String"))
                words = [w for w in PUNCT.sub(" ", txt).split() if w]
                if len(words) >= 3:
                    lines.append(words)
            if len(lines) >= 8:
                yield lines


def lexicon(limit=400):
    """Latin word forms from the Latin Library (plain letters, u for v, i for j)."""
    lex = Counter()
    files = sorted(LATIN.rglob("*.txt"))[:limit]
    for f in files:
        for w in re.findall(r"[a-z]+", f.read_text(encoding="utf-8", errors="ignore").lower()):
            lex[w.replace("v", "u").replace("j", "i")] += 1
    return {w for w, c in lex.items() if c >= 2}


def measure(ms_dir: Path, lex, rng):
    P = list(pages(ms_dir))
    if sum(map(len, P)) < 80:
        return None
    last, other, per_line, null_sd, splits, breaks, openers = [], [], [], [], 0, 0, []
    glyphs = []
    for lines in P:
        s = [[nsigns(w) for w in ln] for ln in lines]
        for ln in s:
            last.append(ln[-1] > 0); other.extend(x > 0 for x in ln[:-1])
        m = np.array([np.mean(ln) for ln in s]); per_line.append(m - m.mean())
        flat = np.array([x for ln in s for x in ln]); sizes = [len(ln) for ln in s]
        sd = []
        for _ in range(50):
            f = rng.permutation(flat); k = 0; mm = []
            for z in sizes:
                mm.append(f[k:k + z].mean()); k += z
            sd.append(np.std(mm))
        null_sd.append(np.mean(sd) ** 2 * len(lines))
        for a, b in zip(lines[:-1], lines[1:]):
            x, y = plain(a[-1]).replace("v", "u").replace("j", "i"), plain(b[0]).replace("v", "u").replace("j", "i")
            breaks += 1
            if len(x) >= 2 and len(y) >= 2 and nsigns(a[-1]) == 0 and nsigns(b[0]) == 0:
                splits += (x + y) in lex and not (x in lex and y in lex)
        seq = [unicodedata.normalize("NFD", ln[0])[0].lower() for ln in lines]
        glyphs.append(seq)
        openers.extend(ln[0][0] in SIGNS for ln in lines)
    allw = np.array(last + other)
    signs = np.mean([nsigns(w) for lines in P for ln in lines for w in ln])
    chars = np.mean([len(plain(w)) + sum(c in SIGNS for c in w) for lines in P for ln in lines for w in ln])
    obs_var = np.mean(np.concatenate(per_line) ** 2)
    n_lines = sum(map(len, P))
    coherent = np.sqrt(max(0.0, obs_var - sum(null_sd) / n_lines))
    same = sum(a == b for s in glyphs for a, b in zip(s[:-1], s[1:])); n = sum(len(s) - 1 for s in glyphs)
    exp = sum((len(s) - 1) * sum(c * (c - 1) for c in Counter(s).values()) / (len(s) * (len(s) - 1)) for s in glyphs)
    return dict(lines=n_lines, words_per_line=float(np.mean([len(ln) for lines in P for ln in lines])),
                abbreviated_words=float(allw.mean()), signs_per_word=float(signs), glyphs_per_word=float(chars),
                letters_saved_share=float(1.5 * signs / (chars + 1.5 * signs)),
                abbreviated_last_word=float(np.mean(last)), abbreviated_other_words=float(np.mean(other)),
                line_swing_signs_per_word=float(coherent),
                line_swing_width=float(1.5 * coherent / (chars + 1.5 * signs)),
                split_breaks_lower_bound=splits / max(breaks, 1),
                lines_opening_with_a_sign=float(np.mean(openers)),
                r1_same=same / n, r1_shuffled=exp / n)


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT
    lex = lexicon()
    rng = np.random.default_rng(0)
    reg = {}
    for row in (root / "data-registry.csv").read_text(encoding="utf-8").splitlines()[1:]:
        cells = row.split(",")
        reg[cells[0]] = row
    out = {}
    for ms in sorted(p for p in (root / "data").iterdir() if p.is_dir()):
        r = measure(ms, lex, rng)
        if r:
            out[ms.name] = r
    json.dump(out, open(HERE / "results" / "real_mss.json", "w"), indent=1)
    cols = ["lines", "words_per_line", "abbreviated_words", "letters_saved_share", "abbreviated_last_word",
            "abbreviated_other_words", "line_swing_width", "split_breaks_lower_bound",
            "lines_opening_with_a_sign", "r1_same", "r1_shuffled"]
    print("%-26s" % "manuscript" + "".join("%10s" % c[:9] for c in cols))
    for ms, r in out.items():
        print("%-26s" % ms[:25] + "".join("%10.3f" % r[c] if isinstance(r[c], float) else "%10d" % r[c] for c in cols))
    tot = sum(r["lines"] for r in out.values())
    avg = lambda k: sum(r[k] * r["lines"] for r in out.values()) / tot
    print("%-26s" % "all (line-weighted)" + "".join("%10.3f" % avg(c) if c != "lines" else "%10d" % tot for c in cols))
    print("median over manuscripts: abbreviated words %.2f, letters saved %.2f, line swing (width) %.3f, splits >= %.2f" % (
        np.median([r["abbreviated_words"] for r in out.values()]), np.median([r["letters_saved_share"] for r in out.values()]),
        np.median([r["line_swing_width"] for r in out.values()]), np.median([r["split_breaks_lower_bound"] for r in out.values()])))


if __name__ == "__main__":
    main()
