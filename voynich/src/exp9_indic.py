"""Experiment 9 - Sanskrit texts, spells (dharani) and Vedic recitation patterns.
Writes results/exp9_indic.json.

Same eight-number fingerprint as exp8, placed in the same standardised space (the eighteen
Latin texts), so distances to the Voynich are comparable with that table.

The recitation patterns are generated from the Rigveda's word-by-word text (padapatha) by
the traditional rules, for consecutive words 1 2 3 ... of a verse:
    krama   1 2 | 2 3 | 3 4 ...
    jata    1 2 2 1 1 2 | 2 3 3 2 2 3 ...
    ghana   1 2 2 1 1 2 3 3 2 1 1 2 3 | 2 3 3 2 2 3 4 4 3 2 2 3 4 ...
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from exp8_genres import FEATS, profile
from vdata import ROOT

G = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/user/tokushige-koyasan/gretil-corpus")
TEXTS = {
    "Rigveda hymns": ("Vedic hymns", ["1_veda/1_sam/sa_Rgveda-edAufrecht.txt"]),
    "Rigveda, word by word": ("Vedic hymns, words separated", ["1_veda/1_sam/sa_RgvedasaMhitApadapATha.txt"]),
    "Samaveda": ("Vedic chant", ["1_veda/1_sam/sa_sAmavedasaMhitA.txt"]),
    "Atharvaveda (Paippalada)": ("charms and spells", ["1_veda/1_sam/sa_paippalAdasaMhitA.txt"]),
    "Mahamayuri": ("Buddhist spell text", ["4_rellit/buddh/sa_mahAmAyUrIvidyArAjJI.txt"]),
    "Buddhist dharanis, pooled": ("spell texts with their prose", ["4_rellit/buddh/*dhAraNI*.txt", "4_rellit/buddh/sa_mahAmantrAnusAriNI.txt"]),
    "Buddhist hymns of praise": ("stotras", ["4_rellit/buddh/sa_108-buddhist-stotras.txt"]),
    "Ashtangahridaya": ("medicine", ["6_sastra/7_ayur/sa_vAgbhaTa-aSTAGgahRdayasUtra.txt"]),
    "Brihajjataka": ("astrology", ["6_sastra/8_jyot/sa_varAhamihira-bRhajjAtaka.txt"]),
    "Agni Purana": ("religious encyclopedia", ["3_purana/sa_agnipurANa.txt"]),
}


def words_of(text):
    text = re.sub(r"//[^/\n]*//|\([^)\n]*\)|-RV_[^-\n]*-|[A-Za-z]+_[\d,.:/]+", " ", text)
    text = text.replace("-", "")
    return [w for w in re.findall(r"[^\W\d_]+", text.lower()) if w]


def read(globs):
    out = []
    for g in globs:
        for f in sorted(G.glob(g)):
            out += words_of(f.read_text(errors="replace"))
    return out


def verses_padapatha(limit=4000):
    raw = (G / "1_veda/1_sam/sa_RgvedasaMhitApadapATha.txt").read_text(errors="replace")
    vs = []
    for chunk in re.split(r"//[^/\n]*//", raw):
        w = [x for x in (words_of(p) for p in chunk.split("|")) if x]
        w = ["".join(x) for x in w]
        if 4 <= len(w) <= 40:
            vs.append(w)
        if len(vs) >= limit:
            break
    return vs


def recite(verses, mode):
    out = []
    for w in verses:
        line = []
        for i in range(len(w) - 1):
            a, b = w[i], w[i + 1]
            if mode == "krama":
                line += [a, b]
            elif mode == "jata":
                line += [a, b, b, a, a, b]
            elif mode == "ghana" and i + 2 < len(w):
                c = w[i + 2]
                line += [a, b, b, a, a, b, c, c, b, a, a, b, c]
        out += line
    return out


def main():
    rng = np.random.default_rng(0)
    base = json.loads((ROOT / "results" / "exp8_genres.json").read_text())
    M = np.array([[base["real"][t][f] for f in FEATS] for t in base["real"]])
    mu, sd = M.mean(0), M.std(0)
    six = [i for i, f in enumerate(FEATS) if f not in ("repeat", "near")]
    z = lambda p: (np.array([p[f] for f in FEATS]) - mu) / sd
    zv = z(base["other"]["Voynich"])
    rows = {}
    for label, (kind, globs) in TEXTS.items():
        w = read(globs)
        if len(w) < 7500:
            print("skipped (too short):", label, len(w)); continue
        rows[label] = profile([w[i:i + 9] for i in range(0, len(w), 9)], rng) | {"kind": kind, "words": len(w)}
    vs = verses_padapatha()
    for mode, name in (("krama", "Rigveda recited: krama"), ("jata", "Rigveda recited: jata"),
                       ("ghana", "Rigveda recited: ghana")):
        w = recite(vs, mode)
        rows[name] = profile([w[i:i + 9] for i in range(0, len(w), 9)], rng) | {
            "kind": "recitation pattern, generated", "words": len(w)}
    for r in rows.values():
        d = z(r) - zv
        r["dist6"] = float(np.linalg.norm(d[six])); r["dist8"] = float(np.linalg.norm(d))
    rows["Voynich"] = base["other"]["Voynich"] | {"kind": "", "dist6": 0.0, "dist8": 0.0}
    (ROOT / "results" / "exp9_indic.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False))
    print("%-28s %-30s %7s %7s %6s %8s %6s %6s" % ("text", "kind", "repeat", "near", "order", "phrases", "dist6", "dist8"))
    for t, r in sorted(rows.items(), key=lambda kv: kv[1]["dist8"]):
        print("%-28s %-30s %6.1f%% %6.1f%% %6.2f %8.1f %6.1f %6.1f" % (
            t, r["kind"], 100 * r["repeat"], 100 * r["near"], r["order"], r["phrases"], r["dist6"], r["dist8"]))


if __name__ == "__main__":
    main()


def mantra_strings():
    """The spell proper: text between the standard opening 'tadyatha' and closing 'svaha'
    in the Buddhist spell texts.  Compared with Voynich windows of the same length."""
    from exp6_llm_context import load
    from exp8_genres import features
    rng = np.random.default_rng(0)
    raw = ""
    for g in ("4_rellit/buddh/*dhAraNI*.txt", "4_rellit/buddh/sa_mahAmantrAnusAriNI.txt",
              "4_rellit/buddh/sa_mahAmAyUrIvidyArAjJI.txt"):
        for f in sorted(G.glob(g)):
            raw += "\n" + f.read_text(errors="replace")
    spans = re.findall(r"tadyath[āa]\b(.{20,2500}?)sv[āa]h[āa]", raw, flags=re.S | re.I)
    w = [x for s in spans for x in words_of(s)]
    lines = [w[i:i + 9] for i in range(0, len(w), 9)]
    m = features(lines, rng)
    vl = [ln for page in load("voynich") for ln in page]
    vw, cur, cnt = [], [], 0
    for ln in vl:
        cur.append(ln); cnt += len(ln)
        if cnt >= len(w):
            vw.append(features(cur, rng)); cur, cnt = [], 0
    v = {k: float(np.nanmean([f[k] for f in vw])) for k in FEATS}
    out = dict(spans=len(spans), words=len(w), mantra=m, voynich_same_length=v, sample=" ".join(w[:60]))
    p = ROOT / "results" / "exp9_indic.json"
    d = json.loads(p.read_text()); d["mantra_strings"] = out
    p.write_text(json.dumps(d, indent=1, ensure_ascii=False))
    print("mantra strings: %d spans, %d words" % (len(spans), len(w)))
    print("sample:", out["sample"][:300])
    print("%-22s " % "" + " ".join("%8s" % f for f in FEATS))
    print("%-22s " % "spell strings" + " ".join("%8.3f" % m[f] for f in FEATS))
    print("%-22s " % "Voynich, same length" + " ".join("%8.3f" % v[f] for f in FEATS))


if __name__ == "__main__":
    mantra_strings()
