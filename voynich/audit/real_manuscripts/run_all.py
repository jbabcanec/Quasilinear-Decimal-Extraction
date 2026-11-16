"""Three Voynich statistics on real manuscripts (HTR ground truth cloned from GitHub).

Usage:  python3 -I run_all.py [--shuffles 1000] [--only <dataset substring>] [--out results]
Writes  results.csv / results.json   one row per manuscript
        pages.csv                    every XML file read, and whether it was used or why not
        pooled.json                  pooled figures per language and form

De-duplication (three passes, logged in pages.csv):
  1. same page key inside a manuscript: the file name without the suffixes .chocomufin and
     .mufichecker (re-exports of one page with normalised characters), and without the
     directory when a repository ships one page as both ALTO and PAGE;
  2. same line coordinates anywhere in the run (the baselines of all main-text lines);
  3. same main-text transcription anywhere in the run.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mslib as M                                                            # noqa: E402

H = Path("/home/user")
SKIP_FILE = re.compile(r"^(mets|metadata|METS)\.xml$|teiMinimal|Masterfile", re.I)


def D(name, root, **kw):
    kw.setdefault("ms", "dir")          # 'dir': one manuscript per directory holding XML files
    kw.setdefault("unlabelled", False)  # accept blocks without a zone label as main text
    kw.setdefault("rtl", False)
    kw.setdefault("prefer", ())         # path substrings preferred when one page has several files
    kw.setdefault("exclude", ())        # path substrings to leave out entirely
    kw.setdefault("paragraphs", False)  # every block is one paragraph (charters, entries)
    return dict(name=name, root=H / root, **kw)


DATASETS = [
    D("Gallicorpora-MSS-15", "Gallicorpora/HTR-MSS-15e-Siecle/data"),
    D("CREMMA-Medieval", "HTR-United/cremma-medieval/data"),
    D("CREMMA-Medieval-LAT", "HTR-United/CREMMA-Medieval-LAT/data"),
    D("HTRomance-French", "HTRomance-Project/medieval-french/data"),
    D("HTRomance-Latin", "HTRomance-Project/medieval-latin/data"),
    D("HTRomance-Italian", "HTRomance-Project/medieval-italian/data"),
    D("HTRomance-Spain", "HTRomance-Project/middle-ages-in-spain/data"),
    D("HTRomance-Occitan", "HTRomance-Project/medieval-occitan/data", unlabelled=True),
    D("HTRogene-French", "HTRogene/french/data", unlabelled=True),
    D("HTRogene-Latin", "HTRogene/latin/data"),
    D("HTRogene-Italian", "HTRogene/italian/data"),
    D("HTRogene-Spanish", "HTRogene/spanish/data"),
    D("HTRogene-Occitan", "HTRogene/occitan/data", unlabelled=True),
    D("Gatmuza-Occitan", "LostMa-ERC/gatmuza", unlabelled=True),
    D("CIHAM-Fabliaux", "CIHAM-HTR/Fabliaux/data"),
    D("CIHAM-Liber", "CIHAM-HTR/Liber/data"),
    D("DecameronFR", "PSL-Chartes-HTR-Students/TNAH-2021-DecameronFR/veriteTerrain", ms="one", unlabelled=True),
    D("GrandesChroniques", "malamatenia/palaeographic-variability-analysis-grandes-chroniques-fr-2813/data/raw_ground_truth/annotations", ms="one"),
    D("Torino-L-II-14", "RESCAPE-Biblissima/Torino_L-II-14/data/transcriptions", ms=r"^(IT-TO\d+_MSL_II_14_Tome_\d+)", prefer=("/alto/",)),
    D("FoNDUE-LA-MSS-MA", "FoNDUE-HTR/FONDUE-LA-MSS-MA/data", unlabelled=True),
    D("ParisBible", "parisbible/ground_truth/PBP 1.0/ALTO", ms="one", unlabelled=True),
    D("Vienna-2023-Latin", "HTR-School-Vienna/2023--late-medieval-latin", ms="one", unlabelled=True),
    D("Vienna-2024-Latin", "HTR-School-Vienna/2024--late-medieval-latin", ms="top", unlabelled=True, prefer=("/final/",)),
    D("Vienna-2025-Latin", "HTR-School-Vienna/2025--late-medieval-latin", ms="top", unlabelled=True),
    D("Vienna-2024-German", "HTR-School-Vienna/2024--German-Group", ms=r"^(.*?)[-_ ]?\d+\.xml$", unlabelled=True),
    D("Vienna-2025-German", "HTR-School-Vienna/2025--early-modern-german", ms="top", unlabelled=True),
    D("Vienna-2023-Czech", "HTR-School-Vienna/2023--medieval-czech", ms="one", unlabelled=True),
    D("Vienna-2025-Czech", "HTR-School-Vienna/2025--medieval-czech", ms="one", unlabelled=True),
    D("Vienna-2025-Hebrew", "HTR-School-Vienna/2025-hebrew", ms="one", unlabelled=False, rtl=True),
    D("Vienna-2024-Syriac", "HTR-School-Vienna/2024--Syriac", ms="one", rtl=True),
    D("Vienna-2025-Syriac", "HTR-School-Vienna/2025-syriac", ms="one", rtl=True, prefer=("/alto/",)),
    D("LiDi-Ligorio", "Giorgiaagostini/LiDi1.0-project/Data/Ground Truth", ms="one", unlabelled=True, prefer=("/page/",), paragraphs=True),
    D("Eutyches", "malamatenia/Eutyches", ms="top", exclude=("XML-XSLT", "kraken-YALTAi")),
    D("CarolineMinuscule", "rescribe/carolineminuscule-groundtruth", ms="one"),
    D("iForal", "Arch-W/iForal-Dataset/data", ms="one", unlabelled=True, paragraphs=True),
    D("Koenigsfelden", "Staatsarchiv-Aargau/koenigsfelden-data/data/docs", ms="koenigsfelden", paragraphs=True),
    D("CREMMA-MSS-16", "HTR-United/CREMMA-MSS-16/data", ms="one"),
    D("GasparoSardi", "PaulineJac/GasparoSardiToponomasia/HTR/Data", ms="one"),
    D("BibliaArabica", "biblia-arabica/academies/htr/ground-truth", ms="top", rtl=True, unlabelled=True, prefer=("/alto/",)),
]


REJECTED = {      # looked at and not used
    "michaelscho/bdd-segmentation-data": "layout only: 1 of 3,538 pages carries a transcription",
    "calfa-co/tarima": "Maghrebi Arabic, 16th-19th c., items undated, all but one under 200 lines",
    "calfa-co/datalab-dulaurier": "19th-century scholar's copies (Armenian)",
    "biblia-arabica/academies (htr/transcriptions)": "model output, not ground truth; only htr/ground-truth is read",
    "DenisaFB/ERC-Philiumm-data-annotation-and-HTR": "17th-century scientific drafts with formula zones",
    "Ioannes-Sapidus-1490-1561/HTR": "printed books",
    "PSL-Chartes-HTR-Students/HN2021-ChateauChavigny": "253 lines of a 16th-century inventory, no zone labels",
    "HTR-United/CREMMA-MSS-15, HTR-United/cremma-medieval-occ, vivianpl/hpgtr, crlsrmrlsz/paleo-ocr, "
    "CIHAM-HTR/SEHM_HTR_2024, Didip-eu/*": "no ALTO/PAGE files in the repository",
    "HTR-School-Vienna/2023--byzantine-greek, *--carolingian-latin": "TEI or text only, no line coordinates",
    "HTR-School-Vienna/2024--medieval-czech, Gallicorpora/HTR-incunable-*, HTR-imprime-*": "printed books",
}
FORMAT_DIR = {"alto", "page", "pagexml", "xml"}


def page_key(p: Path, root: Path) -> str:
    """Path below the dataset root, without format directories (alto/, page/) and without
    the suffixes of normalised re-exports: two files with one key are one page."""
    rel = p.relative_to(root)
    parts = [x.lower() for x in rel.parts[:-1] if x.lower() not in FORMAT_DIR]
    name = re.sub(r"\.(chocomufin|mufichecker)", "", rel.name, flags=re.I).rsplit(".xml", 1)[0].lower()
    return "/".join(parts + [name])


def manuscripts(ds):
    """-> {manuscript id: [xml files]} (before de-duplication)."""
    files = [p for p in sorted(ds["root"].rglob("*.xml")) if not SKIP_FILE.search(p.name)
             and not any(x in str(p) for x in ds["exclude"])]
    out = defaultdict(list)
    for p in files:
        rel = p.relative_to(ds["root"])
        if ds["ms"] == "one":
            key = ds["name"]
        elif ds["ms"] == "top":
            key = rel.parts[0] if len(rel.parts) > 1 else ds["name"]
        elif ds["ms"] == "koenigsfelden":                # pooled by language, century and kind of record
            h = M.tei_header(p)
            if h is None or h[1] is None or not h[0]:
                continue
            lang = {"Frühneuhochdeutsch": "German", "Latein": "Latin"}.get(h[0], h[0])
            key = "%s %dth c. %s" % (lang, (h[1] - 1) // 100 + 1, "charters" if p.name.startswith("U-") else "records")
        elif ds["ms"] == "dir":
            key = str(rel.parent) if len(rel.parts) > 1 else ds["name"]
        else:
            m = re.match(ds["ms"], p.name)
            key = m.group(1).strip() if m else ds["name"]
        out[re.sub(r"^Under Review - ", "", key)].append(p)   # same codex, pages not yet signed off
    return out


def rank(p: Path, prefer) -> tuple:
    s = str(p)
    return (0 if any(x in s for x in prefer) else 1, 0 if ".chocomufin" in p.name else 1, len(p.name), s)


def read_manuscript(ds, files, seen_geom, seen_text, log):
    """-> list of pages; a page is a list of columns; a column is a list of [label, xl, xr, y, text, sig]."""
    zone_ok, line_ok = M.default_zone_ok(ds["unlabelled"]), M.default_line_ok(True)
    by_key = defaultdict(list)
    for p in files:
        by_key[page_key(p, ds["root"])].append(p)
    pages = []
    for key in sorted(by_key):
        cands = sorted(by_key[key], key=lambda p: rank(p, ds["prefer"]))
        for p in cands[1:]:
            log.append((ds["name"], str(p), "duplicate: same page key as " + cands[0].name, 0))
        p = cands[0]
        got = M.parse(p)
        if got is None:
            log.append((ds["name"], str(p), "not ALTO/PAGE or unreadable", 0))
            continue
        kind, image, zones = got
        cols = M.columns_of(zones, zone_ok, line_ok)
        if ds["unlabelled"] and len(cols) > 1:           # no zone labels: drop narrow side blocks
            wid = [float(np.median([r[2] - r[1] for r in c])) for _, c in cols]
            cols = [zc for zc, w in zip(cols, wid) if zc[0] or w >= 0.5 * max(wid)]
        n = sum(len(c) for _, c in cols)
        if n == 0:
            log.append((ds["name"], str(p), "no main-text column of 4+ lines", 0))
            continue
        g, t = M.geom_signature(cols), M.text_signature(cols)
        if g in seen_geom:
            log.append((ds["name"], str(p), "duplicate: same line coordinates as " + seen_geom[g], 0))
            continue
        if t in seen_text:
            log.append((ds["name"], str(p), "duplicate: same transcription as " + seen_text[t], 0))
            continue
        seen_geom[g] = seen_text[t] = str(p)
        if ds["rtl"]:                                    # mirror, so that 'right end' is where lines end
            cols = [(z, [[r[0], -r[2], -r[1], r[3], r[4], r[5]] for r in c]) for z, c in cols]
        log.append((ds["name"], str(p), "used", n))
        pages.append((p, image, [c for _, c in cols]))
    return pages


def last2(text):
    w = M.words(text)
    s = re.sub(r"[^a-zà-ÿ]", "", w[-1]) if w else ""
    return s[-2:]


def overfull(cols, factor=1.8, min_chars=10):
    """Transcription check: share of lines carrying more than `factor` times the letters per
    unit of baseline length that is normal for their column.  Text that was predicted by a
    model and never corrected, or that belongs to another line, shows up here."""
    bad = n = 0
    for c in cols:
        d = [(len(re.sub(r"\s", "", r[4])), r[2] - r[1]) for r in c]
        dens = np.array([k / w for k, w in d if k >= min_chars and w > 0])
        if len(dens) < 4:
            continue
        n += len(dens)
        bad += int((dens > factor * np.median(dens)).sum())
    return bad / n if n else None


def par_final(cols, half=3):
    """Datasets whose blocks are paragraphs (one charter, one entry): the last line of each."""
    n = one = short = 0
    for c in cols:
        if len(c) < 4:
            continue
        xl, xr, nw = c[-1][1], c[-1][2], len(M.words(c[-1][4]))
        m1 = max(r[2] for r in c[-1 - half:-1])
        if m1 <= xl:
            continue
        n += 1
        one += nw == 1
        short += (xr - xl) / (m1 - xl) < 0.3
    return dict(n=n, one_word=one / n if n else None, under_30=short / n if n else None)


def measure(pages, shuffles, seed=0, paragraphs=False):
    cols = [c for _, _, pc in pages for c in pc]
    lines = [r for c in cols for r in c]
    wl = [M.words(r[4]) for r in lines]
    res = dict(pages=len(pages), columns=len(cols), lines=len(lines),
               words_per_line=float(np.mean([len(w) for w in wl])) if wl else None,
               chars_per_line=float(np.mean([sum(map(len, w)) for w in wl])) if wl else None)
    res["s1"] = M.s1([[M.first_letter(r[4]) for r in c] for c in cols], shuffles=shuffles, seed=seed)
    res["s2"] = M.s2([[(r[1], r[2], len(M.words(r[4]))) for r in c] for c in cols])
    res["s3"] = M.s3(wl)
    res["par_final"] = par_final(cols) if paragraphs else dict(n=None, one_word=None, under_30=None)
    res["rhyme"] = M.rhyme([[last2(r[4]) for r in c] for c in cols], seed=seed)
    res["capital_start"] = float(np.mean([(M.first_letter(r[4]) != "" and next(
        (ch for ch in r[4] if M.is_sign(ch)), "a").isupper()) for r in lines])) if lines else None
    res["sample"] = [r[4] for r in lines[len(lines) // 2: len(lines) // 2 + 4]]
    res["overfull_share"] = overfull(cols)
    return res


def auto_form(r) -> str:
    """Layout class from the lines themselves (a check on the catalogue labels, not a label):
    verse = one verse per line (rhyme between lines 1 or 2 apart, or capitals at line starts,
    together with a ragged right edge)."""
    rh = r.get("rhyme") or {}
    r1, s1_ = rh.get("lag1", (0, 0))
    r2, s2_ = rh.get("lag2", (0, 0))
    early = r["s2"]["early_share"] or 0
    cap = r["capital_start"] or 0
    rhymed = (r1 - s1_ >= 0.10) or (r2 - s2_ >= 0.10) or r1 >= 0.25
    if early >= 0.15 and (rhymed or cap >= 0.6):
        return "verse"
    if rhymed:
        return "verse?"
    if early >= 0.25:
        return "ragged"
    return "prose"


def flat(ds, ms, meta, r):
    def g(d, *ks):
        for k in ks:
            d = d.get(k) if isinstance(d, dict) else None
        return d
    row = dict(dataset=ds, manuscript=ms, shelfmark=meta.get("shelfmark", ""), language=meta.get("language", ""),
               century=meta.get("century", ""), form=meta.get("form", ""), hand=meta.get("hand", ""),
               pages=r["pages"], columns=r["columns"], lines=r["lines"],
               words_per_line=r["words_per_line"])
    for lag in ("lag1", "lag2"):
        for k in ("n", "same", "shuffled", "z"):
            row[f"s1_{lag}_{k}"] = g(r, "s1", lag, k)
    for k in ("lines", "early", "early_share", "one_word", "one_word_n", "under_30", "under_30_n"):
        row[f"s2_{k}"] = r["s2"][k]
    row["s2_one_word_per_1000_lines"] = 1000 * r["s2"]["one_word_n"] / r["s2"]["lines"] if r["s2"]["lines"] else None
    for k in ("n", "one_word", "under_30"):
        row[f"par_final_{k}"] = r["par_final"][k]
    for k in ("pairs", "identical_n", "expected_n", "ratio"):
        row[f"s3_{k}"] = r["s3"][k]
    row["rhyme_lag1"], row["rhyme_lag1_shuffled"] = g(r, "rhyme", "lag1") or (None, None)
    row["rhyme_lag2"], _ = g(r, "rhyme", "lag2") or (None, None)
    row["capital_start"] = r["capital_start"]
    row["auto_form"] = auto_form(r)
    row["overfull_share"] = r["overfull_share"]
    row["sample"] = " / ".join(r["sample"][:2])[:120]
    row["note"] = meta.get("note", "")
    return row


def load_manifest():
    p = HERE / "manifest.csv"
    if not p.exists():
        return {}
    with open(p, newline="", encoding="utf-8") as f:
        return {(r["dataset"], r["manuscript"]): r for r in csv.DictReader(f)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shuffles", type=int, default=1000)
    ap.add_argument("--only", default="")
    ap.add_argument("--out", default="results")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    manifest = load_manifest()
    seen_geom, seen_text, log, rows, full = {}, {}, [], [], {}
    store = {}
    for ds in DATASETS:
        if a.only and a.only.lower() not in ds["name"].lower():
            continue
        if not ds["root"].exists():
            print("missing", ds["root"], file=sys.stderr)
            continue
        for ms, files in sorted(manuscripts(ds).items()):
            meta = manifest.get((ds["name"], ms), {})
            if meta.get("form") == "exclude":
                for p in files:
                    log.append((ds["name"], str(p), "manuscript excluded: " + meta.get("note", ""), 0))
                continue
            pages = read_manuscript(ds, files, seen_geom, seen_text, log)
            if not pages:
                continue
            r = measure(pages, a.shuffles, a.seed, ds["paragraphs"])
            full[f"{ds['name']} | {ms}"] = r | {"meta": meta}
            rows.append(flat(ds["name"], ms, meta, r))
            store[(ds["name"], ms)] = pages
            x = rows[-1]
            f = lambda v, p="%.3f": "  -  " if v is None else p % v  # noqa: E731
            print("%-20s %-34s pg %3d col %3d ln %5d w/l %4.1f | S1 %s/%s z %s  lag2 %s/%s z %s | S2 early %s (n=%s) 1w %s <30 %s | S3 %s (%s/%s) | rhyme %s/%s cap %s" % (
                ds["name"][:20], ms[:34], x["pages"], x["columns"], x["lines"], x["words_per_line"] or 0,
                f(x["s1_lag1_same"]), f(x["s1_lag1_shuffled"]), f(x["s1_lag1_z"], "%+.1f"),
                f(x["s1_lag2_same"]), f(x["s1_lag2_shuffled"]), f(x["s1_lag2_z"], "%+.1f"),
                f(x["s2_early_share"], "%.2f"), x["s2_early"], f(x["s2_one_word"]), f(x["s2_under_30"]),
                f(x["s3_ratio"], "%.2f"), x["s3_identical_n"], f(x["s3_expected_n"], "%.1f"),
                f(x["rhyme_lag1"], "%.2f"), f(x["rhyme_lag2"], "%.2f"), f(x["capital_start"], "%.2f")), flush=True)
    out = HERE / a.out
    if rows:
        with open(f"{out}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        json.dump(full, open(f"{out}.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    with open(HERE / ("pages.csv" if a.out == "results" else a.out + "_pages.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["dataset", "file", "status", "main_text_lines"])
        w.writerows(log)
    c = Counter(re.sub(r" as .*", "", s) for _, _, s, _ in log)
    print("\nfiles:", dict(c), file=sys.stderr)
    return store, rows


if __name__ == "__main__":
    main()
