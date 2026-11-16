"""Writes manifest.csv: language, century and form (prose / verse / ...) of every manuscript.

Sources, in order of precedence:
  1. OVERRIDE below: entries checked by hand against the transcriptions themselves;
  2. the tables in each repository's README (columns Type, Century, Content) and
     CREMMA-Medieval-LAT/data-registry.csv, read as text;
  3. the dataset defaults.
`form` describes the layout of the lines that are measured, because that is what the three
statistics depend on:
  prose        continuous text, lines broken where the column ends
  verse        one verse per line
  verse-run-on verse copied continuously like prose (line breaks fall anywhere)
  mixed        prose and verse, or text with a commentary, on the same pages
  list         lists, tables of contents, price lists, glossaries, inscriptions
  documentary  charters, registers, accounts, letters (not book hands)
  exclude      not used (reason in `note`)
`hand`: book / documentary / autograph.

Usage: python3 -I build_manifest.py   (reads a previous results.csv, if there is one, only to
list manuscripts whose automatic layout class disagrees with the label)
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_all as R                                                         # noqa: E402

H = Path("/home/user")

DEFAULTS = {                       # dataset: (language, century, form, hand)
    "Gallicorpora-MSS-15": ("French", "15", "", "book"),
    "CREMMA-Medieval": ("French", "", "", "book"),
    "CREMMA-Medieval-LAT": ("Latin", "", "", "book"),
    "HTRomance-French": ("French", "", "", "book"),
    "HTRomance-Latin": ("Latin", "", "", "book"),
    "HTRomance-Italian": ("Italian", "", "", "book"),
    "HTRomance-Spain": ("Spanish", "", "", "book"),
    "HTRomance-Occitan": ("Occitan", "", "", "book"),
    "HTRogene-French": ("French", "", "", "book"),
    "HTRogene-Latin": ("Latin", "", "", "book"),
    "HTRogene-Italian": ("Italian", "", "", "book"),
    "HTRogene-Spanish": ("Spanish", "", "", "book"),
    "HTRogene-Occitan": ("Occitan", "", "", "book"),
    "Gatmuza-Occitan": ("Occitan", "", "", "book"),
    "CIHAM-Fabliaux": ("French", "", "verse", "book"),
    "CIHAM-Liber": ("French", "", "prose", "book"),
    "DecameronFR": ("French", "15", "prose", "book"),
    "GrandesChroniques": ("French", "14", "prose", "book"),
    "Torino-L-II-14": ("French", "14", "verse", "book"),
    "FoNDUE-LA-MSS-MA": ("Latin", "", "", "book"),
    "ParisBible": ("Latin", "13", "prose", "book"),
    "Vienna-2023-Latin": ("Latin", "15", "prose", "book"),
    "Vienna-2024-Latin": ("Latin", "15", "prose", "book"),
    "Vienna-2025-Latin": ("Latin", "", "prose", "book"),
    "Vienna-2024-German": ("German", "15", "", "book"),
    "Vienna-2025-German": ("German", "17", "documentary", "documentary"),
    "Vienna-2023-Czech": ("Czech", "15", "prose", "book"),
    "Vienna-2025-Czech": ("Czech", "15", "prose", "book"),
    "Vienna-2025-Hebrew": ("Hebrew", "", "prose", "book"),
    "Vienna-2024-Syriac": ("Syriac", "16", "prose", "book"),
    "Vienna-2025-Syriac": ("Syriac", "", "prose", "book"),
    "LiDi-Ligorio": ("Italian", "16", "list", "autograph"),
    "Eutyches": ("Latin", "9", "mixed", "book"),
    "CarolineMinuscule": ("Latin", "", "prose", "book"),
    "BDD-Burchard": ("Latin", "11", "exclude", "book"),
    "iForal": ("Latin/Portuguese", "", "documentary", "documentary"),
    "GasparoSardi": ("Latin", "16", "list", "autograph"),
    "BibliaArabica": ("Arabic", "", "prose", "book"),
    "Tarima-Arabic": ("Arabic", "", "prose", ""),
}

README = {                          # dataset: README with a table whose rows name data/<folder>
    "HTRomance-French": "HTRomance-Project/medieval-french/README.md",
    "HTRomance-Latin": "HTRomance-Project/medieval-latin/README.md",
    "HTRomance-Italian": "HTRomance-Project/medieval-italian/README.md",
    "HTRomance-Spain": "HTRomance-Project/middle-ages-in-spain/README.md",
    "HTRogene-French": "HTRogene/french/README.md",
    "HTRogene-Latin": "HTRogene/latin/README.md",
    "HTRogene-Italian": "HTRogene/italian/README.md",
    "HTRogene-Spanish": "HTRogene/spanish/README.md",
    "HTRogene-Occitan": "HTRogene/occitan/README.md",
}

TYPE = {"prose": "prose", "vers": "verse", "verse": "verse", "mixed": "mixed"}


def slug(s: str) -> list[str]:
    s = s.lower().replace("é", "e").replace("è", "e").replace("ü", "u").replace("ë", "e")
    return re.findall(r"[a-z]+|\d+", s)


def readme_tables():
    """-> {(dataset, folder or shelfmark tokens): dict(shelfmark, type, century, content)}"""
    out = {}
    for ds, rel in README.items():
        p = H / rel
        if not p.exists():
            continue
        head = None
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.startswith("|"):
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if "Shelfmark" in cells[0]:
                head = cells
                continue
            if head is None or set(cells[0]) <= set("-: "):
                continue
            row = dict(zip(head, cells))
            shelf = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", row.get("Shelfmark", "")).strip("[] ")
            m = re.search(r"\(data/([^)]+)\)", row.get("Links", ""))
            rec = dict(shelfmark=shelf, type=TYPE.get(row.get("Type", "").lower(), ""),
                       century=row.get("Century", "").strip(),
                       content=(row.get("Content") or row.get("Genre") or "").strip(" _"))
            if m:
                out.setdefault((ds, m.group(1)), rec)
            else:
                out.setdefault((ds, tuple(slug(shelf))), rec)
    return out


def from_readme(tables, ds, ms):
    if (ds, ms) in tables:
        return tables[(ds, ms)]
    want = set(slug(ms))
    best = None
    for (d, key), rec in tables.items():
        if d != ds or not isinstance(key, tuple):
            continue
        have = set(key)
        nums = {t for t in want if t.isdigit()}
        if nums and nums <= have and len(want & have) >= max(2, len(want) - 2):
            if best is None or len(want & have) > best[0]:
                best = (len(want & have), rec)
    return best[1] if best else {}


def cremma_lat():
    out = {}
    p = H / "HTR-United/CREMMA-Medieval-LAT/data-registry.csv"
    if p.exists():
        with open(p, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                folder = r["Folder"].split("/")[-1]
                try:
                    a, b = int(r["notBefore"]), int(r["notAfter"])
                    cent = str((a + b) // 2 // 100 + 1)
                except ValueError:
                    cent = ""
                out[folder] = dict(shelfmark=r["Longer Shelfmark"], century=cent, content=r["Content"][:80])
    return out


try:                                # entries checked by hand, see manifest_overrides.py
    from manifest_overrides import OVERRIDE
except ImportError:
    OVERRIDE = {}


def main():
    tables, lat = readme_tables(), cremma_lat()
    prev = {}
    if (HERE / "results.csv").exists():
        with open(HERE / "results.csv", newline="", encoding="utf-8") as f:
            prev = {(r["dataset"], r["manuscript"]): r for r in csv.DictReader(f)}
    rows = []
    for ds in R.DATASETS:
        if not ds["root"].exists():
            continue
        lang, cent, form, hand = DEFAULTS[ds["name"]]
        for ms in sorted(R.manuscripts(ds)):
            rec = dict(dataset=ds["name"], manuscript=ms, shelfmark="", language=lang, century=cent,
                       form=form, hand=hand, note="", label_source="dataset default")
            t = from_readme(tables, ds["name"], ms)
            if ds["name"] == "CREMMA-Medieval-LAT" and ms in lat:
                t = lat[ms]
            if t:
                rec["shelfmark"] = t.get("shelfmark", "")
                rec["century"] = t.get("century") or rec["century"]
                rec["form"] = t.get("type") or rec["form"]
                rec["note"] = t.get("content", "")
                rec["label_source"] = "repository README"
            o = OVERRIDE.get((ds["name"], ms))
            if o:
                rec.update(o)
                rec["label_source"] = "checked by hand"
            rows.append(rec)
    with open(HERE / "manifest.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    n = 0
    for r in rows:
        p = prev.get((r["dataset"], r["manuscript"]))
        if p and r["form"] != "exclude":
            auto = p.get("auto_form", "")
            if not r["form"] or (auto.startswith("verse") != (r["form"] == "verse") and r["form"] in ("prose", "verse", "")):
                n += 1
                print("CHECK %-20s %-38s form=%-8s auto=%-7s lines %5s | %s" % (
                    r["dataset"], r["manuscript"][:38], r["form"] or "?", auto, p["lines"], p.get("sample", "")[:70]))
    print(len(rows), "manuscripts;", n, "to check", file=sys.stderr)


if __name__ == "__main__":
    main()
