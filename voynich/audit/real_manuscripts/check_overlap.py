"""Do two 'manuscripts' share pages?  Compares every pair of manuscripts read by run_all.py
by (a) the name of the page image and (b) the main-text lines, reduced to their base letters
(so that different transcription conventions for the same page still match in part)
and (c) runs of three words.

Usage: python3 -I check_overlap.py     Writes overlap.csv
"""
from __future__ import annotations

import csv
import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_all as R                                                         # noqa: E402


def base(text: str) -> str:
    t = unicodedata.normalize("NFKD", text.lower())
    return re.sub(r"[^a-z]", "", t.replace("u", "v").replace("j", "i"))


def main():
    manifest = R.load_manifest()
    lines, images, grams, size = defaultdict(set), defaultdict(set), defaultdict(set), defaultdict(int)
    for ds in R.DATASETS:
        if not ds["root"].exists():
            continue
        for ms, files in sorted(R.manuscripts(ds).items()):
            if manifest.get((ds["name"], ms), {}).get("form") == "exclude":
                continue
            key = f"{ds['name']} | {ms}"
            for p, image, cols in R.read_manuscript(ds, files, {}, {}, []):
                stem = re.sub(r"\.(jpe?g|png|tiff?)$", "", (image or p.name).lower().split("/")[-1])
                stem = re.sub(r"\.(chocomufin|mufichecker)|\.xml$", "", stem)
                if len(stem) >= 8 and not re.fullmatch(r"[\d_\-a-z]{0,7}", stem):
                    images[stem].add(key)
                for c in cols:
                    for r in c:
                        b = base(r[4])
                        if len(b) >= 18:
                            lines[b].add(key)
                        w = [x for x in (base(t) for t in r[4].split()) if x]
                        for i in range(len(w) - 2):
                            g = " ".join(w[i:i + 3])
                            if len(g) >= 14:
                                grams[g].add(key)
                                size[key] += 1
    pair = defaultdict(lambda: [0, 0, 0])
    for store, slot in ((lines, 0), (images, 1), (grams, 2)):
        for owners in store.values():
            owners = sorted(owners)
            if len(owners) > 6:                          # formulas shared by many manuscripts
                continue
            for i in range(len(owners)):
                for j in range(i + 1, len(owners)):
                    pair[(owners[i], owners[j])][slot] += 1
    rows = []
    for (a, b), (n, m, g) in pair.items():
        share = g / max(1, min(size[a], size[b]))        # share of the smaller sample's word triples
        if n >= 3 or m >= 1 or share >= 0.05:
            rows.append((a, b, n, m, g, round(share, 3)))
    rows.sort(key=lambda r: (-r[2], -r[5]))
    with open(HERE / "overlap.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["manuscript_a", "manuscript_b", "identical_lines", "same_image_names",
                    "shared_word_triples", "share_of_smaller_sample"])
        w.writerows(rows)
    for r in rows:
        print("%-50s %-50s lines %4d images %d triples %5d (%.0f%%)" % (r[0][:50], r[1][:50], r[2], r[3], r[4], 100 * r[5]))
    print(len(rows), "pairs with 3+ identical lines, a shared image name, or 5%+ shared word triples")


if __name__ == "__main__":
    main()
