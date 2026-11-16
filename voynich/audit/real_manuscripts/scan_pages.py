"""Page-by-page search for line-start avoidance (S1) in every dataset.

A manuscript-level figure can hide a few pages that behave differently (a table of
chapters at the front of a prose book).  This applies S1 to each page on its own.

Usage: python3 -I scan_pages.py [--shuffles 1000]      Writes pages_s1.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mslib as M                                                            # noqa: E402
import run_all as R                                                          # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shuffles", type=int, default=1000)
    ap.add_argument("--min-pairs", type=int, default=40)
    a = ap.parse_args()
    manifest = R.load_manifest()
    seen_g, seen_t, rows = {}, {}, []
    for ds in R.DATASETS:
        if not ds["root"].exists():
            continue
        for ms, files in sorted(R.manuscripts(ds).items()):
            meta = manifest.get((ds["name"], ms), {})
            if meta.get("form") == "exclude":
                continue
            for p, _, cols in R.read_manuscript(ds, files, seen_g, seen_t, []):
                seqs = [[M.first_letter(r[4]) for r in c] for c in cols]
                r = M.s1(seqs, [[x[6] for x in c] for c in cols], shuffles=a.shuffles)
                if "lag1" not in r or r["lag1"]["n"] < a.min_pairs:
                    continue
                first = Counter(x for s in seqs for x in s if x)
                top, k = first.most_common(1)[0]
                lines = [x for c in cols for x in c]
                v, w = r["lag1"], r.get("lag2", {})
                rows.append(dict(dataset=ds["name"], manuscript=ms, form=meta.get("form", ""),
                                 language=meta.get("language", ""), page=p.name, lines=len(lines),
                                 pairs=v["n"], same=round(v["same"], 4), shuffled=round(v["shuffled"], 4),
                                 z=round(v["z"], 2), lag2_same=round(w.get("same", float("nan")), 4),
                                 lag2_shuffled=round(w.get("shuffled", float("nan")), 4),
                                 lag2_z=round(w.get("z", float("nan")), 2),
                                 top_letter=top, top_share=round(k / sum(first.values()), 3),
                                 sample=" / ".join(x[4] for x in lines[2:5])[:140]))
    rows.sort(key=lambda r: r["z"])
    with open(HERE / "pages_s1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    zs = [r["z"] for r in rows]
    print(len(rows), "pages with %d+ line pairs; z below -3: %d, below -2: %d, above +2: %d, above +3: %d" % (
        a.min_pairs, sum(z < -3 for z in zs), sum(z < -2 for z in zs), sum(z > 2 for z in zs), sum(z > 3 for z in zs)))
    for r in rows[:30]:
        print("%-18s %-26s %-11s %-22s ln %3d  %.3f vs %.3f z %+.1f | lag2 %.3f vs %.3f z %+.1f | top '%s' %.0f%% | %s" % (
            r["dataset"][:18], r["manuscript"][:26], r["form"], r["page"][-22:], r["lines"], r["same"], r["shuffled"], r["z"],
            r["lag2_same"], r["lag2_shuffled"], r["lag2_z"], r["top_letter"], 100 * r["top_share"], r["sample"][:60]))


if __name__ == "__main__":
    main()
