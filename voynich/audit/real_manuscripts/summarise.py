"""Pools results.csv by language and form, and lists the manuscripts closest to the Voynich.

Usage: python3 -I summarise.py [results.csv]      Writes pooled.csv and flagged.csv

Pooling: observed and shuffled counts are summed over manuscripts; the pooled z is
(sum observed - sum expected) / sqrt(sum of the shuffle variances), manuscripts being
shuffled independently.  'core' = book or autograph hands, 13th-16th century, 200+ lines.

Voynich reference (from results/exp13, exp15, exp16 of this study):
  S1 lag 1  9.6% observed, 16.0% shuffled, z -11 (2,584 pairs); lag 2  20.9% vs 16.2%, z +6
  S2 956 of 3,968 lines stop early (24%); of those 0.3% hold one word, 3.0% are under 30%
  S3 1.07
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
VOY = dict(s1_ratio=0.096 / 0.160, s2_one=0.003, s2_early_share=0.24)


def num(v):
    return None if v in ("", None, "None") else float(v)


def centuries(s: str) -> list[int]:
    return [int(x) for x in re.findall(r"\d+", s or "")]


def core(r) -> bool:
    c = centuries(r["century"])
    return bool(c) and max(c) >= 13 and min(c) <= 16 and r["hand"] in ("book", "autograph") and int(r["lines"]) >= 200


def pool(rows):
    out = dict(manuscripts=len(rows), lines=sum(int(r["lines"]) for r in rows))
    for lag in ("lag1", "lag2"):
        n = obs = exp = var = 0.0
        zs = []
        for r in rows:
            k = num(r[f"s1_{lag}_n"])
            if not k:
                continue
            o, e, z = num(r[f"s1_{lag}_same"]), num(r[f"s1_{lag}_shuffled"]), num(r[f"s1_{lag}_z"])
            n += k
            obs += o * k
            exp += e * k
            if z and not math.isnan(z) and abs(o - e) > 0:
                var += ((o - e) * k / z) ** 2
                zs.append(z)
        out[f"s1_{lag}_pairs"] = int(n)
        out[f"s1_{lag}_same"] = obs / n if n else None
        out[f"s1_{lag}_shuffled"] = exp / n if n else None
        out[f"s1_{lag}_z"] = (obs - exp) / math.sqrt(var) if var else None
        if lag == "lag1":
            out["s1_z_min"] = min(zs) if zs else None
            out["s1_z_max"] = max(zs) if zs else None
            out["s1_z_below_minus2"] = sum(z < -2 for z in zs)
            out["s1_z_above_plus2"] = sum(z > 2 for z in zs)
            out["s1_z_sd"] = (sum((z - sum(zs) / len(zs)) ** 2 for z in zs) / (len(zs) - 1)) ** 0.5 if len(zs) > 2 else None
    L = sum(int(num(r["s2_lines"]) or 0) for r in rows)
    E = sum(int(num(r["s2_early"]) or 0) for r in rows)
    one = sum(int(num(r["s2_one_word_n"]) or 0) for r in rows)
    und = sum(int(num(r["s2_under_30_n"]) or 0) for r in rows)
    out.update(s2_lines=L, s2_early=E, s2_early_share=E / L if L else None, s2_one_word=one / E if E else None,
               s2_under_30=und / E if E else None, s2_one_word_per_1000_lines=1000 * one / L if L else None)
    pf = [(int(num(r["par_final_n"])), num(r["par_final_one_word"])) for r in rows if num(r.get("par_final_n"))]
    out["par_final_n"] = sum(n for n, _ in pf) or None
    out["par_final_one_word"] = sum(n * s for n, s in pf) / sum(n for n, _ in pf) if pf else None
    i = sum(num(r["s3_identical_n"]) or 0 for r in rows)
    e = sum(num(r["s3_expected_n"]) or 0 for r in rows)
    out.update(s3_identical_n=int(i), s3_expected_n=e, s3_ratio=i / e if e else None)
    return out


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "results.csv"
    with open(path, newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["form"] != "exclude"]
    groups = defaultdict(list)
    for r in rows:
        groups[("all", r["language"], r["form"])].append(r)
        groups[("all", "ALL", r["form"])].append(r)
        if core(r):
            groups[("core", r["language"], r["form"])].append(r)
            groups[("core", "ALL", r["form"])].append(r)
    table = []
    for (scope, lang, form), g in sorted(groups.items()):
        table.append(dict(scope=scope, language=lang, form=form) | pool(g))
    with open(HERE / "pooled.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    f3 = lambda v, p="%.3f": "  -  " if v is None else p % v  # noqa: E731
    print("%-5s %-17s %-12s %4s %7s | S1 lag1 obs/shuf  z   [min,max] <-2 >+2 | lag2 obs/shuf  z | S2 early%% n   1w%%  <30%% 1w/1000 | S3" % (
        "scope", "language", "form", "mss", "lines"))
    for t in table:
        print("%-5s %-17s %-12s %4d %7d | %s/%s %s [%s,%s] %2d %2d | %s/%s %s | %s %5d %s %s %s | %s (%d/%.0f)" % (
            t["scope"], t["language"][:17], t["form"][:12], t["manuscripts"], t["lines"],
            f3(t["s1_lag1_same"]), f3(t["s1_lag1_shuffled"]), f3(t["s1_lag1_z"], "%+5.1f"),
            f3(t["s1_z_min"], "%+.1f"), f3(t["s1_z_max"], "%+.1f"), t["s1_z_below_minus2"], t["s1_z_above_plus2"],
            f3(t["s1_lag2_same"]), f3(t["s1_lag2_shuffled"]), f3(t["s1_lag2_z"], "%+5.1f"),
            f3(t["s2_early_share"], "%.2f"), t["s2_early"], f3(t["s2_one_word"]), f3(t["s2_under_30"]),
            f3(t["s2_one_word_per_1000_lines"], "%5.1f"), f3(t["s3_ratio"], "%.2f"), t["s3_identical_n"], t["s3_expected_n"]))

    # manuscripts closest to the Voynich on each statistic
    flagged = []
    for r in rows:
        if int(r["lines"]) < 200:
            continue
        o, e, z = num(r["s1_lag1_same"]), num(r["s1_lag1_shuffled"]), num(r["s1_lag1_z"])
        why = []
        if z is not None and z <= -2.5:
            why.append("S1: lag-1 z %.1f (%.3f vs %.3f, ratio %.2f); lag-2 z %+.1f" % (z, o, e, o / e, num(r["s1_lag2_z"]) or 0))
        ne, ow = int(num(r["s2_early"]) or 0), num(r["s2_one_word"])
        if ne >= 50 and ow is not None and ow <= 0.03:
            why.append("S2: %d lines stop early (%.0f%% of lines), %.1f%% of them one word, %.1f%% under 30%%" % (
                ne, 100 * num(r["s2_early_share"]), 100 * ow, 100 * num(r["s2_under_30"])))
        ex, ra = num(r["s3_expected_n"]) or 0, num(r["s3_ratio"])
        if ex >= 5 and ra is not None and ra >= 0.5:
            why.append("S3: ratio %.2f (%d identical neighbours, %.1f expected)" % (ra, num(r["s3_identical_n"]), ex))
        for w in why:
            flagged.append(dict(dataset=r["dataset"], manuscript=r["manuscript"], shelfmark=r["shelfmark"],
                                language=r["language"], century=r["century"], form=r["form"], hand=r["hand"],
                                lines=r["lines"], core=core(r), flag=w))
    with open(HERE / "flagged.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(flagged[0]))
        w.writeheader()
        w.writerows(flagged)
    print()
    for x in sorted(flagged, key=lambda x: (x["flag"][:2], x["form"], x["dataset"])):
        print("%-4s %-12s %-18s %-32s c%-5s %-6s ln%6s | %s" % (
            "core" if x["core"] else "", x["form"], x["dataset"][:18], (x["shelfmark"] or x["manuscript"])[:32],
            x["century"], x["language"][:6], x["lines"], x["flag"]))


if __name__ == "__main__":
    main()
