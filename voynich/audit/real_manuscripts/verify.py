"""Checks on the pipeline (run after run_all.py).

1. S1 again with another random seed: results_seed1.csv against results.csv.
2. S1 (true adjacency included), S2, S2' and S3 recomputed for a sample of manuscripts with
   plain loops that share no code with mslib.s1/s2/s2_after_full/s3 (S3 by actually
   shuffling the words, not by the formula).
3. S1 and S3 of the Voynich itself through mslib, against the figures of exp13 and exp15
   (only if the Voynich data of this study are on disk).

Usage: python3 -I verify.py
"""
from __future__ import annotations

import csv
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mslib as M                                                            # noqa: E402
import run_all as R                                                          # noqa: E402

SAMPLE = [("Gallicorpora-MSS-15", "btv1b90069505"), ("CREMMA-Medieval", "bnf_fr_22549-septSages"),
          ("HTRomance-Italian", "bnf-it-70"), ("Vienna-2024-German", "DTL_7978199"),
          ("Koenigsfelden", "U-17_0954 (German, 1557)"), ("HTRomance-Latin", "bnf-lat-16204")]


def slow(pages, shuffles=2000, seed=7):
    rnd = random.Random(seed)
    cols = [c for _, _, pc in pages for c in pc]
    firsts = [[M.first_letter(r[4]) or ("#%d" % id(r)) for r in c] for c in cols]   # no letter: matches nothing
    out = {}
    for lag in (1, 2):
        # places where line i and line i+lag stand directly under one another
        places = [[i for i in range(len(c) - lag) if all(c[i + j][6] for j in range(1, lag + 1))] for c in cols]
        n = sum(len(p) for p in places)
        obs = sum(s[i] == s[i + lag] for s, p in zip(firsts, places) for i in p)
        null = []
        for _ in range(shuffles):
            tot = 0
            for s, p in zip(firsts, places):
                t = rnd.sample(s, len(s))
                tot += sum(t[i] == t[i + lag] for i in p)
            null.append(tot)
        mu = sum(null) / len(null)
        sd = (sum((x - mu) ** 2 for x in null) / len(null)) ** 0.5
        out[f"s1_lag{lag}"] = (n, obs / n, mu / n, (obs - mu) / sd)
    early = one = short = 0
    b_n = b_one = 0
    for c in cols:
        for i, r in enumerate(c):
            nb = [c[j][2] for j in range(i - 3, i + 4) if 0 <= j < len(c) and j != i]
            if nb and max(nb) > r[1]:
                fill = (r[2] - r[1]) / (max(nb) - r[1])
                if 0 < fill < 0.9:
                    early += 1
                    one += len(M.words(r[4])) == 1
                    short += fill < 0.3
            if i >= 3 and c[i][6] and c[i - 1][6] and c[i - 2][6]:        # three lines directly above
                m = max(c[j][2] for j in (i - 3, i - 2, i - 1))
                if all(m > c[j][1] and (c[j][2] - c[j][1]) / (m - c[j][1]) >= 0.9 for j in (i - 3, i - 2, i - 1)) \
                        and m > r[1] and 0 < (r[2] - r[1]) / (m - r[1]) < 0.9:
                    b_n += 1
                    b_one += len(M.words(r[4])) == 1
    out["s2"] = (early, one / early if early else None, short / early if early else None)
    out["s2b"] = (b_n, b_one)
    lines = [M.words(r[4]) for c in cols for r in c]
    lines = [ln for ln in lines if len(ln) > 1]
    obs = sum(a == b for ln in lines for a, b in zip(ln[:-1], ln[1:]))
    tot = 0
    for _ in range(300):
        for ln in lines:
            t = rnd.sample(ln, len(ln))
            tot += sum(a == b for a, b in zip(t[:-1], t[1:]))
    out["s3"] = (obs, tot / 300, obs / (tot / 300) if tot else None)
    return out


def main():
    # 1. another seed
    subprocess.run([sys.executable, "-I", str(HERE / "run_all.py"), "--shuffles", "2000", "--seed", "1",
                    "--out", "results_seed1"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for extra in ("results_seed1.json", "results_seed1_pages.csv"):       # by-products of the second run
        (HERE / extra).unlink(missing_ok=True)
    a = {(r["dataset"], r["manuscript"]): r for r in csv.DictReader(open(HERE / "results.csv", encoding="utf-8"))}
    b = {(r["dataset"], r["manuscript"]): r for r in csv.DictReader(open(HERE / "results_seed1.csv", encoding="utf-8"))}
    assert a.keys() == b.keys()
    dz, dsh, fixed = [], [], 0
    for k in a:
        for col in ("lines", "s1_lag1_n", "s1_lag1_same", "s2_early", "s2_one_word_n", "s2b_n", "s2b_one_word_n",
                    "s3_identical_n", "s3_expected_n"):
            assert a[k][col] == b[k][col], (k, col)
        fixed += 1
        if a[k]["s1_lag1_z"] and int(a[k]["lines"]) >= 200:
            dz.append(abs(float(a[k]["s1_lag1_z"]) - float(b[k]["s1_lag1_z"])))
            dsh.append(abs(float(a[k]["s1_lag1_shuffled"]) - float(b[k]["s1_lag1_shuffled"])))
    print("1. seed 0 against seed 1: %d manuscripts, all counts identical; for those of 200+ lines the S1 lag-1 z "
          "differs by at most %.2f (median %.3f), the shuffled rate by at most %.4f" % (
              fixed, max(dz), sorted(dz)[len(dz) // 2], max(dsh)))

    # 2. plain-loop recomputation
    print("2. plain-loop recomputation (pipeline value in brackets)")
    for ds in R.DATASETS:
        for ms, files in R.manuscripts(ds).items():
            if (ds["name"], ms) not in SAMPLE:
                continue
            pages = R.read_manuscript(ds, files, {}, {}, [])
            s, r = slow(pages), a[(ds["name"], ms)]
            print("   %-22s %-26s S1 n=%d %.4f/%.4f z %+.2f [n=%s %.4f/%.4f z %+.2f] lag2 z %+.2f [%+.2f] | S2 n=%d 1w %.3f "
                  "<30 %.3f [n=%s %.3f %.3f] | S2' n=%d one-word %d [n=%s %s] | S3 %d/%.1f = %.3f [%s/%.1f = %.3f]" % (
                      ds["name"], ms[:26], s["s1_lag1"][0], s["s1_lag1"][1], s["s1_lag1"][2], s["s1_lag1"][3],
                      r["s1_lag1_n"], float(r["s1_lag1_same"]), float(r["s1_lag1_shuffled"]), float(r["s1_lag1_z"]),
                      s["s1_lag2"][3], float(r["s1_lag2_z"]), s["s2"][0], s["s2"][1], s["s2"][2],
                      r["s2_early"], float(r["s2_one_word"]), float(r["s2_under_30"]),
                      s["s2b"][0], s["s2b"][1], r["s2b_n"], r["s2b_one_word_n"],
                      s["s3"][0], s["s3"][1], s["s3"][2] or 0, r["s3_identical_n"], float(r["s3_expected_n"]),
                      float(r["s3_ratio"] or 0)))

    # 3. the Voynich through the same functions
    src = HERE.parents[1] / "src"
    try:
        sys.path.insert(0, str(src))
        sys.path.insert(0, str(HERE.parents[1] / "control_real"))
        import exp13_line_starts as E
        from exp6_llm_context import load
        Z = E.zl_paragraphs()
        r = M.s1([[ln[0][0] for ln in lines[1:]] for _, lines in Z], None, shuffles=2000)
        v = M.s3([ln for p in load("voynich") for ln in p])
        print("3. Voynich through mslib: S1 lag 1 %.4f vs %.4f, z %+.1f (n=%d); lag 2 %.4f vs %.4f, z %+.1f; S3 %.3f "
              "(exp13: 0.0960 vs 0.1601, z -11.4; lag 2 0.2092 vs 0.1617, z +6.2; exp15: 1.075)" % (
                  r["lag1"]["same"], r["lag1"]["shuffled"], r["lag1"]["z"], r["lag1"]["n"],
                  r["lag2"]["same"], r["lag2"]["shuffled"], r["lag2"]["z"], v["ratio"]))
    except Exception as e:                                # data of the main study not available
        print("3. skipped:", type(e).__name__, e)


if __name__ == "__main__":
    main()
