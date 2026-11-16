"""Checks on the pipeline (run after run_all.py).

1. S1 again with another random seed: results_seed1.csv against results.csv.
2. S1, S2 and S3 recomputed for a sample of manuscripts with plain loops that share no code
   with mslib.s1/s2/s3 (S3 by actually shuffling the words, not by the formula).
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
    firsts = [[M.first_letter(r[4]) for r in c] for c in cols]
    firsts = [[x for x in s if x] for s in firsts]

    def same(seqs, lag):
        return sum(a == b for s in seqs for a, b in zip(s[:-lag], s[lag:]))
    out = {}
    for lag in (1, 2):
        n = sum(max(0, len(s) - lag) for s in firsts)
        obs = same(firsts, lag)
        null = []
        for _ in range(shuffles):
            sh = [rnd.sample(s, len(s)) for s in firsts]
            null.append(same(sh, lag))
        mu = sum(null) / len(null)
        sd = (sum((x - mu) ** 2 for x in null) / len(null)) ** 0.5
        out[f"s1_lag{lag}"] = (n, obs / n, mu / n, (obs - mu) / sd)
    early = one = short = 0
    for c in cols:
        for i, r in enumerate(c):
            nb = [c[j][2] for j in range(i - 3, i + 4) if 0 <= j < len(c) and j != i]
            if not nb or max(nb) <= r[1]:
                continue
            fill = (r[2] - r[1]) / (max(nb) - r[1])
            if 0 < fill < 0.9:
                early += 1
                one += len(M.words(r[4])) == 1
                short += fill < 0.3
    out["s2"] = (early, one / early if early else None, short / early if early else None)
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
    a = {(r["dataset"], r["manuscript"]): r for r in csv.DictReader(open(HERE / "results.csv", encoding="utf-8"))}
    b = {(r["dataset"], r["manuscript"]): r for r in csv.DictReader(open(HERE / "results_seed1.csv", encoding="utf-8"))}
    assert a.keys() == b.keys()
    dz, dsh, fixed = [], [], 0
    for k in a:
        for col in ("lines", "s1_lag1_n", "s1_lag1_same", "s2_early", "s2_one_word_n", "s3_identical_n", "s3_expected_n"):
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
            print("   %-22s %-26s S1 %.4f/%.4f z %+.2f [%.4f/%.4f z %+.2f] lag2 z %+.2f [%+.2f] | S2 n=%d 1w %.3f <30 %.3f "
                  "[n=%s %.3f %.3f] | S3 %d/%.1f = %.3f [%s/%.1f = %.3f]" % (
                      ds["name"], ms[:26], s["s1_lag1"][1], s["s1_lag1"][2], s["s1_lag1"][3],
                      float(r["s1_lag1_same"]), float(r["s1_lag1_shuffled"]), float(r["s1_lag1_z"]),
                      s["s1_lag2"][3], float(r["s1_lag2_z"]), s["s2"][0], s["s2"][1], s["s2"][2],
                      r["s2_early"], float(r["s2_one_word"]), float(r["s2_under_30"]),
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
        r = M.s1([[ln[0][0] for ln in lines[1:]] for _, lines in Z], shuffles=2000)
        v = M.s3([ln for p in load("voynich") for ln in p])
        print("3. Voynich through mslib: S1 lag 1 %.4f vs %.4f, z %+.1f (n=%d); lag 2 %.4f vs %.4f, z %+.1f; S3 %.3f "
              "(exp13: 0.0960 vs 0.1601, z -11.4; lag 2 0.2092 vs 0.1617, z +6.2; exp15: 1.075)" % (
                  r["lag1"]["same"], r["lag1"]["shuffled"], r["lag1"]["z"], r["lag1"]["n"],
                  r["lag2"]["same"], r["lag2"]["shuffled"], r["lag2"]["z"], v["ratio"]))
    except Exception as e:                                # data of the main study not available
        print("3. skipped:", type(e).__name__, e)


if __name__ == "__main__":
    main()
