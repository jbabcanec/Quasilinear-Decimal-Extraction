"""Short tables from results/sweeps.json (run_experiments.py), results/voynich_traces.json
and results/real_mss.json: freedom given -> value reached, for R3, R1 and R2.

    python3 headline.py
"""
from __future__ import annotations

import json
from pathlib import Path

R = Path(__file__).parent / "results"


def m(A, fam, lab, key):
    try:
        return A[fam][lab]["mean"][key]
    except KeyError:
        return float("nan")


def main():
    A = json.load(open(R / "sweeps.json"))
    print("R3  line ends, no other goal            (Voynich: median blank 0.13, coupling +0.05, 95% -0.03..+0.14)")
    print("    %-48s %8s %9s %14s %11s" % ("freedom", "blank", "coupling", "cannot fit", "density sd"))
    for lab in A["R3"]:
        print("    %-48s %8.2f %+9.2f %13.0f%% %11.3f" % (lab, m(A, "R3", lab, "r3_median"), m(A, "R3", lab, "r3_coupling"),
                                                         100 * m(A, "R3", lab, "r3_cannot_reach_margin"), m(A, "R3", lab, "t_adj_sd")))

    looks = ["1", "2", "3", "par"]
    for title, pre in (("Latin letters as they are", "R1 re-flow alone, look-ahead "),
                       ("letters merged to the Voynich's eight shares", "R1 re-flow alone, coarse glyphs, look-ahead ")):
        print("\nR1  same opening glyph as the line above / shuffled, re-flow alone, %s   (Voynich 0.60)" % title)
        print("    %-44s" % "freedom \\ lines planned ahead" + "".join("%8s" % l for l in looks) + "   | blank, coupling at look-ahead 3")
        for lab in A[pre + "1"]:
            print("    %-44s" % lab + "".join("%8.2f" % m(A, pre + l, lab, "r1_ratio") for l in looks) +
                  "   | %.2f %+.2f" % (m(A, pre + "3", lab, "r3_median"), m(A, pre + "3", lab, "r3_coupling")))
    pre = "R1 re-flow alone, coarse glyphs, look-ahead "
    print("\n    lag 2 / lag 3 against shuffled (Voynich x1.29 / x1.07), merged letters")
    for lab in ("elastic last 3 words 20%", "elastic whole line 5%", "elastic whole line 8%"):
        print("    %-44s" % lab + "".join("  x%.2f/x%.2f" % (m(A, pre + l, lab, "r1_lag2_ratio"), m(A, pre + l, lab, "r1_lag3_ratio")) for l in looks))
    for fam in ("R1 spelling variants alone, coarse glyphs", "R1 word-splitting alone", "R1 extra: splitting combined, coarse glyphs",
                "R1 combinations, coarse glyphs", "R1 extra: blank cap 0.6, coarse glyphs", "R1 extra: blank cap 1.0, coarse glyphs",
                "R1 re-flow, coarse glyphs, casting-off error sd 3%", "R1 re-flow, coarse glyphs, casting-off error sd 6%"):
        if fam in A:
            print("\n    " + fam)
            for lab in A[fam]:
                print("      %-52s ratio %.2f  lag2 x%.2f lag3 x%.2f | blank %.2f coupling %+.2f | lines opening with a fragment %2.0f%%" % (
                    lab, m(A, fam, lab, "r1_ratio"), m(A, fam, lab, "r1_lag2_ratio"), m(A, fam, lab, "r1_lag3_ratio"),
                    m(A, fam, lab, "r3_median"), m(A, fam, lab, "r3_coupling"), 100 * m(A, fam, lab, "t_fragment_starts")))

    print("\nR2  paragraph's last line: one word / under 30% / 40-80% / over-full, in %   (Voynich 0.4-2.2 / 5.0-7.2 / 66 / 4.6)")
    hz = [("no goal", "R2 baseline (no goal)")] + [(f"last {l}", f"R2 shape goal, re-flow last {l} lines") for l in ("2", "3", "4")] + \
         [("paragraph", "R2 shape goal, re-flow last par lines")]
    print("    %-38s" % "freedom \\ lines re-flowed" + "".join("%22s" % h for h, _ in hz))
    for lab in A["R2 baseline (no goal)"]:
        print("    %-38s" % lab + "".join("%22s" % ("%.1f/%.1f/%.0f/%.1f" % (
            100 * m(A, f, lab, "r2_one_word"), 100 * m(A, f, lab, "r2_under_30"), 100 * m(A, f, lab, "r2_mid_40_80"),
            100 * m(A, f, lab, "r2_over_full"))) for _, f in hz))
    for fam in A:
        if fam.startswith("R2 shape goal, casting-off") or fam.startswith("R2 goal 'short', re-flow last par"):
            print("\n    " + fam)
            for lab in A[fam]:
                print("      %-46s %.1f / %.1f / %.0f / %.1f" % (lab, 100 * m(A, fam, lab, "r2_one_word"), 100 * m(A, fam, lab, "r2_under_30"),
                                                               100 * m(A, fam, lab, "r2_mid_40_80"), 100 * m(A, fam, lab, "r2_over_full")))

    for fam in A:
        if fam.startswith("ALL THREE"):
            print("\n" + fam)
            for lab in A[fam]:
                print("    %-44s R1 %.3f vs %.3f (%.2f), lag2 x%.2f | R2 %.1f/%.1f/%.0f/%.1f | R3 %.2f %+.2f" % (
                    lab, m(A, fam, lab, "r1_same"), m(A, fam, lab, "r1_shuffled"), m(A, fam, lab, "r1_ratio"),
                    m(A, fam, lab, "r1_lag2_ratio"), 100 * m(A, fam, lab, "r2_one_word"), 100 * m(A, fam, lab, "r2_under_30"),
                    100 * m(A, fam, lab, "r2_mid_40_80"), 100 * m(A, fam, lab, "r2_over_full"),
                    m(A, fam, lab, "r3_median"), m(A, fam, lab, "r3_coupling")))

    if (R / "voynich_traces.json").exists():
        T = json.load(open(R / "voynich_traces.json"))
        a = T["A_displaced_repeats"]
        print("\nTRACES IN THE MANUSCRIPT")
        print("  A. the dodged word: lift over chance at 'last word of the line' and 'second word of the line below'")
        print("     %-62s %.2f (z %+.1f)   %.2f (z %+.1f)" % ("Voynich", a["voynich"]["same line, last word"]["lift"],
              a["voynich"]["same line, last word"]["z"], a["voynich"]["line below, 2nd word"]["lift"], a["voynich"]["line below, 2nd word"]["z"]))
        for k, v in a.items():
            if k.startswith("sim"):
                print("     %-62s %.2f            %.2f" % ((k[5:] + " [R1 %.2f]" % v["r1"]["ratio"])[:62],
                      v["profile"]["same line, last word"]["lift"], v["profile"]["line below, 2nd word"]["lift"]))
        if "B_line_density_swing" in T:
            print("  B. coherent line-to-line swing of written/predicted width (elastic strategies above need 0.035-0.047)")
            for k, v in T["B_line_density_swing"].items():
                if isinstance(v, dict):
                    print("     Voynich, %-60s %.3f" % (k, v["coherent_sd"]))
            print("  C. slope of paragraph density on last-line fill")
            for k, v in T["C_paragraph_density_vs_last_line"].items():
                print("     %-70s %+.4f (se %.4f), corr %+.2f" % (k[:70], v["slope"], v["slope_se"], v["corr"]))
    if (R / "real_mss.json").exists():
        M = json.load(open(R / "real_mss.json"))
        n = sum(v["lines"] for v in M.values())
        w = lambda k: sum(v[k] * v["lines"] for v in M.values()) / n
        print("\nREAL LATIN MANUSCRIPTS (%d, %d lines): words abbreviated %.0f%% (range %.0f-%.0f%%), letters saved about %.0f%%, "
              "line-to-line width swing from abbreviation %.1f%% (max %.1f%%), breaks inside a word >= %.0f%%, "
              "line openings: same as line above %.1f%% vs %.1f%% shuffled" % (
                  len(M), n, 100 * w("abbreviated_words"), 100 * min(v["abbreviated_words"] for v in M.values()),
                  100 * max(v["abbreviated_words"] for v in M.values()), 100 * w("letters_saved_share"),
                  100 * w("line_swing_width"), 100 * max(v["line_swing_width"] for v in M.values()),
                  100 * w("split_breaks_lower_bound"), 100 * w("r1_same"), 100 * w("r1_shuffled")))


if __name__ == "__main__":
    main()
