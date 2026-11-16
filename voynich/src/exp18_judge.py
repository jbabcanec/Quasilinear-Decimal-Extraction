"""Experiment 18 - a judge that tells real Voynich lines from fakes.
Writes results/exp18_judge.json.

For each candidate way of producing Voynich-like text, a classifier is trained to tell its
lines from real ones and scored on pages it has not seen (area under the ROC curve:
0.5 = cannot tell, 1.0 = always tells).  The glyph patterns it relies on most show what the
candidate gets wrong.  Features: runs of 1-4 characters, with ^ and $ marking the line's
start and end and a space between words.  Two views: whole lines, and lines with their
first and last word removed (so line-edge habits cannot be the giveaway).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).parent))
from exp6_llm_context import load
from vdata import ROOT


def judge(real, fake, view):
    """real, fake: lists of (page id, [words]).  Returns AUC and the top giveaway patterns."""
    def text(ws):
        ws = ws[1:-1] if view == "inside" else ws
        return ("^" + " ".join(ws) + "$") if view == "whole" else " ".join(ws)
    rows = [(f"r{p}", text(w), 1) for p, w in real if len(w) >= 4] + \
           [(f"f{p}", text(w), 0) for p, w in fake if len(w) >= 4]
    g = np.array([r[0] for r in rows]); X = [r[1] for r in rows]; y = np.array([r[2] for r in rows])
    aucs = []
    for tr, te in GroupKFold(5).split(X, y, g):
        vec = TfidfVectorizer(analyzer="char", ngram_range=(1, 4), min_df=3, sublinear_tf=True)
        A = vec.fit_transform([X[i] for i in tr]); B = vec.transform([X[i] for i in te])
        clf = LogisticRegression(C=1.0, max_iter=400).fit(A, y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(B)[:, 1]))
    vec = TfidfVectorizer(analyzer="char", ngram_range=(1, 4), min_df=3, sublinear_tf=True)
    A = vec.fit_transform(X); clf = LogisticRegression(C=1.0, max_iter=400).fit(A, y)
    names = np.array(vec.get_feature_names_out()); w = clf.coef_[0]
    return float(np.mean(aucs)), [repr(n) for n in names[np.argsort(w)[-6:][::-1]]], \
        [repr(n) for n in names[np.argsort(w)[:6]]]


def main():
    rng = np.random.default_rng(0)
    V = load("voynich")
    real = [(i, ln) for i, page in enumerate(V) for ln in page]
    half = len(real) // 2
    fakes = {
        "Real Voynich, other half (sanity: should be 0.5)": None,
        "Voynich, words shuffled inside each line": [(i, list(rng.permutation(ln))) for i, ln in real],
        "Voynich, words shuffled across the page": None,
        "Timm-Schinner generator": [(i, ln) for i, page in enumerate(load("timm")) for ln in page],
        "Naibbe cipher of Latin": [(i, ln) for i, page in enumerate(load("naibbe")) for ln in page],
    }
    pg = []
    for i, page in enumerate(V):
        words = [w for ln in page for w in ln]; words = list(rng.permutation(words)); k = 0
        for ln in page:
            pg.append((i, words[k:k + len(ln)])); k += len(ln)
    fakes["Voynich, words shuffled across the page"] = pg
    out = {}
    for name, fake in fakes.items():
        if fake is None:      # split real pages into two arbitrary groups
            pages = rng.permutation(len(V)); a = set(pages[: len(pages) // 2])
            r = [(p, w) for p, w in real if p in a]; fk = [(p, w) for p, w in real if p not in a]
        else:
            r, fk = real, fake
        res = {}
        for view in ("whole", "inside"):
            auc, real_marks, fake_marks = judge(r, fk, view)
            res[view] = dict(auc=auc, marks_of_real=real_marks, marks_of_fake=fake_marks)
        out[name] = res
        print("%-50s whole line %.3f | first and last word removed %.3f" % (name, res["whole"]["auc"], res["inside"]["auc"]))
        print("      real gives itself away by:", " ".join(res["whole"]["marks_of_real"]))
        print("      fake gives itself away by:", " ".join(res["whole"]["marks_of_fake"]), flush=True)
    json.dump(out, open(ROOT / "results" / "exp18_judge.json", "w"), indent=1)


if __name__ == "__main__":
    main()
