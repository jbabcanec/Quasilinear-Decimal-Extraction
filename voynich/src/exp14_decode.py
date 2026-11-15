"""Experiment 14 - a first decipherment attack: are the words letters in disguise?
Writes results/exp14_decode.json.

Hypothesis class: each Voynich word stands for one letter of a hidden text, with many
different words per letter (a homophonic cipher; the Naibbe proposal is of this family).

Engine: a hidden Markov model whose hidden states are the letters of a known language with
that language's letter-to-letter statistics held fixed; which word stands for which letter
is learned by expectation-maximisation (Knight et al. 2006).  Each line is a separate
sequence.  The score is

    gain = bits per word by which the true word order is explained better than the same
           lines with their words shuffled

A text whose words really are letters of that language gives a large gain, and the learned
mapping can be checked against the truth where the truth is known.
"""
from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from exp6_llm_context import load
from vdata import EXT, ROOT, paragraph_text, tokens

K = 800
ALPHA = "abcdefghilmnopqrstuvx"


def letters_of(path_glob, root):
    txt = ""
    for f in sorted(root.glob(path_glob)):
        txt += re.sub(r"<[^>]*>", " ", f.read_text(errors="replace").lower())
    txt = txt.replace("j", "i").replace("y", "i").replace("k", "c").replace("w", "v").replace("z", "s")
    return [c for c in txt if c in ALPHA]


def language_model(seq):
    ix = {c: i for i, c in enumerate(ALPHA)}
    A = np.ones((len(ALPHA), len(ALPHA))) * 0.5
    for a, b in zip(seq[:-1], seq[1:]):
        A[ix[a], ix[b]] += 1
    A /= A.sum(1, keepdims=True)
    pi = np.bincount([ix[c] for c in seq], minlength=len(ALPHA)) / len(seq)
    return A, pi


def encode(lines):
    cnt = Counter(w for ln in lines for w in ln)
    ix = {w: i for i, (w, _) in enumerate(cnt.most_common(K))}
    V = len(ix) + 1
    T = max(len(ln) for ln in lines)
    X = np.zeros((len(lines), T), dtype=np.int64); M = np.zeros((len(lines), T), dtype=bool)
    for i, ln in enumerate(lines):
        X[i, : len(ln)] = [ix.get(w, V - 1) for w in ln]; M[i, : len(ln)] = True
    return X, M, V, ix


def em(X, M, V, A, pi, rng, restarts=6, iters=40):
    """Best log-likelihood (bits per word) and emission table over random restarts."""
    S = len(pi); B, T = X.shape
    n = M.sum()
    best = (-1e18, None)
    for _ in range(restarts):
        E = rng.dirichlet(np.ones(V) * 0.3, size=S)
        for it in range(iters):
            alpha = np.zeros((B, T, S)); c = np.ones((B, T))
            a = pi[None, :] * E[:, X[:, 0]].T
            c[:, 0] = a.sum(1); alpha[:, 0] = a / c[:, 0:1]
            for t in range(1, T):
                a = (alpha[:, t - 1] @ A) * E[:, X[:, t]].T
                s = a.sum(1); s[~M[:, t]] = 1.0
                c[:, t] = s; alpha[:, t] = a / s[:, None]
            ll = np.log2(c[M]).sum() / n
            beta = np.ones((B, S)); cnt = np.zeros((S, V))
            for t in range(T - 1, -1, -1):
                g = alpha[:, t] * beta
                g = g / np.maximum(g.sum(1, keepdims=True), 1e-300)
                g[~M[:, t]] = 0
                np.add.at(cnt.T, X[:, t], g)
                if t:
                    nb = (beta * E[:, X[:, t]].T) @ A.T / c[:, t:t + 1]
                    beta = np.where(M[:, t:t + 1], nb, beta)
            E = (cnt + 0.01) / (cnt + 0.01).sum(1, keepdims=True)
        if ll > best[0]:
            best = (ll, E.copy())
    return best


def attack(lines, A, pi, rng, truth=None):
    X, M, V, ix = encode(lines)
    ll, E = em(X, M, V, A, pi, rng)
    sh = [list(rng.permutation(ln)) for ln in lines]
    Xs, Ms, Vs, _ = encode(sh)
    ll0, _ = em(Xs, Ms, Vs, A, pi, rng, restarts=3)
    out = dict(words=int(M.sum()), bits_true=float(ll), bits_shuffled=float(ll0), gain=float(ll - ll0))
    post = E * pi[:, None]
    guess = post.argmax(0)
    if truth is not None:
        ok = tot = 0
        cnt = Counter(w for ln in lines for w in ln)
        for w, i in ix.items():
            if w in truth:
                tot += cnt[w]; ok += cnt[w] * (ALPHA[guess[i]] == truth[w])
        out["accuracy"] = ok / max(tot, 1)
    out["sample"] = " ".join("".join(ALPHA[guess[ix[w]]] if w in ix else "?" for w in ln) for ln in lines[:6])
    return out


def homophonic(letters, rng, per_letter=30, n=34000):
    """Latin letters enciphered as words: each letter has its own stock of invented words."""
    stock = {c: [f"{c}{j}" for j in range(per_letter)] for c in ALPHA}
    p = 1 / np.arange(1, per_letter + 1); p /= p.sum()
    words = [stock[c][rng.choice(per_letter, p=p)] for c in letters[:n]]
    truth = {w: w[0] for c in ALPHA for w in stock[c]}
    return [words[i:i + 9] for i in range(0, len(words), 9)], truth


def voynich_clean():
    """Layout removed: the paragraph-opening gallows and the y/d/s line-start prefix are
    stripped where what remains is a known word; the line-opening word's turn-taking is
    left alone (it is a choice among real words)."""
    pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
    seen = Counter(pt.word[pt.k > 0])
    lines = []
    for _, g in pt.groupby(["folio", "line"], sort=False):
        ws = list(g.word)
        w = ws[0]
        m = re.match(r"^(cph|cth|ckh|cfh|[ptkf])(.+)$", w) if g.par_first.iloc[0] else None
        if m and seen.get(m.group(2), 0) >= 3:
            ws[0] = m.group(2)
        elif len(w) > 2 and w[0] in "yds" and seen.get(w[1:], 0) >= 3:
            ws[0] = w[1:]
        lines.append(ws)
    return lines


def main():
    rng = np.random.default_rng(0)
    lat = letters_of("voynich-units/herbal/pliny*.tess", EXT)
    ita = letters_of("voynich-units/italian/*.txt", EXT)
    models = {"Latin": language_model(lat[200000:]), "Italian": language_model(ita[100000:])}
    out = {}

    def run(name, lines, lang="Latin", truth=None):
        A, pi = models[lang]
        r = attack(lines, A, pi, rng, truth)
        out[f"{name} | read as {lang}"] = r
        print("%-46s read as %-7s gain %+.3f bits/word  (true %.3f, shuffled %.3f)%s" % (
            name, lang, r["gain"], r["bits_true"], r["bits_shuffled"],
            "  mapping accuracy %.0f%%" % (100 * r["accuracy"]) if "accuracy" in r else ""), flush=True)

    hl, truth = homophonic(lat[:40000], rng)
    run("Control: Latin letters as invented words", hl, "Latin", truth)
    hi, truth_i = homophonic(ita[:40000], rng)
    run("Control: Italian letters as invented words", hi, "Latin", truth_i)
    run("Control: Latin words (not a letter cipher)", [ln for p in load("latin") for ln in p])
    run("Naibbe cipher of Latin", [ln for p in load("naibbe") for ln in p])
    run("Timm-Schinner generator", [ln for p in load("timm") for ln in p])
    v_raw = [ln for p in load("voynich") for ln in p]
    v_clean = voynich_clean()
    for lang in ("Latin", "Italian"):
        run("Voynich as transliterated", v_raw, lang)
        run("Voynich with layout stripped", v_clean, lang)
    json.dump(out, open(ROOT / "results" / "exp14_decode.json", "w"), indent=1)
    print("sample decode, Latin control:", out["Control: Latin letters as invented words | read as Latin"]["sample"][:150])
    print("sample decode, Voynich clean:", out["Voynich with layout stripped | read as Latin"]["sample"][:150])


if __name__ == "__main__":
    main()
