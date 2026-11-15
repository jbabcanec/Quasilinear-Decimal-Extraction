"""Experiment 6 - context ablation with a small transformer language model.

One character-level transformer is trained per corpus to predict each line of text given
up to three lines of context, where the context is one of

    A  the true preceding lines
    S  the same lines with the words inside each line shuffled   (same vocabulary, no order)
    R  random other lines from the same page                     (page vocabulary only)
    N  nothing

A condition token tells the model which it is getting, so one model serves all four and the
comparisons are paired on identical held-out lines.  What the differences mean:

    S - A   what the true ORDER of the preceding text is worth   (does the text run on?)
    R - S   what NEARBY lines add over the rest of the page      (local vocabulary)
    N - R   what the page adds                                   (page vocabulary)

all in bits per character of the predicted line.  Usage:
    python3 src/exp6_llm_context.py <corpus> [seed]
corpus in: voynich french latin naibbe timm
"""
from __future__ import annotations

import json
import math
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "control_real"))
from vdata import EXT, ROOT

K = 3                  # context lines
L = 320                # max sequence length (characters)
CONDS = "ASRN"
FRENCH = ["btv1b100737746", "btv1b55008562q", "btv1b90069505", "btv1b9009670t",
          "btv1b90595325", "btv1b90616967"]


# ------------------------------------------------------------------------------ corpora
def chunk(words, per_line=9, per_page=30, limit=34000):
    words = words[:limit]
    lines = [words[i:i + per_line] for i in range(0, len(words), per_line)]
    return [lines[i:i + per_page] for i in range(0, len(lines), per_page)]


def load(name):
    """-> list of pages; a page is a list of lines; a line is a list of words."""
    if name == "voynich":
        from vdata import paragraph_text, tokens
        pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
        return [[list(h.word) for _, h in g.groupby("line", sort=True)]
                for _, g in pt.groupby("folio", sort=False)]
    if name == "french":
        from control_real import token_table
        root = Path("/home/user/Gallicorpora/HTR-MSS-15e-Siecle/data")
        pages = []
        for ms in FRENCH:
            t = token_table(root / ms)
            t["word"] = t.word.map(lambda s: re.sub(r"[^\w]", "", s.lower()) or "_")
            pages += [[list(h.word) for _, h in g.groupby("line", sort=True)]
                      for _, g in t.groupby("folio", sort=False)]
        return pages
    if name == "latin":
        words = []
        for f in sorted((EXT / "voynich-units" / "herbal").glob("pliny*.tess")):
            for ln in f.read_text(errors="replace").splitlines():
                words += re.findall(r"[a-z]+", re.sub(r"<[^>]*>", " ", ln.lower()))
        return chunk(words)
    if name == "naibbe":
        return chunk((EXT / "controls" / "naibbe" / "greshko_nathist_output_ciphertext.txt").read_text().split())
    if name == "timm":
        lines = [ln.split() for ln in (EXT / "controls" / "selfcitation" / "selfcite_seed1.txt")
                 .read_text().splitlines() if ln.strip()]
        return [lines[i:i + 30] for i in range(0, len(lines), 30)]
    raise SystemExit("unknown corpus")


# ------------------------------------------------------------------------------- examples
class Data:
    def __init__(self, pages, seed):
        self.pages = [p for p in pages if len(p) >= 2]
        counts = Counter(c for p in self.pages for ln in p for w in ln for c in w)
        self.chars = sorted(c for c, n in counts.items() if n >= 20)
        self.ix = {c: i + 8 for i, c in enumerate(self.chars)}       # 0 pad, 1-4 cond, 5 nl, 6 space, 7 unk
        self.V = len(self.chars) + 8
        order = np.random.default_rng(seed).permutation(len(self.pages))
        n = len(order)
        self.split = dict(train=order[: int(.7 * n)], val=order[int(.7 * n): int(.8 * n)],
                          test=order[int(.8 * n):])

    def enc(self, words):
        out = []
        for i, w in enumerate(words):
            if i:
                out.append(6)
            out += [self.ix.get(c, 7) for c in w]
        return out

    def example(self, p, t, cond, rng):
        """tokens, target mask, first-word mask for target line t of page p."""
        page = self.pages[p]
        if cond == "N":
            ctx = []
        elif cond == "R":
            far = [i for i in range(len(page)) if abs(i - t) > K] or [i for i in range(len(page)) if i != t]
            pick = rng.choice(far, size=min(K, len(far), max(t, 1)), replace=False)
            ctx = [page[i] for i in pick]
        else:
            ctx = [list(page[i]) for i in range(max(0, t - K), t)]
            if cond == "S":
                ctx = [list(rng.permutation(ln)) for ln in ctx]
        body = []
        for ln in ctx:
            body += self.enc(ln) + [5]
        tgt = self.enc(page[t]) + [5]
        body = body[-(L - 1 - len(tgt)):] if len(tgt) < L - 1 else []
        tgt = tgt[: L - 1]
        toks = [1 + CONDS.index(cond)] + body + tgt
        mask = [0] * (1 + len(body)) + [1] * len(tgt)
        nfirst = len(page[t][0]) + 1
        first = [0] * (1 + len(body)) + [1] * min(nfirst, len(tgt)) + [0] * max(0, len(tgt) - nfirst)
        return toks, mask, first

    def batch(self, items, rng):
        ex = [self.example(p, t, c, rng) for p, t, c in items]
        n = max(len(e[0]) for e in ex)
        X = np.zeros((len(ex), n), dtype=np.int64); M = np.zeros((len(ex), n), dtype=np.float32)
        Fm = np.zeros((len(ex), n), dtype=np.float32)
        for i, (tk, mk, fk) in enumerate(ex):
            X[i, : len(tk)] = tk; M[i, : len(mk)] = mk; Fm[i, : len(fk)] = fk
        return torch.from_numpy(X), torch.from_numpy(M), torch.from_numpy(Fm)


# --------------------------------------------------------------------------------- model
class Block(nn.Module):
    def __init__(self, d, h, drop):
        super().__init__()
        self.h = h
        self.ln1, self.ln2 = nn.LayerNorm(d), nn.LayerNorm(d)
        self.qkv, self.proj = nn.Linear(d, 3 * d), nn.Linear(d, d)
        self.mlp = nn.Sequential(nn.Linear(d, 4 * d), nn.GELU(), nn.Linear(4 * d, d))
        self.drop = drop

    def forward(self, x):
        B, T, D = x.shape
        q, k, v = (t.contiguous() for t in
                   self.qkv(self.ln1(x)).view(B, T, 3, self.h, D // self.h).permute(2, 0, 3, 1, 4))
        a = F.scaled_dot_product_attention(q, k, v, is_causal=True,
                                           dropout_p=self.drop if self.training else 0.0)
        x = x + F.dropout(self.proj(a.transpose(1, 2).reshape(B, T, D)), self.drop, self.training)
        return x + F.dropout(self.mlp(self.ln2(x)), self.drop, self.training)


class GPT(nn.Module):
    def __init__(self, V, d=96, layers=3, heads=4, drop=0.0):     # dropout off: slow on CPU; early stopping instead
        super().__init__()
        self.tok, self.pos = nn.Embedding(V, d), nn.Embedding(L, d)
        self.blocks = nn.ModuleList(Block(d, heads, drop) for _ in range(layers))
        self.ln, self.out = nn.LayerNorm(d), nn.Linear(d, V)

    def forward(self, x):
        h = self.tok(x) + self.pos(torch.arange(x.shape[1]))
        for b in self.blocks:
            h = b(h)
        return self.out(self.ln(h))


def nll(model, X):
    """per-position negative log-likelihood (nats) of token t+1 given tokens <= t."""
    logits = model(X[:, :-1])
    return F.cross_entropy(logits.reshape(-1, logits.shape[-1]), X[:, 1:].reshape(-1),
                           reduction="none").view(X.shape[0], -1)


# ---------------------------------------------------------------------------------- main
def evaluate(model, D, which, rng, draws=2):
    """per page: nats and character counts for each condition, all chars and first word."""
    model.eval()
    res = {}
    with torch.no_grad():
        for p in D.split[which]:
            n = len(D.pages[p])
            acc = {c: np.zeros(4) for c in CONDS}           # nats, chars, nats_first, chars_first
            for c in CONDS:
                reps = draws if c in "SR" else 1
                for _ in range(reps):
                    X, M, Fm = D.batch([(p, t, c) for t in range(1, n)], rng)
                    loss = nll(model, X)
                    m, f = M[:, 1:], Fm[:, 1:]
                    acc[c] += np.array([float((loss * m).sum()), float(m.sum()),
                                        float((loss * f).sum()), float(f.sum())]) / reps
            res[int(p)] = {c: acc[c].tolist() for c in CONDS}
    model.train()
    return res


def bits(res, cond, first=False):
    a = np.array([r[cond] for r in res.values()])
    return float(a[:, 2 if first else 0].sum() / a[:, 3 if first else 1].sum() / math.log(2))


def main():
    name = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    torch.manual_seed(seed); torch.set_num_threads(2)
    rng = np.random.default_rng(seed)
    D = Data(load(name), seed)
    items = [(p, t) for p in D.split["train"] for t in range(1, len(D.pages[p]))]
    model = GPT(D.V)
    steps, bs = 1500, 32
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=0.01)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=2e-3, total_steps=steps, pct_start=0.05)
    best, best_state, t0 = 1e9, None, time.time()
    for step in range(1, steps + 1):
        pick = rng.integers(0, len(items), bs)
        X, M, _ = D.batch([(*items[i], CONDS[rng.integers(4)]) for i in pick], rng)
        loss = (nll(model, X) * M[:, 1:]).sum() / M[:, 1:].sum()
        opt.zero_grad(); loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
        if step % 250 == 0 or step == steps:
            v = evaluate(model, D, "val", np.random.default_rng(1), draws=1)
            va = np.mean([bits(v, c) for c in CONDS])
            if va < best:
                best, best_state = va, {k: t.clone() for k, t in model.state_dict().items()}
            print(f"{name} step {step} train {loss.item() / math.log(2):.3f} val {va:.3f} "
                  f"A {bits(v, 'A'):.3f} S {bits(v, 'S'):.3f} R {bits(v, 'R'):.3f} N {bits(v, 'N'):.3f} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
    model.load_state_dict(best_state)
    res = evaluate(model, D, "test", np.random.default_rng(2), draws=3)

    # paired differences with a page bootstrap
    pages = list(res)
    arr = {c: np.array([res[p][c] for p in pages]) for c in CONDS}

    def diff(c1, c2, first=False, w=None):
        w = np.ones(len(pages)) if w is None else w
        i, j = (2, 3) if first else (0, 1)
        return float(((w @ arr[c1][:, i]) / (w @ arr[c1][:, j]) - (w @ arr[c2][:, i]) / (w @ arr[c2][:, j]))
                     / math.log(2))
    brng = np.random.default_rng(3)
    boots = [np.bincount(brng.integers(0, len(pages), len(pages)), minlength=len(pages)).astype(float)
             for _ in range(2000)]
    out = dict(corpus=name, seed=seed, pages=len(D.pages), test_pages=len(pages),
               test_chars=float(arr["A"][:, 1].sum()),
               bits={c: bits(res, c) for c in CONDS}, bits_first={c: bits(res, c, True) for c in CONDS})
    for key, (c1, c2) in {"order": ("S", "A"), "nearby": ("R", "S"), "page": ("N", "R")}.items():
        for first in (False, True):
            d = diff(c1, c2, first)
            bs_ = [diff(c1, c2, first, w) for w in boots]
            out[key + ("_first" if first else "")] = [d, float(np.percentile(bs_, 2.5)), float(np.percentile(bs_, 97.5))]
    (ROOT / "results" / f"exp6_{name}_s{seed}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
