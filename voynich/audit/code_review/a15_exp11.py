"""exp11 key contrast: a 250-word sample is ~3 plant pages but ~0.7 of one bathing/text page.
Does the 'same pictures > same scribe' gap survive when every sample mixes the same number of pages?"""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import numpy as np
from collections import Counter
from vdata import paragraph_text, tokens
pt = paragraph_text(tokens()).sort_values(["folio", "line", "k"])
pages = lambda m: {f: list(h.word) for f, h in pt[m].groupby("folio")}
A = pages((pt.hand == "2") & (pt.illust == "H"))
B = pages(pt.hand.isin(["3", "5"]) & (pt.illust == "H") & (pt.lang == "B"))
C = pages((pt.hand == "2") & pt.illust.isin(["B", "T"]))
for n, X in (("A scribe2 plants", A), ("B scribes3,5 plants", B), ("C scribe2 bathing/text", C)):
    print(n, "pages", len(X), "words/page median", np.median([len(v) for v in X.values()]))
SAMPLE = 250
def overlap(a, b):
    ca, cb = Counter(a), Counter(b)
    return sum(min(ca[w], cb[w]) for w in ca) / SAMPLE
def as_coded(X, rng):
    words = []; npg = 0
    for p in rng.permutation(list(X)):
        words += X[p]; npg += 1
        if len(words) >= SAMPLE: return words[:SAMPLE], npg
    return None, npg
def k_pages(X, rng, k=3):
    """SAMPLE words taken evenly from k different pages (random contiguous stretch of each)."""
    keys = [p for p in X if len(X[p]) >= SAMPLE // k + 1]
    pick = rng.choice(len(keys), k, replace=False); words = []
    for j, i in enumerate(pick):
        w = X[keys[i]]; need = SAMPLE // k + (1 if j < SAMPLE % k else 0)
        s = rng.integers(0, len(w) - need + 1); words += w[s:s + need]
    return words, k
rng = np.random.default_rng(1)
def mean_ov(fa, fb, X, Y, n=3000):
    v, pa, pb = [], [], []
    for _ in range(n):
        a, na = fa(X, rng); b, nb = fb(Y, rng)
        if a is None or b is None: continue
        v.append(overlap(a, b)); pa.append(na); pb.append(nb)
    return np.mean(v), np.std(v) / np.sqrt(len(v)), np.mean(pa), np.mean(pb)
print("as coded:  A vs B (same pictures, other scribes) %.3f ±%.3f  pages per sample %.1f / %.1f" % mean_ov(as_coded, as_coded, A, B))
print("as coded:  A vs C (same scribe, other pictures)  %.3f ±%.3f  pages per sample %.1f / %.1f" % mean_ov(as_coded, as_coded, A, C))
print("3 pages each: A vs B %.3f ±%.3f (%.0f/%.0f pages)" % mean_ov(k_pages, k_pages, A, B))
print("3 pages each: A vs C %.3f ±%.3f (%.0f/%.0f pages)" % mean_ov(k_pages, k_pages, A, C))
# page-bootstrap of the matched difference
def resample(X):
    keys = list(X); pick = rng.integers(0, len(keys), len(keys))
    return {f"{keys[i]}#{j}": X[keys[i]] for j, i in enumerate(pick)}
d = []
for _ in range(300):
    a, b, c = resample(A), resample(B), resample(C)
    d.append(mean_ov(k_pages, k_pages, a, b, n=60)[0] - mean_ov(k_pages, k_pages, a, c, n=60)[0])
print("matched difference: mean %.3f  95%% [%.3f, %.3f]  share>0 %.2f" % (np.mean(d), np.percentile(d, 2.5), np.percentile(d, 97.5), np.mean(np.array(d) > 0)))
