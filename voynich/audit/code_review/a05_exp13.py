"""Independent re-derivation of exp13 (same first glyph as the line above) from ZL3b."""
import sys, json
sys.dont_write_bytecode = True
import numpy as np
from collections import Counter
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/audit/code_review")
import ivtff
pages, loci = ivtff.read()
print("pages", len(pages), "loci", len(loci))

def stat(seqs, lag=1, key=lambda w: w[0], shuffles=2000, seed=1):
    seqs = [[key(w) for w in s] for s in seqs if len(s) >= 2]
    n = sum(max(0, len(s) - lag) for s in seqs)
    obs = sum(a == b for s in seqs for a, b in zip(s[:-lag], s[lag:]))
    # exact expectation under within-sequence permutation
    exp = 0.0
    for s in seqs:
        m = len(s)
        if m - lag <= 0: continue
        c = Counter(s)
        exp += (m - lag) * sum(v * (v - 1) for v in c.values()) / (m * (m - 1))
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(shuffles):
        t = 0
        for s in seqs:
            p = rng.permutation(len(s))
            t += sum(s[p[i]] == s[p[i + lag]] for i in range(len(s) - lag))
        null.append(t)
    null = np.array(null)
    return dict(n=n, same=obs / n, exp=exp / n, null=null.mean() / n, z=(obs - null.mean()) / null.std())

for name, kw in (("strict <%>..<$>, labels split paragraphs", dict(break_on_other=True)),
                 ("labels between lines do not split", dict(break_on_other=False))):
    P = ivtff.paragraphs(loci, **kw)
    firsts = [[L["words"][0] for L in p[1:]] for p in P]
    print(name, "paragraphs", len(P))
    for lag in (1, 2, 3):
        print("  lag", lag, {k: round(v, 4) for k, v in stat(firsts, lag).items()})
    # only strictly consecutive '+' loci
    seqs = []
    for p in P:
        run = []
        for a in p[1:]:
            if a["pos"] != "+": 
                if run: seqs.append(run); run = []
                continue
            run.append(a["words"][0])
        if run: seqs.append(run)
    print("  only '+' continuation loci:", {k: round(v, 4) for k, v in stat(seqs, 1).items()})
    # glyph-level key: treat ch, sh, cth, ckh, cph, cfh as single glyphs
    import re
    def g1(w):
        m = re.match(r"c[tkpf]h|ch|sh", w)
        return m.group(0) if m else w[0]
    print("  first GLYPH (ch/sh/cXh as units):", {k: round(v, 4) for k, v in stat(firsts, 1, key=g1).items()})
    # exclude '?' openers
    f2 = [[w for w in s] for s in firsts]
    print("  share of line openers that are '?':", np.mean([w[0] == "?" for s in firsts for w in s]))
