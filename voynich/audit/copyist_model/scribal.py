"""What a medieval Latin scribe could vary in the written form of a word.

Three tables, all deliberately generous (each is an upper bound on the freedom; the
simulator scales them down):

  abbr_saving(word)       letters of width saved if the word is written with every common
                          abbreviation of the late-medieval system (Cappelli): nasal
                          strokes, -us -ur -rum -bus -que signs, per/pro/prae, con-, the
                          er/re hook, superscript vowels, and the usual whole-word forms
                          (et est non quod enim dominus ...).
  initial_variants(word)  other glyphs with which the word could legitimately BEGIN.
        'modest'    u/v, i/j, et/& (tironian sign)
        'generous'  + h dropped or added, ae->e, i/y, ca-/ka-, ph/f, chr->xp, con-/com-/cum
                    as the 9-shaped sign, and whole-word signs that start differently
                    (enim .n., ergo/igitur g, vel l-stroke, est as a division sign)
        'max'       + p-with-stroke signs for per/pro/prae and q-with-stroke signs for
                    quod/quam/quae/qui, counted as glyphs distinct from p and q
  split_points(word)      offsets at which the word may be divided across a line break
        'syllable'  Latin syllable boundaries, at least two letters on each side
        'any'       anywhere with at least two letters on each side
"""
from __future__ import annotations

import re
from functools import lru_cache

V = "aeiouy"

SPECIAL = {
    "et": 1, "est": 2, "non": 2, "in": 1, "cum": 2, "quod": 2, "quam": 3, "quae": 3, "qui": 1, "quo": 1,
    "per": 2, "pro": 2, "prae": 3, "pre": 2, "post": 3, "sunt": 2, "esse": 2, "enim": 3, "autem": 3,
    "ergo": 3, "igitur": 5, "sed": 2, "sicut": 2, "secundum": 4, "super": 2, "inter": 3, "propter": 4,
    "quoniam": 5, "quia": 2, "quidem": 3, "tamen": 3, "tantum": 4, "tunc": 2, "vel": 2, "uel": 2,
    "unde": 2, "nisi": 3, "nihil": 3, "mihi": 3, "tibi": 3, "sibi": 3, "omnis": 2, "omnes": 2,
    "omnia": 2, "omnium": 3, "omnibus": 3, "omni": 2, "omnem": 3, "deus": 2, "dei": 1, "deo": 1,
    "deum": 2, "dominus": 4, "domini": 4, "domino": 4, "dominum": 4, "domine": 4, "sanctus": 4,
    "sancti": 4, "spiritus": 5, "christus": 5, "iesus": 2, "ecclesia": 2, "gratia": 3, "gloria": 3,
    "noster": 3, "nostra": 3, "nostri": 3, "nostrum": 4, "vester": 3, "frater": 3, "fratres": 3,
    "pater": 3, "mater": 3, "homo": 2, "hominem": 2, "homines": 2, "habet": 2, "habent": 2,
    "esset": 2, "erat": 1, "item": 1, "contra": 4, "nunc": 1, "ante": 2, "apud": 2, "atque": 2,
    "neque": 2, "quibus": 4, "quid": 2, "quis": 1, "hoc": 1, "haec": 3, "hic": 1, "huius": 2,
    "eius": 2, "cuius": 2, "modo": 2, "sine": 1, "ipse": 1, "ipsa": 1, "ipsum": 2, "ille": 1, "illa": 1,
    "aliquid": 3, "aliquis": 2, "semper": 3, "similiter": 5, "videlicet": 5, "scilicet": 5,
}

_NASAL = re.compile(r"[aeiouy][mn](?=[^aeiouy]|$)")
_PER = re.compile(r"p[ae]r|por")
_ER = re.compile(r"(?<=[bcdfglmnstv])(er|re)")
_SUP = re.compile(r"(?<=[qtpgc])(ra|ri|ua|ui)")


@lru_cache(maxsize=None)
def abbr_saving(word: str) -> float:
    if word in SPECIAL:
        return float(min(SPECIAL[word], len(word) - 1))
    w, s = word, 0
    for end, k in (("que", 2), ("rum", 2), ("bus", 2), ("tur", 2), ("us", 1), ("ur", 1)):
        if w.endswith(end) and len(w) > len(end) + 1:
            s += k; w = w[:-len(end)]
            break
    if w.startswith(("con", "com")) and len(w) > 4:
        s += 2; w = w[3:]
    for pat, k in ((re.compile(r"prae|pre|pro"), 2), (_PER, 2)):
        s += k * len(pat.findall(w)); w = pat.sub("p", w)
    s += len(_NASAL.findall(w))
    w = _NASAL.sub(lambda m: m.group(0)[0], w)
    s += len(_ER.findall(w)); w = _ER.sub("", w)
    s += len(_SUP.findall(w))
    return float(max(0, min(s, len(word) - 1)))


H_ADD = re.compile(r"^(abund|ost|onus|oner|ii$|is$|eis$|er[aeiou]|aren|umer|um[oi]|olus|oler|eder|arusp|"
                   r"ort|erem|ebdom|ymn|edus|ed[aeiou]|ac$|ora|ordeu|ar[uiu]nd|elleb|ibern|iem|omel|"
                   r"ysop|eli[ao]|erb|ospit|umil|onor|onest|abit|abe)")
_KA = re.compile(r"^ca(r|l|p|t|n)")
P_SIGN = re.compile(r"^(per|par|por|pro|prae|pre)")
Q_SIGN = {"quod", "quam", "quae", "qui", "quia", "quoniam", "quoque", "quibus", "quid", "quem", "quo"}


@lru_cache(maxsize=None)
def initial_variants(word: str, level: str) -> tuple:
    """Alternative first glyphs for the word (the default is word[0]); () if none."""
    if level in (None, "none") or not word:
        return ()
    c, out = word[0], []
    nxt = word[1] if len(word) > 1 else ""
    if c == "u":
        out.append("v")
    elif c == "v":
        out.append("u")
    if c == "i":
        out.append("j")
    elif c == "j":
        out.append("i")
    if word == "et" or word == "etiam":
        out.append("&")
    if level in ("generous", "max"):
        if c == "h" and nxt and nxt in V:
            out.append(nxt)                           # habeo -> abeo
        if c in V and H_ADD.match(word):
            out.append("h")                           # abundo -> habundo
        if word.startswith(("ae", "oe")) and len(word) > 2:
            out.append("e")                           # aetas -> etas
        if c == "i":
            out.append("y")                           # imago -> ymago
        if _KA.match(word):
            out.append("k")                           # caput -> kaput
        if word.startswith("ph"):
            out.append("f")
        if word.startswith("chr"):
            out.append("x")                           # christus -> xps
        if word.startswith(("con", "com")) or word == "cum":
            out.append("9")                           # the con/cum sign
        if word == "enim":
            out.append("n")                           # .n.
        if word in ("ergo", "igitur"):
            out.append("g")
        if word in ("vel", "uel"):
            out.append("l")                           # l with a stroke
        if word == "est":
            out.append("%")                           # the est sign
    if level == "max":
        if P_SIGN.match(word):
            out.append("P")                           # p with a stroke, counted as its own glyph
        if word in Q_SIGN:
            out.append("Q")
    return tuple(dict.fromkeys(g for g in out if g != c))


_ONSET2 = {"bl", "br", "cl", "cr", "dr", "fl", "fr", "gl", "gr", "pl", "pr", "tr", "ch", "ph", "th",
           "qu", "gn", "sc", "sp", "st"}


@lru_cache(maxsize=None)
def split_points(word: str, mode: str, min_frag: int = 2) -> tuple:
    n = len(word)
    if mode in (None, "none", "off") or n < 2 * min_frag:
        return ()
    if mode == "any":
        return tuple(range(min_frag, n - min_frag + 1))
    # vowel nuclei; 'u' after q is part of the consonant
    isv = [ch in V and not (ch == "u" and i > 0 and word[i - 1] == "q") for i, ch in enumerate(word)]
    nuclei, i = [], 0
    while i < n:
        if isv[i]:
            j = i + 1
            if j < n and isv[j] and word[i:j + 1] in ("ae", "au", "oe", "eu"):
                j += 1
            nuclei.append((i, j)); i = j
        else:
            i += 1
    pts = []
    for (a0, a1), (b0, b1) in zip(nuclei[:-1], nuclei[1:]):
        cons = b0 - a1
        if cons <= 1:
            p = a1
        elif word[b0 - 2:b0] in _ONSET2:
            p = b0 - 2
        else:
            p = b0 - 1
        if min_frag <= p <= n - min_frag:
            pts.append(p)
    return tuple(dict.fromkeys(pts))


if __name__ == "__main__":
    import numpy as np
    from texts import TEXTS, load
    for nm in TEXTS:
        W = [w for p in load(nm) for w in p]
        L = np.array([len(w) for w in W]); S = np.array([abbr_saving(w) for w in W])
        line = "%-15s letters/word %.2f  abbreviable letters %.1f%%  words with an abbreviation %.0f%%" % (
            nm, L.mean(), 100 * S.sum() / L.sum(), 100 * (S > 0).mean())
        for lev in ("modest", "generous", "max"):
            line += "  | %s: %.1f%% of words have another first glyph" % (
                lev, 100 * np.mean([bool(initial_variants(w, lev)) for w in W]))
        print(line)
    for w in ("consuetudine", "appellantur", "quibuscum", "provinciae", "fortissimi", "effeminandos",
              "aqua", "institutis", "christus", "habet"):
        print(w, [w[:p] + "-" + w[p:] for p in split_points(w, "syllable")], abbr_saving(w),
              initial_variants(w, "max"))
