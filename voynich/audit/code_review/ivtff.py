"""Independent minimal IVTFF reader for the audit (does not import the project's parser)."""
import re
from pathlib import Path
ZL = Path("/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/voynich-units/ZL3b.txt")

def clean(body):
    s = body
    s = re.sub(r"<!.*?>", "", s)                 # inline comments
    s = re.sub(r"<@\w=\w+>", "", s)              # in-line variable changes (hand etc.)
    pstart, pend = "<%>" in s, "<$>" in s
    s = s.replace("<%>", "").replace("<$>", "")
    s = s.replace("<->", " ").replace("<~>", " ")
    assert "<" not in s and ">" not in s, body
    # alternative readings: take the first
    while True:
        t = re.sub(r"\[([^\[\]:]*):[^\[\]]*\]", r"\1", s)
        if t == s: break
        s = t
    s = re.sub(r"@\d+;", "?", s)
    s = s.replace("{", "").replace("}", "").replace("'", "")
    assert not re.search(r"[\[\]{}@;:<>]", s), (body, s)
    words = [w for w in re.split(r"[.,\s]+", s) if w]
    hard = [w for w in re.split(r"[.\s]+", s.replace(",", "")) if w]   # uncertain spaces joined
    return words, hard, pstart, pend

def read():
    pages, loci, cur = {}, [], None
    for n, raw in enumerate(ZL.read_text(encoding="latin-1").splitlines(), 1):
        if not raw.strip() or raw.startswith("#"): continue
        m = re.match(r"^<(f\w+)>\s*(?:<!(.*)>)?\s*$", raw)
        if m:
            cur = m.group(1)
            pages[cur] = dict(re.findall(r"\$(\w)=(\w+)", m.group(2) or ""))
            continue
        m = re.match(r"^<(f\w+)\.(\w+),(.)(\w+)(?:;\w+)?>\s*(.*)$", raw)
        assert m, (n, raw)
        folio, num, pos, ltype, body = m.groups()
        assert folio == cur, (n, raw)
        words, hard, ps, pe = clean(body)
        loci.append(dict(folio=folio, num=num, pos=pos, ltype=ltype, words=words, hard=hard,
                         pstart=ps, pend=pe, fileline=n, raw=body))
    return pages, loci

def paragraphs(loci, include_pt=False, break_on_other=True):
    """Paragraphs strictly from <%> ... <$>, P-type loci only.
    break_on_other=False: a label locus between two lines of a paragraph does not split it."""
    out, cur = [], None
    for L in loci:
        isP = L["ltype"].startswith("P") and (include_pt or L["ltype"] != "Pt")
        if not isP or not L["words"]:
            if break_on_other and cur: out.append(cur); cur = None
            if cur is not None and L["folio"] != cur[0]["folio"]: out.append(cur); cur = None
            continue
        if cur is not None and (L["pstart"] or L["folio"] != cur[0]["folio"]):
            out.append(cur); cur = None
        if cur is None: cur = []
        cur.append(L)
        if L["pend"]: out.append(cur); cur = None
    if cur: out.append(cur)
    return out
