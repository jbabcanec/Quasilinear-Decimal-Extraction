import sys, re
from collections import Counter
sys.dont_write_bytecode = True
sys.path.insert(0, "/home/user/Quasilinear-Decimal-Extraction/voynich/src")
import vdata
from vdata import ZL_PATH, _PAGE_RE, _LOCUS_RE, load_zl, clean_zl_text
raw = ZL_PATH.read_text(encoding="utf-8", errors="replace").splitlines()
un = []
nloc = 0
for r in raw:
    if r.startswith("#") or not r.strip(): continue
    m = _PAGE_RE.match(r)
    if m and "." not in r.split(">")[0]: continue
    if _LOCUS_RE.match(r): nloc += 1; continue
    un.append(r)
print("unmatched non-comment lines:", len(un), "matched loci", nloc)
c = Counter(re.match(r"^<([^>]*)>", u).group(1)[:6] if u.startswith("<") else "CONT:" + u[:10] for u in un)
print(c.most_common(40))
for u in un[:25]: print("  ", u[:150])
pages, loci = load_zl()
print("pages", len(pages), "loci", len(loci))
print("ltype", Counter(L["ltype"] for L in loci).most_common())
print("pos", Counter(L["pos"] for L in loci).most_common())
print("pos x Ptype", Counter((L["pos"], L["ltype"]) for L in loci if L["ltype"].startswith("P")).most_common())
# all raw tags inside locus bodies
tags = Counter()
chars = Counter()
for r in raw:
    m = _LOCUS_RE.match(r)
    if not m: continue
    body = m.group(5)
    for t in re.findall(r"<[^>]*>", body):
        tags[t if not t.startswith("<!") else "<!...>"] += 1
    w, _, _ = clean_zl_text(body)
    for x in w: chars.update(x)
print("tags", tags.most_common())
print("chars in cleaned words", sorted(chars.items(), key=lambda x: -x[1]))
