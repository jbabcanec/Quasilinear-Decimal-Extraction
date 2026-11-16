"""Map voynichese word boxes onto ZL3b loci by word-sequence alignment (own code)."""
import json, os, re, difflib
import numpy as np
from collections import defaultdict, Counter
import zl

D = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/voynich-units/morphometry_voynichese/voynichese_boxes/'


def norm(w):
    return re.sub(r'[^a-z]', '', w.lower())


def load_boxes(page):
    fn = D + page + '.js'
    if not os.path.exists(fn):
        return None
    d = json.load(open(fn))
    words = d[0]
    return [(norm(words[b[0]][0]), b[1], b[2], b[3], b[4]) for b in d[1]]


def box_lines(boxes, tol=30):
    """segment reading-order boxes into visual lines by x-reset"""
    lines = []; cur = []
    for b in boxes:
        if cur and b[1] < cur[-1][1] - tol:
            lines.append(cur); cur = []
        cur.append(b)
    if cur:
        lines.append(cur)
    return lines


def align_page(page, rows):
    """rows: ZL loci of that page (all types, file order). returns dict locus_num -> list of boxes,
    plus match-rate"""
    boxes = load_boxes(page)
    if boxes is None:
        return None, 0
    zw = []; zloc = []
    for r in rows:
        for w in zl.words(r['raw']):
            zw.append(norm(w)); zloc.append(r['num'])
    bw = [b[0] for b in boxes]
    sm = difflib.SequenceMatcher(None, zw, bw, autojunk=False)
    assign = [None] * len(bw)
    nm = 0
    for a, b, n in sm.get_matching_blocks():
        for k in range(n):
            assign[b + k] = zloc[a + k]; nm += 1
    # fill gaps between anchor blocks: positional / proportional assignment of the
    # unmatched box run to the unmatched ZL run that sits between the same anchors
    blocks = sm.get_matching_blocks()
    pa, pb = 0, 0
    for a, b, n in blocks:
        za = list(range(pa, a)); bb = list(range(pb, b))
        if bb and za:
            for k, bi in enumerate(bb):
                zi = za[min(len(za) - 1, int(k * len(za) / len(bb)))]
                assign[bi] = zloc[zi]
        elif bb and not za:
            # extra boxes with no ZL counterpart: attach to neighbouring anchor in same visual line
            pass
        pa, pb = a + n, b + n
    vl = []; lid = 0
    for i, b in enumerate(boxes):
        if i and b[1] < boxes[i - 1][1] - 30:
            lid += 1
        vl.append(lid)
    for i in range(len(bw)):
        if assign[i] is None:
            best = None
            for dlt in range(1, 12):
                for j in (i - dlt, i + dlt):
                    if 0 <= j < len(bw) and assign[j] is not None and vl[j] == vl[i]:
                        best = assign[j]; break
                if best:
                    break
            assign[i] = best
    out = defaultdict(list)
    for a, b in zip(assign, boxes):
        if a is not None:
            out[a].append(b)
    return out, nm / max(1, len(bw))


if __name__ == '__main__':
    rows = zl.read()
    bypage = defaultdict(list)
    for r in rows:
        bypage[r['page']].append(r)
    rates = []
    missing = []
    for pg, rr in bypage.items():
        out, rate = align_page(pg, rr)
        if out is None:
            missing.append(pg); continue
        rates.append((rate, pg))
    print('pages in ZL', len(bypage), 'missing boxes', missing)
    rates.sort()
    print('lowest match rates', rates[:15])
    print('median match rate', np.median([r for r, _ in rates]))
