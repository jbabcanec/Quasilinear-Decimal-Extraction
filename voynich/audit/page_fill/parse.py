"""Independent IVTFF parser for the page-fill audit (written from scratch)."""
import re, json, os, collections
import numpy as np

BASE = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/voynich-units'
ZL = os.path.join(BASE, 'ZL3b.txt')
BOX = os.path.join(BASE, 'morphometry_voynichese', 'voynichese_boxes')

HDR = re.compile(r'^<(f\w+)>\s+<!(.*)>')
LOC = re.compile(r'^<(f\w+)\.(\w+),(.)(\w\w)(?:;\w)?>\s*(.*)$')


def clean(text):
    """Return list of words of a locus text (rough EVA tokenisation)."""
    t = re.sub(r'<![^>]*>', '', text)          # inline comments
    t = re.sub(r'<[%$]>', '', t)
    t = re.sub(r'<->', '.', t)                 # drawing intrusion = word break
    t = re.sub(r'<[^>]*>', '', t)
    t = re.sub(r'\[([^:\]]*):[^\]]*\]', r'\1', t)   # alternative readings: first
    t = re.sub(r'@\d+;', '?', t)
    t = t.replace('{', '').replace('}', '').replace("'", '')
    t = t.replace(',', '.')                    # uncertain spaces counted as spaces
    return [w for w in t.split('.') if w.strip()]


def load():
    pages = collections.OrderedDict()
    cur = None
    for raw in open(ZL, encoding='latin-1'):
        ln = raw.rstrip('\n')
        if not ln.strip():
            continue
        if ln.startswith('#'):
            if cur is not None:
                pages[cur]['comments'].append(ln)
            continue
        m = HDR.match(ln)
        if m:
            cur = m.group(1)
            pages[cur] = dict(meta=dict(re.findall(r'\$(\w)=(\w+)', m.group(2))),
                              loci=[], comments=[])
            continue
        m = LOC.match(ln)
        if not m:
            raise ValueError(ln)
        pg, num, mark, typ, text = m.groups()
        pages[pg]['loci'].append(dict(num=num, mark=mark, typ=typ, text=text,
                                      pstart='<%>' in text, pend='<$>' in text,
                                      words=clean(text)))
    return pages


def paragraphs(page):
    """Split the P-type loci of a page into paragraphs using <%> / <$>.
    Returns list of dicts {lines:[locus,...], closed:bool, opened:bool}."""
    out = []
    cur = None
    for L in page['loci']:
        if L['typ'][0] != 'P':
            continue
        if L['pstart'] or cur is None:
            if cur is not None:
                out.append(cur)
            cur = dict(lines=[], closed=False, opened=L['pstart'])
        cur['lines'].append(L)
        if L['pend']:
            cur['closed'] = True
            out.append(cur)
            cur = None
    if cur is not None:
        out.append(cur)
    return out


def boxes(pg):
    f = os.path.join(BOX, pg + '.js')
    if not os.path.exists(f):
        return None
    d = json.load(open(f))
    words = [w[0] for w in d[0]]
    b = [(words[i], x, y, w, h) for i, x, y, w, h in d[1]]
    return b


if __name__ == '__main__':
    P = load()
    print(len(P))
