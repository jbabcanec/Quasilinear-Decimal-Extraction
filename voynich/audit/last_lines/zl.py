"""Independent ZL3b (IVTFF) parser for the last-line audit. Written from scratch."""
import re
from collections import defaultdict

ZL = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/voynich-units/ZL3b.txt'

HDR = re.compile(r'^<(f[^.,>]+)>\s*<!(.*)>')
LOC = re.compile(r'^<(f[^.,>]+)\.([^,>]+),(.)(..)>\s*(.*)$')


def clean(text):
    """Return text with comments/alternatives/ligatures normalised; keeps . , <-> <~>"""
    t = re.sub(r'<![^>]*>', '', text)
    t = re.sub(r'<@[^>]*>', '', t)
    t = t.replace('<%>', '').replace('<$>', '')
    t = re.sub(r'\[([^:\]]*):[^\]]*\]', r'\1', t)   # first alternative
    t = t.replace('{', '').replace('}', '')
    t = re.sub(r'@\d+;?', '*', t)
    t = t.replace("'", '')
    return t.strip()


def words(text, comma_breaks=False, gap='<->'):
    """split cleaned text into words. <-> and <~> always break words."""
    t = clean(text)
    t = t.replace('<->', '.').replace('<~>', '.')
    if comma_breaks:
        t = t.replace(',', '.')
    else:
        t = t.replace(',', '')
    return [w for w in t.split('.') if w.strip()]


def segments(text, comma_breaks=False):
    """list of word-lists, split at drawing interruptions <->"""
    t = clean(text)
    out = []
    for seg in re.split(r'<->|<~>', t):
        if comma_breaks:
            seg = seg.replace(',', '.')
        else:
            seg = seg.replace(',', '')
        out.append([w for w in seg.split('.') if w.strip()])
    return out


def read():
    """yield dict per locus, in file order, with page header vars"""
    page = None
    hv = {}
    rows = []
    for ln, raw in enumerate(open(ZL, encoding='latin-1'), 1):
        raw = raw.rstrip('\n').rstrip('\r')
        if raw.startswith('#') or not raw.strip():
            continue
        m = HDR.match(raw)
        if m:
            page = m.group(1)
            hv = dict(re.findall(r'\$(\w)=(\S)', m.group(2)))
            continue
        m = LOC.match(raw)
        if not m:
            print('UNPARSED', ln, raw[:60])
            continue
        pg, num, pos, typ, text = m.groups()
        rows.append(dict(ln=ln, page=pg, num=num, pos=pos, typ=typ, raw=text,
                         pstart='<%>' in text, pend='<$>' in text,
                         I=hv.get('I'), H=hv.get('H'), L=hv.get('L'), Q=hv.get('Q'),
                         X=hv.get('X')))
    return rows


def paragraphs(rows, ptypes=('P0', 'P1', 'Pb', 'Pt', 'Pc', 'Pr')):
    """Group loci into paragraphs using <%> ... <$>.  Only loci whose type is in ptypes
    are considered.  Returns list of lists of rows.  Also returns anomalies."""
    paras = []
    cur = None
    anomalies = []
    for r in rows:
        if r['typ'] not in ptypes:
            if cur is not None and r['typ'][0] == 'P':
                # a P-type line that we are skipping inside an open paragraph
                cur_skipped.append(r)
                if r['pend']:
                    # paragraph end sits on a skipped line: close paragraph here
                    paras.append(dict(lines=cur, skipped=cur_skipped, closed_by_skipped=True))
                    cur = None
            continue
        if r['pstart']:
            if cur is not None:
                anomalies.append(('start_while_open', r['page'], r['num']))
                paras.append(dict(lines=cur, skipped=cur_skipped, closed_by_skipped=False, unclosed=True))
            cur = []
            cur_skipped = []
        if cur is None:
            anomalies.append(('line_outside_para', r['page'], r['num'], r['typ']))
            continue
        cur.append(r)
        if r['pend']:
            paras.append(dict(lines=cur, skipped=cur_skipped, closed_by_skipped=False))
            cur = None
    if cur is not None:
        anomalies.append(('unclosed_at_eof',))
    return paras, anomalies


if __name__ == '__main__':
    rows = read()
    print(len(rows), 'loci')
    paras, an = paragraphs(rows)
    print(len(paras), 'paragraphs', len(an), 'anomalies')
    for a in an[:40]:
        print(a)
    from collections import Counter
    print(Counter(tuple(sorted(set(l['typ'] for l in p['lines']))) for p in paras).most_common())
    print(Counter(len(p['lines']) for p in paras).most_common())
