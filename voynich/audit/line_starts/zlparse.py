"""Independent IVTFF (ZL3b) parser written from scratch for the line-start audit.

Produces a list of line records for paragraph-type loci (and optionally all loci).
"""
import re

ZL = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/voynich-units/ZL3b.txt'

HDR = re.compile(r'^<(f\w+)>\s*<!(.*)>\s*$')
LOC = re.compile(r'^<(f\w+)\.(\w+),([@+*=&~/!])(\w+)>\s*(.*)$')

BENCH = ['cth', 'ckh', 'cph', 'cfh', 'ch', 'sh']


def clean(raw, alt=0):
    """Return (cleaned text with only glyphs and '.' ',' separators, flags)."""
    flags = {
        'pstart': '<%>' in raw,
        'pend': '<$>' in raw,
        'ndraw': raw.count('<->'),
        'ntilde': raw.count('<~>'),
    }
    t = raw
    t = re.sub(r'<![^>]*>', '', t)          # inline comments
    t = t.replace('<%>', '').replace('<$>', '')
    t = t.replace('<->', '\x01').replace('<~>', '\x02')   # keep position of interruptions
    t = re.sub(r'<@[^>]*>', '', t)           # inline variable changes e.g. <@H=3>
    t = re.sub(r'<[^>]*>', '', t)            # any other tag
    # alternative readings [a:b:c] -> choose alt-th (default first)
    def pick(m):
        opts = m.group(1).split(':')
        return opts[alt] if alt < len(opts) else opts[0]
    # ligature braces may be nested inside alternatives; strip braces first
    t = t.replace('{', '').replace('}', '')
    t = re.sub(r'\[([^\]]*)\]', pick, t)
    t = re.sub(r'@(\d+);', lambda m: chr(0x2460 + (int(m.group(1)) % 200)), t)  # special glyph -> one char
    t = t.replace("'", '')
    t = t.strip()
    return t, flags


def words_of(t, comma_break=True):
    s = t.replace('\x01', '.').replace('\x02', '.')
    if comma_break:
        s = s.replace(',', '.')
    else:
        s = s.replace(',', '')
    return [w for w in s.split('.') if w != '']


def first_sign(word, mode='eva'):
    if not word:
        return None
    if mode == 'eva':
        return word[0]
    if mode == 'bench':
        for b in BENCH:
            if word.startswith(b):
                return b
        return word[0]
    if mode == 'bench_qo':
        for b in BENCH + ['qo']:
            if word.startswith(b):
                return b
        return word[0]
    raise ValueError(mode)


def last_sign(word, mode='eva'):
    if not word:
        return None
    if mode == 'bench':
        for b in BENCH:
            if word.endswith(b):
                return b
    return word[-1]


def parse(path=ZL, alt=0):
    pages = {}
    recs = []
    order = []
    cur = None
    for ln, line in enumerate(open(path, encoding='latin-1'), 1):
        line = line.rstrip('\n').rstrip('\r')
        if not line or line.startswith('#'):
            continue
        m = HDR.match(line)
        if m:
            cur = m.group(1)
            meta = dict(re.findall(r'\$(\w)=(\S+)', m.group(2)))
            pages[cur] = meta
            order.append(cur)
            continue
        m = LOC.match(line)
        if not m:
            raise RuntimeError('unparsed line %d: %r' % (ln, line[:80]))
        folio, num, pre, typ, raw = m.groups()
        t, flags = clean(raw, alt=alt)
        meta = pages[folio]
        # in-line hand change
        hm = re.search(r'<@H=(\w+)>', raw)
        if hm:
            meta = dict(meta)
            meta['H'] = hm.group(1)
            pages[folio] = meta   # applies to the rest of the page
        recs.append(dict(folio=folio, num=num, pre=pre, typ=typ, raw=raw, text=t,
                         srcline=ln, meta=meta, **flags))
    return pages, order, recs


def paragraphs_by_markers(recs, types=('P',), strict_types=None):
    """Group P-type loci into paragraphs using <%> ... <$> markers.
    Returns list of paragraphs (each a list of recs) and list of orphan recs."""
    paras, orphans = [], []
    cur = None
    cur_folio = None
    for r in recs:
        if not r['typ'].startswith(types):
            continue
        if strict_types is not None and r['typ'] not in strict_types:
            # line of a non-wanted subtype: still passes markers through
            pass
        if cur is not None and r['folio'] != cur_folio:
            # paragraph ran over the page end without <$>
            paras.append(cur)
            cur = None
        if r['pstart']:
            if cur is not None:
                paras.append(cur)
            cur = []
            cur_folio = r['folio']
        if cur is None:
            orphans.append(r)
            if r['pend']:
                pass
            continue
        cur.append(r)
        if r['pend']:
            paras.append(cur)
            cur = None
    if cur is not None:
        paras.append(cur)
    return paras, orphans


if __name__ == '__main__':
    pages, order, recs = parse()
    print('pages', len(pages), 'loci', len(recs))
    P = [r for r in recs if r['typ'].startswith('P')]
    print('P loci', len(P))
    paras, orphans = paragraphs_by_markers(recs)
    print('paragraphs', len(paras), 'orphans', len(orphans))
    import collections
    print(collections.Counter(len(p) for p in paras).most_common(40))
    for o in orphans[:20]:
        print('ORPHAN', o['folio'], o['num'], o['pre'], o['typ'], o['raw'][:60])
    # markers on non-P loci?
    nonP = [r for r in recs if not r['typ'].startswith('P') and (r['pstart'] or r['pend'])]
    print('non-P loci with para markers', len(nonP))
    for o in nonP[:10]:
        print('  ', o['folio'], o['num'], o['pre'], o['typ'], o['raw'][:60])
    # paragraph start vs prefix
    print(collections.Counter((p[0]['pre'], p[0]['typ']) for p in paras))
    print(collections.Counter((r['pre'], r['typ']) for p in paras for r in p[1:]))
    # paragraphs with pstart in the middle or missing pend
    bad = [p for p in paras if not p[-1]['pend']]
    print('paras without end marker', len(bad))
    for p in bad[:10]:
        print('  ', p[0]['folio'], p[0]['num'], '->', p[-1]['num'])
    # residual characters
    cs = collections.Counter(c for r in P for c in r['text'])
    print(sorted(cs.items(), key=lambda x: -x[1]))
    print(collections.Counter(tuple(sorted(pages[f].items()))[0:0] or (pages[f].get('I'), pages[f].get('L'), pages[f].get('H')) for f in order).most_common(50))
