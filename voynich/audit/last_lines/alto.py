"""Own ALTO reader for Gallicorpora 15th c. MSS. One copy per page (skip *.chocomufin.xml)."""
import glob, os, re
import xml.etree.ElementTree as ET
import numpy as np

ROOT = '/home/user/Gallicorpora/HTR-MSS-15e-Siecle/data/'
NS = {'a': 'http://www.loc.gov/standards/alto/ns-v4#'}


def pages():
    for ms in sorted(os.listdir(ROOT)):
        d = ROOT + ms
        if not os.path.isdir(d):
            continue
        for fn in sorted(glob.glob(d + '/*.xml')):
            if fn.endswith('.chocomufin.xml'):
                continue
            yield ms, fn


def read_page(fn):
    """returns list of blocks: dict(label, lines=[dict(x0,x1,y,text,ltype)])"""
    root = ET.parse(fn).getroot()
    tags = {t.get('ID'): t.get('LABEL') for t in root.iter('{%s}OtherTag' % NS['a'])}
    blocks = []
    for tb in root.iter('{%s}TextBlock' % NS['a']):
        label = tags.get(tb.get('TAGREFS'), None)
        lines = []
        for tl in tb.findall('a:TextLine', NS):
            bl = tl.get('BASELINE')
            if not bl:
                continue
            pts = list(map(float, bl.split()))
            xs = pts[0::2]; ys = pts[1::2]
            text = ' '.join(s.get('CONTENT') or '' for s in tl.findall('a:String', NS)).strip()
            lines.append(dict(x0=min(xs), x1=max(xs), y=float(np.mean(ys)), text=text,
                              ltype=tags.get(tl.get('TAGREFS'), None)))
        blocks.append(dict(label=label, lines=lines, id=tb.get('ID')))
    return blocks


if __name__ == '__main__':
    from collections import Counter, defaultdict
    per = defaultdict(lambda: dict(pages=0, lines=0, words=0, blocks=0))
    lab = Counter(); lt = Counter()
    for ms, fn in pages():
        per[ms]['pages'] += 1
        for b in read_page(fn):
            lab[b['label']] += 1
            if b['label'] and b['label'].startswith('MainZone'):
                per[ms]['blocks'] += 1
                for l in b['lines']:
                    lt[l['ltype']] += 1
                    per[ms]['lines'] += 1
                    per[ms]['words'] += len(l['text'].split())
    print(lab); print(lt)
    for ms, v in per.items():
        print(ms, v, 'w/line %.1f' % (v['words'] / max(1, v['lines'])))
