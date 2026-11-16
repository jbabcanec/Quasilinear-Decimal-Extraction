"""Independent parser for the v101 (Glen Claston) transliteration."""
import re
import pandas as pd

V101 = '/home/user/Quasilinear-Decimal-Extraction/voynich/external/voynich-units/data/v101/voyn_101.txt'
LOC = re.compile(r'^<([^>]*)>(.*)$')
NUMLOC = re.compile(r'^(\d+[rv]\d*(?:-[rv]\d+)?)\.(\d+)$')


def parse_v101(path=V101):
    rows = []
    for ln, line in enumerate(open(path, encoding='latin-1'), 1):
        line = line.rstrip('\r\n')
        if not line:
            continue
        m = LOC.match(line)
        if not m:
            raise RuntimeError(line[:50])
        loc, body = m.groups()
        end = body[-1] if body and body[-1] in '-=' else ''
        text = body[:-1] if end else body
        mm = NUMLOC.match(loc)
        rows.append(dict(loc=loc, page=(mm.group(1) if mm else loc.split('.')[0]),
                         num=(int(mm.group(2)) if mm else None), numeric=bool(mm),
                         text=text, end=end, srcline=ln))
    return pd.DataFrame(rows)


def v101_table():
    d = parse_v101()
    d = d[d.numeric].reset_index(drop=True)
    # paragraphs: consecutive numeric-locus lines on a page up to a line ending '='
    pid = []
    cur = 0
    prev_page = None
    newpar = True
    li = []
    k = 0
    for page, end in zip(d.page, d.end):
        if page != prev_page and not newpar:
            cur += 1
            k = 0
        prev_page = page
        k += 1
        pid.append(cur)
        li.append(k)
        newpar = False
        if end == '=':
            cur += 1
            k = 0
            newpar = True
    d['pid'] = pid
    d['li'] = li
    d['npar'] = d.groupby('pid').li.transform('max')

    def words(t):
        return [w for w in re.split(r'[.,]', t) if w]
    d['words'] = d.text.map(words)
    d['nw'] = d.words.map(len)
    d['g'] = d.words.map(lambda w: w[0][0] if w else None)
    d['w2g'] = d.words.map(lambda w: w[1][0] if len(w) > 1 else None)
    d['lastg'] = d.words.map(lambda w: w[-1][-1] if w else None)
    for c in ['g', 'w2g', 'lastg']:
        d[c] = d[c].map(lambda x: None if x in ('?', '*', None) else x)
    d['folio'] = 'f' + d.page
    return d


if __name__ == '__main__':
    d = v101_table()
    print(len(d), d.pid.nunique())
    print(d.npar.groupby(d.pid).first().value_counts().sort_index().to_string())
    print(d.g.value_counts(dropna=False).head(40).to_string())
    print(d.head(30)[['loc', 'page', 'num', 'pid', 'li', 'npar', 'g', 'text']].to_string())
