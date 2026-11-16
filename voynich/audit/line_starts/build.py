"""Build a per-line table of paragraph text from ZL3b with my own parser."""
import numpy as np
import pandas as pd
from zlparse import parse, paragraphs_by_markers, words_of, first_sign, last_sign


def strip_lead(t):
    # drop any leading separators / interruption markers
    return t.lstrip('.,\x01\x02')


def build(alt=0):
    pages, order, recs = parse(alt=alt)
    paras, orphans = paragraphs_by_markers(recs)
    rows = []
    for pid, p in enumerate(paras):
        n = len(p)
        for i, r in enumerate(p):
            t = strip_lead(r['text'])
            w = words_of(t, comma_break=True)
            w2 = words_of(t, comma_break=False)
            w0 = w[0] if w else ''
            m = r['meta']
            # text up to first drawing interruption / after it
            segs = [s for s in r['text'].split('\x01')]
            rows.append(dict(
                pid=pid, folio=r['folio'], num=r['num'], pre=r['pre'], typ=r['typ'],
                li=i + 1, npar=n, last=(i == n - 1), text=t,
                nw=len(w), nw_nocomma=len(w2), ndraw=r['ndraw'],
                g_eva=first_sign(w0, 'eva'), g_bench=first_sign(w0, 'bench'),
                g_bq=first_sign(w0, 'bench_qo'),
                w1=w0,
                w2g=first_sign(w[1], 'eva') if len(w) > 1 else None,
                w2g_nc=first_sign(w2[1], 'eva') if len(w2) > 1 else None,
                w3g=first_sign(w[2], 'eva') if len(w) > 2 else None,
                lastg=last_sign(w[-1], 'eva') if w else None,
                lastw_first=first_sign(w[-1], 'eva') if w else None,
                g2=w0[1] if len(w0) > 1 else None,
                Q=m.get('Q'), I=m.get('I'), L=m.get('L', '-'), H=m.get('H', '-'),
                C=m.get('C', '-'),
                seg2g=first_sign(words_of(segs[1])[0], 'eva') if len(segs) > 1 and words_of(segs[1]) else None,
                nchar=len(t.replace('.', '').replace(',', '').replace('\x01', '').replace('\x02', '')),
            ))
    df = pd.DataFrame(rows)
    for c in ['g_eva', 'g_bench', 'g_bq', 'w2g', 'w2g_nc', 'w3g', 'lastg', 'lastw_first', 'g2', 'seg2g']:
        df[c] = df[c].where(df[c] != '?', None)
        df[c] = df[c].where(df[c].notna(), None)
    return df, pages, order, recs, paras, orphans


if __name__ == '__main__':
    df, *_ = build()
    print(df.shape)
    print(df.head(12).to_string())
    print(df.g_eva.value_counts(dropna=False).to_string())
    print(df.typ.value_counts())
