import numpy as np, json
from parse import load, paragraphs, boxes
from geom import page_geometry

def build(P=None, types=('S','B','T')):
    P = P or load()
    rows = {}
    for pg, d in P.items():
        if d['meta'].get('I') not in types: continue
        r = page_geometry(pg, d)
        if r is None: continue
        out, bl = r
        if not out: continue
        b = boxes(pg)
        # interpolate missing line geometry from neighbours (same physical line as previous)
        ys = [o['g']['yc'] if o['g'] else None for o in out]
        for i, y in enumerate(ys):
            if y is None:
                ys[i] = ys[i-1] if i > 0 else ys[i+1]
        ys = np.array(ys, float)
        dy = np.diff(ys)
        dyp = dy[(dy > 5) & (dy < 40)]
        n = len(ys)
        third = max(3, n // 3)
        cw = sum(bb[3] for bb in b) / max(1, sum(len(bb[0]) for bb in b))
        # paragraph structure: map loci to paragraphs
        ps = paragraphs(d)
        plen = [len(p['lines']) for p in ps]
        # paragraph start line indices
        starts = np.cumsum([0] + plen[:-1])
        rows[pg] = dict(meta=d['meta'], n=n, ys=ys.tolist(), plen=plen, starts=starts.tolist(),
                        typs=[o['L']['typ'] for o in out],
                        y_first=float(ys[0]), y_last=float(ys[-1]),
                        box_ymin=min(bb[2] for bb in b), box_ymax=max(bb[2]+bb[4] for bb in b),
                        box_xmax=max(bb[1]+bb[3] for bb in b),
                        pitch=float(np.median(dyp)), pitch_top=float(np.median(dy[:third])),
                        pitch_bot=float(np.median(dy[-third:])), cw=cw,
                        hmed=float(np.median([bb[4] for bb in b])),
                        closed_last=ps[-1]['closed'])
    return rows

if __name__ == '__main__':
    R = build()
    json.dump(R, open('pages_SBT.json','w'))
    print('page I Q B  n npar yfirst ylast boxymax pitch ptop pbot cw hmed  extent')
    for pg, r in R.items():
        m = r['meta']
        print(f"{pg:6s} {m['I']} {m['Q']} {m['B']} {r['n']:3d} {len(r['plen']):3d} {r['y_first']:6.1f} {r['y_last']:6.1f} {r['box_ymax']:4d} "
              f"{r['pitch']:5.1f} {r['pitch_top']:5.1f} {r['pitch_bot']:5.1f} {r['cw']:5.2f} {r['hmed']:4.1f} {r['y_last']-r['y_first']:6.1f} {(r['y_last']-r['y_first'])/(r['n']-1):5.2f}")
