"""Nulls for 'text simply runs out'.
N1: word-count renewal: paragraph keeps its total word count; line capacities are
    bootstrapped from its own non-last lines (or page pool); last line = remainder.
N2: phase null: actual capacities, but total word count jittered by U{-J..J}.
N3: geometric re-wrap of word boxes at width*f and with random first-line indent (phase).
"""
import numpy as np, pickle, sys
from collections import defaultdict, Counter
import zl

rng = np.random.default_rng(12345)
rows = zl.read()
MODE = sys.argv[1] if len(sys.argv) > 1 else 'A'
ptypes = ('P0', 'P1', 'Pb', 'Pt', 'Pc', 'Pr') if MODE == 'A' else ('P0', 'P1', 'Pb')
paras, an = zl.paragraphs(rows, ptypes)

P = []
pagepool = defaultdict(list)
for p in paras:
    wc = [len(zl.words(l['raw'])) for l in p['lines']]
    P.append(dict(wc=wc, page=p['lines'][0]['page'], I=p['lines'][0]['I'], H=p['lines'][0]['H']))
    pagepool[p['lines'][0]['page']] += wc[:-1]


def sim_once(P, minl, pool='own', cond='actual'):
    """returns (n, singles, le2, ratios)"""
    n = one = two = 0; ratios = []; lasts = []
    for p in P:
        wc = p['wc']
        if cond == 'actual' and len(wc) < minl:
            continue
        N = sum(wc)
        caps_src = np.array(wc[:-1]) if pool == 'own' else np.array(pagepool[p['page']])
        if len(caps_src) == 0:
            continue
        rem = N; nl = 0; last = None
        while rem > 0:
            c = int(rng.choice(caps_src))
            nl += 1
            if rem <= c:
                last = rem; rem = 0
            else:
                rem -= c
        if cond == 'sim' and nl < minl:
            continue
        n += 1; one += last == 1; two += last <= 2
        ratios.append(last / np.median(caps_src)); lasts.append(last)
    return n, one, two, ratios, lasts


def run(minl, pool, cond, R=1000, label=''):
    ones = []; twos = []; ns = []; allr = []; alll = []
    for _ in range(R):
        n, o, t, r, l = sim_once(P, minl, pool, cond)
        ones.append(o); twos.append(t); ns.append(n); allr += r; alll += l
    ones = np.array(ones); twos = np.array(twos); ns = np.array(ns)
    h, _ = np.histogram(allr, bins=[0, .2, .4, .6, .8, 1.0, 1.2, 10])
    print(f'{label} min{minl} pool={pool} cond={cond}: n~{ns.mean():.0f}  singles exp={ones.mean():.1f} '
          f'[{np.percentile(ones,2.5):.0f},{np.percentile(ones,97.5):.0f}] ({100*ones.mean()/ns.mean():.1f}%)  '
          f'<=2 exp={twos.mean():.1f} [{np.percentile(twos,2.5):.0f},{np.percentile(twos,97.5):.0f}] ({100*twos.mean()/ns.mean():.1f}%) '
          f' min singles over sims={ones.min()}  ratio hist share={np.round(h/h.sum(),3)}')
    c = Counter(alll); tot = sum(c.values())
    print('     last-wc share (1..12):', [round(c[k] / tot, 3) for k in range(1, 13)])
    return ones, twos, ns


# observed
for minl in (4, 3):
    obs = [(p['wc'][-1]) for p in P if len(p['wc']) >= minl]
    c = Counter(obs)
    print(f'OBS mode={MODE} min{minl}: n={len(obs)} singles={c[1]} <=2={c[1]+c[2]}  share(1..12)={[round(c[k]/len(obs),3) for k in range(1,13)]}')
    r = [p['wc'][-1] / np.median(p['wc'][:-1]) for p in P if len(p['wc']) >= minl]
    h, _ = np.histogram(r, bins=[0, .2, .4, .6, .8, 1.0, 1.2, 10]); print('     obs ratio hist share', np.round(h / h.sum(), 3))
    for pool in ('own', 'page'):
        for cond in ('actual', 'sim'):
            run(minl, pool, cond, R=500, label='N1')

# N2 phase null: jitter total words
def phase(minl, J, R=1000):
    ones = []; twos = []
    for _ in range(R):
        o = t = 0
        for p in P:
            wc = p['wc']
            if len(wc) < minl:
                continue
            N = sum(wc) + int(rng.integers(-J, J + 1))
            caps = wc[:-1]
            rem = N; i = 0; last = None
            while rem > 0:
                c = caps[i] if i < len(caps) else int(rng.choice(caps))
                i += 1
                if rem <= c:
                    last = rem; rem = 0
                else:
                    rem -= c
            o += last == 1; t += last <= 2
        ones.append(o); twos.append(t)
    ones = np.array(ones); twos = np.array(twos)
    print(f'N2 phase jitter +-{J} min{minl}: singles exp={ones.mean():.1f} [{np.percentile(ones,2.5):.0f},{np.percentile(ones,97.5):.0f}]  '
          f'<=2 exp={twos.mean():.1f} [{np.percentile(twos,2.5):.0f},{np.percentile(twos,97.5):.0f}]')


for J in (2, 4, 8):
    phase(4, J)
