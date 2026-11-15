"""Vertical-similarity kernel in physical (pixel) space, with page-bootstrap intervals.

For a pair table P (see pairs.py) this estimates

    lift(d, dx) = sum_pairs sim(a, b) / sum_pairs p_a

where b is a word d lines above a, horizontally offset by dx page-median word widths, and
p_a is the chance that a is similar to a random other word on its own page.  lift == 1
means "no more similar than any two words on this page".

Two things other than local copying raise the lift and must be removed before talking
about copying:
  * position classes  - words at the same place in their lines resemble each other at ANY
                        line distance (line-initial / line-final / drift along the line);
  * slow drift        - nearby lines share vocabulary whatever the horizontal offset.
The far-line rows (d >= 7) of the same map measure the first; the far-dx columns of the
near rows measure the second.  What is left - a bump at dx ~ 0 that is present at d = 1-2
and absent far away - is the signature of copying from the words physically above.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

D_BANDS = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 9), (10, 14), (15, 20)]
DX_EDGES = np.array([-6, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6])
ADX_EDGES = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 6])


def _cells(P, sim, col, edges, bands=D_BANDS):
    """Per-folio sums of observed and expected similar pairs in each (d band, dx bin)."""
    base = {"ex": "p_ex", "e1": "p_e1", "e2": "p_e2"}[sim]
    fol = P["folio"].astype("category")
    f = fol.cat.codes.to_numpy()
    nf = len(fol.cat.categories)
    d = P["d"].to_numpy()
    band = np.full(len(P), -1)
    for i, (lo, hi) in enumerate(bands):
        band[(d >= lo) & (d <= hi)] = i
    x = np.abs(P["dx"].to_numpy()) if col == "adx" else P[col].to_numpy()
    xb = np.digitize(x, edges) - 1
    ok = (band >= 0) & (xb >= 0) & (xb < len(edges) - 1)
    nb, nx = len(bands), len(edges) - 1
    flat = (f[ok] * nb + band[ok]) * nx + xb[ok]
    obs = np.bincount(flat, weights=P[sim].to_numpy()[ok].astype(float), minlength=nf * nb * nx)
    exp = np.bincount(flat, weights=P[base].to_numpy()[ok], minlength=nf * nb * nx)
    cnt = np.bincount(flat, minlength=nf * nb * nx)
    sh = (nf, nb, nx)
    return obs.reshape(sh), exp.reshape(sh), cnt.reshape(sh)


def lift_map(P, sim="e1", col="dx", edges=DX_EDGES, bands=D_BANDS, nboot=1000, seed=0):
    """Returns dict(lift, lo, hi, n, obs) arrays of shape (len(bands), len(edges)-1).
    Intervals: 95% percentile bootstrap resampling folios."""
    obs, exp, cnt = _cells(P, sim, col, edges, bands)
    L = obs.sum(0) / exp.sum(0)
    rng = np.random.default_rng(seed)
    nf = obs.shape[0]
    boots = np.empty((nboot,) + L.shape)
    for i in range(nboot):
        w = np.bincount(rng.integers(0, nf, nf), minlength=nf).astype(float)
        boots[i] = np.tensordot(w, obs, 1) / np.tensordot(w, exp, 1)
    return dict(lift=L, lo=np.nanpercentile(boots, 2.5, 0), hi=np.nanpercentile(boots, 97.5, 0),
                n=cnt.sum(0), obs=obs.sum(0), boots=boots)


def local_contrast(P, sim="e1", near=(0.0, 0.75), far=(2.0, 6.0), d_near=(1, 2), d_far=(7, 20),
                   nboot=2000, seed=0):
    """The headline number.

        C = [lift(d near, |dx| near) / lift(d near, |dx| far)]
            / [lift(d far, |dx| near) / lift(d far, |dx| far)]

    The numerator asks: on the lines just above, is the word straight overhead more similar
    than words further along the same line?  The denominator asks the same of lines far
    above, where copying from overhead cannot be the cause, so it absorbs position-class
    effects.  C = 1: no physically local copying.  Returns C, its bootstrap CI and the four
    lifts it is built from.
    """
    base = {"ex": "p_ex", "e1": "p_e1", "e2": "p_e2"}[sim]
    adx = P["dx"].abs().to_numpy(); d = P["d"].to_numpy()
    fol = P["folio"].astype("category"); f = fol.cat.codes.to_numpy(); nf = len(fol.cat.categories)
    s = P[sim].to_numpy().astype(float); p = P[base].to_numpy()
    cells = {}
    for dn, (dlo, dhi) in (("n", d_near), ("f", d_far)):
        for xn, (xlo, xhi) in (("n", near), ("f", far)):
            m = (d >= dlo) & (d <= dhi) & (adx >= xlo) & (adx < xhi)
            cells[dn + xn] = (np.bincount(f[m], weights=s[m], minlength=nf),
                              np.bincount(f[m], weights=p[m], minlength=nf))

    def stat(w):
        r = {key: (w @ o) / (w @ e) for key, (o, e) in cells.items()}
        return (r["nn"] / r["nf"]) / (r["fn"] / r["ff"]), r

    C, parts = stat(np.ones(nf))
    rng = np.random.default_rng(seed)
    bs = np.array([stat(np.bincount(rng.integers(0, nf, nf), minlength=nf).astype(float))[0]
                   for _ in range(nboot)])
    return dict(C=float(C), lo=float(np.nanpercentile(bs, 2.5)), hi=float(np.nanpercentile(bs, 97.5)),
                near_line_overhead=float(parts["nn"]), near_line_aside=float(parts["nf"]),
                far_line_overhead=float(parts["fn"]), far_line_aside=float(parts["ff"]),
                events_overhead=int(sum(cells["nn"][0])))


def line_decay(P, sim="e1", nboot=1000, seed=0, dmax=20):
    """Lift as a function of line distance alone (all horizontal offsets pooled)."""
    edges = np.array([0, 1e9])
    bands = [(d, d) for d in range(1, dmax + 1)]
    r = lift_map(P, sim, "adx", edges, bands, nboot, seed)
    return {k: r[k][:, 0] for k in ("lift", "lo", "hi", "n")}
