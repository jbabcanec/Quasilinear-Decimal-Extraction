"""Exact binomial/Euler transform for the actual weighted BBP V.

Sum_{j>=0} 16^j/(A-a*j)
 = -1/(15*A) Sum_{r>=0} r!*(-16*a)^r / (15^r Prod_{h=1}^r(A-a*h)).

All denominators are odd. Equality is 2-adic; convergence gains the
factorial valuation. This changes constants, not the proven log exponent.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

from gmpy2 import mpz

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from decimal_bbp_extract import bbp_pi_terms, cutoff
from decimal_extract_v5 import transformed_sum, build_v_euler


def band_data(N, G=192):
    X = mpz(5)**(N-1)
    limit = X.bit_length()+G+8
    K = cutoff(N, G, 'pi')
    band = [(sg, -E, mu) for sg,E,mu in bbp_pi_terms(N,K)
            if E < 0 and -E < limit]
    return X, K, limit, band


def v_euler(N, G=192):
    X = mpz(5)**(N-1)
    K = cutoff(N,G,'pi')
    return build_v_euler(N,K,X.bit_length()+G+8)


def v_original(band):
    Jmax = max(j for _,j,_ in band)
    nodes = [(mpz(sg),mpz(mu),Jmax-j) for sg,j,mu in sorted(band,key=lambda x:x[1])]
    while len(nodes)>1:
        nxt=[]
        for i in range(0,len(nodes)-1,2):
            nl,dl,sl = nodes[i]
            nr,dr,sr = nodes[i+1]
            s=min(sl,sr)
            cap=Jmax-s
            nu=((nl*dr)<<(sl-s))+((nr*dl)<<(sr-s))
            de=dl*dr
            if nu.bit_length()>cap:
                nu &= (mpz(1)<<cap)-1
            if de.bit_length()>cap:
                de &= (mpz(1)<<cap)-1
            nxt.append((nu,de,s))
        if len(nodes)&1:
            nxt.append(nodes[-1])
        nodes=nxt
    nu,de,s=nodes[0]
    mod=mpz(1)<<(Jmax-s)
    return ((nu*pow(de,-1,mod))%mod)<<s


def selftest():
    rng=random.Random(20260906)
    for _ in range(200):
        bits=rng.randrange(1,512)
        A=2*rng.randrange(1,1000)+1
        a=1<<rng.randrange(1,6)
        mod=1<<bits
        ref=sum(pow(16,j,mod)*pow(A-a*j,-1,mod)
                for j in range((bits+3)//4))%mod
        got,_=transformed_sum(A,a,bits)
        assert got==ref,(A,a,bits,int(got),ref)
    for N in list(range(1,61))+[100,997,2500,10000]:
        _,_,_,band=band_data(N)
        V,m,_,_=v_euler(N)
        assert V==v_original(band),N
    print(json.dumps(dict(selftest='PASS',random_series_cases=200,actual_V_cases=64)))


def benchmark(sizes):
    for N in sizes:
        _,_,_,band=band_data(N)
        t=time.perf_counter()
        original=v_original(band)
        orig_time=time.perf_counter()-t
        t=time.perf_counter()
        actual,bits,terms,retained=v_euler(N)
        euler_time=time.perf_counter()-t
        assert actual==original
        print(json.dumps(dict(N=N,bits=bits,original_terms=terms,transformed_terms=retained,
                              original_seconds=orig_time,euler_seconds=euler_time,match=True)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selftest',action='store_true')
    ap.add_argument('--bench',nargs='+',type=int)
    args=ap.parse_args()
    if args.selftest: selftest()
    if args.bench: benchmark(args.bench)
