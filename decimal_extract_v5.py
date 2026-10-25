"""Experimental decimal pi extractor with an Euler-transformed 2-adic core.

This uses the same BBP word phase and shared 5**(N-1) as v3. The weighted
core uses about 56% as many series terms. The proved bit complexity is
still O(N log^3 N); this is a practical acceleration, not an exponent claim.

Run: python decimal_extract_v5.py --pos 1000000
Verify: python decimal_extract_v5.py --selftest
Compare: python decimal_extract_v5.py --bench 10000 100000 1000000
"""
from __future__ import annotations

import argparse
import json
import time

from gmpy2 import mpz

from decimal_bbp_extract import bbp_pi_terms, cutoff, extract_v3, pi_reference


def transformed_rational(A, a, bits):
    """Return numerator, odd denominator, term count for the transformed sum.

    For positive power-of-two a>=2 and odd A, the binomial transform is
      -1/(15*A) sum_{r>=0} r!*(-16*a)^r /
                   (15^r product_{h=1}^r (A-a*h)).
    The r-th term has valuation (5+v2(a))*r - popcount(r).
    All denominators are units. See notes/prime_crt_and_euler_2026-09-06.md.
    """
    if bits < 1 or not A & 1 or a < 2 or a & (a-1):
        raise ValueError("require bits>=1, odd A, power-of-two a>=2")
    rate = 5 + (a.bit_length()-1)
    R = max(1, bits//rate)
    while rate*R - R.bit_count() < bits:
        R += 1
    while R > 1 and rate*(R-1) - (R-1).bit_count() >= bits:
        R -= 1
    mask = (mpz(1) << bits)-1

    def cap(x, precision=bits):
        # Keep small negatives signed: reducing every negative at a leaf
        # would allocate a full-precision residue per term.
        return x & ((mpz(1)<<precision)-1) if x.bit_length() > precision else x

    def bs(lo, hi):
        prefix=max(0,lo-1)
        precision=bits-(rate*prefix-prefix.bit_count())
        if hi-lo == 1:
            if lo == 0:
                return mpz(1), mpz(1), mpz(1)
            p = mpz(-16*a*lo)
            q = mpz(15*(A-a*lo))
            return cap(p,precision), cap(q,precision), cap(p,precision)
        mid = (lo+hi)//2
        pl, ql, tl = bs(lo, mid)
        pr, qr, tr = bs(mid, hi)
        return cap(pl*pr,precision), cap(ql*qr,precision), cap(tl*qr + pl*tr,precision)

    _, q, t = bs(0, R)
    return cap(-t), cap(15*A*q), R


def inverse_2exp(d, bits):
    inv=mpz(1)
    precision=1
    while precision<bits:
        precision=min(2*precision,bits)
        mask=(mpz(1)<<precision)-1
        inv=(inv*(2-(d&mask)*inv))&mask
    return inv


def transformed_sum(A, a, bits):
    """Evaluate the transformed series modulo 2^bits (validation API)."""
    num,den,R=transformed_rational(A,a,bits)
    return (num*inverse_2exp(den,bits))&((mpz(1)<<bits)-1),R


def build_v_euler(N, K, limit):
    """Construct exactly the same finite-band V as v3, without a band array."""
    families = []
    original_terms = 0
    for a,b,c,sign in [(8,1,2,1),(2,1,-1,-1),(8,5,0,-1),(4,3,-1,-1)]:
        first = max(0, (N-1+c)//4+1)
        last = min(K, (limit-1+N-1+c)//4)
        if last < first:
            continue
        bits = 4*last-(N-1)-c
        families.append((a, a*last+b, sign, bits))
        original_terms += last-first+1
    if not families:
        return mpz(0), 1, 0, 0
    Jmax = max(bits for _,_,_,bits in families)
    num,den=mpz(0),mpz(1)
    retained = 0
    mask=(mpz(1)<<Jmax)-1
    for a,A,sign,bits in families:
        n,d,R = transformed_rational(A,a,bits)
        num=(num*d+sign*(n*den << (Jmax-bits)))&mask
        den=(den*d)&mask
        retained += R
    V=(num*inverse_2exp(den,Jmax))&mask
    return V, Jmax, original_terms, retained


def extract_v5(N, G=192, window=12):
    """Experimental pi-only replacement for v3's V phase.

    Uses the existing fixed-guard output convention. It does not add the
    adaptive carry certification discussed in the main paper.
    """
    if N < 1 or G < 8 or window < 1:
        raise ValueError("require N>=1, G>=8, window>=1")
    if (10**window).bit_length()+32 > G:
        raise ValueError("output window needs more guard bits; increase G")
    start = time.perf_counter()
    K = cutoff(N,G,'pi')
    X = mpz(5)**(N-1)
    limit = X.bit_length()+G+8
    t_pre = time.perf_counter()-start
    MOD = 1<<G
    acc = 0
    start = time.perf_counter()
    for sign,E,mu in bbp_pi_terms(N,K):
        if E >= 0:
            if mu == 1:
                continue
            r = (pow(2,E,mu)*pow(5,N-1,mu))%mu
            if r:
                acc = (acc + sign*((r<<G)//mu))%MOD
        elif -E < limit and mu > 1:
            R = pow(5,N-1,mu)
            u = pow((mu+1)>>1,-E,mu)
            ro = (R*u)%mu
            if ro:
                acc = (acc + sign*((ro<<G)//mu))%MOD
    t_word = time.perf_counter()-start
    start = time.perf_counter()
    V,Jmax,old_terms,new_terms = build_v_euler(N,K,limit)
    W = (X*V)&((mpz(1)<<Jmax)-1)
    dy = int(W>>(Jmax-G)) if Jmax>=G else int(W<<(G-Jmax))
    acc = (acc+dy)%MOD
    t_v = time.perf_counter()-start
    stats = dict(N=N,pre=t_pre,word=t_word,vphase=t_v,total=t_pre+t_word+t_v,
                 original_terms=old_terms,transformed_terms=new_terms,
                 Jmax=Jmax,gmp=True,engine='v5-euler')
    return format((acc*10**window)>>G, '0%dd'%window),stats


def selftest():
    positions=list(range(1,61))+[100,997,2500,10000]
    reference=pi_reference(max(positions)+30)
    for N in positions:
        got,_=extract_v5(N)
        want=reference[N-1:N+11]
        assert got==want,(N,got,want)
        old,_=extract_v3(N)
        assert got==old,(N,got,old)
    print(json.dumps(dict(selftest='PASS',positions=len(positions),
                         references=['integer Chudnovsky','v3'])))


def benchmark(sizes, repeats=1):
    # Alternate order across repeats to reduce ordering bias.
    for N in sizes:
        runs=[]
        for repeat in range(repeats):
            current={}
            functions=[('v3',extract_v3),('v5',extract_v5)]
            if repeat & 1:
                functions.reverse()
            for name,fn in functions:
                digits,stats=fn(N)
                current[name]=stats
                current[name]['digits']=digits
            assert current['v3']['digits']==current['v5']['digits'],N
            runs.append(current)
        print(json.dumps(dict(N=N,match=True,runs=runs)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pos',type=int)
    ap.add_argument('--window',type=int,default=12)
    ap.add_argument('--guard',type=int,default=192)
    ap.add_argument('--selftest',action='store_true')
    ap.add_argument('--bench',nargs='+',type=int)
    ap.add_argument('--repeats',type=int,default=1)
    args=ap.parse_args()
    if args.selftest: selftest()
    if args.bench: benchmark(args.bench,args.repeats)
    if args.pos is not None:
        digits,stats=extract_v5(args.pos,args.guard,args.window)
        print(json.dumps(dict(digits=digits,stats=stats)))
    if not args.selftest and not args.bench and args.pos is None:
        ap.print_help()
