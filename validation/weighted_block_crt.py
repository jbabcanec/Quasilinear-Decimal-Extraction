"""Experimental exact weighted BBP sums via block CRT and real Euler tails.

Bit bound: O(N log^(5/2) N), with proof in
notes/weighted_block_crt_2026-09-06.md. This file is for mathematical
validation; it does not claim competitive prototype timings.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import random
import sys
import time
from pathlib import Path

from gmpy2 import mpz,invert

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from decimal_bbp_extract import bbp_pi_terms,cutoff,pi_reference
from decimal_extract_v5 import inverse_2exp,build_v_euler
from validation.harmonic_prime_crt import primes_to,rational_tree


def factor_ap(A,a,length,small_primes):
    """Segmented factorization; no sieve of length Theta(N) words."""
    remainders=[A+a*j for j in range(length)]
    factors=[[] for _ in range(length)]
    maximums={}
    for p in small_primes:
        if p*p>A+a*(length-1): break
        first=(-A*pow(a,-1,p))%p
        for j in range(first,length,p):
            pv=1
            while remainders[j]%p==0:
                remainders[j]//=p
                pv*=p
            if pv>1:
                factors[j].append((p,pv))
                maximums[p]=max(maximums.get(p,1),pv)
    for j,left in enumerate(remainders):
        if left>1:
            factors[j].append((left,left))
            maximums[left]=max(maximums.get(left,1),left)
    return factors,maximums


def real_euler_approx(A,a,length,guard):
    """Return floor(2^guard * approximation to the local weighted sum).

    The exact alternating-series remainder is certified below 1/64.
    Extending the finite length to infinity adds <1/[15(A+a*length)].
    """
    # An integer lower bound for log2(15*A/(a*length)). The subsequent
    # exact tail check is authoritative, including outside the asymptotic
    # short-block regime. No floating-point estimate controls correctness.
    log_ratio=max(1,(15*A//(a*length)).bit_length()-1)
    count=max(1,(4*length+16+log_ratio-1)//log_ratio)

    def bs(lo,hi):
        if hi-lo==1:
            if lo==0: return mpz(1),mpz(1),mpz(1)
            p=mpz(-a*lo)
            q=mpz(15*(A+a*lo))
            return p,q,p
        mid=(lo+hi)//2
        pl,ql,tl=bs(lo,mid)
        pr,qr,tr=bs(mid,hi)
        return pl*pr,ql*qr,tl*qr+pl*tr

    while True:
        p,q,t=bs(0,count)
        # First omitted term, including the 16^length/(15A) prefactor.
        error_num=abs(p)*a*count
        error_den=q*225*A*(A+a*count)
        if (error_num<<(4*length+6))<error_den:
            break
        count*=2
    approximation=(t<<(4*length+guard))//(15*A*q)
    return approximation,count,dict(real_product_bits=q.bit_length(),
                                     real_count=count,tail_bound_below='1/64')


def exact_weighted_block(A,a,length,small_primes=None):
    """Exact sum_{j=0}^{length-1} 16^(length-1-j)/(A+a*j)."""
    if A<1 or not A&1 or a<2 or a&(a-1) or length<1:
        raise ValueError('positive odd A, power-of-two a>=2, positive length required')
    maximum=A+a*(length-1)
    if small_primes is None:
        small_primes=list(primes_to(math.isqrt(maximum)))
    factors,powers=factor_ap(A,a,length,small_primes)
    coefficients=dict.fromkeys(powers,0)
    incidences=0
    for j,row in enumerate(factors):
        mu=A+a*j
        for p,pv in row:
            q=powers[p]
            term=(q//pv)*pow(16,length-1-j,q)*int(invert(mu//pv,q))
            coefficients[p]=(coefficients[p]+term)%q
            incidences+=1
    guard=(maximum+length+1).bit_length()+8
    unit=1<<guard
    rational_parts=[]
    t_lower=0
    lcm_bits=0
    for p,q in powers.items():
        c=coefficients[p]
        lcm_bits+=q.bit_length()
        if c:
            rational_parts.append((mpz(c),mpz(q)))
            t_lower+=(c*unit)//q
    approx,R,stats=real_euler_approx(A,a,length,guard)
    correction=(approx-t_lower+unit//2)//unit
    num,den=rational_tree(rational_parts)
    num+=correction*den
    stats.update(A=A,a=a,length=length,local_denominator_bits=den.bit_length(),
                 lcm_leaf_bits=lcm_bits,correction_bits=abs(correction).bit_length(),
                 incidences=incidences)
    return num,den,stats


def merge_nodes(nodes,bits):
    while len(nodes)>1:
        nxt=[]
        for i in range(0,len(nodes)-1,2):
            nl,dl,sl=nodes[i]
            nr,dr,sr=nodes[i+1]
            s=min(sl,sr)
            precision=bits-s
            num=((nl*dr)<<(sl-s))+((nr*dl)<<(sr-s))
            den=dl*dr
            if num.bit_length()>precision: num&=(mpz(1)<<precision)-1
            if den.bit_length()>precision: den&=(mpz(1)<<precision)-1
            nxt.append((num,den,s))
        if len(nodes)&1: nxt.append(nodes[-1])
        nodes=nxt
    num,den,shift=nodes[0]
    precision=bits-shift
    return ((num*inverse_2exp(den,precision))&((mpz(1)<<precision)-1))<<shift


def block_core(N,G=192,block_length=None,*,x_bits=None):
    if N<1 or G<8 or (block_length is not None and block_length<1):
        raise ValueError('require N>=1, G>=8, and positive block length')
    if x_bits is None:
        x_bits=(mpz(5)**(N-1)).bit_length()
    K=cutoff(N,G,'pi')
    limit=x_bits+G+8
    families=[]
    Jmax=1
    max_mu=1
    count=0
    for a,b,c,sign in [(8,1,2,1),(2,1,-1,-1),(8,5,0,-1),(4,3,-1,-1)]:
        first=max(0,(N-1+c)//4+1)
        last=min(K,(limit-1+N-1+c)//4)
        if last<first: continue
        Jmax=max(Jmax,4*last-(N-1)-c)
        max_mu=max(max_mu,a*last+b)
        families.append((a,b,c,sign,first,last))
        count+=last-first+1
    L=max(1,(N-1).bit_length())
    u=math.isqrt(L)
    if u*u<L: u+=1
    if block_length is None:
        block_length=max(1,(N+(1<<u)-1)//(1<<u))
    small_primes=list(primes_to(math.isqrt(max_mu)))
    pending=[]
    pending_mass=0
    V=mpz(0)
    mask=(mpz(1)<<Jmax)-1
    stats=dict(N=N,Jmax=Jmax,terms=count,block_length=block_length,blocks=0,
               real_terms=0,lcm_leaf_bits=0,real_product_bits=0,forests=0,
               maximum_pending_bits=0)
    for a,b,c,sign,first,last in families:
        for lo in range(first,last+1,block_length):
            hi=min(last+1,lo+block_length)
            num,den,st=exact_weighted_block(a*lo+b,a,hi-lo,small_primes)
            shift=Jmax-(4*(hi-1)-(N-1)-c)
            pending.append((sign*num,den,shift))
            pending_mass+=num.bit_length()+den.bit_length()+4*(hi-lo)
            stats['maximum_pending_bits']=max(stats['maximum_pending_bits'],pending_mass)
            stats['blocks']+=1
            stats['real_terms']+=st['real_count']
            stats['lcm_leaf_bits']+=st['lcm_leaf_bits']
            stats['real_product_bits']+=st['real_product_bits']
            # Forests contain consecutive, disjoint intervals of one family.
            # Thus the SUM of shift spans at any tree level is O(N), as is
            # the sum of denominator bit lengths. A span bound alone would
            # be insufficient for arbitrary overlapping/orderless intervals.
            if pending_mass>=Jmax:
                V=(V+merge_nodes(pending,Jmax))&mask
                pending=[]
                pending_mass=0
                stats['forests']+=1
        if pending:
            V=(V+merge_nodes(pending,Jmax))&mask
            pending=[]
            pending_mass=0
            stats['forests']+=1
    return V,Jmax,stats


def extract_block_crt(N,G=192,window=12):
    """Complete pi experiment, retaining the existing fixed-guard convention.

    The weighted core and each block reconstruction are exact. Output carry
    certification is a separate issue, as in the existing v3/v5 engines.
    """
    if N<1 or G<8 or window<1 or (10**window).bit_length()+32>G:
        raise ValueError('require positive position/window and sufficient guard bits')
    start=time.perf_counter()
    X=mpz(5)**(N-1)
    K=cutoff(N,G,'pi')
    limit=X.bit_length()+G+8
    MOD=1<<G
    acc=0
    for sign,E,mu in bbp_pi_terms(N,K):
        if E>=0:
            if mu==1: continue
            r=pow(2,E,mu)*pow(5,N-1,mu)%mu
            acc=(acc+sign*((r<<G)//mu))%MOD
        elif -E<limit and mu>1:
            r=pow(5,N-1,mu)*pow((mu+1)>>1,-E,mu)%mu
            acc=(acc+sign*((r<<G)//mu))%MOD
    V,m,stats=block_core(N,G,x_bits=X.bit_length())
    W=(X*V)&((mpz(1)<<m)-1)
    dy=int(W>>(m-G)) if m>=G else int(W<<(G-m))
    acc=(acc+dy)%MOD
    stats['elapsed_seconds']=time.perf_counter()-start
    stats['engine']='weighted-block-crt'
    return format((acc*10**window)>>G,'0%dd'%window),stats


def check_lcm_lemma(A,a,length):
    product=mpz(1)
    actual_lcm=mpz(1)
    for j in range(length):
        mu=A+a*j
        product*=mu
        actual_lcm=math.lcm(actual_lcm,mu)
    factorial=math.factorial(length)
    odd_factorial=factorial//(factorial&-factorial)
    small_lcm=math.lcm(*range(1,length+1))
    odd_lcm=small_lcm//(small_lcm&-small_lcm)
    assert product%odd_factorial==0,(A,a,length)
    multiple=(product//odd_factorial)*odd_lcm
    assert multiple%actual_lcm==0,(A,a,length)


def selftest():
    rng=random.Random(20260906)
    lemma_cases=0
    for A in [1,3,9,101,1001,2**64+1,2**127+1]:
        for a in [2,4,8,16]:
            for length in range(1,61):
                check_lcm_lemma(A,a,length)
                lemma_cases+=1
    for _ in range(300):
        A=2*rng.randrange(0,10000)+1
        a=1<<rng.randrange(1,5)
        length=rng.randrange(1,60)
        num,den,_=exact_weighted_block(A,a,length)
        ref=sum((Fraction(1<<(4*(length-1-j)),A+a*j) for j in range(length)),Fraction())
        assert num*ref.denominator==den*ref.numerator,(A,a,length)
    positions=list(range(1,41))+[100,997,2500,10000]
    digits_reference=pi_reference(max(positions)+30)
    for N in positions:
        X=mpz(5)**(N-1)
        reference,_,_,_=build_v_euler(N,cutoff(N,192,'pi'),X.bit_length()+200)
        for length in ([1,7,None] if N<101 else [None]):
            V,m,_=block_core(N,block_length=length)
            assert V==reference,(N,length)
        digits,_=extract_block_crt(N)
        assert digits==digits_reference[N-1:N+11],(N,digits)
    print(json.dumps(dict(selftest='PASS',exact_weighted_blocks=300,
                         lcm_divisibility_cases=lemma_cases,
                         complete_V_positions=len(positions),
                         independent_digit_positions=len(positions))))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selftest',action='store_true')
    ap.add_argument('--positions',nargs='+',type=int)
    ap.add_argument('--pos',type=int)
    args=ap.parse_args()
    if args.selftest: selftest()
    if args.pos is not None:
        digits,stats=extract_block_crt(args.pos)
        print(json.dumps(dict(digits=digits,stats=stats)),flush=True)
    for N in args.positions or []:
        V,m,stats=block_core(N)
        X=mpz(5)**(N-1)
        ref,_,_,_=build_v_euler(N,cutoff(N,192,'pi'),X.bit_length()+200)
        assert V==ref,N
        stats['exact_V_match']=True
        print(json.dumps(stats),flush=True)
