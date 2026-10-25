"""Exact small witnesses for the integer correction in weighted prime CRT.

This validator uses slow rational arithmetic deliberately. It is not the
proposed evaluator; it checks why the unweighted algorithm does not solve V.
"""
from fractions import Fraction
import json
import random

from harmonic_prime_crt import primes_to


def check(weights):
    K=len(weights)
    B=8*K+1
    S=sum((Fraction(w,8*k+1) for k,w in enumerate(weights,1)),Fraction())
    T=Fraction()
    for p in primes_to(B):
        q=p
        while q<=B//p: q*=p
        c=0
        for k,w in enumerate(weights,1):
            u=8*k+1
            if u%p: continue
            pv=1
            while u%p==0:
                pv*=p
                u//=p
            c=(c+(q//pv)*(w%q)*pow(u,-1,q))%q
        T+=Fraction(c,q)
    correction=S-T
    assert correction.denominator==1
    m=4*K+8
    mod=1<<m
    direct=(S.numerator*pow(S.denominator,-1,mod))%mod
    compressed=(T.numerator*pow(T.denominator,-1,mod))%mod
    assert (direct-compressed)%mod==correction.numerator%mod
    return dict(K=K,maximum_weight_bits=max(abs(w).bit_length() for w in weights),
                correction_bits=abs(correction.numerator).bit_length(),
                correction_is_integer=True)


if __name__=='__main__':
    rng=random.Random(20260906)
    for _ in range(100):
        K=rng.randrange(1,35)
        check([rng.randrange(-10000,10001) for _ in range(K)])
    print(json.dumps(dict(random_weighted_identities='PASS',cases=100)))
    for K in [10,20,40,80,160]:
        print(json.dumps(check([1<<(4*(K-k)) for k in range(1,K+1)])))
