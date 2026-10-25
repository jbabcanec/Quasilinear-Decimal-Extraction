// Shared BBP word residues. Each modulus is checked to fit uint32_t.
// The returned fixed-point accumulator agrees exactly with the original
// per-family convention: floor each positive residue, then apply its sign.
#pragma once
#include <gmpxx.h>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

namespace weighted_words {
using word = std::uint64_t;
using signed_word = std::int64_t;
static word multiply(word a, word b, word m) { return a*b%m; }
static word joint_power(word n, signed_word exponent, word m) {
    const word power=exponent<0 ? word(-exponent) : word(exponent);
    const word mask=n|power;
    if(!mask) return 1%m;
    word z=1;
    for(int bit=63-__builtin_clzll(mask);bit>=0;--bit) {
        z=multiply(z,z,m);
        if((n>>bit)&1) z=(5*z)%m;
        if((power>>bit)&1) {
            if(exponent<0) { if(z&1) z+=m; z>>=1; }
            else { z+=z; if(z>=m) z-=m; }
        }
    }
    return z;
}
static mpz_class accumulate(word N, signed_word K, signed_word limit,
                            word guard, word* rounding_terms=nullptr) {
    const signed_word n=static_cast<signed_word>(N-1);
    struct Family { word a,b; signed_word e0; };
    const Family families[]={{8,1,n+2},{2,1,n-1},{8,5,n},{4,3,n-1}};
    word ends[4],maximum=0;
    for(unsigned i=0;i<4;++i) {
        const signed_word last=std::min(K,(families[i].e0+limit-1)/4);
        ends[i]=last<0 ? 0 : families[i].a*last+families[i].b;
        maximum=std::max(maximum,ends[i]);
    }
    if(maximum>UINT32_MAX) throw std::overflow_error("word modulus exceeds uint32_t");
    mpz_class acc=0,t;
    word count=0;
    for(word mu=3;mu<=maximum;mu+=2) {
        const bool first=mu%8==1 && mu<=ends[0];
        const bool second=mu<=ends[1];
        const bool third=mu%8==5 && mu<=ends[2];
        const bool fourth=mu%4==3 && mu<=ends[3];
        if(!(first||second||third||fourth)) continue;
        auto add=[&](word r,bool positive) {
            if(!r) return;
            // All residues and denominators fit even a 32-bit unsigned long.
            mpz_set_ui(t.get_mpz_t(),static_cast<unsigned long>(r));
            mpz_mul_2exp(t.get_mpz_t(),t.get_mpz_t(),guard);
            mpz_fdiv_q_ui(t.get_mpz_t(),t.get_mpz_t(),static_cast<unsigned long>(mu));
            if(positive) acc+=t; else acc-=t;
            ++count;
        };
        if(unsigned(first)+unsigned(second)+unsigned(third)+unsigned(fourth)==1) {
            const signed_word q=static_cast<signed_word>((mu+1)/2);
            add(joint_power(n,n+3-(second ? 4*q : fourth ? 2*q : q),mu),first);
        } else {
            const word q=(mu+1)/2;
            const word T=joint_power(n,n,mu),Y=joint_power(0,-static_cast<signed_word>(q),mu);
            const word factor=multiply(8,T,mu);
            if(first||third) add(multiply(factor,Y,mu),first);
            if(second||fourth) {
                const word y2=multiply(Y,Y,mu);
                if(fourth) add(multiply(factor,y2,mu),false);
                if(second) add(multiply(factor,multiply(y2,y2,mu),mu),false);
            }
        }
    }
    mpz_fdiv_r_2exp(acc.get_mpz_t(),acc.get_mpz_t(),guard);
    if(rounding_terms) *rounding_terms=count;
    return acc;
}
} // namespace weighted_words
