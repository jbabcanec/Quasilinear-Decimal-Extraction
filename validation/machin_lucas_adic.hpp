#pragma once
// Experimental exact fixed-prime projection of the positive Euler sum.
// Global modulus caps are valid; no per-subtree precision bound is assumed.
#include "machin_lucas_common.hpp"

namespace machin_lucas {

class AdicEuler {
    struct Node {
        mpz_class p,q,t;
        u64 vp=0,vq=0,vt=0;
    };
    u64 A,a,c,mult,prime,exponent;
    mpz_class modulus;
    std::map<u64,mpz_class> powers;

    u64 strip(u64& value) const {
        u64 valuation=0;
        while(value%prime==0) { value/=prime; ++valuation; }
        return valuation;
    }
    void reduce_unit(mpz_class& value) const {
        if(value>=modulus)
            mpz_fdiv_r(value.get_mpz_t(),value.get_mpz_t(),modulus.get_mpz_t());
    }
    const mpz_class& power(u64 n) {
        auto it=powers.find(n);
        if(it==powers.end()) it=powers.emplace(n,powerz(prime,n)).first;
        return it->second;
    }
    void shift(mpz_class& value,u64 n) {
        if(!n || value==0) return;
        if(n>=exponent) { value=0; return; }
        value*=power(n); reduce_unit(value);
    }
    mpz_class inverse_unit(const mpz_class& value) const {
        std::vector<u64> precision;
        for(u64 n=exponent;;n=(n+1)/2) {
            precision.push_back(n);
            if(n==1) break;
        }
        std::reverse(precision.begin(),precision.end());
        std::vector<mpz_class> moduli;
        moduli.reserve(precision.size()); moduli.push_back(integer(prime));
        for(std::size_t i=1;i<precision.size();++i) {
            mpz_class next=moduli.back()*moduli.back();
            if(precision[i]&1)
                mpz_divexact_ui(next.get_mpz_t(),next.get_mpz_t(),static_cast<unsigned long>(prime));
            moduli.push_back(std::move(next));
        }
        std::vector<mpz_class> residues(moduli.size());
        residues.back()=modz(value,moduli.back());
        for(std::size_t i=residues.size()-1;i>0;--i)
            residues[i-1]=modz(residues[i],moduli[i-1]);
        mpz_class result=invz(residues.front(),moduli.front());
        for(std::size_t i=1;i<moduli.size();++i)
            result=modz(mpz_class(result*(2-residues[i]*result)),moduli[i]);
        return result;
    }
    Node build(u64 lo,u64 hi) {
        const u128 largest=static_cast<u128>(c)*(A+static_cast<u128>(a)*(hi-1));
        if((hi-lo<=2 && largest<(static_cast<u128>(1)<<62)) ||
           (hi-lo<=3 && largest<(static_cast<u128>(1)<<41)) ||
           (hi-lo<=4 && largest<(1ULL<<31))) {
            u128 p=1,q=1,t=0;
            for(u64 j=hi;j-->lo;) {
                const u128 pj=j ? static_cast<u128>(mult)*j : 1;
                const u128 qj=j ? static_cast<u128>(c)*(A+static_cast<u128>(a)*j) : 1;
                t=pj*(q+t); p*=pj; q*=qj;
            }
            Node result{integer128(p),integer128(q),integer128(t)};
            const mpz_class factor=integer(prime);
            result.vp=mpz_remove(result.p.get_mpz_t(),result.p.get_mpz_t(),factor.get_mpz_t());
            result.vq=mpz_remove(result.q.get_mpz_t(),result.q.get_mpz_t(),factor.get_mpz_t());
            result.vt=mpz_remove(result.t.get_mpz_t(),result.t.get_mpz_t(),factor.get_mpz_t());
            reduce_unit(result.p); reduce_unit(result.q); reduce_unit(result.t);
            return result;
        }
        if(hi-lo==1) {
            if(!lo) return {mpz_class(1),mpz_class(1),mpz_class(1),0,0,0};
            u64 pu=mult*lo,qu=A+a*lo,cu=c;
            const u64 vp=strip(pu),vq=strip(qu)+strip(cu);
            const mpz_class p=integer(pu),q=integer(qu)*integer(cu);
            return {p,q,p,vp,vq,vp};
        }
        const u64 mid=lo+(hi-lo)/2;
        auto left=build(lo,mid),right=build(mid,hi);
        const u64 lv=left.vt+right.vq,rv=left.vp+right.vt;
        const u64 vt=std::min(lv,rv);
        mpz_class t=left.t*right.q,other=left.p*right.t;
        shift(t,lv-vt); shift(other,rv-vt);
        t+=other; reduce_unit(t);
        left.p*=right.p; reduce_unit(left.p);
        left.q*=right.q; reduce_unit(left.q);
        return {std::move(left.p),std::move(left.q),std::move(t),
                left.vp+right.vp,left.vq+right.vq,vt};
    }
public:
    AdicEuler(u64 origin,u64 step,u64 C,u64 p,u64 precision)
        :A(origin),a(step),c(C/std::gcd(step,C)),mult(step/std::gcd(step,C)),
         prime(p),exponent(precision),modulus(powerz(p,precision)) {}

    Rational fraction(u64 count) {
        if(!exponent) return {0,1};
        auto root=build(0,count);
        if(root.vq<root.vt) return {0,1};
        const u64 pole=root.vq-root.vt;
        if(pole>exponent)
            throw std::runtime_error("fixed-prime Euler precision bound failed");
        mpz_class numerator=root.t*inverse_unit(root.q);
        reduce_unit(numerator);
        shift(numerator,exponent-pole);
        // Removing only powers of the selected prime is exact and avoids a
        // general big-integer gcd. This also gives the oracle canonical form.
        if(numerator==0) return {0,1};
        mpz_class unit;
        const mpz_class factor=integer(prime);
        const u64 cancellation=mpz_remove(unit.get_mpz_t(),numerator.get_mpz_t(),factor.get_mpz_t());
        return {std::move(unit),powerz(prime,exponent-cancellation)};
    }
};

static Rational fixed_prime_fraction(u64 A,u64 a,u64 C,u64 r,u64 prime=13) {
    if(!A || !r || C<2 || prime<3 || a<2 || (a&(a-1)) || !(A&1))
        throw std::invalid_argument("invalid fixed-prime Euler input");
    if(r==1) return {0,1};
    if(static_cast<u128>(A)+static_cast<u128>(a)*(r-1)>UINT32_MAX ||
       static_cast<u128>(a)*r>UINT32_MAX)
        throw std::overflow_error("fixed-prime Euler input exceeds audited range");
    u64 c=C/std::gcd(a,C),e=0;
    while(c%prime==0) { c/=prime; ++e; }
    if(!e || c!=1)
        throw std::invalid_argument("Euler scale is not a selected-prime power");
    // Every prefix valuation satisfies
    // v_p(product_{i=1}^j(A+a*i))-v_p(j!) <= floor(log_p(maximum)).
    // At each p-power the two counting functions differ by at most one.
    u64 maximum=A+a*(r-1),padding=0;
    while(maximum>=prime) { maximum/=prime; ++padding; }
    const u64 precision=e*(r-1)+padding;
    return AdicEuler(A,a,C,prime,precision).fraction(r);
}

} // namespace machin_lucas
