// BENCHMARK CONTROL ONLY: common arithmetic optimizations; root source preserved.
// Euler decimal extraction with shared residues across BBP families.
// The guard and 2-adic core match the archived Euler implementation.
// Build: g++ -O3 -march=native -std=c++17 decimal_extract_shared_residues.cpp -lgmpxx -lgmp -o extract_shared_residues
// Run: extract_shared_residues 1000000 --timing
// Exact word-accumulator checks: extract_shared_residues --check-shared

#include <gmpxx.h>
#include "../../decimal_extract_weighted_words.hpp"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

using u64=std::uint64_t;
using i64=std::int64_t;
using u128=unsigned __int128;
using Clock=std::chrono::steady_clock;

static double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now()-start).count();
}
static mpz_class from_u64(u64 x) {
    mpz_class z;
    mpz_import(z.get_mpz_t(),1,1,sizeof(x),0,0,&x);
    return z;
}
static mpz_class from_i64(i64 x) {
    if (x>=0) return from_u64(static_cast<u64>(x));
    return -from_u64(static_cast<u64>(-(x+1))+1);
}
static u64 mulmod(u64 a,u64 b,u64 m) {
    // Current callers pass residues (or the small constant 8), so inputs are
    // at most UINT32_MAX whenever the modulus is. Keep a general fallback.
    if(m<=UINT32_MAX && a<=UINT32_MAX && b<=UINT32_MAX) return a*b%m;
    return static_cast<u128>(a)*b%m;
}
static u64 powmod(u64 b,u64 e,u64 m) {
    u64 r=1%m;
    b%=m;
    while(e) { if(e&1) r=mulmod(r,b,m); b=mulmod(b,b,m); e>>=1; }
    return r;
}
// Simultaneous left-to-right powering. Each bit requires one square;
// multiplication by 5 and signed powers of 2 use small-word operations.
// The supported extraction range keeps 5*(m-1) within uint64_t.
static u64 pow5_times_pow2(u64 n,i64 e,u64 m) {
    const u64 power=e<0 ? static_cast<u64>(-e) : static_cast<u64>(e);
    const u64 combined=n|power;
    if(!combined) return 1%m;
    u64 z=1;
    for(int bit=63-__builtin_clzll(combined);bit>=0;--bit) {
        z=mulmod(z,z,m);
        if((n>>bit)&1) z=(5*z)%m;
        if((power>>bit)&1) {
            if(e<0) { if(z&1) z+=m; z>>=1; }
            else { z+=z; if(z>=m) z-=m; }
        }
    }
    return z;
}
static void cap(mpz_class& z,mp_bitcnt_t bits) {
    // Preserve small signed integers to avoid full-width negative leaves.
    if(mpz_sizeinbase(z.get_mpz_t(),2)>bits)
        mpz_fdiv_r_2exp(z.get_mpz_t(),z.get_mpz_t(),bits);
}
static void reduce(mpz_class& z,mp_bitcnt_t bits) {
    mpz_fdiv_r_2exp(z.get_mpz_t(),z.get_mpz_t(),bits);
}
static mpz_class inverse_2exp(const mpz_class& d,mp_bitcnt_t bits) {
    if(!mpz_odd_p(d.get_mpz_t())) throw std::runtime_error("even denominator");
    mpz_class inv=1;
    for(mp_bitcnt_t precision=1;precision<bits;) {
        const mp_bitcnt_t previous=precision;
        precision=precision>=bits-precision ? bits : 2*precision;
        mpz_class low;
        mpz_fdiv_r_2exp(low.get_mpz_t(),d.get_mpz_t(),precision);
        // Only the newly exposed half of the Newton correction is needed.
        mpz_class correction=mpz_class(1)-low*inv;
        mpz_fdiv_q_2exp(correction.get_mpz_t(),correction.get_mpz_t(),previous);
        reduce(correction,precision-previous);
        correction*=inv; reduce(correction,precision-previous);
        correction<<=previous;
        inv+=correction;
    }
    return inv;
}

struct PQT { mpz_class p,q,t; };

static PQT euler_bs(i64 A,u64 a,u64 lo,u64 hi,mp_bitcnt_t bits) {
    const u64 prefix=lo==0 ? 0 : lo-1;
    const unsigned rate=5+__builtin_ctzll(a);
    const mp_bitcnt_t precision=bits-(rate*prefix-__builtin_popcountll(prefix));
    if(hi-lo==1) {
        if(lo==0) return {mpz_class(1),mpz_class(1),mpz_class(1)};
        mpz_class p=-from_u64(16*a*lo);
        mpz_class q=from_i64(15*(A-static_cast<i64>(a*lo)));
        cap(p,precision); cap(q,precision);
        return {p,q,p};
    }
    const u64 mid=lo+(hi-lo)/2;
    PQT l=euler_bs(A,a,lo,mid,bits);
    PQT r=euler_bs(A,a,mid,hi,bits);
    PQT out{l.p*r.p,l.q*r.q,l.t*r.q+l.p*r.t};
    cap(out.p,precision); cap(out.q,precision); cap(out.t,precision);
    return out;
}

struct Rational { mpz_class n,d; u64 terms; };
static Rational euler_rational(i64 A,u64 a,mp_bitcnt_t bits) {
    const unsigned rate=5+__builtin_ctzll(a);
    u64 R=std::max<u64>(1,bits/rate);
    auto valuation=[rate](u64 r) { return rate*r-__builtin_popcountll(r); };
    while(valuation(R)<bits) ++R;
    while(R>1 && valuation(R-1)>=bits) --R;
    PQT root=euler_bs(A,a,0,R,bits);
    Rational result{-root.t,from_i64(15*A)*root.q,R};
    cap(result.n,bits); cap(result.d,bits);
    return result;
}

struct Family { int sign; u64 a,b; i64 e0; };

static mpz_class original_words(i64 n,i64 K,i64 cutoff,unsigned guard,bool joint=false) {
    const std::vector<Family> families{{1,8,1,n+2},{-1,2,1,n-1},
                                       {-1,8,5,n},{-1,4,3,n-1}};
    mpz_class acc=0;
    for(const auto& f:families) for(i64 k=0;k<=K;++k) {
        const i64 e=f.e0-4*k;
        const u64 mu=f.a*k+f.b;
        if(mu==1 || e<=-cutoff) continue;
        const u64 base=e>=0 ? 2 : (mu+1)/2;
        const u64 power=e>=0 ? e : -e;
        const u64 r=joint ? pow5_times_pow2(n,e,mu)
                          : mulmod(powmod(5,n,mu),powmod(base,power,mu),mu);
        mpz_class t=from_u64(r)<<guard;
        mpz_class d=from_u64(mu);
        mpz_fdiv_q(t.get_mpz_t(),t.get_mpz_t(),d.get_mpz_t());
        if(f.sign>0) acc+=t; else acc-=t;
    }
    reduce(acc,guard);
    return acc;
}

static mpz_class shared_words(i64 n,i64 K,i64 cutoff,unsigned guard,bool joint=true) {
    // Same denominator sharing, joint powering, and per-term floors.
    // The header checks maximum<=UINT32_MAX before its 64-bit products.
    if(joint) return weighted_words::accumulate(static_cast<u64>(n)+1,K,cutoff,guard);
    const std::vector<Family> families{{1,8,1,n+2},{-1,2,1,n-1},
                                       {-1,8,5,n},{-1,4,3,n-1}};
    u64 ends[4],maximum=0;
    for(unsigned i=0;i<4;++i) {
        const i64 last=std::min(K,(families[i].e0+cutoff-1)/4);
        ends[i]=last<0 ? 0 : families[i].a*last+families[i].b;
        maximum=std::max(maximum,ends[i]);
    }
    mpz_class acc=0,t,d;
    for(u64 mu=3;mu<=maximum;mu+=2) {
        const bool first=(mu%8==1 && mu<=ends[0]);
        const bool third=(mu%8==5 && mu<=ends[2]);
        const bool fourth=(mu%4==3 && mu<=ends[3]);
        const bool second=mu<=ends[1];
        if(!(first||second||third||fourth)) continue;
        d=from_u64(mu);
        auto add=[&](u64 r,bool positive) {
            if(!r) return;
            t=from_u64(r)<<guard;
            mpz_fdiv_q(t.get_mpz_t(),t.get_mpz_t(),d.get_mpz_t());
            if(positive) acc+=t; else acc-=t;
        };
        if(joint && (unsigned(first)+unsigned(second)+unsigned(third)+unsigned(fourth)==1)) {
            const i64 q=static_cast<i64>((mu+1)/2);
            const i64 exponent=n+3-(second ? 4*q : fourth ? 2*q : q);
            add(pow5_times_pow2(n,exponent,mu),first);
            continue;
        }
        const u64 q=(mu+1)/2;
        const u64 T=joint ? pow5_times_pow2(n,n,mu) : powmod(10,n,mu);
        const u64 Y=joint ? pow5_times_pow2(0,-static_cast<i64>(q),mu)
                          : powmod(q,q,mu);
        const u64 factor=mulmod(8,T,mu);
        if(first||third) add(mulmod(factor,Y,mu),first);
        if(second||fourth) {
            const u64 y2=mulmod(Y,Y,mu);
            if(fourth) add(mulmod(factor,y2,mu),false);
            if(second) add(mulmod(factor,mulmod(y2,y2,mu),mu),false);
        }
    }
    reduce(acc,guard);
    return acc;
}

static void check_shared() {
    u64 checked=0,residues=0;
    for(u64 n:{0ULL,1ULL,7ULL,32ULL,9999ULL,999999ULL})
        for(i64 e:{0LL,1LL,-1LL,6LL,-6LL,127LL,-127LL,1000000000LL,-1000000000LL})
            for(u64 mu:{3ULL,5ULL,7ULL,9ULL,15ULL,25ULL,45ULL,97ULL,127ULL,
                        257ULL,65535ULL,1000001ULL,3321928221ULL}) {
                const u64 base=e<0 ? (mu+1)/2 : 2;
                const u64 power=e<0 ? -e : e;
                const u64 expected=mulmod(powmod(5,n,mu),powmod(base,power,mu),mu);
                if(pow5_times_pow2(n,e,mu)!=expected)
                    throw std::runtime_error("joint exponentiation mismatch");
                ++residues;
            }
    for(i64 n=0;n<=128;++n) for(unsigned g:{8u,31u,64u,192u}) {
        const i64 K=n+17;
        const i64 cutoff=3*n+g+8;
        const mpz_class original=original_words(n,K,cutoff,g);
        if(original!=shared_words(n,K,cutoff,g) || original!=shared_words(n,K,cutoff,g,false)
           || original!=original_words(n,K,cutoff,g,true))
            throw std::runtime_error("shared accumulator mismatch");
        ++checked;
    }
    for(i64 n:{996,9999,99999}) {
        const i64 K=n+71, cutoff=3*n+200;
        const mpz_class original=original_words(n,K,cutoff,192);
        if(original!=shared_words(n,K,cutoff,192) || original!=shared_words(n,K,cutoff,192,false)
           || original!=original_words(n,K,cutoff,192,true))
            throw std::runtime_error("large shared accumulator mismatch");
        ++checked;
    }
    std::printf("{\"shared_word_selftest\":\"PASS\",\"exact_accumulators\":%llu,"
                "\"modes_per_accumulator\":3,\"individual_residues\":%llu}\n",
                static_cast<unsigned long long>(checked),
                static_cast<unsigned long long>(residues));
}

struct Core { mpz_class V; mp_bitcnt_t bits; u64 old_terms,new_terms; };
static Core build_core(const std::vector<Family>& families,i64 K,i64 Jcut) {
    struct Range { Family f; i64 last; mp_bitcnt_t bits; };
    std::vector<Range> ranges;
    mp_bitcnt_t Jmax=1;
    u64 old_terms=0,new_terms=0;
    for(const auto& f:families) {
        // e0 may be -1 at N=1. C++ integer division truncates toward zero.
        const i64 first=std::max<i64>(0,f.e0>=0 ? f.e0/4+1 : 0);
        const i64 last=std::min(K,(Jcut-1+f.e0)/4);
        if(last<first) continue;
        const auto bits=static_cast<mp_bitcnt_t>(4*last-f.e0);
        ranges.push_back({f,last,bits});
        old_terms+=last-first+1;
        Jmax=std::max(Jmax,bits);
    }
    mpz_class num=0,den=1;
    for(const auto& r:ranges) {
        const i64 A=static_cast<i64>(r.f.a*r.last+r.f.b);
        Rational value=euler_rational(A,r.f.a,r.bits);
        mpz_class contribution=value.n*den;
        contribution<<=Jmax-r.bits;
        if(r.f.sign<0) contribution=-contribution;
        num=num*value.d+contribution;
        den*=value.d;
        reduce(num,Jmax); reduce(den,Jmax);
        new_terms+=value.terms;
    }
    mpz_class V=num*inverse_2exp(den,Jmax);
    reduce(V,Jmax);
    return {std::move(V),Jmax,old_terms,new_terms};
}

static void check_core() {
    u64 checks=0;
    for(u64 a:{2u,4u,8u,16u})
        for(i64 A:{1,9,73,769})
            for(mp_bitcnt_t bits:{1u,7u,32u,97u,256u,511u}) {
                Rational r=euler_rational(A,a,bits);
                mpz_class got=r.n*inverse_2exp(r.d,bits);
                reduce(got,bits);
                mpz_class ref=0;
                for(u64 j=0;4*j<bits;++j) {
                    mpz_class d=from_i64(A-static_cast<i64>(a*j));
                    ref+=inverse_2exp(d,bits)<<(4*j);
                }
                reduce(ref,bits);
                if(got!=ref) throw std::runtime_error("Euler core mismatch");
                ++checks;
            }
    std::printf("{\"core_selftest\":\"PASS\",\"cases\":%llu}\n",
                static_cast<unsigned long long>(checks));
}

int main(int argc,char** argv) {
    try {
        if(argc>1 && std::strcmp(argv[1],"--check-shared")==0) {
            check_shared(); return 0;
        }
        if(argc>1 && std::strcmp(argv[1],"--check-core")==0) {
            check_core(); return 0;
        }
        u64 N=1000000;
        if(argc>1) {
            std::size_t used=0;
            const std::string arg=argv[1];
            if(arg.empty() || arg[0]=='-') throw std::invalid_argument("invalid position");
            N=std::stoull(arg,&used);
            if(used!=arg.size()) throw std::invalid_argument("invalid position");
        }
        if(N<1 || N-1>ULONG_MAX || N>500000000ULL)
            throw std::invalid_argument("position must be in 1..500000000");
        bool timing=false,joint=true,joint_only=false;
        for(int i=2;i<argc;++i) {
            if(std::strcmp(argv[i],"--timing")==0) timing=true;
            else if(std::strcmp(argv[i],"--shared-only")==0) joint=false;
            else if(std::strcmp(argv[i],"--joint-only")==0) joint_only=true;
            else throw std::invalid_argument("unknown option");
        }
        const int G=192,window=12;
        const i64 n5=static_cast<i64>(N-1);
        const i64 K=static_cast<i64>((n5*std::log2(10.0)+G+24)/4)+2;
        auto start=Clock::now();
        mpz_class X;
        mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(n5));
        const i64 Jcut=mpz_sizeinbase(X.get_mpz_t(),2)+G+8;
        const double pre=elapsed(start);
        std::vector<Family> families{{1,8,1,n5+2},{-1,2,1,n5-1},
                                     {-1,8,5,n5},{-1,4,3,n5-1}};
        mpz_class acc=0;
        start=Clock::now();
        acc=joint_only ? original_words(n5,K,Jcut,G,true) : shared_words(n5,K,Jcut,G,joint);
        const double word=elapsed(start);
        start=Clock::now();
        Core core=build_core(families,K,Jcut);
        mpz_class W=X*core.V;
        reduce(W,core.bits);
        if(core.bits>=G) W>>=core.bits-G; else W<<=G-core.bits;
        acc+=W; reduce(acc,G);
        const double vphase=elapsed(start);
        mpz_class tp;
        mpz_ui_pow_ui(tp.get_mpz_t(),10,window);
        acc*=tp; acc>>=G;
        gmp_printf("%012Zd\n",acc.get_mpz_t());
        if(timing) std::fprintf(stderr,
            "{\"engine\":\"%s\",\"N\":%llu,\"pre\":%.6f,\"word\":%.6f,"
            "\"vphase\":%.6f,\"total\":%.6f,\"original_terms\":%llu,\"transformed_terms\":%llu}\n",
            joint_only ? "euler-joint" : joint ? "euler-shared-joint" : "euler-shared-denominators",
            static_cast<unsigned long long>(N),pre,word,vphase,pre+word+vphase,
            static_cast<unsigned long long>(core.old_terms),static_cast<unsigned long long>(core.new_terms));
        return 0;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"error: %s\n",e.what()); return 1;
    }
}
