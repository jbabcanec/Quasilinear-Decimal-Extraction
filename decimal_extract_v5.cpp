// Experimental Euler-transformed BBP decimal extraction of pi.
// Same O(N log^3 N) bit bound as v3; fewer terms and no materialized band.
// Build: g++ -O3 -march=native -std=c++17 decimal_extract_v5.cpp -lgmpxx -lgmp -o extract_v5
// Run: extract_v5 1000000 --timing
// Exact core checks: extract_v5 --check-core
// Derivation: notes/prime_crt_and_euler_2026-09-06.md

#include <gmpxx.h>
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
static u64 mulmod(u64 a,u64 b,u64 m) { return static_cast<u128>(a)*b%m; }
static u64 powmod(u64 b,u64 e,u64 m) {
    u64 r=1%m;
    b%=m;
    while(e) { if(e&1) r=mulmod(r,b,m); b=mulmod(b,b,m); e>>=1; }
    return r;
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
    mpz_class inv=1;
    for(mp_bitcnt_t precision=1;precision<bits;) {
        precision=std::min(2*precision,bits);
        mpz_class low;
        mpz_fdiv_r_2exp(low.get_mpz_t(),d.get_mpz_t(),precision);
        inv*=mpz_class(2)-low*inv;
        reduce(inv,precision);
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
        if(N<1 || N-1>ULONG_MAX || N>4294967295ULL)
            throw std::invalid_argument("position must be in 1..4294967295");
        const bool timing=argc>2 && std::strcmp(argv[2],"--timing")==0;
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
        mpz_class acc=0,t;
        start=Clock::now();
        // Same arithmetic and accumulation as decimal_extract_v3.cpp.
        // The band is implicit in four endpoint ranges.
        for(const auto& f:families) {
            for(i64 k=0;k<=K;++k) {
                const i64 E=f.e0-4*k;
                const u64 mu=f.a*k+f.b;
                if(mu==1) continue;
                u64 r=0;
                if(E>=0) r=mulmod(powmod(2,E,mu),powmod(5,n5,mu),mu);
                else if(-E<Jcut)
                    r=mulmod(powmod(5,n5,mu),powmod((mu+1)>>1,-E,mu),mu);
                if(r) {
                    t=from_u64(r)<<G;
                    mpz_class mz=from_u64(mu);
                    mpz_fdiv_q(t.get_mpz_t(),t.get_mpz_t(),mz.get_mpz_t());
                    if(f.sign>0) acc+=t; else acc-=t;
                }
            }
        }
        reduce(acc,G);
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
            "{\"engine\":\"v5-euler-cpp\",\"N\":%llu,\"pre\":%.6f,\"word\":%.6f,"
            "\"vphase\":%.6f,\"total\":%.6f,\"original_terms\":%llu,\"transformed_terms\":%llu}\n",
            static_cast<unsigned long long>(N),pre,word,vphase,pre+word+vphase,
            static_cast<unsigned long long>(core.old_terms),static_cast<unsigned long long>(core.new_terms));
        return 0;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"error: %s\n",e.what()); return 1;
    }
}
