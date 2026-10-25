// Exact weighted-block CRT core with real Euler correction.
// Native implementation of validation/weighted_block_crt.py.
// The decimal output retains the experimental G=192 convention; the
// complete certified word aggregation/output procedure is in the paper.
// Build: g++ -O3 -march=native -std=c++17 decimal_extract_weighted.cpp -lgmpxx -lgmp -o extract_weighted
#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = unsigned __int128;
using Clock = std::chrono::steady_clock;
static double elapsed(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now()-t).count();
}
static mpz_class integer(u64 x) {
    mpz_class z;
    mpz_import(z.get_mpz_t(),1,1,sizeof(x),0,0,&x);
    return z;
}
static u64 bits(u64 x) { return x ? 64-__builtin_clzll(x) : 0; }
static u64 bits(const mpz_class& x) {
    return x==0 ? 0 : mpz_sizeinbase(x.get_mpz_t(),2);
}
static mp_bitcnt_t bitcount(u64 n) {
    if(n>std::numeric_limits<mp_bitcnt_t>::max())
        throw std::overflow_error("GMP bit-count limit exceeded");
    return static_cast<mp_bitcnt_t>(n);
}
static void reduce(mpz_class& x,u64 b) {
    mpz_fdiv_r_2exp(x.get_mpz_t(),x.get_mpz_t(),bitcount(b));
}
static void cap(mpz_class& x,u64 b) {
    if(bits(x)>b) reduce(x,b);
}
static mpz_class floor_div(const mpz_class& x,const mpz_class& y) {
    mpz_class z;
    mpz_fdiv_q(z.get_mpz_t(),x.get_mpz_t(),y.get_mpz_t());
    return z;
}
static u64 mulmod(u64 x,u64 y,u64 m) { return static_cast<u128>(x)*y%m; }
static u64 powmod(u64 x,u64 e,u64 m) {
    u64 r=1%m; x%=m;
    while(e) { if(e&1) r=mulmod(r,x,m); x=mulmod(x,x,m); e>>=1; }
    return r;
}
static u64 invmod(u64 x,u64 m) {
    // Inputs in this executable are below 2^32. Signed 128-bit
    // coefficients make the intermediate Euclidean updates explicit.
    __int128 t=0,nt=1;
    u64 r=m,nr=x%m;
    while(nr) {
        u64 q=r/nr;
        u64 next=r-q*nr; r=nr; nr=next;
        __int128 next_t=t-static_cast<__int128>(q)*nt;
        t=nt; nt=next_t;
    }
    if(r!=1) throw std::runtime_error("noninvertible local denominator");
    t%=m; if(t<0) t+=m;
    return static_cast<u64>(t);
}
static mpz_class inverse_2exp(const mpz_class& d,u64 b) {
    if(!mpz_odd_p(d.get_mpz_t())) throw std::runtime_error("even denominator");
    mpz_class inv=1;
    for(u64 precision=1;precision<b;) {
        // Saturating doubling avoids the 32-bit mp_bitcnt_t overflow.
        precision=precision>=b-precision ? b : 2*precision;
        mpz_class low=d; reduce(low,precision);
        inv*=mpz_class(2)-low*inv; reduce(inv,precision);
    }
    return inv;
}
static std::vector<u64> primes_to(u64 bound) {
    std::vector<bool> prime(bound+1,true);
    prime[0]=false;
    if(bound>=1) prime[1]=false;
    for(u64 p=2;p<=bound/p;++p)
        if(prime[p]) for(u64 j=p*p;j<=bound;j+=p) prime[j]=false;
    std::vector<u64> out;
    for(u64 p=3;p<=bound;p+=2) if(prime[p]) out.push_back(p);
    return out;
}
static u64 isqrt(u64 n) {
    u64 r=static_cast<u64>(std::sqrt(static_cast<long double>(n)));
    while(r+1<=n/(r+1)) ++r;
    while(r && r>n/r) --r;
    return r;
}
struct Rational { mpz_class n,d; };
struct PQT { mpz_class p,q,t; };
struct Stats {
    u64 N=0,Jmax=0,terms=0,block_length=0,blocks=0,real_terms=0;
    u64 lcm_leaf_bits=0,real_product_bits=0,forests=0,maximum_pending_bits=0;
    double factor=0,coefficients=0,real=0,rational=0,assembly=0;
};
static PQT real_bs(u64 A,u64 a,u64 lo,u64 hi) {
    if(hi-lo==1) {
        if(lo==0) return {mpz_class(1),mpz_class(1),mpz_class(1)};
        mpz_class p=-integer(a*lo),q=integer(15*(A+a*lo));
        return {p,q,p};
    }
    u64 mid=lo+(hi-lo)/2;
    PQT l=real_bs(A,a,lo,mid),r=real_bs(A,a,mid,hi);
    return {l.p*r.p,l.q*r.q,l.t*r.q+l.p*r.t};
}
static mpz_class real_approx(u64 A,u64 a,u64 length,u64 guard,Stats& stats) {
    const u64 ratio=15*A/(a*length);
    const u64 log_ratio=std::max<u64>(1,bits(ratio)>0 ? bits(ratio)-1 : 0);
    u64 count=std::max<u64>(1,(4*length+16+log_ratio-1)/log_ratio);
    for(;;) {
        PQT root=real_bs(A,a,0,count);
        mpz_class error_num=abs(root.p)*integer(a*count);
        mpz_class error_den=root.q*integer(225*A)*integer(A+a*count);
        error_num<<=bitcount(4*length+6);
        if(error_num<error_den) {
            stats.real_terms+=count;
            stats.real_product_bits+=bits(root.q);
            root.t<<=bitcount(4*length+guard);
            return floor_div(root.t,integer(15*A)*root.q);
        }
        if(count>std::numeric_limits<u64>::max()/2)
            throw std::overflow_error("real-series count overflow");
        count*=2;
    }
}
static Rational rational_tree(std::vector<Rational> nodes) {
    if(nodes.empty()) return {mpz_class(0),mpz_class(1)};
    while(nodes.size()>1) {
        std::vector<Rational> next;
        next.reserve((nodes.size()+1)/2);
        for(std::size_t i=0;i+1<nodes.size();i+=2) {
            auto& l=nodes[i]; auto& r=nodes[i+1];
            next.push_back({l.n*r.d+r.n*l.d,l.d*r.d});
        }
        if(nodes.size()&1) next.push_back(std::move(nodes.back()));
        nodes=std::move(next);
    }
    return std::move(nodes.front());
}
static Rational weighted_block(u64 A,u64 a,u64 length,
                               const std::vector<u64>& primes,Stats& stats) {
    if(!A || !(A&1) || a<2 || (a&(a-1)) || !length)
        throw std::invalid_argument("require odd A, power-of-two a>=2 and positive length");
    const u64 maximum=A+a*(length-1);
    if(maximum>UINT32_MAX || length>500000000ULL)
        throw std::invalid_argument("block exceeds the audited word range");
    auto start=Clock::now();
    std::vector<u64> remainder(length);
    std::vector<std::vector<std::pair<u64,u64>>> factors(length);
    struct Prime {
        u64 q=1,c=0,weight=0,step=0,last=0;
        bool seen=false,cached=false;
    };
    std::unordered_map<u64,Prime> powers;
    for(u64 j=0;j<length;++j) remainder[j]=A+a*j;
    for(u64 p:primes) {
        if(p>maximum/p) break;
        u64 first=mulmod((p-A%p)%p,invmod(a,p),p);
        for(u64 j=first;j<length;j+=p) {
            u64 pv=1;
            while(remainder[j]%p==0) { remainder[j]/=p; pv*=p; }
            if(pv>1) {
                factors[j].push_back({p,pv});
                auto& item=powers[p]; item.q=std::max(item.q,pv);
            }
        }
    }
    for(u64 j=0;j<length;++j) if(remainder[j]>1) {
        u64 p=remainder[j];
        factors[j].push_back({p,p});
        auto& item=powers[p]; item.q=std::max(item.q,p);
    }
    stats.factor+=elapsed(start); start=Clock::now();
    for(u64 j=0;j<length;++j) {
        const u64 mu=A+a*j;
        for(auto [p,pv]:factors[j]) {
            auto& item=powers.at(p); const u64 q=item.q;
            if(!item.seen) {
                item.weight=powmod(16,length-1-j,q);
                item.cached=(length-1-j)/p>=2;
                if(item.cached) item.step=powmod(invmod(16,q),p,q);
                item.seen=true;
            } else if(item.cached) {
                // Consecutive denominators divisible by p are p indices apart:
                // 16^(length-1-(j+p)) = 16^(-p)*16^(length-1-j) mod q.
                if(j-item.last!=p) throw std::runtime_error("nonconsecutive prime incidence");
                item.weight=mulmod(item.weight,item.step,q);
            } else {
                item.weight=powmod(16,length-1-j,q);
            }
            item.last=j;
            u64 term=mulmod(q/pv,item.weight,q);
            term=mulmod(term,invmod(mu/pv,q),q);
            item.c=(item.c+term)%q;
        }
    }
    const u64 guard=bits(maximum+length+1)+8;
    if(guard>=64) throw std::overflow_error("local guard exceeds word accumulator");
    mpz_class lower=0;
    std::vector<Rational> parts;
    parts.reserve(powers.size());
    // Ordered leaves make assembly and statistics reproducible across runs.
    std::vector<u64> keys;
    keys.reserve(powers.size());
    for(const auto& entry:powers) keys.push_back(entry.first);
    std::sort(keys.begin(),keys.end());
    for(u64 p:keys) {
        const auto item=powers.at(p);
        stats.lcm_leaf_bits+=bits(item.q);
        if(item.c) {
            parts.push_back({integer(item.c),integer(item.q)});
            lower+=integer(static_cast<u64>((static_cast<u128>(item.c)<<guard)/item.q));
        }
    }
    stats.coefficients+=elapsed(start); start=Clock::now();
    mpz_class approx=real_approx(A,a,length,guard,stats);
    mpz_class correction=approx-lower+(mpz_class(1)<<bitcount(guard-1));
    mpz_fdiv_q_2exp(correction.get_mpz_t(),correction.get_mpz_t(),bitcount(guard));
    stats.real+=elapsed(start); start=Clock::now();
    Rational result=rational_tree(std::move(parts));
    result.n+=correction*result.d;
    stats.rational+=elapsed(start);
    return result;
}
struct Node { mpz_class n,d; u64 shift; };
static mpz_class merge_nodes(std::vector<Node> nodes,u64 b) {
    while(nodes.size()>1) {
        std::vector<Node> next; next.reserve((nodes.size()+1)/2);
        for(std::size_t i=0;i+1<nodes.size();i+=2) {
            auto& l=nodes[i]; auto& r=nodes[i+1];
            u64 shift=std::min(l.shift,r.shift),precision=b-shift;
            mpz_class nl=l.n*r.d,nr=r.n*l.d;
            nl<<=bitcount(l.shift-shift); nr<<=bitcount(r.shift-shift);
            mpz_class n=nl+nr,d=l.d*r.d;
            cap(n,precision); cap(d,precision);
            next.push_back({std::move(n),std::move(d),shift});
        }
        if(nodes.size()&1) next.push_back(std::move(nodes.back()));
        nodes=std::move(next);
    }
    auto& root=nodes.front();
    mpz_class value=root.n*inverse_2exp(root.d,b-root.shift);
    reduce(value,b-root.shift); value<<=bitcount(root.shift);
    return value;
}
struct Family { u64 a,b; i64 c; int sign; };
static const Family families[]={{8,1,2,1},{2,1,-1,-1},{8,5,0,-1},{4,3,-1,-1}};
static i64 cutoff(u64 N,u64 G) {
    return static_cast<i64>(((N-1)*std::log2(10.0)+G+24)/4)+2;
}
static mpz_class block_core(u64 N,u64 G,u64 x_bits,Stats& stats,u64 length=0) {
    struct Range { Family f; i64 first,last; };
    std::vector<Range> ranges;
    const i64 K=cutoff(N,G),limit=static_cast<i64>(x_bits+G+8);
    u64 Jmax=1,max_mu=1;
    stats.N=N;
    for(auto f:families) {
        const i64 e0=static_cast<i64>(N-1)+f.c;
        const i64 first=std::max<i64>(0,e0>=0 ? e0/4+1 : 0);
        const i64 last=std::min(K,(limit-1+e0)/4);
        if(last<first) continue;
        Jmax=std::max<u64>(Jmax,4*last-e0);
        max_mu=std::max<u64>(max_mu,f.a*last+f.b);
        stats.terms+=last-first+1;
        ranges.push_back({f,first,last});
    }
    if(!length) {
        u64 L=std::max<u64>(1,bits(N-1)),u=isqrt(L);
        if(u*u<L) ++u;
        length=std::max<u64>(1,(N+(1ULL<<u)-1)/(1ULL<<u));
    }
    stats.Jmax=Jmax; stats.block_length=length;
    auto primes=primes_to(isqrt(max_mu));
    std::vector<Node> pending;
    u64 mass=0;
    mpz_class V=0;
    auto flush=[&]() {
        if(pending.empty()) return;
        auto t=Clock::now();
        V+=merge_nodes(std::move(pending),Jmax); reduce(V,Jmax);
        pending.clear(); mass=0; ++stats.forests;
        stats.assembly+=elapsed(t);
    };
    for(const auto& range:ranges) {
        const auto f=range.f;
        for(i64 lo=range.first;lo<=range.last;lo+=length) {
            const i64 hi=std::min<i64>(range.last+1,lo+length);
            auto value=weighted_block(f.a*lo+f.b,f.a,hi-lo,primes,stats);
            const u64 shift=Jmax-(4*(hi-1)-static_cast<i64>(N-1)-f.c);
            mass+=bits(value.n)+bits(value.d)+4*(hi-lo);
            stats.maximum_pending_bits=std::max(stats.maximum_pending_bits,mass);
            if(f.sign<0) value.n=-value.n;
            pending.push_back({std::move(value.n),std::move(value.d),shift});
            ++stats.blocks;
            if(mass>=Jmax) flush();
        }
        flush();
    }
    return V;
}
static void print_stats(const Stats& s,double pre,double word,double core,double total) {
    std::fprintf(stderr,
        "{\"engine\":\"weighted-block-cpp\",\"N\":%llu,\"Jmax\":%llu,\"terms\":%llu,"
        "\"block_length\":%llu,\"blocks\":%llu,\"forests\":%llu,\"real_terms\":%llu,"
        "\"lcm_leaf_bits\":%llu,\"real_product_bits\":%llu,\"maximum_pending_bits\":%llu,"
        "\"pre\":%.9f,\"word\":%.9f,\"core\":%.9f,\"total\":%.9f,"
        "\"factor\":%.9f,\"coefficients\":%.9f,\"real\":%.9f,\"rational\":%.9f,\"assembly\":%.9f}\n",
        (unsigned long long)s.N,(unsigned long long)s.Jmax,(unsigned long long)s.terms,
        (unsigned long long)s.block_length,(unsigned long long)s.blocks,
        (unsigned long long)s.forests,(unsigned long long)s.real_terms,
        (unsigned long long)s.lcm_leaf_bits,(unsigned long long)s.real_product_bits,
        (unsigned long long)s.maximum_pending_bits,pre,word,core,total,
        s.factor,s.coefficients,s.real,s.rational,s.assembly);
}
static void selftest() {
    std::mt19937_64 rng(20260907);
    for(unsigned i=0;i<300;++i) {
        u64 A=2*(rng()%10000)+1,a=1ULL<<(1+rng()%4),length=1+rng()%59;
        Stats st;
        auto r=weighted_block(A,a,length,primes_to(isqrt(A+a*(length-1))),st);
        mpq_class ref=0;
        for(u64 j=0;j<length;++j)
            ref+=mpq_class(mpz_class(1)<<bitcount(4*(length-1-j)),integer(A+a*j));
        if(r.n*ref.get_den()!=r.d*ref.get_num())
            throw std::runtime_error("exact weighted block mismatch");
    }
    unsigned core_cases=0;
    for(u64 N=1;N<=40;++N) {
        mpz_class X; mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(N-1));
        for(u64 length:{0ULL,1ULL,7ULL}) {
            Stats st; mpz_class got=block_core(N,192,bits(X),st,length),ref=0;
            for(auto f:families) {
                for(i64 k=0;k<=cutoff(N,192);++k) {
                    const i64 E=static_cast<i64>(N-1)+f.c-4*k;
                    if(E>=0 || -E>=static_cast<i64>(bits(X)+200)) continue;
                    mpz_class term=inverse_2exp(integer(f.a*k+f.b),st.Jmax);
                    term<<=bitcount(st.Jmax+E);
                    if(f.sign>0) ref+=term; else ref-=term;
                }
            }
            reduce(ref,st.Jmax);
            if(got!=ref) throw std::runtime_error("full modular core mismatch");
            ++core_cases;
        }
    }
    std::printf("{\"selftest\":\"PASS\",\"exact_rational_blocks\":300,\"full_core_cases\":%u}\n",core_cases);
}
int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--selftest") { selftest(); return 0; }
        if((argc==3 || argc==4) && std::string(argv[1])=="--core") {
            u64 N=std::stoull(argv[2]);
            if(N<1 || N>500000000ULL) throw std::invalid_argument("position outside 1..500000000");
            u64 length=argc>=4 ? std::stoull(argv[3]) : 0;
            if(length>500000000ULL) throw std::invalid_argument("block length exceeds supported range");
            mpz_class X; mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(N-1));
            Stats st; mpz_class V=block_core(N,192,bits(X),st,length);
            std::cout<<"{\"bits\":"<<st.Jmax<<",\"V_hex\":\""<<V.get_str(16)<<"\"}\n";
            return 0;
        }
        if(argc<2 || argc>3) throw std::invalid_argument("usage: extract_weighted POSITION [--timing]");
        if(argc==3 && std::string(argv[2])!="--timing") throw std::invalid_argument("unknown option");
        std::size_t used=0; std::string arg=argv[1];
        u64 N=std::stoull(arg,&used);
        if(arg.empty() || arg[0]=='-' || used!=arg.size() || N<1 || N>500000000ULL)
            throw std::invalid_argument("position outside 1..500000000");
        const u64 G=192;
        auto all=Clock::now(),start=all;
        mpz_class X; mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(N-1));
        const i64 K=cutoff(N,G),limit=bits(X)+G+8;
        double pre=elapsed(start); start=Clock::now();
        mpz_class acc=0,t;
        for(auto f:families) for(i64 k=0;k<=K;++k) {
            const i64 E=static_cast<i64>(N-1)+f.c-4*k;
            const u64 mu=f.a*k+f.b;
            if(mu==1) continue;
            u64 r=0;
            if(E>=0) r=mulmod(powmod(2,E,mu),powmod(5,N-1,mu),mu);
            else if(-E<limit) r=mulmod(powmod(5,N-1,mu),powmod((mu+1)>>1,-E,mu),mu);
            if(r) {
                t=integer(r)<<bitcount(G);
                t=floor_div(t,integer(mu));
                if(f.sign>0) acc+=t; else acc-=t;
            }
        }
        reduce(acc,G); double word=elapsed(start); start=Clock::now();
        Stats st; mpz_class V=block_core(N,G,bits(X),st),W=X*V;
        reduce(W,st.Jmax);
        if(st.Jmax>=G) W>>=bitcount(st.Jmax-G); else W<<=bitcount(G-st.Jmax);
        acc+=W; reduce(acc,G);
        double core=elapsed(start);
        mpz_class scale; mpz_ui_pow_ui(scale.get_mpz_t(),10,12);
        acc*=scale; mpz_fdiv_q_2exp(acc.get_mpz_t(),acc.get_mpz_t(),bitcount(G));
        double total=elapsed(all);
        gmp_printf("%012Zd\n",acc.get_mpz_t());
        if(argc==3) {
            print_stats(st,pre,word,core,total);
        }
        return 0;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"error: %s\n",e.what()); return 1;
    }
}
