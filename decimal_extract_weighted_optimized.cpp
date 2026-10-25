// Exact weighted-block CRT core with real Euler correction.
// Native implementation of validation/weighted_block_crt.py.
// Optimized 24 October 2026; original source and measurements are preserved.
// Successful decimal outputs carry an integer interval certificate. A failed
// fixed-guard certificate returns an error; this is not the paper's complete
// adaptive exact-word procedure or a new end-to-end complexity proof.
// Build: g++ -O3 -march=native -std=c++17 decimal_extract_weighted_optimized.cpp -lgmpxx -lgmp -o extract_weighted_optimized
#include <gmpxx.h>
#include "decimal_extract_weighted_words.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
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
static mpz_class integer128(u128 x) {
    mpz_class z;
    mpz_import(z.get_mpz_t(),1,1,sizeof(x),0,0,&x);
    return z;
}
static mpz_class signed_integer128(__int128 x) {
    const bool negative=x<0;
    u128 magnitude=negative ? static_cast<u128>(-x) : static_cast<u128>(x);
    mpz_class z;
    mpz_import(z.get_mpz_t(),1,1,sizeof(magnitude),0,0,&magnitude);
    if(negative) mpz_neg(z.get_mpz_t(),z.get_mpz_t());
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
// Every modulus in the supported position range is <= UINT32_MAX, and
// operands are <= UINT32_MAX. Their product fits exactly in uint64_t.
static u64 mulmod(u64 x,u64 y,u64 m) { return x*y%m; }
static u64 invmod(u64 x,u64 m) {
    // Inputs are <= UINT32_MAX. Extended-Euclid coefficients (including
    // the terminal coefficient) have magnitude <= m, so int64_t suffices.
    i64 t=0,nt=1;
    std::uint32_t r=m,nr=x%m;
    while(nr) {
        std::uint32_t q=r/nr;
        std::uint32_t next=r-q*nr; r=nr; nr=next;
        i64 next_t=t-static_cast<i64>(q)*nt;
        t=nt; nt=next_t;
    }
    if(r!=1) throw std::runtime_error("noninvertible local denominator");
    t%=static_cast<i64>(m); if(t<0) t+=static_cast<i64>(m);
    return static_cast<u64>(t);
}
static mpz_class inverse_2exp(const mpz_class& d,u64 b) {
    if(!mpz_odd_p(d.get_mpz_t())) throw std::runtime_error("even denominator");
    mpz_class inv=1;
    for(u64 precision=1;precision<b;) {
        // Saturating doubling avoids the 32-bit mp_bitcnt_t overflow.
        const u64 previous=precision;
        precision=precision>=b-precision ? b : 2*precision;
        mpz_class low;
        mpz_fdiv_r_2exp(low.get_mpz_t(),d.get_mpz_t(),bitcount(precision));
        // 1-d*inv is divisible by 2^previous. Only the next bits of the
        // Newton correction are needed, so its second multiply is half-width.
        mpz_class correction=mpz_class(1)-low*inv;
        mpz_fdiv_q_2exp(correction.get_mpz_t(),correction.get_mpz_t(),bitcount(previous));
        reduce(correction,precision-previous);
        correction*=inv; reduce(correction,precision-previous);
        correction<<=bitcount(previous);
        inv+=correction;
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
static void sort_residual_primes(std::vector<u64>& residual,u64 maximum) {
    if(residual.size()<128) { std::sort(residual.begin(),residual.end()); return; }
    std::vector<u64> scratch(residual.size());
    // Residual records enter in index order. Stable radix passes on just
    // their high 32-bit prime retain increasing indices within each group.
    for(u64 shift=32;shift<32+bits(maximum);shift+=8) {
        std::size_t offsets[256]={};
        for(auto record:residual) ++offsets[(record>>shift)&255];
        std::size_t total=0;
        for(auto& offset:offsets) { const auto count=offset; offset=total; total+=count; }
        for(auto record:residual) scratch[offsets[(record>>shift)&255]++]=record;
        residual.swap(scratch);
    }
}
struct Rational { mpz_class n,d; };
struct PQT { mpz_class p,q,t; };
struct Stats {
    u64 N=0,Jmax=0,terms=0,block_length=0,blocks=0,real_terms=0;
    u64 lcm_leaf_bits=0,real_product_bits=0,forests=0,maximum_pending_bits=0;
    double factor=0,coefficients=0,real=0,rational=0,assembly=0;
};
static PQT real_bs(u64 A,u64 a,u64 lo,u64 hi) {
    // A few ratios fit in 128 bits: perform their recurrence before creating
    // GMP objects. q<2^31 bounds four-factor Q by 2^124, and |T|<2Q.
    const u128 largest_q=15*(static_cast<u128>(A)+static_cast<u128>(a)*(hi-1));
    if((hi-lo<=2 && largest_q<(static_cast<u128>(1)<<62)) ||
       (hi-lo<=4 && largest_q<(1ULL<<31))) {
        __int128 p=1,q=1,t=0;
        for(u64 r=hi;r-->lo;) {
            const __int128 pr=r ? -static_cast<__int128>(a*r) : 1;
            const __int128 qr=r ? static_cast<__int128>(15*(A+a*r)) : 1;
            t=pr*(q+t); p*=pr; q*=qr;
        }
        return {signed_integer128(p),signed_integer128(q),signed_integer128(t)};
    }
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
    // Select a tighter starting count from |h_R| <= R! (a/(15A))^R.
    // Floating arithmetic is only a scheduling hint: the exact integer
    // first-omitted-term certificate below remains the acceptance test.
    const long double ratio_log=std::log2(static_cast<long double>(a)/(15*A));
    const long double scale_log=4*static_cast<long double>(length)+6
        -std::log2(static_cast<long double>(15*A));
    auto error_log=[&](u64 r) {
        return scale_log+std::lgamma(static_cast<long double>(r)+1)/std::log(2.0L)
            +static_cast<long double>(r)*ratio_log;
    };
    if(error_log(count)<-8) {
        u64 low=1,high=count;
        while(low<high) {
            const u64 mid=low+(high-low)/2;
            if(error_log(mid)<-8) high=mid; else low=mid+1;
        }
        count=low;
    }
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
        const std::size_t count=nodes.size(); std::size_t out=0;
        for(std::size_t i=0;i+1<count;i+=2) {
            auto& l=nodes[i]; auto& r=nodes[i+1];
            Rational value{l.n*r.d+r.n*l.d,l.d*r.d};
            nodes[out++]=std::move(value);
        }
        if(count&1) nodes[out++]=std::move(nodes.back());
        nodes.resize(out);
    }
    return std::move(nodes.front());
}
struct WordModulus {
    u64 modulus,reciprocal;
    explicit WordModulus(u64 m):modulus(m),reciprocal(UINT64_MAX/m) {}
    u64 multiply(u64 x,u64 y) const {
        // m is odd and >=3, so reciprocal=floor(2^64/m). The high
        // product underestimates floor(x*y/m) by at most one.
        const u64 product=x*y;
        const u64 quotient=static_cast<u64>((static_cast<u128>(product)*reciprocal)>>64);
        u64 remainder=product-quotient*modulus;
        if(remainder>=modulus) remainder-=modulus;
        return remainder;
    }
    u64 power(u64 x,u64 e) const {
        u64 r=1; x%=modulus;
        while(e) { if(e&1) r=multiply(r,x); x=multiply(x,x); e>>=1; }
        return r;
    }
};
static Rational weighted_block(u64 A,u64 a,u64 length,
                               const std::vector<u64>& primes,Stats& stats) {
    if(!A || !(A&1) || a<2 || (a&(a-1)) || !length)
        throw std::invalid_argument("require odd A, power-of-two a>=2 and positive length");
    const u64 maximum=A+a*(length-1);
    if(maximum>UINT32_MAX || length>500000000ULL)
        throw std::invalid_argument("block exceeds the audited word range");
    auto start=Clock::now();
    std::vector<std::uint32_t> remainder(length);
    struct Incidence { std::uint32_t index,power; };
    std::vector<Incidence> factors;
    factors.reserve(3*length);
    struct Prime {
        std::uint32_t p,q;
        std::size_t begin,end;
    };
    std::vector<Prime> powers;
    for(u64 j=0;j<length;++j) remainder[j]=A+a*j;
    for(u64 p:primes) {
        if(p>maximum/p) break;
        u64 first=mulmod((p-A%p)%p,invmod(a,p),p);
        u64 q=1;
        const auto begin=factors.size();
        for(u64 j=first;j<length;j+=p) {
            u64 pv=1;
            while(remainder[j]%p==0) { remainder[j]/=p; pv*=p; }
            if(pv>1) {
                factors.push_back({static_cast<std::uint32_t>(j),static_cast<std::uint32_t>(pv)});
                q=std::max(q,pv);
            }
        }
        if(q>1) powers.push_back({static_cast<std::uint32_t>(p),static_cast<std::uint32_t>(q),begin,factors.size()});
    }
    // The residual factors exceed sqrt(maximum), hence are prime. Sorting
    // (prime,index) groups repeated large primes and preserves increasing
    // index within each group. Small primes above are already ordered.
    std::vector<u64> residual;
    residual.reserve(length);
    for(u64 j=0;j<length;++j) if(remainder[j]>1) {
        residual.push_back((static_cast<u64>(remainder[j])<<32)|j);
    }
    sort_residual_primes(residual,maximum);
    for(std::size_t i=0;i<residual.size();) {
        const auto p=static_cast<std::uint32_t>(residual[i]>>32);
        const auto begin=factors.size();
        do {
            factors.push_back({static_cast<std::uint32_t>(residual[i]),p});
            ++i;
        } while(i<residual.size() && residual[i]>>32==p);
        powers.push_back({p,p,begin,factors.size()});
    }
    stats.factor+=elapsed(start); start=Clock::now();
    const u64 guard=bits(maximum+length+1)+8;
    if(guard>=64) throw std::overflow_error("local guard exceeds word accumulator");
    // There are fewer than 2^32 distinct primes, and guard<=41 in the
    // audited range. The sum of the per-fraction floors fits in 128 bits.
    u128 lower_word=0;
    std::vector<Rational> parts;
    // If k*b+ceil(log2(k)) <=128 for b=bits(maximum), the product of
    // k denominators is <2^(k*b), and the numerator is less than k times
    // that product. This packs 5 leaves at typical million-digit positions,
    // while reducing to 3 at the full audited 32-bit denominator bound.
    unsigned packed_leaves=1;
    const auto denominator_bits=bits(maximum);
    while((packed_leaves+1)*denominator_bits+bits(packed_leaves)<=128) ++packed_leaves;
    parts.reserve((powers.size()+packed_leaves-1)/packed_leaves);
    u128 pending_c=0,pending_q=1;
    unsigned pending_count=0;
    std::vector<std::uint32_t> prefix;
    for(const auto& item:powers) {
        const u64 p=item.p,q=item.q;
        const WordModulus modulus(q);
        const auto count=item.end-item.begin;
        // All mu/p^v are units modulo q=p^e. Prefix products recover all
        // inverses using one extended-Euclid call per sufficiently large group.
        const bool batch=count>=4;
        u64 inverse=0;
        if(batch) {
            prefix.resize(count);
            u64 product=1;
            for(auto i=item.begin;i<item.end;++i) {
                prefix[i-item.begin]=product;
                const auto incidence=factors[i];
                product=modulus.multiply(product,(A+a*incidence.index)/incidence.power);
            }
            inverse=invmod(product,q);
        }
        u64 j=factors[item.end-1].index;
        u64 weight=modulus.power(16,length-1-j),c=0;
        const bool cached=count>=3;
        const u64 step=cached ? modulus.power(16,p) : 0;
        for(auto i=item.end;i-- >item.begin;) {
            const u64 next=factors[i].index,pv=factors[i].power;
            if(i!=item.end-1) {
                // Consecutive denominators divisible by p are p indices apart.
                if(cached) {
                    if(j-next!=p) throw std::runtime_error("nonconsecutive prime incidence");
                    weight=modulus.multiply(weight,step);
                } else weight=modulus.power(16,length-1-next);
            }
            j=next;
            const u64 mu=A+a*j;
            const u64 inverse_den=batch ? modulus.multiply(inverse,prefix[i-item.begin]) : invmod(mu/pv,q);
            if(batch && i!=item.begin) inverse=modulus.multiply(inverse,mu/pv);
            u64 term=q==pv ? weight : modulus.multiply(q/pv,weight);
            term=modulus.multiply(term,inverse_den);
            c+=term; if(c>=q) c-=q;
        }
        stats.lcm_leaf_bits+=bits(q);
        if(c) {
            if(c<=(UINT64_MAX>>guard)) lower_word+=(c<<guard)/q;
            else lower_word+=(static_cast<u128>(c)<<guard)/q;
            // The block-specific packing bound above keeps both products
            // and the sum exact in 128-bit arithmetic.
            pending_c=pending_c*q+static_cast<u128>(c)*pending_q;
            pending_q*=q;
            if(++pending_count==packed_leaves) {
                parts.push_back({integer128(pending_c),integer128(pending_q)});
                pending_c=0; pending_q=1; pending_count=0;
            }
        }
    }
    if(pending_count) parts.push_back({integer128(pending_c),integer128(pending_q)});
    const mpz_class lower=integer128(lower_word);
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
static Node merge_rational_nodes(std::vector<Node> nodes,u64 b) {
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
    return std::move(nodes.front());
}
static mpz_class node_residue(const Node& root,u64 b) {
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
    // Keep forest roots as odd-denominator rationals. A root truncated at
    // Jmax-shift is a valid lift after multiplication by 2^shift. Merging
    // into one bounded accumulator preserves that congruence and needs
    // just one final inverse, without retaining all forest roots in memory.
    Node accumulated{0,1,Jmax};
    bool have_accumulated=false;
    auto flush=[&]() {
        if(pending.empty()) return;
        auto t=Clock::now();
        Node root=merge_rational_nodes(std::move(pending),Jmax);
        if(!have_accumulated) { accumulated=std::move(root); have_accumulated=true; }
        else {
            std::vector<Node> pair; pair.push_back(std::move(accumulated)); pair.push_back(std::move(root));
            accumulated=merge_rational_nodes(std::move(pair),Jmax);
        }
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
    auto t=Clock::now();
    mpz_class V=have_accumulated ? node_residue(accumulated,Jmax) : mpz_class(0);
    stats.assembly+=elapsed(t);
    return V;
}
static void print_stats(const Stats& s,double pre,double word,double core,double total,u64 guard,u64 word_errors) {
    std::fprintf(stderr,
        "{\"engine\":\"weighted-block-optimized\",\"certified\":true,\"N\":%llu,\"Jmax\":%llu,\"terms\":%llu,"
        "\"block_length\":%llu,\"blocks\":%llu,\"forests\":%llu,\"real_terms\":%llu,"
        "\"lcm_leaf_bits\":%llu,\"real_product_bits\":%llu,\"maximum_pending_bits\":%llu,"
        "\"pre\":%.9f,\"word\":%.9f,\"core\":%.9f,\"total\":%.9f,"
        "\"factor\":%.9f,\"coefficients\":%.9f,\"real\":%.9f,\"rational\":%.9f,\"assembly\":%.9f,\"guard\":%llu,\"word_errors\":%llu,\"passes\":1}\n",
        (unsigned long long)s.N,(unsigned long long)s.Jmax,(unsigned long long)s.terms,
        (unsigned long long)s.block_length,(unsigned long long)s.blocks,
        (unsigned long long)s.forests,(unsigned long long)s.real_terms,
        (unsigned long long)s.lcm_leaf_bits,(unsigned long long)s.real_product_bits,
        (unsigned long long)s.maximum_pending_bits,pre,word,core,total,
        s.factor,s.coefficients,s.real,s.rational,s.assembly,
        (unsigned long long)guard,(unsigned long long)word_errors);
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
static u64 parse_number(const std::string& arg,u64 minimum,u64 maximum,const char* name) {
    std::size_t used=0;
    if(arg.empty() || arg[0]=='-' || arg[0]=='+') throw std::invalid_argument(name);
    u64 value=std::stoull(arg,&used);
    if(used!=arg.size() || value<minimum || value>maximum) throw std::invalid_argument(name);
    return value;
}
int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--selftest") { selftest(); return 0; }
        if((argc==3 || argc==4) && std::string(argv[1])=="--core") {
            const u64 N=parse_number(argv[2],1,500000000,"position outside 1..500000000");
            const u64 length=argc==4 ? parse_number(argv[3],0,500000000,"block range") : 0;
            mpz_class X; mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(N-1));
            Stats st; mpz_class V=block_core(N,192,bits(X),st,length);
            std::cout<<"{\"bits\":"<<st.Jmax<<",\"V_hex\":\""<<V.get_str(16)<<"\"}\n";
            return 0;
        }
        if(argc<2) throw std::invalid_argument("usage: extract_weighted_optimized POSITION [--timing] [--guard=N] [--block=N]");
        const u64 N=parse_number(argv[1],1,500000000,"position outside 1..500000000");
        u64 G=192,requested_block=0; bool timing=false;
        for(int i=2;i<argc;++i) {
            std::string option=argv[i];
            if(option=="--timing") timing=true;
            else if(option.rfind("--block=",0)==0) requested_block=parse_number(option.substr(8),1,500000000,"block range");
            else if(option.rfind("--guard=",0)==0) G=parse_number(option.substr(8),1,1000000,"guard outside 1..1000000");
            else throw std::invalid_argument("unknown option");
        }
        auto all=Clock::now(),start=all;
        mpz_class X; mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(N-1));
        const i64 K=cutoff(N,G),limit=bits(X)+G+8;
        // Absolute BBP bracket < 8 and the 1/16 geometric tail give
        // |10^(N-1) tail| < (128/15)*10^(N-1)*16^(-K-1).
        // Floating cutoff is only a proposal; this integer test proves <2^-G.
        if((X<<bitcount(N-1+G+7)) >= (mpz_class(15)<<bitcount(4*(K+1))))
            throw std::runtime_error("BBP tail cutoff certificate failed");
        double pre=elapsed(start); start=Clock::now();
        u64 word_errors=0;
        mpz_class acc=weighted_words::accumulate(N,K,limit,G,&word_errors);
        double word=elapsed(start); start=Clock::now();
        Stats st; mpz_class V=block_core(N,G,bits(X),st,requested_block),W=X*V;
        reduce(W,st.Jmax);
        if(st.Jmax>=G) W>>=bitcount(st.Jmax-G); else W<<=bitcount(G-st.Jmax);
        acc+=W; reduce(acc,G);
        double core=elapsed(start);
        mpz_class scale; mpz_ui_pow_ui(scale.get_mpz_t(),10,12);
        const mpz_class unit=mpz_class(1)<<bitcount(G),numerator=acc*scale;
        mpz_class digit,remainder;
        mpz_fdiv_qr(digit.get_mpz_t(),remainder.get_mpz_t(),numerator.get_mpz_t(),unit.get_mpz_t());
        // <1 ulp per word floor, <1 for core projection, <1 for the BBP
        // tail, and <1/60 for omitted J>=bits(X)+G+8 terms. The symmetric
        // bound also handles signed floors and wraparound modulo one.
        const mpz_class margin=scale*integer(word_errors+3);
        if(remainder<=margin || unit-remainder<=margin)
            throw std::runtime_error("decimal certificate failed; increase --guard");
        double total=elapsed(all);
        gmp_printf("%012Zd\n",digit.get_mpz_t());
        if(timing) print_stats(st,pre,word,core,total,G,word_errors);
        return 0;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"error: %s\n",e.what()); return 1;
    }
}
