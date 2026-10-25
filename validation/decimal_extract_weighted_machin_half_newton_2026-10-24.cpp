// Certified decimal windows from Machin's identity with weighted CRT blocks.
// Experimental implementation, 24 October 2026. The weighted prime-power cores
// use exact local CRT fractions and certified signed real correction, then
// bounded-memory rational forests. Successful output retains Machin's integer
// interval certificate; guard retries are adaptive. The proposed half-log core
// bound requires separate analysis and does not bound arbitrary guard retries.
// Build: g++ -O3 -march=native -std=c++17 decimal_extract_weighted_machin.cpp -lgmpxx -lgmp -o extract_weighted_machin.exe
#include <gmpxx.h>
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

static mpz_class prime_power(u64 p,u64 exponent) {
    mpz_class out;
    mpz_ui_pow_ui(out.get_mpz_t(),static_cast<unsigned long>(p),static_cast<unsigned long>(exponent));
    return out;
}
static u64 powmod(u64 x,u64 exponent,u64 modulus) {
    u64 result=1%modulus; x%=modulus;
    while(exponent) {
        if(exponent&1) result=mulmod(result,x,modulus);
        x=mulmod(x,x,modulus); exponent>>=1;
    }
    return result;
}
static PQT machin_real_bs(u64 A,u64 a,u64 C,u64 lo,u64 hi) {
    const u128 largest=static_cast<u128>(C)*(A+static_cast<u128>(a)*(hi-1));
    if((hi-lo<=2 && largest<(static_cast<u128>(1)<<62)) ||
       (hi-lo<=3 && largest<(static_cast<u128>(1)<<41)) ||
       (hi-lo<=4 && largest<(1ULL<<31))) {
        u128 p=1,q=1,t=0;
        for(u64 r=hi;r-->lo;) {
            const u128 pr=r ? static_cast<u128>(a)*r : 1;
            const u128 qr=r ? static_cast<u128>(C)*(A+static_cast<u128>(a)*r) : 1;
            t=pr*(q+t); p*=pr; q*=qr;
        }
        return {integer128(p),integer128(q),integer128(t)};
    }
    if(hi-lo==1) {
        if(!lo) return {mpz_class(1),mpz_class(1),mpz_class(1)};
        mpz_class p=integer(a*lo),q=integer(C*(A+a*lo));
        return {p,q,p};
    }
    const u64 mid=lo+(hi-lo)/2;
    auto left=machin_real_bs(A,a,C,lo,mid),right=machin_real_bs(A,a,C,mid,hi);
    return {left.p*right.p,left.q*right.q,left.t*right.q+left.p*right.t};
}
static mpz_class machin_real_approx(u64 prime,u64 local_power,u64 A,u64 a,u64 length,u64 guard,Stats& stats) {
    const u64 C=prime*prime+1;
    const long double logp=std::log2(static_cast<long double>(prime));
    const long double logC=std::log2(static_cast<long double>(C));
    const long double target=std::log2(static_cast<long double>(local_power))+2*length*logp+12;
    const u64 target_bits=bits(local_power)+2*length*bits(prime)+12;
    const u64 geometric_bits=bits(C)-1;
    u64 count=(target_bits+geometric_bits-1)/geometric_bits;
    const u64 ratio=C*A/(a*count);
    const u64 ratio_bits=ratio ? bits(ratio)-1 : 0;
    if(ratio_bits) count=std::min(count,(target_bits+ratio_bits-1)/ratio_bits);
    const long double alpha=static_cast<long double>(A)/a;
    const long double scale=std::log2(static_cast<long double>(local_power))+2*length*logp
        -std::log2(static_cast<long double>(C)*A)+std::log2(static_cast<long double>(C)/(C-1));
    auto error_log=[&](u64 r) {
        return scale+(std::lgamma(static_cast<long double>(r)+1)+std::lgamma(alpha+1)
            -std::lgamma(alpha+r+1))/std::log(2.0L)-r*logC;
    };
    if(error_log(count)<-8) {
        u64 lo=1,hi=count;
        while(lo<hi) { const u64 mid=lo+(hi-lo)/2; if(error_log(mid)<-8) hi=mid; else lo=mid+1; }
        count=lo;
    }
    const mpz_class scale_integer=integer(local_power)*prime_power(prime,2*length);
    for(;;) {
        auto root=machin_real_bs(A,a,C,0,count);
        // H_R <= (P/Q)*aR/[C(A+aR)], remaining positive terms <=H_R*C/(C-1).
        // After multiplication by local_power*|B|^length/(CA), prove error<1/64.
        const mpz_class error_num=(scale_integer*root.p*integer(a*count))<<6;
        const mpz_class error_den=root.q*integer(C*A)*integer(A+a*count)*integer(C-1);
        if(error_num<error_den) {
            stats.real_terms+=count; stats.real_product_bits+=bits(root.q);
            mpz_class numerator=scale_integer*root.t;
            if(!(length&1)) numerator=-numerator; // B^L S_B has sign (-1)^(L+1).
            numerator<<=bitcount(guard);
            return floor_div(numerator,integer(C*A)*root.q);
        }
        count*=2;
    }
}
struct WeightedMachinBlock { Rational value; u64 valuation; };

static WeightedMachinBlock weighted_machin_block(u64 prime,u64 A,u64 a,u64 length,
                               const std::vector<u64>& primes,Stats& stats) {
    if(!A || !(A&1) || a<2 || (a&(a-1)) || !length)
        throw std::invalid_argument("require odd A, power-of-two a>=2 and positive length");
    const u64 maximum=A+a*(length-1);
    if(maximum>UINT32_MAX || length>500000000ULL)
        throw std::invalid_argument("block exceeds the audited word range");
    u64 local_power=1,local_valuation=0;
    while(local_power<=maximum/prime) { local_power*=prime; ++local_valuation; }
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
        if(p==prime) continue;
        const u64 base=(q-(prime*prime)%q)%q;
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
        u64 weight=modulus.multiply(local_power,modulus.power(base,length-1-j)),c=0;
        const bool cached=count>=3;
        const u64 step=cached ? modulus.power(base,p) : 0;
        for(auto i=item.end;i-- >item.begin;) {
            const u64 next=factors[i].index,pv=factors[i].power;
            if(i!=item.end-1) {
                // Consecutive denominators divisible by p are p indices apart.
                if(cached) {
                    if(j-next!=p) throw std::runtime_error("nonconsecutive prime incidence");
                    weight=modulus.multiply(weight,step);
                } else weight=modulus.multiply(local_power,modulus.power(base,length-1-next));
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
    mpz_class approx=machin_real_approx(prime,local_power,A,a,length,guard,stats);
    mpz_class correction=approx-lower+(mpz_class(1)<<bitcount(guard-1));
    mpz_fdiv_q_2exp(correction.get_mpz_t(),correction.get_mpz_t(),bitcount(guard));
    stats.real+=elapsed(start); start=Clock::now();
    Rational result=rational_tree(std::move(parts));
    result.n+=correction*result.d;
    stats.rational+=elapsed(start);
    return {std::move(result),local_valuation};
}
#include <map>
#include <sstream>
static mpz_class from_u64(u64 x) { return integer(x); }
namespace machin {
static mpz_class power(u64 p,u64 e) {
    if(p>ULONG_MAX || e>ULONG_MAX) throw std::overflow_error("power range");
    mpz_class v; mpz_ui_pow_ui(v.get_mpz_t(),static_cast<unsigned long>(p),static_cast<unsigned long>(e)); return v;
}
static mpz_class modulo(mpz_class a,const mpz_class& m) {
    mpz_fdiv_r(a.get_mpz_t(),a.get_mpz_t(),m.get_mpz_t()); return a;
}
static mpz_class inverse(const mpz_class& a,const mpz_class& m) {
    mpz_class b;
    if(!mpz_invert(b.get_mpz_t(),a.get_mpz_t(),m.get_mpz_t())) throw std::runtime_error("nonunit inverse");
    return b;
}
static mpz_class lift_inverse(const mpz_class& a,u64 p,u64 exponent) {
    mpz_class modulus=from_u64(p);
    mpz_class b=inverse(modulo(a,modulus),modulus);
    for(u64 precision=1;precision<exponent;) {
        const u64 previous=precision;
        precision=precision>=exponent-precision ? exponent : 2*precision;
        const u64 added=precision-previous;
        const mpz_class step=added==previous ? modulus : power(p,added);
        const mpz_class next_modulus=modulus*step;
        mpz_class correction=1-modulo(a,next_modulus)*b;
        // a*b == 1 mod p^previous, so this signed division is exact.
        mpz_divexact(correction.get_mpz_t(),correction.get_mpz_t(),modulus.get_mpz_t());
        correction=modulo(std::move(correction),step);
        correction=modulo(mpz_class(correction*b),step);
        b+=modulus*correction;
        modulus=next_modulus;
    }
    return b;
}

struct Stats {
    ::Stats weighted;
    u64 block_length=0;
    double core5=0,core239=0;
    u64 terms=0,core_terms=0,passes=0,guard=0,word_errors=0,core_errors=0;
    double pre=0,word=0,core=0,total=0;
};

// b^n p^e modulo a p-coprime word, allowing a negative exponent e.
// Division by p uses the exact representative (z+t*m)/p.
static u64 joint(u64 b,u64 n,u64 p,i64 e,u64 m,const std::vector<u64>& inverse_mod_p) {
    if(m==1) return 0;
    u64 exponent=e<0 ? -e : e,mask=n|exponent,z=1%m;
    if(!mask) return z;
    const u64 inv=inverse_mod_p[m%p];
    for(int bit=63-__builtin_clzll(mask);bit>=0;--bit) {
        z=mulmod(z,z,m);
        if((n>>bit)&1) z=(b*z)%m;
        if((exponent>>bit)&1) {
            if(e>=0) z=(p*z)%m;
            else {
                u64 t=((p-z%p)*inv)%p;
                z=(z+t*m)/p;
            }
        }
    }
    return z;
}
static std::vector<u64> small_inverses(u64 p) {
    std::vector<u64> inv(p);
    for(u64 j=1;j<p;++j) inv[j]=powmod(j,p-2,p);
    return inv;
}
struct Fraction { mpz_class n,d; };
class GeometricTree {
    u64 p,A,h;
    const mpz_class& modulus;
    std::size_t modulus_bits;
    std::map<u64,mpz_class> shifts;
    std::vector<mpz_class> scales;
    void trim(mpz_class& a) const {
        if(mpz_sizeinbase(a.get_mpz_t(),2)>=modulus_bits) a=modulo(std::move(a),modulus);
    }
    const mpz_class& shift(u64 length) {
        auto found=shifts.find(length);
        if(found!=shifts.end()) return found->second;
        mpz_class value=power(p,2*length);
        if(length&1) value=-value;
        trim(value);
        return shifts.emplace(length,std::move(value)).first->second;
    }
public:
    GeometricTree(u64 prime,u64 endpoint,u64 valuation,const mpz_class& m):
        p(prime),A(endpoint),h(valuation),modulus(m),modulus_bits(mpz_sizeinbase(m.get_mpz_t(),2)) {
        scales.push_back(1);
        for(u64 j=1;j<=h;++j) scales.push_back(scales.back()*from_u64(p));
    }
    Fraction sum(u64 lo,u64 hi) {
        if(hi-lo==1) {
            u64 u=A-2*lo,v=0;
            while(u%p==0) { u/=p; ++v; }
            return {scales[h-v],from_u64(u)};
        }
        u64 mid=lo+(hi-lo)/2;
        auto left=sum(lo,mid),right=sum(mid,hi);
        Fraction out{left.n*right.d+shift(mid-lo)*left.d*right.n,left.d*right.d};
        trim(out.n); trim(out.d); return out;
    }
};
struct PrimeCore { mpz_class residue,modulus; u64 exponent,terms; };
static PrimeCore reference_prime_core(u64 p,u64 K,u64 cancelled,const mpz_class& Y,int sign) {
    u64 A=2*K+1,h=0,q=1;
    while(q<=A/p) { q*=p; ++h; }
    i64 signed_m=static_cast<i64>(A+h)-static_cast<i64>(cancelled);
    if(signed_m<=0) return {0,1,0,0};
    u64 m=signed_m,R=std::min(K+1,(m+1)/2);
    mpz_class modulus=power(p,m);
    GeometricTree tree(p,A,h,modulus); auto value=tree.sum(0,R);
    auto inv=lift_inverse(value.d,p,m);
    mpz_class W=modulo(mpz_class(value.n*inv),modulus);
    W=modulo(mpz_class(W*modulo(Y,modulus)),modulus);
    if(((K&1) ? -sign : sign)<0 && W!=0) W=modulus-W;
    return {std::move(W),std::move(modulus),m,R};
}
struct WeightedNode { mpz_class n,d; u64 shift; };
static PrimeCore prime_core(u64 p,u64 K,u64 cancelled,const mpz_class& Y,int sign,Stats* profile=nullptr) {
    const u64 A=2*K+1;
    u64 h=0,q=1;
    while(q<=A/p) { q*=p; ++h; }
    const i64 signed_m=static_cast<i64>(A+h)-static_cast<i64>(cancelled);
    if(signed_m<=0) return {0,1,0,0};
    const u64 m=signed_m,R=std::min(K+1,(m+1)/2);
    const u64 logarithm=std::max<u64>(1,bits(R));
    u64 root=isqrt(logarithm); if(root*root<logarithm) ++root;
    const u64 block_scale=1ULL<<root;
    const u64 length=profile && profile->block_length ? profile->block_length : (R+block_scale-1)/block_scale;
    const u64 floor_log_p=bits(p)-1;
    const u64 mass_limit=m*bits(p);
    ::Stats scratch;
    auto& weighted=profile ? profile->weighted : scratch;
    weighted.block_length=length;
    auto trim_pair=[&](mpz_class& n,mpz_class& d,u64 precision) {
        // This exact integer lower bound avoids forming a full-precision
        // modulus for small leaf values. No modulus cache survives this call.
        const u64 lower=precision*floor_log_p;
        if(bits(n)<=lower && bits(d)<=lower) return;
        const mpz_class modulus=power(p,precision);
        if(bits(n)>=bits(modulus)) n=modulo(std::move(n),modulus);
        if(bits(d)>=bits(modulus)) d=modulo(std::move(d),modulus);
    };
    auto merge_pair=[&](WeightedNode left,WeightedNode right) {
        const u64 shift=std::min(left.shift,right.shift);
        mpz_class ln=left.n*right.d,rn=right.n*left.d;
        if(left.shift>shift) ln*=power(p,left.shift-shift);
        if(right.shift>shift) rn*=power(p,right.shift-shift);
        mpz_class n=ln+rn,d=left.d*right.d;
        trim_pair(n,d,m-shift);
        return WeightedNode{std::move(n),std::move(d),shift};
    };
    auto merge_tree=[&](std::vector<WeightedNode> nodes) {
        while(nodes.size()>1) {
            std::vector<WeightedNode> next;
            next.reserve((nodes.size()+1)/2);
            for(std::size_t i=0;i+1<nodes.size();i+=2)
                next.push_back(merge_pair(std::move(nodes[i]),std::move(nodes[i+1])));
            if(nodes.size()&1) next.push_back(std::move(nodes.back()));
            nodes=std::move(next);
        }
        return std::move(nodes.front());
    };
    std::vector<WeightedNode> pending;
    WeightedNode accumulated{0,1,m};
    bool have_accumulated=false;
    u64 mass=0;
    auto flush=[&]() {
        if(pending.empty()) return;
        const auto start=Clock::now();
        auto value=merge_tree(std::move(pending));
        pending.clear();
        if(have_accumulated) accumulated=merge_pair(std::move(accumulated),std::move(value));
        else { accumulated=std::move(value); have_accumulated=true; }
        mass=0; ++weighted.forests;
        weighted.assembly+=elapsed(start);
    };
    auto append=[&](WeightedNode node,u64 span) {
        if(node.shift>=m) return;
        const auto start=Clock::now();
        trim_pair(node.n,node.d,m-node.shift);
        mass+=bits(node.n)+bits(node.d)+2*span*bits(p);
        weighted.maximum_pending_bits=std::max(weighted.maximum_pending_bits,mass);
        pending.push_back(std::move(node));
        weighted.assembly+=elapsed(start);
        if(mass>=mass_limit) flush();
    };
    const auto primes=primes_to(isqrt(A));
    // Descending r intervals yield increasing denominators and decreasing
    // p-adic shifts. Accumulated precision therefore grows only as needed.
    for(u64 hi=R;hi;) {
        const u64 low_mu=A-2*(hi-1);
        const u64 count=std::min(hi,length),lo=hi-count;
        auto block=weighted_machin_block(p,low_mu,2,count,primes,weighted);
        const u64 shift=2*lo+h-block.valuation;
        if(shift<m) {
            if(lo&1) block.value.n=-block.value.n;
            append({std::move(block.value.n),std::move(block.value.d),shift},count);
            ++weighted.blocks;
        }
        hi=lo;
    }
    flush();
    const auto start=Clock::now();
    const mpz_class modulus=power(p,m);
    mpz_class W=0;
    if(have_accumulated) {
        const auto& value=accumulated;
        const mpz_class local_modulus=value.shift ? power(p,m-value.shift) : modulus;
        W=modulo(mpz_class(value.n*lift_inverse(value.d,p,m-value.shift)),local_modulus);
        if(value.shift) W*=power(p,value.shift);
        W=modulo(mpz_class(W*modulo(Y,modulus)),modulus);
        if(((K&1) ? -sign : sign)<0 && W!=0) W=modulus-W;
    }
    weighted.assembly+=elapsed(start);
    return {std::move(W),modulus,m,R};
}

static u64 cutoff(u64 p,u64 coefficient,u64 n,u64 G,const mpz_class& ten) {
    // Floating arithmetic only proposes a candidate. Integer comparisons
    // establish the alternating-tail bound and determine the final cutoff.
    u64 K=static_cast<u64>((n*std::log2(10.0)+G+10)/(2*std::log2(double(p))));
    mpz_class threshold=(ten*from_u64(coefficient))<<(G+1);
    auto denominator=power(p,2*K+3);
    while(denominator<=threshold) { ++K; denominator*=from_u64(p*p); }
    while(K && denominator/from_u64(p*p)>threshold) { --K; denominator/=from_u64(p*p); }
    return K;
}
static void add_word(mpz_class& acc,u64 residue,u64 divisor,int sign,u64 G) {
    if(!residue) return;
    mpz_class t;
    if(residue<=ULONG_MAX && divisor<=ULONG_MAX) {
        mpz_set_ui(t.get_mpz_t(),static_cast<unsigned long>(residue));
        mpz_mul_2exp(t.get_mpz_t(),t.get_mpz_t(),G);
        mpz_fdiv_q_ui(t.get_mpz_t(),t.get_mpz_t(),static_cast<unsigned long>(divisor));
    } else {
        t=from_u64(residue)<<G;
        mpz_class d=from_u64(divisor);
        mpz_fdiv_q(t.get_mpz_t(),t.get_mpz_t(),d.get_mpz_t());
    }
    if(sign>0) acc+=t; else acc-=t;
}
static mpz_class attempt(u64 N,u64 G,Stats& stats) {
    u64 n=N-1; auto start=Clock::now();
    auto ten=power(10,n); mpz_class acc=0;
    stats.pre+=elapsed(start);
    for(u64 p:{5ULL,239ULL}) {
        start=Clock::now();
        u64 c=p==5 ? 16 : 4,s=p==5 ? n : 0;
        int sign=p==5 ? 1 : -1;
        u64 K=cutoff(p,c,n,G,ten);
        mpz_class Y=p==5 ? mpz_class(mpz_class(1)<<(n+4)) : mpz_class(4*ten);
        auto inverses=small_inverses(p); stats.pre+=elapsed(start);
        start=Clock::now();
        for(u64 k=0;k<=K;++k) {
            u64 mu=2*k+1,u=mu,v=0;
            while(u%p==0) { u/=p; ++v; }
            ++stats.terms;
            if(u==1) continue;
            i64 e=static_cast<i64>(s)-static_cast<i64>(2*k+1+v);
            u64 residue=joint(p==5 ? 2 : 10,n,p,e,u,inverses);
            residue=(c*residue)%u;
            add_word(acc,residue,u,(k&1) ? -sign : sign,G);
            ++stats.word_errors;
        }
        stats.word+=elapsed(start); start=Clock::now();
        auto core=prime_core(p,K,s,Y,sign,&stats); stats.core_terms+=core.terms;
        if(core.exponent) {
            mpz_class t=core.residue<<G;
            mpz_fdiv_q(t.get_mpz_t(),t.get_mpz_t(),core.modulus.get_mpz_t());
            acc+=t; ++stats.core_errors;
        }
        const double core_elapsed=elapsed(start);
        stats.core+=core_elapsed;
        if(p==5) stats.core5+=core_elapsed; else stats.core239+=core_elapsed;
    }
    reduce(acc,G); return acc;
}
static std::string extract(u64 N,u64 initial_guard,bool euler,Stats& stats) {
    auto start=Clock::now(); u64 G=initial_guard;
    for(;;) {
        u64 prior_words=stats.word_errors,prior_cores=stats.core_errors;
        ++stats.passes; stats.guard=G;
        auto A=attempt(N,G,stats);
        mpz_class scale=power(10,12),unit=mpz_class(1)<<G,digit,remainder;
        mpz_class numerator=A*scale;
        mpz_fdiv_qr(digit.get_mpz_t(),remainder.get_mpz_t(),numerator.get_mpz_t(),unit.get_mpz_t());
        mpz_class error=from_u64(stats.word_errors-prior_words+stats.core_errors-prior_cores+1);
        mpz_class margin=scale*error;
        if(remainder>margin && unit-remainder>margin) {
            stats.total=elapsed(start); auto result=digit.get_str();
            return std::string(12-result.size(),'0')+result;
        }
        if(G>100000000) throw std::overflow_error("certificate guard range");
        G*=2;
    }
}
static void selftest() {
    u64 residue_checks=0,core_checks=0;
    for(u64 p:{3ULL,5ULL,239ULL}) {
        auto inv=small_inverses(p);
        for(u64 m=2;m<200;++m) if(m%p) for(i64 e:{-101LL,-1LL,0LL,1LL,87LL}) {
            u64 n=57,b=10;
            auto factor=e>=0 ? powmod(p,e,m) : powmod(mpz_get_ui(inverse(from_u64(p),from_u64(m)).get_mpz_t()),-e,m);
            if(joint(b,n,p,e,m,inv)!=mulmod(powmod(b,n,m),factor,m)) throw std::runtime_error("joint residue mismatch");
            ++residue_checks;
        }
        for(u64 K:{0ULL,1ULL,4ULL,7ULL,31ULL,128ULL}) for(u64 s:{0ULL,K,2*K}) for(int sign:{-1,1}) {
            mpz_class Y=1234567; auto core=prime_core(p,K,s,Y,sign); mpz_class expected=0;
            for(u64 k=0;k<=K && core.exponent;++k) {
                u64 u=2*k+1,v=0; while(u%p==0) {u/=p; ++v;}
                i64 J=static_cast<i64>(2*k+1+v)-static_cast<i64>(s);
                if(J<=0) continue;
                mpz_class term=Y*power(p,core.exponent-J)*inverse(from_u64(u),core.modulus);
                if((k&1) ? sign>0 : sign<0) expected-=term; else expected+=term;
            }
            if(core.residue!=modulo(std::move(expected),core.modulus)) throw std::runtime_error("prime core mismatch");
            ++core_checks;
        }
    }
    std::cout<<"{\"selftest\":\"PASS\",\"residue_checks\":"<<residue_checks<<",\"prime_core_checks\":"<<core_checks<<"}\n";
}

} // namespace machin

static void weighted_selftest() {
    machin::selftest();
    std::mt19937_64 rng(20261024);
    unsigned blocks=0,cores=0;
    for(u64 p:{3ULL,5ULL,239ULL}) {
        for(unsigned i=0;i<60;++i) {
            const u64 A=i<8 ? 1 : 2*(rng()%10000)+1;
            const u64 length=1+rng()%80;
            Stats stats;
            auto actual=weighted_machin_block(p,A,2,length,primes_to(isqrt(A+2*(length-1))),stats);
            mpq_class expected=0;
            for(u64 j=0;j<length;++j) {
                mpz_class weight=prime_power(p,2*(length-1-j));
                if((length-1-j)&1) weight=-weight;
                expected+=mpq_class(weight,integer(A+2*j));
            }
            expected*=prime_power(p,actual.valuation);
            if(actual.value.n*expected.get_den()!=actual.value.d*expected.get_num())
                throw std::runtime_error("exact weighted Machin block mismatch");
            ++blocks;
        }
        for(u64 K:{0ULL,1ULL,7ULL,31ULL,128ULL,257ULL}) for(u64 s:{0ULL,K,2*K}) for(u64 length:{1ULL,7ULL,32ULL,128ULL}) {
            machin::Stats stats; stats.block_length=length;
            auto actual=machin::prime_core(p,K,s,integer(1234567),-1,&stats);
            auto expected=machin::reference_prime_core(p,K,s,integer(1234567),-1);
            if(actual.exponent!=expected.exponent || actual.modulus!=expected.modulus || actual.residue!=expected.residue)
                throw std::runtime_error("weighted Machin core differs from original geometric core");
            ++cores;
        }
    }
    std::cout<<"{\"weighted_selftest\":\"PASS\",\"exact_blocks\":"<<blocks<<",\"reference_core_cases\":"<<cores<<"}\n";
}

static void weighted_core_check(u64 K,u64 length) {
    unsigned checks=0;
    for(u64 p:{3ULL,5ULL,239ULL}) for(u64 cancelled:{0ULL,K,2*K}) for(int sign:{-1,1}) {
        machin::Stats stats; stats.block_length=length;
        auto actual=machin::prime_core(p,K,cancelled,integer(1234567),sign,&stats);
        auto expected=machin::reference_prime_core(p,K,cancelled,integer(1234567),sign);
        if(actual.residue!=expected.residue || actual.modulus!=expected.modulus || actual.exponent!=expected.exponent)
            throw std::runtime_error("weighted-vs-geometric exact core mismatch");
        ++checks;
    }
    std::cout<<"{\"core_check\":\"PASS\",\"K\":"<<K<<",\"cases\":"<<checks<<"}\n";
}

static u64 parse_cli_number(const std::string& text,u64 minimum,u64 maximum,const char* message) {
    if(text.empty()) throw std::invalid_argument(message);
    for(char digit:text) if(digit<'0' || digit>'9') throw std::invalid_argument(message);
    const u64 value=std::stoull(text);
    if(value<minimum || value>maximum) throw std::invalid_argument(message);
    return value;
}
int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--selftest") { weighted_selftest(); return 0; }
        if((argc==3 || argc==4) && std::string(argv[1])=="--core-check") {
            const u64 K=parse_cli_number(argv[2],0,500000000,"core K outside 0..500000000");
            const u64 block=argc==4 ? parse_cli_number(argv[3],1,500000000,"block outside 1..500000000") : 0;
            weighted_core_check(K,block); return 0;
        }
        if(argc<2) throw std::invalid_argument("usage: weighted_machin POSITION [--timing] [--guard=N] [--block=N]");
        const u64 N=parse_cli_number(argv[1],1,500000000,"position outside 1..500000000");
        u64 guard=192,requested_block=0;
        bool timing=false;
        for(int j=2;j<argc;++j) {
            const std::string option=argv[j];
            if(option=="--timing") timing=true;
            else if(option.rfind("--block=",0)==0)
                requested_block=parse_cli_number(option.substr(8),1,500000000,"block outside 1..500000000");
            else if(option.rfind("--guard=",0)==0)
                guard=parse_cli_number(option.substr(8),1,100000000,"guard outside 1..100000000");
            else throw std::invalid_argument("unknown option");
        }
        machin::Stats stats; stats.block_length=requested_block;
        const auto result=machin::extract(N,guard,false,stats);
        std::cout<<result<<'\n';
        if(timing) {
            std::ostringstream out; out.precision(10);
            out<<"{\"engine\":\"weighted-machin\",\"N\":"<<N<<",\"certified\":true";
#define FIELD(name) out<<",\"" #name "\":"<<stats.name
            FIELD(terms); FIELD(core_terms); FIELD(passes); FIELD(guard); FIELD(word_errors); FIELD(core_errors);
            FIELD(pre); FIELD(word); FIELD(core); FIELD(core5); FIELD(core239); FIELD(total); FIELD(block_length);
#undef FIELD
            out<<",\"factor\":"<<stats.weighted.factor<<",\"coefficients\":"<<stats.weighted.coefficients
               <<",\"real\":"<<stats.weighted.real<<",\"rational\":"<<stats.weighted.rational
               <<",\"assembly\":"<<stats.weighted.assembly<<",\"blocks\":"<<stats.weighted.blocks
               <<",\"real_terms\":"<<stats.weighted.real_terms<<",\"forests\":"<<stats.weighted.forests
               <<",\"maximum_pending_bits\":"<<stats.weighted.maximum_pending_bits<<"}";
            std::cerr<<out.str()<<'\n';
        }
        return 0;
    } catch(const std::exception& error) { std::cerr<<"error: "<<error.what()<<'\n'; return 1; }
}
