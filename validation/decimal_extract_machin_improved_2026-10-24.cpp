// BENCHMARK CONTROL ONLY: common arithmetic optimizations; root source preserved.
#ifndef COMMON_EULER_DEFAULT
#define COMMON_EULER_DEFAULT 0
#endif
// Certified decimal windows from Machin's identity and prime-power cores.
// With G=O(log N), this prototype has an O(N log^3 N) bound, not a new exponent.
// Unrestricted guard retries use a wider word accumulator and need separate
// worst-case accounting; the native program does not inherit the Lucas bound.
// Build: g++ -O3 -march=native -std=c++17 decimal_extract_machin.cpp -lgmpxx -lgmp -o extract_machin.exe
#define main shared_reference_main
#include "weighted_controls_common/decimal_extract_shared_residues.cpp"
#undef main
#include <map>
#include <sstream>
#include <iostream>

namespace machin {
static mpz_class power(u64 p,u64 e) {
    if(p>ULONG_MAX || e>ULONG_MAX) throw std::overflow_error("power range");
    mpz_class v; mpz_ui_pow_ui(v.get_mpz_t(),static_cast<unsigned long>(p),static_cast<unsigned long>(e)); return v;
}
static mpz_class modulo(mpz_class a,const mpz_class& m) {
    mpz_fdiv_r(a.get_mpz_t(),a.get_mpz_t(),m.get_mpz_t()); return a;
}
// Remove binary zero limbs before the final large multiplication.
// This is exact for arbitrary signed Y, including Y == 0 modulo modulus.
static mpz_class scale_residue(mpz_class value,const mpz_class& Y,const mpz_class& modulus) {
    mpz_class odd=modulo(Y,modulus);
    if(odd==0) return mpz_class(0);
    const mp_bitcnt_t shift=mpz_scan1(odd.get_mpz_t(),0);
    if(shift) mpz_fdiv_q_2exp(odd.get_mpz_t(),odd.get_mpz_t(),shift);
    if(odd!=1) value*=odd;
    if(shift) mpz_mul_2exp(value.get_mpz_t(),value.get_mpz_t(),shift);
    return modulo(std::move(value),modulus);
}
static mpz_class inverse(const mpz_class& a,const mpz_class& m) {
    mpz_class b;
    if(!mpz_invert(b.get_mpz_t(),a.get_mpz_t(),m.get_mpz_t())) throw std::runtime_error("nonunit inverse");
    return b;
}
// Balanced Newton lifting with nested operand residues and reused moduli.
// Moduli and residues have geometrically decreasing sizes: O(N) live bits.
static mpz_class lift_inverse(const mpz_class& a,u64 p,u64 exponent) {
    if(p<2 || p>ULONG_MAX || !exponent || exponent>ULONG_MAX)
        throw std::invalid_argument("prime-power inverse range");
    std::vector<u64> precision;
    for(u64 j=exponent;;j=j/2+j%2) {
        precision.push_back(j);
        if(j==1) break;
    }
    std::reverse(precision.begin(),precision.end());
    std::vector<mpz_class> moduli;
    moduli.reserve(precision.size());
    moduli.push_back(from_u64(p));
    for(std::size_t i=1;i<precision.size();++i) {
        mpz_class next=moduli.back()*moduli.back();
        // j is either twice the preceding precision or one less.
        if(precision[i]&1)
            mpz_divexact_ui(next.get_mpz_t(),next.get_mpz_t(),static_cast<unsigned long>(p));
        moduli.push_back(std::move(next));
    }
    std::vector<mpz_class> residues(moduli.size());
    residues.back()=modulo(a,moduli.back());
    for(std::size_t i=residues.size()-1;i>0;--i)
        residues[i-1]=modulo(residues[i],moduli[i-1]);
    mpz_class b=inverse(residues.front(),moduli.front());
    for(std::size_t i=1;i<moduli.size();++i)
        b=modulo(mpz_class(b*(2-residues[i]*b)),moduli[i]);
    return b;
}

struct Stats {
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
static PrimeCore prime_core(u64 p,u64 K,u64 cancelled,const mpz_class& Y,int sign) {
    u64 A=2*K+1,h=0,q=1;
    while(q<=A/p) { q*=p; ++h; }
    i64 signed_m=static_cast<i64>(A+h)-static_cast<i64>(cancelled);
    if(signed_m<=0) return {0,1,0,0};
    u64 m=signed_m,R=std::min(K+1,(m+1)/2);
    mpz_class modulus=power(p,m);
    GeometricTree tree(p,A,h,modulus); auto value=tree.sum(0,R);
    auto inv=lift_inverse(value.d,p,m);
    mpz_class W=modulo(mpz_class(value.n*inv),modulus);
    W=scale_residue(std::move(W),Y,modulus);
    if(((K&1) ? -sign : sign)<0 && W!=0) W=modulus-W;
    return {std::move(W),std::move(modulus),m,R};
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
        auto core=prime_core(p,K,s,Y,sign); stats.core_terms+=core.terms;
        if(core.exponent) {
            mpz_class t=core.residue<<G;
            mpz_fdiv_q(t.get_mpz_t(),t.get_mpz_t(),core.modulus.get_mpz_t());
            acc+=t; ++stats.core_errors;
        }
        stats.core+=elapsed(start);
    }
    reduce(acc,G); return acc;
}
static mpz_class euler_attempt(u64 N,u64 G,Stats& stats) {
    u64 n=N-1; auto start=Clock::now(); auto X=power(5,n);
    // 4(K+1) >= n+bits(X)+G+8, so the BBP absolute tail is <2^-G.
    i64 K=(n+mpz_sizeinbase(X.get_mpz_t(),2)+G+11)/4-1;
    i64 Jcut=4*K-static_cast<i64>(n)+2;
    std::vector<Family> families{{1,8,1,static_cast<i64>(n)+2},{-1,2,1,static_cast<i64>(n)-1},
        {-1,8,5,static_cast<i64>(n)},{-1,4,3,static_cast<i64>(n)-1}};
    stats.pre+=elapsed(start); start=Clock::now();
    auto acc=shared_words(n,K,Jcut,static_cast<unsigned>(G));
    stats.word_errors+=4*(K+1); stats.terms+=4*(K+1);
    stats.word+=elapsed(start); start=Clock::now();
    auto core=build_core(families,K,Jcut);
    mpz_class W=X*core.V; reduce(W,core.bits);
    if(core.bits>=G) W>>=core.bits-G; else W<<=G-core.bits;
    acc+=W; reduce(acc,G);
    ++stats.core_errors; stats.core_terms+=core.new_terms; stats.core+=elapsed(start);
    return acc;
}
static std::string extract(u64 N,u64 initial_guard,bool euler,Stats& stats) {
    auto start=Clock::now(); u64 G=initial_guard;
    for(;;) {
        u64 prior_words=stats.word_errors,prior_cores=stats.core_errors;
        ++stats.passes; stats.guard=G;
        auto A=euler ? euler_attempt(N,G,stats) : attempt(N,G,stats);
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
int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--selftest") {machin::selftest();return 0;}
        if(argc<2) throw std::invalid_argument("usage: extract_machin POSITION [--timing] [--euler] [--guard=N]");
        std::string arg=argv[1]; std::size_t used=0; u64 N=std::stoull(arg,&used),guard=192;
        if(arg.empty() || arg[0]=='-' || used!=arg.size() || N<1 || N>500000000) throw std::invalid_argument("position range");
        bool euler=COMMON_EULER_DEFAULT,timing=false;
        for(int j=2;j<argc;++j) {
            std::string option=argv[j];
            if(option=="--euler") euler=true; else if(option=="--timing") timing=true;
            else if(option.rfind("--guard=",0)==0) {
                auto value=option.substr(8); guard=std::stoull(value,&used);
                if(value.empty() || value[0]=='-' || used!=value.size() || !guard || guard>100000000) throw std::invalid_argument("guard range");
            } else throw std::invalid_argument("unknown option");
        }
        machin::Stats stats; auto result=machin::extract(N,guard,euler,stats); std::cout<<result<<'\n';
        if(timing) {
            std::ostringstream out; out.precision(10);
            out<<"{\"engine\":\""<<(euler ? "euler-certified-interval" : "machin-certified-prime-cores")<<"\",\"N\":"<<N<<",\"certified\":true";
#define FIELD(name) out<<",\"" #name "\":"<<stats.name
            FIELD(terms); FIELD(core_terms); FIELD(passes); FIELD(guard); FIELD(word_errors); FIELD(core_errors);
            FIELD(pre); FIELD(word); FIELD(core); FIELD(total);
#undef FIELD
            out<<"}"; std::cerr<<out.str()<<'\n';
        }
        return 0;
    } catch(const std::exception& e) {std::cerr<<"error: "<<e.what()<<'\n';return 1;}
}
