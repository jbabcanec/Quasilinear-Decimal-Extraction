#pragma once
// Experimental shared real-correction reconstruction. Included after the
// weighted Machin integer, rational-tree, and exact real-PQT primitives.
#include <map>
#include <sstream>
namespace machin_lucas {
struct Block { u64 A,a; };
struct Query { u64 k,p,mod,value=0; };
struct Counters {
    u64 levels=0,base_sums=0,prime_queries=0,power_queries=0;
    u64 ordinary_batches=0,power_batches=0,nonunit_divisions=0,incidences=0;
    u64 max_root_bits=0,root_terms=0,batch_blocks=0;
    double collect=0,trees=0,lucas=0,total=0;
};
static mpz_class modz(mpz_class x,const mpz_class& m) {
    mpz_fdiv_r(x.get_mpz_t(),x.get_mpz_t(),m.get_mpz_t()); return x;
}
static mpz_class gcdz(const mpz_class& x,const mpz_class& y) {
    mpz_class g; mpz_gcd(g.get_mpz_t(),x.get_mpz_t(),y.get_mpz_t()); return g;
}
static mpz_class invz(const mpz_class& x,const mpz_class& m) {
    if(m==1) return 0;
    mpz_class z;
    if(!mpz_invert(z.get_mpz_t(),x.get_mpz_t(),m.get_mpz_t()))
        throw std::runtime_error("nonunit in Machin Lucas reconstruction");
    return z;
}
static mpz_class powerz(u64 p,u64 e) { return prime_power(p,e); }
static Rational canonical(Rational x) {
    if(x.d<=0) throw std::runtime_error("nonpositive Machin Lucas denominator");
    const auto g=gcdz(x.n,x.d);
    mpz_divexact(x.n.get_mpz_t(),x.n.get_mpz_t(),g.get_mpz_t());
    mpz_divexact(x.d.get_mpz_t(),x.d.get_mpz_t(),g.get_mpz_t());
    return x;
}
static bool equal(const Rational& x,const Rational& y) { return x.n*y.d==y.n*x.d; }
static Rational direct_sum(u64 A,u64 a,u64 C,u64 r) {
    auto root=machin_real_bs(A,a,C,0,r);
    return canonical({std::move(root.t),std::move(root.q)});
}
static mpz_class nearest_difference(const Rational& x,const Rational& y) {
    const mpz_class denominator=x.d*y.d,numerator=x.n*y.d-y.n*x.d;
    return floor_div(mpz_class(2*numerator+denominator),mpz_class(2*denominator));
}
}
