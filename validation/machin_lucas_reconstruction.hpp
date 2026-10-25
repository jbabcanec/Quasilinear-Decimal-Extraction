#pragma once
// Exact recursive reconstruction of the positive Machin Euler sums.
// Requires the weighted Machin integer/rational primitives before inclusion.
// This experimental implementation is not a claim of a full extractor bound.
#include "machin_lucas_coefficients.hpp"
#include "machin_lucas_adic.hpp"
#include <map>

namespace machin_lucas {

class Reconstruction {
    struct Factor { u64 p,pv; };
    struct Exceptional { u64 k,pv; std::size_t query; };
    struct Group {
        std::size_t block=0;
        u64 p=0,q=1,k0=0,low=0;
        std::vector<Exceptional> exceptional;
    };
    u64 machin_prime,D,C,base_terms,small_cutoff,quota;
    const std::vector<Block>& blocks;
    Counters& stats;
    bool audit,primary;
    u64 prime_bound=0;
    std::vector<u64> primes;

    u64 maximum_denominator(const Block& b,u64 r) const {
        if(!r || !b.A || !(b.A&1) || b.a<2 || (b.a&(b.a-1)))
            throw std::invalid_argument("invalid positive Machin progression");
        const u128 maximum=static_cast<u128>(b.A)+static_cast<u128>(b.a)*(r-1);
        if(maximum>UINT32_MAX)
            throw std::overflow_error("Lucas denominator exceeds audited 32-bit range");
        return static_cast<u64>(maximum);
    }
    std::vector<std::vector<Factor>> factor_segment(u64 A,u64 a,u64 n) const {
        std::vector<u64> residual(n);
        std::vector<std::vector<Factor>> factors(n);
        if(!n) return factors;
        const u64 maximum=A+a*(n-1);
        for(u64 k=0;k<n;++k) residual[k]=A+a*k;
        for(u64 p:primes) {
            if(p>maximum/p) break;
            const u64 first=mulmod((p-A%p)%p,invmod(a%p,p),p);
            for(u64 k=first;k<n;k+=p) {
                u64 pv=1;
                while(residual[k]%p==0) { residual[k]/=p; pv*=p; }
                if(pv>1) factors[k].push_back({p,pv});
            }
        }
        for(u64 k=0;k<n;++k)
            if(residual[k]>1) factors[k].push_back({residual[k],residual[k]});
        return factors;
    }
    void audit_sums(const std::vector<Rational>& values,u64 r) const {
        if(!audit) return;
        for(std::size_t i=0;i<blocks.size();++i)
            if(!equal(values[i],direct_sum(blocks[i].A,blocks[i].a,C,r)))
                throw std::runtime_error("recursive positive Machin sum mismatch");
    }

public:
    Reconstruction(u64 prime,const std::vector<Block>& input,u64 base,u64 small,
                   u64 forest_quota,Counters& counters,bool exact_audit=false,
                   bool direct_primary=true)
        :machin_prime(prime),D(prime*prime),C(D+1),base_terms(std::max<u64>(1,base)),
         small_cutoff(small),quota(std::max<u64>(1,forest_quota)),blocks(input),
         stats(counters),audit(exact_audit),primary(direct_primary) {
        if(prime<3 || !(prime&1) || static_cast<u128>(prime)*prime+1>UINT32_MAX)
            throw std::invalid_argument("unsupported Machin prime");
        u64 odd=C/2;
        while(odd%13==0) odd/=13;
        primary=primary && odd==1;
    }

    // Every positive Euler ratio up to r is at most ar/(CA).  With
    // ds >= ceil(log2 C)*(r-1)+18, the positive tail after scaling by
    // C^(r-1) is <2^-18*C/(C-1)<2^-17. No floating proposal is used.
    u64 child_length(u64 r) const {
        if(!r) throw std::invalid_argument("zero Euler truncation length");
        u64 d=UINT64_MAX;
        for(const auto& b:blocks) {
            maximum_denominator(b,r);
            const u128 ratio=static_cast<u128>(C)*b.A/(static_cast<u128>(b.a)*r);
            const u64 local=ratio ? bits(static_cast<u64>(ratio))-1 : 0;
            d=std::min(d,local);
        }
        if(blocks.empty() || !d) return r;
        const u64 target=bits(C-1)*(r-1)+18;
        return target/d+(target%d!=0);
    }

    std::vector<Rational> solve(u64 r) {
        if(!r) throw std::invalid_argument("zero Euler truncation length");
        if(blocks.empty()) return {};
        // This global-cap fixed-prime kernel has no improved asymptotic
        // bound yet. Keep it a fixed finite backend; larger problems use
        // the recursive shared-coefficient construction.
        const bool use_primary=primary && r<=16384;
        u64 maximum=0;
        for(const auto& b:blocks) maximum=std::max(maximum,maximum_denominator(b,r));
        const u64 smaller=child_length(r);
        if(r<=base_terms || (!use_primary && smaller>=r)) {
            std::vector<Rational> result; result.reserve(blocks.size());
            for(const auto& b:blocks) {
                result.push_back(direct_sum(b.A,b.a,C,r));
                ++stats.base_sums;
            }
            return result;
        }
        ++stats.levels;
        auto start=Clock::now();
        const u64 bound=isqrt(maximum);
        if(bound>prime_bound) { primes=primes_to(bound); prime_bound=bound; }

        std::vector<Group> groups;
        std::vector<Query> exceptional;
        for(std::size_t bi=0;bi<blocks.size();++bi) {
            const auto& b=blocks[bi];
            const auto factors=factor_segment(b.A+b.a,b.a,r-1);
            std::map<u64,Group> local;
            for(u64 j=0;j+1<r;++j) for(const auto f:factors[j]) {
                if(use_primary && f.p==13) continue;
                auto it=local.find(f.p);
                if(it==local.end()) {
                    Group g; g.block=bi; g.p=f.p; g.q=f.p;
                    g.k0=mulmod((f.p-b.A%f.p)%f.p,invmod(b.a%f.p,f.p),f.p);
                    it=local.emplace(f.p,std::move(g)).first;
                }
                auto& g=it->second; g.q=std::max(g.q,f.pv);
                // D=p0^2 is a nonunit only at the fixed Machin prime.
                // That entire row is computed backwards, never queried.
                if(f.p!=machin_prime && f.p>small_cutoff && f.pv>f.p) {
                    const auto qi=exceptional.size();
                    exceptional.push_back({j+1,f.p,f.pv,0});
                    g.exceptional.push_back({j+1,f.pv,qi});
                }
            }
            for(auto& entry:local) groups.push_back(std::move(entry.second));
        }
        std::sort(groups.begin(),groups.end(),[](const auto& x,const auto& y) {
            return x.p!=y.p ? x.p<y.p : x.block<y.block;
        });
        std::vector<Query> ordinary;
        std::vector<std::size_t> destinations;
        for(std::size_t j=0;j<groups.size();++j) {
            const auto& g=groups[j];
            if(g.p!=machin_prime && g.p>small_cutoff && g.p>r/g.p) {
                ordinary.push_back({g.k0,g.p,g.p,0}); destinations.push_back(j);
            }
        }
        stats.prime_queries+=ordinary.size(); stats.power_queries+=exceptional.size();
        stats.collect+=elapsed(start); start=Clock::now();
        if(!ordinary.empty() || !exceptional.empty()) {
            // A fixed finite cutoff changes constants without replacing the
            // asymptotic shared-tree backend. Generate the exact row once.
            if(r<=16384) {
                std::vector<Query> combined;
                combined.reserve(ordinary.size()+exceptional.size());
                combined.insert(combined.end(),ordinary.begin(),ordinary.end());
                combined.insert(combined.end(),exceptional.begin(),exceptional.end());
                evaluate_streamed(r,machin_prime,combined,quota,stats);
                for(std::size_t j=0;j<ordinary.size();++j) ordinary[j].value=combined[j].value;
                for(std::size_t j=0;j<exceptional.size();++j)
                    exceptional[j].value=combined[ordinary.size()+j].value;
            } else {
                ValueTree tree(r,machin_prime,stats);
                tree.evaluate(ordinary,false,quota);
                tree.evaluate(exceptional,true,quota);
            }
            for(std::size_t j=0;j<ordinary.size();++j)
                groups[destinations[j]].low=ordinary[j].value;
        }
        stats.trees+=elapsed(start); start=Clock::now();

        std::vector<std::vector<Rational>> parts(blocks.size());
        for(std::size_t begin=0;begin<groups.size();) {
            const u64 p=groups[begin].p;
            std::size_t end=begin+1;
            while(end<groups.size() && groups[end].p==p) ++end;
            const bool small=p<=small_cutoff || p==machin_prime;
            const bool medium=p<=r/p;
            const bool leaf_only=!small && p>r;
            std::vector<u64> full;
            WeightedRow high,low;
            const u64 r0=r%p,r1=r/p;
            const u64 powerC=leaf_only || small ? 0 : powmod(C,r,p);
            const u64 full_low=leaf_only || small ? 0 : powmod(C,r0,p);
            const u64 invD=leaf_only || small ? 0 : invmod(D%p,p);
            if(small) {
                u64 q=p; while(q<=maximum/p) q*=p;
                full=coefficient_row(r,machin_prime,p,q);
            } else if(!leaf_only) {
                high=weighted_row(r1,D,p);
                if(medium) low=weighted_row(r0,D,p);
            }
            for(std::size_t gi=begin;gi<end;++gi) {
                const auto& g=groups[gi]; const auto& b=blocks[g.block];
                u64 c=0,clow=0;
                const u64 inverse_scale=use_primary ? invmod(powmod(C,r-1,g.q),g.q) : 1;
                if(!small && !leaf_only) {
                    if(medium) clow=low.prefix[std::min(g.k0,r0)+1];
                    else {
                        u64 term=mulmod(powmod(D,g.k0+1,p),g.low,p);
                        if((g.k0&1) && term) term=p-term;
                        clow=(powerC+p-term)%p;
                    }
                }
                const u64 first=g.k0 ? g.k0 : p;
                u64 inverse_weight=(!small && !leaf_only) ? powmod(invD,first+1,p) : 0;
                if((first&1) && inverse_weight) inverse_weight=p-inverse_weight;
                std::size_t ei=0;
                for(u64 k=first;k<r;k+=p) {
                    ++stats.incidences;
                    const u64 mu=b.A+b.a*k;
                    u64 pv=p;
                    while(mu/pv%p==0) pv*=p;
                    u64 coefficient;
                    if(small) coefficient=full[k]%pv;
                    else if(pv>p) {
                        if(ei>=g.exceptional.size() || g.exceptional[ei].k!=k || g.exceptional[ei].pv!=pv)
                            throw std::runtime_error("missing Machin prime-power coefficient");
                        coefficient=exceptional[g.exceptional[ei++].query].value;
                    } else if(leaf_only) coefficient=g.low;
                    else {
                        const u64 h=k/p;
                        const u64 prefix=(mulmod(full_low,high.prefix[h],p)
                                          +mulmod(high.weight[h],clow,p))%p;
                        coefficient=mulmod((powerC+p-prefix)%p,inverse_weight,p);
                    }
                    u64 numerator=mulmod((b.a*k)%g.q,coefficient,g.q);
                    if(numerator) numerator=g.q-numerator;
                    if(use_primary) numerator=mulmod(numerator,inverse_scale,g.q);
                    u64 term=mulmod(g.q/pv,numerator,g.q);
                    term=mulmod(term,invmod((mu/pv)%g.q,g.q),g.q);
                    c=(c+term)%g.q;
                    if(!small && !leaf_only)
                        inverse_weight=mulmod(inverse_weight,p-invD,p);
                }
                if(c) parts[g.block].push_back({integer(c),integer(g.q)});
            }
            begin=end;
        }
        if(use_primary) {
            for(std::size_t bi=0;bi<blocks.size();++bi) {
                const auto& b=blocks[bi];
                auto value=fixed_prime_fraction(b.A,b.a,C,r,13);
                if(value.n!=0) parts[bi].push_back(std::move(value));
            }
        }
        std::vector<Rational> fractions; fractions.reserve(blocks.size());
        for(auto& leaves:parts) fractions.push_back(rational_tree(std::move(leaves)));
        stats.lucas+=elapsed(start);
        // Retain only rational results across recursion. Free old row/query
        // capacities, rather than clearing objects whose allocations survive.
        std::vector<Group>().swap(groups);
        std::vector<Query>().swap(ordinary);
        std::vector<Query>().swap(exceptional);
        std::vector<std::size_t>().swap(destinations);
        std::vector<std::vector<Rational>>().swap(parts);
        if(use_primary) {
            // H_r lies between1 and C/(C-1), so approximation1 is already
            // within1/2. The missing integer follows without a child sum.
            start=Clock::now();
            for(auto& value:fractions) {
                const mpz_class correction=nearest_difference({mpz_class(1),mpz_class(1)},value);
                value.n+=correction*value.d;
            }
            stats.lucas+=elapsed(start);
            audit_sums(fractions,r);
            return fractions;
        }
        auto child=solve(smaller);
        start=Clock::now();
        const mpz_class scale=powerz(C,r-1);
        for(std::size_t bi=0;bi<blocks.size();++bi) {
            const Rational approximation{scale*child[bi].n,child[bi].d};
            const mpz_class correction=nearest_difference(approximation,fractions[bi]);
            fractions[bi]=canonical({fractions[bi].n+correction*fractions[bi].d,
                                      fractions[bi].d*scale});
        }
        stats.lucas+=elapsed(start);
        audit_sums(fractions,r);
        return fractions;
    }
};

} // namespace machin_lucas
