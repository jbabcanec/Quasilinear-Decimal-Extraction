// Shared Lucas reduction, sparse prime-power lifting, recursive exact Euler
// reconstruction, and exact certified decimal output. See the proof note
// notes/lucas_sparse_prime_power_batch_2026-09-07.md.
// Build with GMP: g++ -O3 -march=native -std=c++17 decimal_extract_lucas.cpp -lgmpxx -lgmp -o extract_lucas
// Reuse the audited integer, rational-tree, and weighted-forest primitives.
#define main weighted_reference_main
#include "decimal_extract_weighted.cpp"
#undef main
#include <map>
#include <sstream>

namespace lucas {
struct Options {
    u64 block=0,small=0,base=0,quota=0,chunk=0,guard=192;
    bool direct=false,audit=false,timing=false;
};
struct Counters {
    u64 blocks=0,levels=0,base_sums=0,root_terms=0,prime_queries=0,power_queries=0;
    u64 ordinary_batches=0,power_batches=0,nonunit_divisions=0,incidences=0;
    u64 registry_passes=0,word_residues=0,certificate_passes=0;
    u64 max_root_bits=0,block_length=0,chunk_length=0,guard=0;
    double word=0,core=0,corrections=0,collect=0,trees=0,lucas=0,assembly=0,total=0;
};
static mpz_class modz(mpz_class x,const mpz_class& m) {
    mpz_fdiv_r(x.get_mpz_t(),x.get_mpz_t(),m.get_mpz_t()); return x;
}
static mpz_class gcdz(const mpz_class& x,const mpz_class& y) {
    mpz_class g; mpz_gcd(g.get_mpz_t(),x.get_mpz_t(),y.get_mpz_t()); return g;
}
static mpz_class invz(const mpz_class& x,const mpz_class& m) {
    if(m==1) return 0;
    mpz_class v;
    if(!mpz_invert(v.get_mpz_t(),x.get_mpz_t(),m.get_mpz_t()))
        throw std::runtime_error("nonunit in exact remainder tree");
    return v;
}
static mpz_class powerz(u64 a,u64 e) {
    if(a>ULONG_MAX || e>ULONG_MAX) throw std::overflow_error("power word range");
    mpz_class z; mpz_ui_pow_ui(z.get_mpz_t(),static_cast<unsigned long>(a),static_cast<unsigned long>(e)); return z;
}
static u64 scaled_power_ceiling(u64 value,u64 power,u64 scale) {
    mpz_class numerator=powerz(value,power),denominator=powerz(scale,power);
    mpz_class quotient=(numerator+denominator-1)/denominator;
    if(bits(quotient)>63) throw std::overflow_error("parameter range exceeded");
    u64 result=0; std::size_t count=0;
    mpz_export(&result,&count,1,sizeof(result),0,0,quotient.get_mpz_t());
    return result;
}
static Options resolve_options(u64 N,Options o) {
    u64 L=std::max<u64>(1,bits(N-1));
    // Fixed rescalings preserve B=Theta(L^32), Y=Theta(L^8), and a
    // Theta(sqrt(N)) terminal cutoff asymptotically, while permitting the
    // genuinely new paths to run at accessible input sizes.
    u64 block_count=std::max<u64>(128,scaled_power_ceiling(L,32,32));
    if(!o.block) o.block=std::max<u64>(1,(N+block_count-1)/block_count);
    if(!o.small) o.small=std::max<u64>(64,scaled_power_ceiling(L+2,8,32));
    if(!o.base) o.base=std::max<u64>(8,isqrt(N)/1048576);
    return o;
}
static Rational canonical(Rational q) {
    if(q.d<=0) throw std::runtime_error("nonpositive rational denominator");
    auto g=gcdz(q.n,q.d); q.n/=g; q.d/=g; return q;
}
static bool equal(const Rational& x,const Rational& y) { return x.n*y.d==y.n*x.d; }
static Rational direct_sum(u64 A,u64 a,u64 r) {
    auto p=real_bs(A,a,0,r); return canonical({std::move(p.t),std::move(p.q)});
}
static mpz_class nearest_difference(const Rational& x,const Rational& y) {
    mpz_class d=x.d*y.d,n=x.n*y.d-y.n*x.d;
    return floor_div(mpz_class(2*n+d),mpz_class(2*d));
}
static u64 joint_power(u64 n,i64 e,u64 m) {
    u64 power=e<0 ? static_cast<u64>(-e) : static_cast<u64>(e),mask=n|power,z=1%m;
    if(!mask) return z;
    for(int b=63-__builtin_clzll(mask);b>=0;--b) {
        z=mulmod(z,z,m);
        if((n>>b)&1) z=5*z%m;
        if((power>>b)&1) {
            if(e<0) { if(z&1) z+=m; z>>=1; }
            else { z+=z; if(z>=m) z-=m; }
        }
    }
    return z;
}
struct Factor { u64 p,pv; };
static std::vector<std::vector<Factor>> factor_segment(u64 A,u64 a,u64 n,const std::vector<u64>& primes) {
    std::vector<u64> rem(n); std::vector<std::vector<Factor>> out(n);
    if(!n) return out;
    u64 maximum=A+a*(n-1);
    if(maximum>UINT32_MAX) throw std::overflow_error("denominator exceeds audited 32-bit word range");
    for(u64 k=0;k<n;++k) rem[k]=A+a*k;
    for(u64 p:primes) {
        if(p>maximum/p) break;
        u64 first=mulmod((p-A%p)%p,invmod(a,p),p);
        for(u64 k=first;k<n;k+=p) {
            u64 pv=1;
            while(rem[k]%p==0) { rem[k]/=p; pv*=p; }
            if(pv>1) out[k].push_back({p,pv});
        }
    }
    for(u64 k=0;k<n;++k) if(rem[k]>1) out[k].push_back({rem[k],rem[k]});
    return out;
}
struct Block { u64 A,a,length,shift; int sign; };
struct Query { u64 k,p,mod,value=0; };
struct Exceptional { u64 k,pv; std::size_t query; };
struct Group {
    std::size_t block;
    u64 p,q,k0,low=0;
    std::vector<Exceptional> exceptional;
};

// Integer numerator entries (a,b;0,c), with a positive common denominator.
struct Transfer { mpz_class a=1,b=0,c=1,d=1; };
static Transfer compose(const Transfer& after,const Transfer& before,bool cancel) {
    Transfer t{after.a*before.a,after.a*before.b+after.b*before.c,
               after.c*before.c,after.d*before.d};
    if(cancel) {
        auto g=gcdz(t.d,t.a); g=gcdz(g,t.b); g=gcdz(g,t.c);
        t.a/=g; t.b/=g; t.c/=g; t.d/=g;
    }
    return t;
}
class ValueTree {
    struct Node {
        u64 lo=0,hi=0,depth=0;
        int left=-1,right=-1;
        Transfer t;
        mpz_class q=1,s=1,m=1;
        std::vector<std::size_t> queries;
        std::vector<mpz_class> query_products;
    };
    std::vector<Node> nodes;
    std::vector<int> leaves;
    u64 r,height=0;
    mpz_class b0,e0;
    Counters& stats;
    int build(u64 lo,u64 hi,u64 depth) {
        int idx=static_cast<int>(nodes.size()); nodes.emplace_back();
        nodes[idx].lo=lo; nodes[idx].hi=hi; nodes[idx].depth=depth;
        height=std::max(height,depth);
        if(hi-lo==1) {
            u64 k=lo+1;
            nodes[idx].t={integer(k),-integer(r-k+1),-integer(16*(r-k+1)),integer(16*k)};
            auto& t=nodes[idx].t;
            auto g=gcdz(t.d,t.a); g=gcdz(g,t.b); g=gcdz(g,t.c);
            t.a/=g; t.b/=g; t.c/=g; t.d/=g;
            leaves[lo]=idx;
        } else {
            u64 mid=lo+(hi-lo)/2;
            int left=build(lo,mid,depth+1),right=build(mid,hi,depth+1);
            nodes[idx].left=left; nodes[idx].right=right;
            nodes[idx].t=compose(nodes[right].t,nodes[left].t,true);
        }
        return idx;
    }
    u64 lambda(u64 p) const {
        u64 power=1,e=0;
        do { ++e; if(power>(r+1)/p) break; power*=p; } while(power<r+1);
        return e;
    }
    void modulus_tree(int index,bool sparse,const std::vector<Query>& queries) {
        auto& v=nodes[index];
        if(v.left>=0) {
            modulus_tree(v.left,sparse,queries); modulus_tree(v.right,sparse,queries);
            v.q=nodes[v.left].q*nodes[v.right].q;
            v.s=nodes[v.left].s*nodes[v.right].s;
        } else if(!v.queries.empty()) {
            const auto count=v.queries.size();
            v.query_products.assign(2*count,mpz_class(1));
            std::vector<mpz_class> precision;
            if(sparse) precision.assign(2*count,mpz_class(1));
            for(std::size_t j=0;j<count;++j) {
                const auto& q=queries[v.queries[j]];
                v.query_products[count+j]=integer(q.mod);
                if(sparse) precision[count+j]=powerz(q.p,lambda(q.p));
            }
            // Balance products within a shared leaf as well as between indices.
            for(std::size_t j=count;j-->1;) {
                v.query_products[j]=v.query_products[2*j]*v.query_products[2*j+1];
                if(sparse) precision[j]=precision[2*j]*precision[2*j+1];
            }
            v.q=v.query_products[1];
            if(sparse) v.s=std::move(precision[1]);
        }
        if(sparse) {
            mpz_class power;
            mpz_pow_ui(power.get_mpz_t(),v.s.get_mpz_t(),static_cast<unsigned long>(height-v.depth));
            v.m=v.q*power;
        } else v.m=v.q;
    }
    void distribute_leaf(const Node& leaf,std::size_t index,mpz_class value,std::vector<Query>& queries) {
        value=modz(std::move(value),leaf.query_products[index]);
        const auto count=leaf.queries.size();
        if(index>=count) {
            // value is already reduced to this query's checked word modulus.
            queries[leaf.queries[index-count]].value=mpz_get_ui(value.get_mpz_t());
        } else {
            distribute_leaf(leaf,2*index,value,queries);
            distribute_leaf(leaf,2*index+1,std::move(value),queries);
        }
    }
    void plain_visit(int idx,mpz_class x,mpz_class y,mpz_class d,std::vector<Query>& queries) {
        auto& v=nodes[idx]; if(v.m==1) return;
        x=modz(std::move(x),v.m); y=modz(std::move(y),v.m); d=modz(std::move(d),v.m);
        if(v.left<0) {
            mpz_class result=modz(mpz_class(x*invz(d,v.m)),v.m);
            distribute_leaf(v,1,std::move(result),queries);
            return;
        }
        plain_visit(v.left,x,y,d,queries);
        auto& target=nodes[v.right]; if(target.m==1) return;
        auto& t=nodes[v.left].t;
        x=modz(std::move(x),target.m); y=modz(std::move(y),target.m); d=modz(std::move(d),target.m);
        mpz_class nx=modz(mpz_class(modz(t.a,target.m)*x+modz(t.b,target.m)*y),target.m);
        mpz_class ny=modz(mpz_class(modz(t.c,target.m)*y),target.m);
        mpz_class nd=modz(mpz_class(modz(t.d,target.m)*d),target.m);
        plain_visit(v.right,std::move(nx),std::move(ny),std::move(nd),queries);
    }
    void sparse_visit(int idx,mpz_class x,mpz_class y,std::vector<Query>& queries) {
        auto& v=nodes[idx]; if(v.m==1) return;
        x=modz(std::move(x),v.m); y=modz(std::move(y),v.m);
        if(v.left<0) {
            distribute_leaf(v,1,std::move(x),queries);
            return;
        }
        sparse_visit(v.left,x,y,queries);
        auto& target=nodes[v.right]; if(target.m==1) return;
        auto& t=nodes[v.left].t;
        mpz_class source=target.m*target.s;
        x=modz(std::move(x),source); y=modz(std::move(y),source);
        mpz_class wx=modz(mpz_class(modz(t.a,source)*x+modz(t.b,source)*y),source);
        mpz_class wy=modz(mpz_class(modz(t.c,source)*y),source);
        auto g=gcdz(t.d,source);
        if(g!=1) ++stats.nonunit_divisions;
        if(!mpz_divisible_p(wx.get_mpz_t(),g.get_mpz_t()) ||
           !mpz_divisible_p(wy.get_mpz_t(),g.get_mpz_t()))
            throw std::runtime_error("precision-loss invariant failed");
        mpz_class remaining=source/g;
        if(!mpz_divisible_p(remaining.get_mpz_t(),target.m.get_mpz_t()))
            throw std::runtime_error("insufficient inflated precision");
        mpz_class unit=t.d/g,inv=invz(unit,target.m);
        wx/=g; wy/=g;
        wx=modz(mpz_class(wx*inv),target.m); wy=modz(mpz_class(wy*inv),target.m);
        sparse_visit(v.right,std::move(wx),std::move(wy),queries);
    }
public:
    ValueTree(u64 length,Counters& st):leaves(length),r(length),stats(st) {
        nodes.reserve(2*r); build(0,r,0);
        b0=(powerz(15,r)-((r&1) ? -1 : 1))/16;
        e0=(r&1) ? 1 : -1;
    }
    void evaluate(std::vector<Query>& queries,bool sparse,u64 quota) {
        std::vector<std::size_t> order(queries.size());
        for(std::size_t j=0;j<order.size();++j) order[j]=j;
        std::sort(order.begin(),order.end(),[&](auto i,auto j){return queries[i].k<queries[j].k;});
        for(std::size_t begin=0;begin<order.size();) {
            for(auto& v:nodes) {
                // Release old GMP capacities when advancing a bounded forest.
                mpz_class(1).swap(v.q); mpz_class(1).swap(v.s); mpz_class(1).swap(v.m);
                v.queries.clear(); v.query_products.clear();
            }
            std::size_t end=begin; u64 mass=0;
            while(end<order.size()) {
                auto qi=order[end]; auto& q=queries[qi];
                if(q.k>=r || q.mod>ULONG_MAX) throw std::runtime_error("invalid remainder query");
                u64 cost=bits(q.mod)+(sparse ? height*lambda(q.p)*bits(q.p) : 0);
                if(end>begin && mass+cost>quota) break;
                auto& leaf=nodes[leaves[q.k]];
                leaf.queries.push_back(qi); mass+=cost; ++end;
            }
            modulus_tree(0,sparse,queries);
            stats.max_root_bits=std::max(stats.max_root_bits,bits(nodes[0].m));
            if(sparse) { ++stats.power_batches; sparse_visit(0,b0,e0,queries); }
            else { ++stats.ordinary_batches; plain_visit(0,b0,e0,mpz_class(1),queries); }
            begin=end;
        }
    }
};

// Complete binomial row modulo a prime power, retaining its unit part even
// when the residue vanishes. This avoids all division by nonunits.
static std::vector<u64> binomial_row(u64 n,u64 p,u64 q) {
    std::vector<u64> row(n+1); row[0]=1%q;
    i64 valuation=0; u64 unit=1,exponent=0,temp=q;
    while(temp>1) { temp/=p; ++exponent; }
    for(u64 j=1;j<=n;++j) {
        u64 num=n-j+1,den=j;
        while(num%p==0) { num/=p; ++valuation; }
        while(den%p==0) { den/=p; --valuation; }
        if(valuation<0) throw std::runtime_error("negative binomial valuation");
        unit=mulmod(unit,num%q,q); unit=mulmod(unit,invmod(den%q,q),q);
        row[j]=static_cast<u64>(valuation)>=exponent ? 0 : mulmod(unit,powmod(p,valuation,q),q);
    }
    return row;
}
static std::vector<u64> coefficient_row(u64 r,u64 p,u64 q) {
    auto choose=binomial_row(r,p,q); std::vector<u64> row(r);
    u64 sign=(r&1) ? q-1 : 1,inv16=invmod(16,q);
    row[0]=mulmod((powmod(15,r,q)+q-sign)%q,inv16,q);
    for(u64 k=1;k<r;++k) {
        u64 e=((r-k+1)&1) && choose[k] ? q-choose[k] : choose[k];
        row[k]=mulmod((row[k-1]+e)%q,inv16,q);
    }
    return row;
}
struct WeightedRow { std::vector<u64> weight,prefix; };
static WeightedRow weighted_row(u64 n,u64 p) {
    auto choose=binomial_row(n,p,p);
    WeightedRow row; row.weight.resize(n+1); row.prefix.resize(n+2);
    u64 power=1,c=(p-16%p)%p;
    for(u64 k=0;k<=n;++k) {
        row.weight[k]=mulmod(choose[k],power,p);
        row.prefix[k+1]=(row.prefix[k]+row.weight[k])%p;
        power=mulmod(power,c,p);
    }
    return row;
}

class Reconstruction {
    const std::vector<Block>& blocks;
    const std::vector<u64>& primes;
    const Options& options;
    Counters& stats;
    u64 quota;
    void audit_sums(const std::vector<Rational>& sums,u64 r) {
        if(!options.audit) return;
        for(std::size_t b=0;b<blocks.size();++b)
            if(!equal(sums[b],direct_sum(blocks[b].A,blocks[b].a,r)))
                throw std::runtime_error("recursive Euler rational mismatch");
    }
public:
    Reconstruction(const std::vector<Block>& b,const std::vector<u64>& p,
                   const Options& o,Counters& s,u64 q):blocks(b),primes(p),options(o),stats(s),quota(q) {}
    std::vector<Rational> solve(u64 r) {
        u64 d=UINT64_MAX;
        for(auto& b:blocks) {
            u64 ratio=15*b.A/(b.a*r);
            d=std::min(d,ratio ? bits(ratio)-1 : 0);
        }
        u64 smaller=d ? (4*r+16+d-1)/d : r;
        if(options.direct || r<=options.base || smaller>=r) {
            std::vector<Rational> out; out.reserve(blocks.size());
            for(auto& b:blocks) { out.push_back(direct_sum(b.A,b.a,r)); ++stats.base_sums; }
            return out;
        }
        ++stats.levels;
        auto start=Clock::now();
        std::vector<Group> groups; std::vector<Query> exceptional;
        u64 maximum=0;
        for(std::size_t bi=0;bi<blocks.size();++bi) {
            auto& b=blocks[bi]; maximum=std::max(maximum,b.A+b.a*(r-1));
            auto factors=factor_segment(b.A+b.a,b.a,r-1,primes);
            std::map<u64,Group> local;
            for(u64 j=0;j+1<r;++j) for(auto f:factors[j]) {
                auto it=local.find(f.p);
                if(it==local.end()) {
                    Group g; g.block=bi; g.p=f.p; g.q=f.p;
                    g.k0=mulmod((f.p-b.A%f.p)%f.p,invmod(b.a,f.p),f.p);
                    it=local.emplace(f.p,std::move(g)).first;
                }
                auto& g=it->second; g.q=std::max(g.q,f.pv);
                if(f.p>options.small && f.pv>f.p) {
                    std::size_t qi=exceptional.size();
                    exceptional.push_back({j+1,f.p,f.pv,0});
                    g.exceptional.push_back({j+1,f.pv,qi});
                }
            }
            for(auto& kv:local) groups.push_back(std::move(kv.second));
        }
        std::sort(groups.begin(),groups.end(),[](const auto& x,const auto& y){
            return x.p!=y.p ? x.p<y.p : x.block<y.block;
        });
        std::vector<Query> ordinary; std::vector<std::size_t> destinations;
        for(std::size_t j=0;j<groups.size();++j) {
            auto& g=groups[j];
            if(g.p>options.small && g.p>r/g.p) {
                ordinary.push_back({g.k0,g.p,g.p,0}); destinations.push_back(j);
            }
        }
        stats.prime_queries+=ordinary.size(); stats.power_queries+=exceptional.size();
        stats.collect+=elapsed(start); start=Clock::now();
        if(!ordinary.empty() || !exceptional.empty()) {
            ValueTree tree(r,stats);
            tree.evaluate(ordinary,false,quota); tree.evaluate(exceptional,true,quota);
            for(std::size_t j=0;j<ordinary.size();++j) groups[destinations[j]].low=ordinary[j].value;
        }
        stats.trees+=elapsed(start); start=Clock::now();
        std::vector<std::vector<Rational>> parts(blocks.size());
        for(std::size_t begin=0;begin<groups.size();) {
            u64 p=groups[begin].p; std::size_t end=begin+1;
            while(end<groups.size() && groups[end].p==p) ++end;
            std::vector<u64> full;
            WeightedRow high,low;
            const bool small=p<=options.small,medium=p<=r/p,leaf_only=!small && p>r;
            u64 power15=leaf_only ? 0 : powmod(15,r,p),r0=r%p,r1=r/p;
            u64 full_low=leaf_only ? 0 : powmod((p-15%p)%p,r0,p);
            u64 inv16=leaf_only || small ? 0 : invmod(16,p);
            if(small) {
                u64 q=p; while(q<=maximum/p) q*=p;
                full=coefficient_row(r,p,q);
            } else if(!leaf_only) {
                high=weighted_row(r1,p);
                if(medium) low=weighted_row(r0,p);
            }
            for(std::size_t gi=begin;gi<end;++gi) {
                auto& g=groups[gi]; auto& b=blocks[g.block]; u64 c=0,clow=0;
                if(!small && !leaf_only) {
                    if(medium) clow=low.prefix[std::min(g.k0,r0)+1];
                    else {
                        clow=(power15+p-mulmod(powmod(16,g.k0+1,p),g.low,p))%p;
                        if((r&1) && clow) clow=p-clow;
                    }
                }
                const u64 first=g.k0 ? g.k0 : p;
                u64 inverse_weight=(!small && !leaf_only) ? powmod(inv16,first+1,p) : 0;
                std::size_t ei=0;
                for(u64 k=g.k0 ? g.k0 : p;k<r;k+=p) {
                    ++stats.incidences;
                    u64 mu=b.A+b.a*k,pv=p;
                    while(mu/pv%p==0) pv*=p;
                    u64 coefficient;
                    if(small) coefficient=full[k]%pv;
                    else if(pv>p) {
                        if(ei>=g.exceptional.size() || g.exceptional[ei].k!=k || g.exceptional[ei].pv!=pv)
                            throw std::runtime_error("missing prime-power coefficient");
                        coefficient=exceptional[g.exceptional[ei++].query].value;
                    } else if(leaf_only) coefficient=g.low;
                    else {
                        u64 h=k/p;
                        u64 prefix=(mulmod(full_low,high.prefix[h],p)+mulmod(high.weight[h],clow,p))%p;
                        if((r&1) && prefix) prefix=p-prefix;
                        coefficient=mulmod((power15+p-prefix)%p,inverse_weight,p);
                    }
                    u64 numerator=mulmod(b.a*k%g.q,coefficient,g.q);
                    if(numerator) numerator=g.q-numerator;
                    u64 term=mulmod(g.q/pv,numerator,g.q);
                    term=mulmod(term,invmod(mu/pv,g.q),g.q);
                    c=(c+term)%g.q;
                    if(!small && !leaf_only) inverse_weight=mulmod(inverse_weight,inv16,p);
                }
                if(c) parts[g.block].push_back({integer(c),integer(g.q)});
            }
            begin=end;
        }
        std::vector<Rational> fractions; fractions.reserve(blocks.size());
        for(auto& p:parts) fractions.push_back(rational_tree(std::move(p)));
        stats.lucas+=elapsed(start);
        // Discard factor/query state before recursion; retain only its rational results.
        groups.clear(); groups.shrink_to_fit(); ordinary.clear(); ordinary.shrink_to_fit();
        exceptional.clear(); exceptional.shrink_to_fit();
        auto child=solve(smaller); mpz_class scale=powerz(15,r-1);
        for(std::size_t bi=0;bi<blocks.size();++bi) {
            Rational approximate{scale*child[bi].n,child[bi].d};
            auto z=nearest_difference(approximate,fractions[bi]);
            fractions[bi]=canonical({fractions[bi].n+z*fractions[bi].d,fractions[bi].d*scale});
        }
        audit_sums(fractions,r); return fractions;
    }
};

static Rational original_block(const Block& b,const Rational& H,const std::vector<u64>& primes) {
    auto factors=factor_segment(b.A,b.a,b.length,primes);
    struct Prime { u64 q=1,c=0; }; std::map<u64,Prime> values;
    for(auto& row:factors) for(auto f:row) values[f.p].q=std::max(values[f.p].q,f.pv);
    for(u64 k=0;k<b.length;++k) {
        u64 mu=b.A+b.a*k,w=powmod(16,b.length-1-k,mu);
        for(auto f:factors[k]) {
            auto& item=values[f.p];
            u64 term=mulmod(item.q/f.pv,w,item.q);
            term=mulmod(term,invmod(mu/f.pv,item.q),item.q);
            item.c=(item.c+term)%item.q;
        }
    }
    std::vector<Rational> parts;
    for(auto& kv:values) if(kv.second.c) parts.push_back({integer(kv.second.c),integer(kv.second.q)});
    auto T=rational_tree(std::move(parts));
    Rational approx{H.n<<bitcount(4*b.length),H.d*integer(15*b.A)};
    auto z=nearest_difference(approx,T);
    return canonical({T.n+z*T.d,T.d});
}
struct CoreResult { mpz_class value; u64 precision; };
static CoreResult exact_core(u64 N,u64 G,const Options& o,Counters& stats) {
    u64 n=N-1,K=n+(G+3)/4+2,m=4*K-n+1;
    const u64 L=std::max<u64>(1,bits(N-1));
    const u64 length=o.block ? o.block : std::max<u64>(1,(N+127)/128);
    const u64 quota=o.quota ? o.quota : std::max<u64>(4096,256*N/(L*L));
    stats.block_length=length;
    std::vector<Block> blocks; u64 ell=0,minimum_log=UINT64_MAX;
    for(auto f:families) {
        i64 e0=static_cast<i64>(n)+f.c;
        u64 first=e0>=0 ? e0/4+1 : 0;
        for(u64 lo=first;lo<=K;lo+=length) {
            u64 hi=std::min(K+1,lo+length),t=hi-lo,A=f.a*lo+f.b;
            blocks.push_back({A,f.a,t,m-(4*(hi-1)-e0),f.sign}); ell=std::max(ell,t);
        }
    }
    for(auto& b:blocks) {
        u64 ratio=15*b.A/(b.a*ell);
        minimum_log=std::min(minimum_log,ratio ? bits(ratio)-1 : 0);
    }
    u64 r=minimum_log ? (4*ell+16+minimum_log-1)/minimum_log : ell+1;
    if(r>ell) r=(4*ell+18)/3; // the uniform ratio |h_(j+1)/h_j| < 1/8
    u64 maximum=1;
    for(auto& b:blocks) maximum=std::max(maximum,b.A+b.a*std::max(r,b.length));
    if(maximum>UINT32_MAX) throw std::overflow_error("core exceeds audited denominator range");
    auto primes=primes_to(isqrt(maximum));
    stats.blocks+=blocks.size(); stats.root_terms=r;
    auto start=Clock::now();
    Reconstruction reconstruction(blocks,primes,o,stats,quota);
    auto sums=reconstruction.solve(r);
    stats.corrections+=elapsed(start);
    std::vector<Node> pending; u64 mass=0; mpz_class V=0;
    auto flush=[&]() {
        if(pending.empty()) return;
        auto tick=Clock::now();
        V+=merge_nodes(std::move(pending),m); reduce(V,m); pending.clear(); mass=0;
        stats.assembly+=elapsed(tick);
    };
    int previous_sign=0; u64 previous_a=0;
    for(std::size_t bi=0;bi<blocks.size();++bi) {
        auto& b=blocks[bi];
        if(b.a!=previous_a || b.sign!=previous_sign) flush();
        auto value=original_block(b,sums[bi],primes);
        if(o.audit) {
            mpq_class ref=0;
            for(u64 k=0;k<b.length;++k) ref+=mpq_class(mpz_class(1)<<bitcount(4*(b.length-1-k)),integer(b.A+b.a*k));
            if(value.n*ref.get_den()!=value.d*ref.get_num()) throw std::runtime_error("weighted block rational mismatch");
        }
        mass+=bits(value.n)+bits(value.d)+4*b.length;
        if(b.sign<0) value.n=-value.n;
        pending.push_back({std::move(value.n),std::move(value.d),b.shift});
        previous_a=b.a; previous_sign=b.sign;
        if(mass>=m) flush();
    }
    flush(); return {std::move(V),m};
}

static Rational exact_words(u64 N,u64 G,const Options& o,Counters& stats) {
    const u64 n=N-1,K=n+(G+3)/4+2,maximum=8*K+5,L=std::max<u64>(1,bits(N-1));
    if(maximum>UINT32_MAX) throw std::overflow_error("word registry exceeds audited denominator range");
    auto all_primes=primes_to(maximum),base=primes_to(isqrt(maximum));
    struct Registry { u64 p,q,c; };
    std::vector<Registry> registry; registry.reserve(all_primes.size());
    for(u64 p:all_primes) { u64 q=p; while(q<=maximum/p) q*=p; registry.push_back({p,q,0}); }
    all_primes.clear(); all_primes.shrink_to_fit();
    u64 chunk=o.chunk ? o.chunk : std::max<u64>(1,(64*N+L*L-1)/(L*L));
    stats.chunk_length=chunk;
    struct Record { u64 p,pv,mu,weight; int sign; };
    const u64 g=bits(4*(K+1)+maximum+1)+5;
    mpz_class real=0;
    for(auto f:families) for(u64 lo=0;lo<=K;lo+=chunk) {
        u64 count=std::min(chunk,K+1-lo);
        auto factors=factor_segment(f.a*lo+f.b,f.a,count,base);
        std::vector<Record> records;
        for(u64 j=0;j<count;++j) {
            u64 k=lo+j,mu=f.a*k+f.b;
            if(mu==1) continue;
            i64 e=static_cast<i64>(n)+f.c-4*static_cast<i64>(k);
            u64 w=joint_power(n,e,mu); ++stats.word_residues;
            mpz_class term=integer(w)<<bitcount(g);
            if(f.sign<0) term=-term;
            real+=floor_div(term,integer(mu));
            for(auto pf:factors[j]) records.push_back({pf.p,pf.pv,mu,w,f.sign});
        }
        std::sort(records.begin(),records.end(),[](auto& a,auto& b){return a.p<b.p;});
        // Explicit full sorted-registry pass, including copies, as in the proof.
        std::vector<Registry> next; next.reserve(registry.size()); std::size_t ri=0;
        for(auto entry:registry) {
            while(ri<records.size() && records[ri].p==entry.p) {
                auto& rec=records[ri++];
                u64 term=mulmod(entry.q/rec.pv,rec.weight,entry.q);
                term=mulmod(term,invmod(rec.mu/rec.pv,entry.q),entry.q);
                if(rec.sign<0 && term) term=entry.q-term;
                entry.c=(entry.c+term)%entry.q;
            }
            next.push_back(entry);
        }
        if(ri!=records.size()) throw std::runtime_error("unrouted word prime");
        registry=std::move(next); ++stats.registry_passes;
    }
    mpz_class lower=0; std::vector<Rational> parts; parts.reserve(registry.size());
    for(auto entry:registry) if(entry.c) {
        lower+=floor_div(mpz_class(integer(entry.c)<<bitcount(g)),integer(entry.q));
        parts.push_back({integer(entry.c),integer(entry.q)});
    }
    mpz_class correction=real-lower+(mpz_class(1)<<bitcount(g-1));
    mpz_fdiv_q_2exp(correction.get_mpz_t(),correction.get_mpz_t(),bitcount(g));
    auto result=rational_tree(std::move(parts)); result.n+=correction*result.d;
    return canonical(std::move(result));
}

static std::string extract(u64 N,const Options& options,Counters& stats) {
    auto all=Clock::now(); u64 G=options.guard;
    for(;;) {
        ++stats.certificate_passes; stats.guard=G;
        auto t=Clock::now(); auto R=exact_words(N,G,options,stats); stats.word+=elapsed(t);
        t=Clock::now(); auto core=exact_core(N,G,options,stats); stats.core+=elapsed(t);
        auto X=powerz(5,N-1); mpz_class W=X*core.value; reduce(W,core.precision);
        mpz_class Q=R.d<<bitcount(core.precision);
        mpz_class U=modz(R.n,R.d); U<<=bitcount(core.precision); U+=W*R.d; U=modz(std::move(U),Q);
        auto scale=powerz(10,12); mpz_class numerator=scale*U,digit,remainder;
        mpz_fdiv_qr(digit.get_mpz_t(),remainder.get_mpz_t(),numerator.get_mpz_t(),Q.get_mpz_t());
        mpz_class left=remainder<<bitcount(G),right=(Q-remainder)<<bitcount(G),margin=scale*Q;
        if(left>margin && right>margin) {
            stats.total=elapsed(all); std::string result=digit.get_str();
            if(result.size()>12) throw std::runtime_error("decimal window width");
            return std::string(12-result.size(),'0')+result;
        }
        if(G>UINT32_MAX/2) throw std::overflow_error("certificate precision range exceeded");
        G*=2;
    }
}
static void print(const Counters& s,const Options& o,u64 N) {
    std::ostringstream out;
    out.precision(10);
    out<<"{\"engine\":\"lucas-sparse-certified\",\"correction_mode\":\""<<(o.direct?"direct":"recursive")
       <<"\",\"N\":"<<N<<",\"certified\":true,\"profile\":\"experimental-scaled-parameters\"";
#define FIELD(name) out<<",\"" #name "\":"<<s.name
    FIELD(blocks); FIELD(levels); FIELD(base_sums); FIELD(root_terms); FIELD(prime_queries); FIELD(power_queries);
    FIELD(ordinary_batches); FIELD(power_batches); FIELD(nonunit_divisions); FIELD(incidences);
    FIELD(registry_passes); FIELD(word_residues); FIELD(certificate_passes); FIELD(max_root_bits);
    FIELD(block_length); FIELD(chunk_length); FIELD(guard); FIELD(word); FIELD(core); FIELD(corrections);
    FIELD(collect); FIELD(trees); FIELD(lucas); FIELD(assembly); FIELD(total);
#undef FIELD
    out<<",\"small_prime_cutoff\":"<<o.small<<",\"base_terms\":"<<o.base<<"}";
    std::cerr<<out.str()<<'\n';
}
static void check() {
    Counters counters; u64 residue_cases=0,reconstruction_cases=0,core_cases=0;
    for(u64 r:{2ULL,7ULL,16ULL,31ULL,65ULL,127ULL}) {
        ValueTree tree(r,counters); std::vector<Query> ordinary,sparse;
        std::vector<mpz_class> exact(r); exact[0]=(powerz(15,r)-((r&1)?-1:1))/16;
        mpz_class choose=1;
        for(u64 k=1;k<r;++k) {
            choose*=integer(r-k+1); choose/=integer(k);
            mpz_class e=((r-k+1)&1) ? mpz_class(-choose) : choose;
            exact[k]=(exact[k-1]+e)/16;
        }
        for(u64 p:{3ULL,5ULL,7ULL,13ULL,67ULL,131ULL}) {
            auto row=coefficient_row(r,p,p*p*p);
            for(u64 k=0;k<r;++k) {
                if(row[k]!=mpz_fdiv_ui(exact[k].get_mpz_t(),static_cast<unsigned long>(p*p*p)))
                    throw std::runtime_error("unit-part coefficient row mismatch");
                ++residue_cases;
                sparse.push_back({k,p,p*p,0});
                if(k<p) ordinary.push_back({k,p,p,0});
            }
        }
        // Repeated primes at the same and different leaves, nonuniform leaf depths,
        // and deliberately tiny forests exercise the precision invariant.
        if(!sparse.empty()) sparse.push_back(sparse.front());
        tree.evaluate(ordinary,false,211); tree.evaluate(sparse,true,307);
        for(auto* list:{&ordinary,&sparse}) for(auto& q:*list) {
            if(q.value!=mpz_fdiv_ui(exact[q.k].get_mpz_t(),static_cast<unsigned long>(q.mod)))
                throw std::runtime_error("batched coefficient mismatch");
            ++residue_cases;
        }
        // Many requests at one index require a second, balanced tree within
        // that leaf. Exercise both a full batch and splits inside that leaf.
        if(r==31) {
            std::vector<Query> shared_ordinary,shared_sparse;
            for(u64 j=0;j<1025;++j) {
                u64 p=(j%3==0) ? 3 : ((j%3==1) ? 5 : 67);
                shared_ordinary.push_back({1,p,p,0});
                shared_sparse.push_back({r-1,p,p*p,0});
            }
            for(u64 quota:{1000000ULL,701ULL}) {
                tree.evaluate(shared_ordinary,false,quota);
                tree.evaluate(shared_sparse,true,quota);
                for(auto* list:{&shared_ordinary,&shared_sparse}) for(auto& q:*list) {
                    if(q.value!=mpz_fdiv_ui(exact[q.k].get_mpz_t(),static_cast<unsigned long>(q.mod)))
                        throw std::runtime_error("shared-leaf remainder mismatch");
                    ++residue_cases;
                }
            }
        }
    }
    std::mt19937_64 rng(20260908);
    for(u64 r:{17ULL,48ULL,113ULL,257ULL}) {
        std::vector<Block> blocks;
        for(int j=0;j<20;++j) blocks.push_back({2*(10000+rng()%300000)+1,1ULL<<(1+rng()%3),128,0,1});
        u64 maximum=1; for(auto& b:blocks) maximum=std::max(maximum,b.A+b.a*r);
        auto primes=primes_to(isqrt(maximum)); Options o; o.small=3; o.base=2; o.audit=true;
        Reconstruction rec(blocks,primes,o,counters,8192); auto result=rec.solve(r);
        for(std::size_t j=0;j<blocks.size();++j) {
            if(!equal(result[j],direct_sum(blocks[j].A,blocks[j].a,r))) throw std::runtime_error("Euler reconstruction mismatch");
            ++reconstruction_cases;
        }
    }
    for(u64 N:{1ULL,2ULL,17ULL,97ULL,1000ULL}) {
        Options o; o.block=17; o.small=3; o.base=2; o.audit=true;
        Counters st; auto core=exact_core(N,192,o,st); mpz_class reference=0;
        u64 n=N-1,K=n+48+2;
        for(auto f:families) for(u64 k=0;k<=K;++k) {
            i64 E=static_cast<i64>(n)+f.c-4*static_cast<i64>(k); if(E>=0) continue;
            mpz_class term=inverse_2exp(integer(f.a*k+f.b),core.precision);
            term<<=bitcount(core.precision+E);
            if(f.sign>0) reference+=term; else reference-=term;
        }
        reduce(reference,core.precision);
        if(reference!=core.value) throw std::runtime_error("complete core mismatch");
        auto words=exact_words(N,192,o,st); mpq_class exact=0;
        for(auto f:families) for(u64 k=0;k<=K;++k) {
            u64 mu=f.a*k+f.b; if(mu==1) continue;
            i64 E=static_cast<i64>(n)+f.c-4*static_cast<i64>(k);
            auto w=joint_power(n,E,mu); mpq_class term(integer(w),integer(mu));
            if(f.sign>0) exact+=term; else exact-=term;
        }
        if(words.n*exact.get_den()!=words.d*exact.get_num()) throw std::runtime_error("exact word aggregation mismatch");
        ++core_cases;
    }
    if(!counters.nonunit_divisions || !counters.power_queries || !counters.prime_queries || !counters.levels)
        throw std::runtime_error("selftest did not exercise all new mechanisms");
    std::cout<<"{\"selftest\":\"PASS\",\"coefficient_checks\":"<<residue_cases
             <<",\"reconstructed_sums\":"<<reconstruction_cases<<",\"exact_core_and_word_sums\":"<<core_cases
             <<",\"nonunit_divisions\":"<<counters.nonunit_divisions<<"}\n";
}
} // namespace lucas

int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--selftest") { lucas::check(); return 0; }
        if(argc<2) throw std::invalid_argument("usage: extract_lucas POSITION [--timing] [--direct] [--block=N] [--small=N] [--base=N] [--quota=N] [--chunk=N] [--guard=N] [--audit]");
        std::size_t used=0; std::string first=argv[1]; u64 N=std::stoull(first,&used);
        if(first.empty() || first[0]=='-' || used!=first.size() || N<1 || N>500000000ULL)
            throw std::invalid_argument("position outside 1..500000000");
        lucas::Options o;
        for(int j=2;j<argc;++j) {
            std::string arg=argv[j];
            if(arg=="--timing") o.timing=true;
            else if(arg=="--direct") o.direct=true;
            else if(arg=="--audit") o.audit=true;
            else {
                auto pos=arg.find('='); if(pos==std::string::npos) throw std::invalid_argument("unknown option");
                auto name=arg.substr(0,pos),value=arg.substr(pos+1); std::size_t consumed=0;
                u64 v=std::stoull(value,&consumed);
                if(value.empty() || value[0]=='-' || consumed!=value.size() || !v || v>500000000ULL)
                    throw std::invalid_argument("invalid option value");
                if(name=="--block") o.block=v; else if(name=="--small") o.small=v;
                else if(name=="--base") o.base=v; else if(name=="--quota") o.quota=v;
                else if(name=="--chunk") o.chunk=v; else if(name=="--guard") o.guard=v;
                else throw std::invalid_argument("unknown option");
            }
        }
        o=lucas::resolve_options(N,o);
        lucas::Counters stats; auto result=lucas::extract(N,o,stats);
        std::cout<<result<<'\n'; if(o.timing) lucas::print(stats,o,N);
        return 0;
    } catch(const std::exception& error) { std::cerr<<"error: "<<error.what()<<'\n'; return 1; }
}
