#pragma once
#include "machin_lucas_common.hpp"
#include <array>
#include <map>

namespace machin_lucas {

// Binomial residues retain their unit part across zeros. Denominator units
// are inverted in bounded native batches, not once per coefficient.
inline std::vector<u64> binomial_row(u64 n,u64 p,u64 q) {
    if(p<2 || q<p || q>UINT32_MAX || n>500000000ULL)
        throw std::invalid_argument("binomial row outside audited word range");
    u64 temporary=q,exponent=0;
    while(temporary%p==0) { temporary/=p; ++exponent; }
    if(temporary!=1) throw std::invalid_argument("binomial modulus is not a prime power");
    std::vector<u64> pp(exponent+1,1);
    for(u64 i=1;i<=exponent;++i) pp[i]=pp[i-1]*p;
    std::vector<u64> row(n+1); row[0]=1;
    i64 valuation=0; u64 unit=1;
    constexpr std::size_t batch=256;
    std::array<u64,batch> denominator,inverses,removed;
    std::array<u64,batch+1> prefix;
    for(u64 begin=1;begin<=n;) {
        const auto count=static_cast<std::size_t>(std::min<u64>(batch,n-begin+1));
        prefix[0]=1;
        for(std::size_t i=0;i<count;++i) {
            u64 d=begin+i,v=0;
            while(d%p==0) { d/=p; ++v; }
            denominator[i]=d%q; removed[i]=v;
            prefix[i+1]=mulmod(prefix[i],denominator[i],q);
        }
        u64 inverse=invmod(prefix[count],q);
        for(std::size_t i=count;i--;) {
            inverses[i]=mulmod(inverse,prefix[i],q);
            inverse=mulmod(inverse,denominator[i],q);
        }
        for(std::size_t i=0;i<count;++i) {
            const u64 j=begin+i;
            u64 numerator=n-j+1;
            while(numerator%p==0) { numerator/=p; ++valuation; }
            valuation-=static_cast<i64>(removed[i]);
            if(valuation<0) throw std::runtime_error("negative binomial valuation");
            unit=mulmod(mulmod(unit,numerator%q,q),inverses[i],q);
            row[j]=static_cast<u64>(valuation)>=exponent ? 0 : mulmod(unit,pp[valuation],q);
        }
        begin+=count;
    }
    return row;
}

inline std::vector<u64> coefficient_row(u64 r,u64 machin_prime,u64 p,u64 q) {
    if(!r || machin_prime<2) throw std::invalid_argument("empty Machin coefficient row");
    auto choose=binomial_row(r,p,q);
    std::vector<u64> row(r);
    const u64 D=mulmod(machin_prime%q,machin_prime%q,q);
    if(p==machin_prime) {
        row[r-1]=(r&1) ? 1 : q-1;
        for(u64 k=r-1;k>0;--k) {
            const u64 b=(k&1) || !choose[k] ? choose[k] : q-choose[k];
            row[k-1]=(b+q-mulmod(D,row[k],q))%q;
        }
    } else {
        const u64 inverse=invmod(D,q);
        row[0]=mulmod((powmod((D+1)%q,r,q)+q-1)%q,inverse,q);
        for(u64 k=1;k<r;++k) {
            const u64 b=(k&1) || !choose[k] ? choose[k] : q-choose[k];
            row[k]=mulmod((b+q-row[k-1])%q,inverse,q);
        }
    }
    return row;
}

struct WeightedRow { std::vector<u64> weight,prefix; };
inline WeightedRow weighted_row(u64 n,u64 weight,u64 p) {
    const auto choose=binomial_row(n,p,p);
    WeightedRow row; row.weight.resize(n+1); row.prefix.resize(n+2);
    u64 power=1; weight%=p;
    for(u64 k=0;k<=n;++k) {
        row.weight[k]=mulmod(choose[k],power,p);
        row.prefix[k+1]=(row.prefix[k]+row.weight[k])%p;
        power=mulmod(power,weight,p);
    }
    return row;
}

// Integer numerator (a,b;0,c) and positive denominator of a transfer.
struct Transfer { mpz_class a=1,b=0,c=1,d=1; };

class TransferNormalizers {
    std::vector<u64> primes;
    std::map<std::pair<u64,u64>,mpz_class> cache;
public:
    explicit TransferNormalizers(u64 r) {
        std::vector<bool> sieve(r+1,true);
        sieve[0]=false; if(r) sieve[1]=false;
        for(u64 p=2;p<=r/p;++p) if(sieve[p])
            for(u64 j=p*p;j<=r;j+=p) sieve[j]=false;
        for(u64 p=2;p<=r;++p) if(sieve[p]) primes.push_back(p);
    }
    const mpz_class& get(u64 left,u64 right) {
        const auto key=std::make_pair(left,right);
        auto found=cache.find(key); if(found!=cache.end()) return found->second;
        const u64 s=left+right;
        std::vector<mpz_class> factors;
        u128 packed=1; unsigned packed_bits=0;
        for(u64 p:primes) {
            if(p>s) break;
            i64 exponent=0;
            for(u64 power=p;power<=s;) {
                exponent+=static_cast<i64>(s/power-left/power-right/power)
                         -1+static_cast<i64>(power<=left)+static_cast<i64>(power<=right);
                if(power>s/p) break;
                power*=p;
            }
            if(exponent<0) throw std::runtime_error("nonintegral transfer normalizer");
            if(!exponent) continue;
            u64 factor=1;
            for(i64 j=0;j<exponent;++j) {
                if(factor>UINT64_MAX/p) throw std::overflow_error("transfer prime-power normalizer range");
                factor*=p;
            }
            const unsigned width=bits(factor);
            if(packed_bits+width>128) {
                factors.push_back(integer128(packed)); packed=1; packed_bits=0;
            }
            packed*=factor; packed_bits+=width;
        }
        if(packed_bits) factors.push_back(integer128(packed));
        while(factors.size()>1) {
            const auto count=factors.size(); std::size_t out=0;
            for(std::size_t j=0;j+1<count;j+=2) {
                mpz_class product=factors[j]*factors[j+1];
                factors[out++]=std::move(product);
            }
            if(count&1) factors[out++]=std::move(factors.back());
            factors.resize(out);
        }
        return cache.emplace(key,factors.empty() ? mpz_class(1) : std::move(factors.front())).first->second;
    }
};

class ValueTree {
    struct Node {
        u64 lo=0,hi=0,depth=0;
        int parent=-1,left=-1,right=-1;
        bool active=false;
        Transfer t;
        mpz_class q=1,s=1,m=1;
        std::vector<std::size_t> queries;
        std::vector<mpz_class> query_products;
    };
    std::vector<Node> nodes;
    std::vector<int> leaves,active;
    u64 r,machin_prime,D,height=0;
    mpz_class e0,b0=-1;
    Counters& stats;
    TransferNormalizers normalizers;

    int build(u64 lo,u64 hi,u64 depth,int parent) {
        const int index=static_cast<int>(nodes.size()); nodes.emplace_back();
        nodes[index].lo=lo; nodes[index].hi=hi; nodes[index].depth=depth; nodes[index].parent=parent;
        height=std::max(height,depth);
        if(hi-lo==1) {
            const u64 k=lo+1,h=r-k+1;
            nodes[index].t={-integer(k),-integer(h),-integer(D)*integer(h),integer(D)*integer(k)};
            leaves[lo]=index;
        } else {
            const u64 mid=lo+(hi-lo)/2;
            const int left=build(lo,mid,depth+1,index),right=build(mid,hi,depth+1,index);
            nodes[index].left=left; nodes[index].right=right;
            const auto& a=nodes[right].t; const auto& b=nodes[left].t;
            Transfer out{a.a*b.a,a.a*b.b+a.b*b.c,a.c*b.c,a.d*b.d};
            const auto& factor=normalizers.get(mid-lo,hi-mid);
            if(factor!=1) {
                mpz_divexact(out.a.get_mpz_t(),out.a.get_mpz_t(),factor.get_mpz_t());
                mpz_divexact(out.b.get_mpz_t(),out.b.get_mpz_t(),factor.get_mpz_t());
                mpz_divexact(out.c.get_mpz_t(),out.c.get_mpz_t(),factor.get_mpz_t());
                mpz_divexact(out.d.get_mpz_t(),out.d.get_mpz_t(),factor.get_mpz_t());
            }
            nodes[index].t=std::move(out);
        }
        return index;
    }
    u64 lambda(u64 p) const {
        u64 t=r,count=0;
        while(t>=p) { t/=p; ++count; }
        return 2*count;
    }
    void mark(int index) {
        while(index>=0 && !nodes[index].active) {
            nodes[index].active=true; active.push_back(index); index=nodes[index].parent;
        }
    }
    void clear_batch() {
        for(int index:active) {
            auto& v=nodes[index]; v.active=false;
            mpz_class(1).swap(v.q); mpz_class(1).swap(v.s); mpz_class(1).swap(v.m);
            std::vector<std::size_t>().swap(v.queries);
            std::vector<mpz_class>().swap(v.query_products);
        }
        active.clear();
    }
    void modulus_tree(int index,bool sparse,const std::vector<Query>& queries) {
        auto& v=nodes[index]; if(!v.active) return;
        if(v.left>=0) {
            modulus_tree(v.left,sparse,queries); modulus_tree(v.right,sparse,queries);
            v.q=nodes[v.left].q*nodes[v.right].q;
            v.s=nodes[v.left].s*nodes[v.right].s;
        } else {
            const auto count=v.queries.size();
            v.query_products.assign(2*count,mpz_class(1));
            std::vector<mpz_class> precision;
            if(sparse) precision.assign(2*count,mpz_class(1));
            for(std::size_t j=0;j<count;++j) {
                const auto& query=queries[v.queries[j]];
                v.query_products[count+j]=integer(query.mod);
                if(sparse) precision[count+j]=powerz(query.p,lambda(query.p));
            }
            for(std::size_t j=count;j-->1;) {
                v.query_products[j]=v.query_products[2*j]*v.query_products[2*j+1];
                if(sparse) precision[j]=precision[2*j]*precision[2*j+1];
            }
            v.q=v.query_products[1];
            if(sparse) v.s=std::move(precision[1]);
        }
        if(sparse) {
            mpz_class padding;
            mpz_pow_ui(padding.get_mpz_t(),v.s.get_mpz_t(),static_cast<unsigned long>(height-v.depth));
            v.m=v.q*padding;
        } else v.m=v.q;
    }
    void distribute_leaf(const Node& leaf,std::size_t index,mpz_class value,std::vector<Query>& queries) {
        value=modz(std::move(value),leaf.query_products[index]);
        const auto count=leaf.queries.size();
        if(index>=count) queries[leaf.queries[index-count]].value=mpz_get_ui(value.get_mpz_t());
        else {
            distribute_leaf(leaf,2*index,value,queries);
            distribute_leaf(leaf,2*index+1,std::move(value),queries);
        }
    }
    void plain_visit(int index,mpz_class x,mpz_class y,mpz_class d,std::vector<Query>& queries) {
        auto& v=nodes[index]; if(!v.active) return;
        x=modz(std::move(x),v.m); y=modz(std::move(y),v.m); d=modz(std::move(d),v.m);
        if(v.left<0) {
            distribute_leaf(v,1,modz(mpz_class(x*invz(d,v.m)),v.m),queries); return;
        }
        plain_visit(v.left,x,y,d,queries);
        const auto& target=nodes[v.right]; if(!target.active) return;
        const auto& t=nodes[v.left].t;
        x=modz(std::move(x),target.m); y=modz(std::move(y),target.m); d=modz(std::move(d),target.m);
        mpz_class nx=modz(mpz_class(modz(t.a,target.m)*x+modz(t.b,target.m)*y),target.m);
        mpz_class ny=modz(mpz_class(modz(t.c,target.m)*y),target.m);
        mpz_class nd=modz(mpz_class(modz(t.d,target.m)*d),target.m);
        plain_visit(v.right,std::move(nx),std::move(ny),std::move(nd),queries);
    }
    void sparse_visit(int index,mpz_class x,mpz_class y,std::vector<Query>& queries) {
        auto& v=nodes[index]; if(!v.active) return;
        x=modz(std::move(x),v.m); y=modz(std::move(y),v.m);
        if(v.left<0) { distribute_leaf(v,1,std::move(x),queries); return; }
        sparse_visit(v.left,x,y,queries);
        const auto& target=nodes[v.right]; if(!target.active) return;
        const auto& t=nodes[v.left].t;
        const mpz_class source=target.m*target.s;
        x=modz(std::move(x),source); y=modz(std::move(y),source);
        mpz_class wx=modz(mpz_class(modz(t.a,source)*x+modz(t.b,source)*y),source);
        mpz_class wy=modz(mpz_class(modz(t.c,source)*y),source);
        // All gcd/inverse operands are bounded by the query modulus. Since
        // target.m divides source/g, (Q mod source)/g is Q/g mod target.m.
        mpz_class denominator=modz(t.d,source),g=gcdz(denominator,source);
        if(g!=1) ++stats.nonunit_divisions;
        if(!mpz_divisible_p(wx.get_mpz_t(),g.get_mpz_t()) ||
           !mpz_divisible_p(wy.get_mpz_t(),g.get_mpz_t()))
            throw std::runtime_error("Machin sparse exact-division invariant failed");
        mpz_class remaining=source/g;
        if(!mpz_divisible_p(remaining.get_mpz_t(),target.m.get_mpz_t()))
            throw std::runtime_error("Machin sparse padded precision exhausted");
        mpz_divexact(denominator.get_mpz_t(),denominator.get_mpz_t(),g.get_mpz_t());
        const mpz_class inverse=invz(modz(std::move(denominator),target.m),target.m);
        mpz_divexact(wx.get_mpz_t(),wx.get_mpz_t(),g.get_mpz_t());
        mpz_divexact(wy.get_mpz_t(),wy.get_mpz_t(),g.get_mpz_t());
        wx=modz(mpz_class(wx*inverse),target.m); wy=modz(mpz_class(wy*inverse),target.m);
        sparse_visit(v.right,std::move(wx),std::move(wy),queries);
    }
public:
    ValueTree(u64 length,u64 prime,Counters& counters):
        leaves(length),r(length),machin_prime(prime),D(prime*prime),stats(counters),normalizers(length) {
        if(!r || r>500000000ULL || prime<2 || prime>65535)
            throw std::invalid_argument("Machin value tree outside audited range");
        nodes.reserve(2*r); build(0,r,0,-1);
        e0=(powerz(D+1,r)-1)/integer(D);
    }
    void evaluate(std::vector<Query>& queries,bool sparse,u64 quota) {
        if(queries.empty()) return;
        quota=std::max<u64>(1,quota);
        std::vector<std::size_t> order(queries.size()),destination(queries.size());
        std::iota(order.begin(),order.end(),0);
        for(const auto& q:queries) {
            if(q.k>=r || q.p<2 || q.p==machin_prime || q.mod<q.p || q.mod>UINT32_MAX || (!sparse && q.k>=q.p))
                throw std::invalid_argument("invalid Machin shared coefficient query");
            u64 remaining=q.mod;
            while(remaining%q.p==0) remaining/=q.p;
            if(remaining!=1) throw std::invalid_argument("coefficient query modulus is not a prime power");
        }
        std::sort(order.begin(),order.end(),[&](auto i,auto j) {
            return queries[i].k!=queries[j].k ? queries[i].k<queries[j].k :
                   queries[i].p!=queries[j].p ? queries[i].p<queries[j].p : queries[i].mod>queries[j].mod;
        });
        // Repeated requests for one coefficient and prime share their largest
        // requested precision. Fan-out afterward is a native word remainder.
        std::vector<Query> unique; unique.reserve(queries.size());
        for(auto index:order) {
            const auto& query=queries[index];
            if(unique.empty() || query.k!=unique.back().k || query.p!=unique.back().p)
                unique.push_back(query);
            destination[index]=unique.size()-1;
        }
        for(std::size_t begin=0;begin<unique.size();) {
            clear_batch();
            std::size_t end=begin; u64 mass=0;
            while(end<unique.size()) {
                const auto& query=unique[end];
                const u64 cost=bits(query.mod)+(sparse ? height*lambda(query.p)*bits(query.p) : 0);
                if(end>begin && mass+cost>quota) break;
                const auto leaf=leaves[query.k]; mark(leaf);
                nodes[leaf].queries.push_back(end); mass+=cost; ++end;
            }
            modulus_tree(0,sparse,unique);
            stats.max_root_bits=std::max(stats.max_root_bits,bits(nodes[0].m));
            if(sparse) { ++stats.power_batches; sparse_visit(0,e0,b0,unique); }
            else { ++stats.ordinary_batches; plain_visit(0,e0,b0,mpz_class(1),unique); }
            begin=end;
        }
        for(std::size_t i=0;i<queries.size();++i) queries[i].value=unique[destination[i]].value%queries[i].mod;
        clear_batch();
    }
};

// Finite-size backend: stream one exact coefficient and signed binomial at a
// time. Select this only below a fixed r cutoff; larger r retains ValueTree's
// shared normalized transfer construction. Storage does not retain the row.
inline void evaluate_streamed(u64 r,u64 machin_prime,std::vector<Query>& queries,
                              u64 quota,Counters& stats) {
    if(queries.empty()) return;
    if(!r || r>500000000ULL || machin_prime<2 || machin_prime>65535)
        throw std::invalid_argument("streamed Machin coefficients outside audited range");
    quota=std::max<u64>(1,quota);
    const u64 D=machin_prime*machin_prime;
    std::vector<std::size_t> order(queries.size()),destination(queries.size());
    std::iota(order.begin(),order.end(),0);
    for(const auto& query:queries) {
        if(query.k>=r || query.p<2 || query.mod<query.p || query.mod>UINT32_MAX)
            throw std::invalid_argument("invalid streamed coefficient query");
        u64 remaining=query.mod;
        while(remaining%query.p==0) remaining/=query.p;
        if(remaining!=1) throw std::invalid_argument("streamed query modulus is not a prime power");
    }
    std::sort(order.begin(),order.end(),[&](auto i,auto j) {
        return queries[i].k!=queries[j].k ? queries[i].k<queries[j].k :
               queries[i].p!=queries[j].p ? queries[i].p<queries[j].p : queries[i].mod>queries[j].mod;
    });
    std::vector<Query> unique; unique.reserve(queries.size());
    for(auto index:order) {
        const auto& query=queries[index];
        if(unique.empty() || query.k!=unique.back().k || query.p!=unique.back().p) unique.push_back(query);
        destination[index]=unique.size()-1;
    }
    mpz_class coefficient=(powerz(D+1,r)-1)/integer(D),binomial=-1;
    u64 current=0;
    for(std::size_t begin=0;begin<unique.size();) {
        const u64 target=unique[begin].k;
        while(current<target) {
            ++current;
            mpz_mul_ui(binomial.get_mpz_t(),binomial.get_mpz_t(),static_cast<unsigned long>(r-current+1));
            mpz_neg(binomial.get_mpz_t(),binomial.get_mpz_t());
            mpz_divexact_ui(binomial.get_mpz_t(),binomial.get_mpz_t(),static_cast<unsigned long>(current));
            mpz_sub(coefficient.get_mpz_t(),binomial.get_mpz_t(),coefficient.get_mpz_t());
            mpz_divexact_ui(coefficient.get_mpz_t(),coefficient.get_mpz_t(),static_cast<unsigned long>(D));
        }
        std::size_t group_end=begin+1;
        while(group_end<unique.size() && unique[group_end].k==target) ++group_end;
        while(begin<group_end) {
            std::size_t end=begin; u64 mass=0;
            while(end<group_end) {
                const auto cost=bits(unique[end].mod);
                if(end>begin && mass+cost>quota) break;
                mass+=cost; ++end;
            }
            const auto count=end-begin;
            std::vector<mpz_class> products(2*count,mpz_class(1));
            for(std::size_t j=0;j<count;++j) products[count+j]=integer(unique[begin+j].mod);
            for(std::size_t j=count;j-->1;) products[j]=products[2*j]*products[2*j+1];
            stats.max_root_bits=std::max(stats.max_root_bits,bits(products[1]));
            ++stats.ordinary_batches;
            auto distribute=[&](auto&& self,std::size_t node,mpz_class value)->void {
                value=modz(std::move(value),products[node]);
                if(node>=count) unique[begin+node-count].value=mpz_get_ui(value.get_mpz_t());
                else {
                    self(self,2*node,value);
                    self(self,2*node+1,std::move(value));
                }
            };
            distribute(distribute,1,coefficient);
            begin=end;
        }
    }
    for(std::size_t i=0;i<queries.size();++i) queries[i].value=unique[destination[i]].value%queries[i].mod;
}

} // namespace machin_lucas
