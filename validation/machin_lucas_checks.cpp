// Independent exact checks for the prospective shared Machin Lucas engine.
// Build after the prototype and its headers are present; no timing assertions.
#define main weighted_machin_lucas_program_main
#include "../decimal_extract_weighted_machin_lucas.cpp"
#undef main

static u64 reference_mod(const mpz_class& value,u64 modulus) {
    return static_cast<u64>(mpz_fdiv_ui(value.get_mpz_t(),static_cast<unsigned long>(modulus)));
}

static std::vector<mpz_class> exact_coefficient_row(u64 r,u64 machin_prime) {
    const u64 D=machin_prime*machin_prime,C=D+1;
    std::vector<mpz_class> row(r);
    mpz_class power;
    mpz_ui_pow_ui(power.get_mpz_t(),static_cast<unsigned long>(C),static_cast<unsigned long>(r));
    row[0]=(power-1)/integer(D);
    mpz_class choose=1;
    for(u64 k=1;k<r;++k) {
        choose*=integer(r-k+1);
        mpz_divexact_ui(choose.get_mpz_t(),choose.get_mpz_t(),static_cast<unsigned long>(k));
        mpz_class signed_choose=(k&1) ? choose : mpz_class(-choose);
        mpz_class numerator=signed_choose-row[k-1];
        if(!mpz_divisible_ui_p(numerator.get_mpz_t(),static_cast<unsigned long>(D)))
            throw std::runtime_error("reference integer coefficient recurrence is not exact");
        mpz_divexact_ui(row[k].get_mpz_t(),numerator.get_mpz_t(),static_cast<unsigned long>(D));
    }
    return row;
}

static std::vector<u64> query_powers(u64 prime) {
    std::vector<u64> powers;
    u64 value=1;
    for(unsigned exponent=1;exponent<=3 && value<=UINT32_MAX/prime;++exponent) {
        value*=prime; powers.push_back(value);
    }
    return powers;
}

int main() {
    try {
        u64 coefficient_cases=0,exceptional_row_cases=0,ordinary_cases=0,sparse_cases=0;
        u64 shared_leaf_cases=0,weighted_row_cases=0,reconstruction_cases=0;
        u64 forced_recursive_H_cases=0,direct_primary_H_cases=0;
        machin_lucas::Counters tree_stats,reconstruction_stats;
        for(u64 p0:{3ULL,5ULL,239ULL}) {
            for(u64 r:{1ULL,2ULL,7ULL,16ULL,31ULL,65ULL,127ULL}) {
                const auto exact=exact_coefficient_row(r,p0);
                machin_lucas::ValueTree tree(r,p0,tree_stats);
                std::vector<machin_lucas::Query> ordinary,sparse;
                for(u64 query_prime:{3ULL,5ULL,7ULL,13ULL,239ULL,65521ULL}) {
                    for(const u64 modulus:query_powers(query_prime)) {
                        const auto row=machin_lucas::coefficient_row(r,p0,query_prime,modulus);
                        if(row.size()!=r) throw std::runtime_error("coefficient row length mismatch");
                        for(u64 k=0;k<r;++k) {
                            if(row[k]!=reference_mod(exact[k],modulus))
                                throw std::runtime_error("modular coefficient row differs from exact integers");
                            ++coefficient_cases;
                            if(query_prime==p0) ++exceptional_row_cases;
                            // Rows check every index; sample larger tree query
                            // sets while preserving boundary and nonuniform paths.
                            if(query_prime!=p0 && (r<=31 || k==r-1 || k%std::max<u64>(1,r/13)==0)) {
                                sparse.push_back({k,query_prime,modulus,0});
                                if(k<query_prime) ordinary.push_back({k,query_prime,modulus,0});
                            }
                        }
                    }
                }
                // Repeated prime powers and repeated indices must remain exact
                // when one leaf spans multiple bounded forests.
                if(!ordinary.empty()) ordinary.push_back(ordinary.front());
                if(!sparse.empty()) sparse.push_back(sparse.back());
                for(u64 quota:{4096ULL,307ULL}) {
                    tree.evaluate(ordinary,false,quota);
                    tree.evaluate(sparse,true,quota);
                    for(const auto& query:ordinary) {
                        if(query.value!=reference_mod(exact[query.k],query.mod))
                            throw std::runtime_error("ordinary coefficient query mismatch");
                        ++ordinary_cases;
                    }
                    for(const auto& query:sparse) {
                        if(query.value!=reference_mod(exact[query.k],query.mod))
                            throw std::runtime_error("sparse coefficient query mismatch");
                        ++sparse_cases;
                    }
                }
                if(r==31) {
                    std::vector<machin_lucas::Query> shared_ordinary,shared_sparse;
                    for(u64 j=0;j<257;++j) {
                        u64 q=j%3==0 ? 3 : j%3==1 ? 5 : 67;
                        if(q==p0) q=13;
                        const auto powers=query_powers(q);
                        const u64 modulus=powers[j%powers.size()];
                        shared_ordinary.push_back({1,q,modulus,0});
                        shared_sparse.push_back({r-1,q,modulus,0});
                    }
                    for(u64 quota:{1000000ULL,71ULL}) {
                        tree.evaluate(shared_ordinary,false,quota);
                        tree.evaluate(shared_sparse,true,quota);
                        for(const auto* list:{&shared_ordinary,&shared_sparse})
                            for(const auto& query:*list) {
                                if(query.value!=reference_mod(exact[query.k],query.mod))
                                    throw std::runtime_error("shared-leaf coefficient query mismatch");
                                ++shared_leaf_cases;
                            }
                    }
                }
                for(u64 q:{3ULL,5ULL,13ULL,239ULL}) for(u64 weight:{0ULL,1ULL,p0*p0,q-1}) {
                    const auto row=machin_lucas::weighted_row(r,weight,q);
                    if(row.weight.size()!=r+1 || row.prefix.size()!=r+2)
                        throw std::runtime_error("weighted prefix row length mismatch");
                    mpz_class choose=1,power=1,total=0;
                    if(row.prefix[0]!=0) throw std::runtime_error("weighted prefix origin mismatch");
                    for(u64 k=0;k<=r;++k) {
                        const mpz_class term=choose*power;
                        total+=term;
                        if(row.weight[k]!=reference_mod(term,q) || row.prefix[k+1]!=reference_mod(total,q))
                            throw std::runtime_error("Lucas weighted prefix differs from exact binomial sum");
                        ++weighted_row_cases;
                        if(k<r) {
                            choose*=integer(r-k);
                            mpz_divexact_ui(choose.get_mpz_t(),choose.get_mpz_t(),static_cast<unsigned long>(k+1));
                            power*=integer(weight);
                        }
                    }
                }
            }
            for(u64 r:{1ULL,2ULL,7ULL,17ULL,48ULL,113ULL,257ULL}) {
                // The first family covers small A and direct-fallback boundaries.
                // The second forces contraction and includes a denominator p0^2
                // at k=1, ensuring the backwards fixed-prime row is required.
                const u64 square=p0*p0;
                u64 multiple=(2000001+square-1)/square;
                if(!(multiple&1)) ++multiple;
                const u64 exceptional_A=multiple*square-2;
                const std::vector<std::vector<machin_lucas::Block>> families={
                    {{1,2},{3,4},{p0,2},{p0*p0-2,2},{101,8}},
                    {{exceptional_A,2},{2000001,4},{10000001,8},{4294960001ULL,16}}
                };
                for(const auto& blocks:families) for(const bool primary:{false,true}) {
                    machin_lucas::Reconstruction engine(p0,blocks,2,3,307,reconstruction_stats,true,primary);
                    const auto actual=engine.solve(r);
                    if(actual.size()!=blocks.size()) throw std::runtime_error("reconstruction batch length mismatch");
                    for(std::size_t b=0;b<blocks.size();++b) {
                        const auto expected=machin_real_bs(blocks[b].A,blocks[b].a,p0*p0+1,0,r);
                        if(actual[b].d<=0 || actual[b].n*expected.q!=actual[b].d*expected.t)
                            throw std::runtime_error("recursive Machin H_r differs from independent direct PQT");
                        ++reconstruction_cases;
                        if(primary && p0!=3) ++direct_primary_H_cases;
                        else ++forced_recursive_H_cases;
                    }
                }
            }
        }
        if(!exceptional_row_cases || !ordinary_cases || !sparse_cases || !tree_stats.nonunit_divisions ||
           !tree_stats.ordinary_batches || !tree_stats.power_batches || !reconstruction_stats.levels)
            throw std::runtime_error("checks did not exercise every required Lucas route");
        std::cout<<"{\"status\":\"PASS\",\"exact_coefficient_cases\":"<<coefficient_cases
                 <<",\"backwards_fixed_prime_row_cases\":"<<exceptional_row_cases
                 <<",\"ordinary_query_cases\":"<<ordinary_cases<<",\"sparse_query_cases\":"<<sparse_cases
                 <<",\"shared_leaf_query_cases\":"<<shared_leaf_cases
                 <<",\"exact_weighted_prefix_cases\":"<<weighted_row_cases
                 <<",\"exact_reconstructed_H_cases\":"<<reconstruction_cases
                 <<",\"forced_recursive_H_cases\":"<<forced_recursive_H_cases
                 <<",\"direct_primary_H_cases\":"<<direct_primary_H_cases
                 <<",\"recursive_levels\":"<<reconstruction_stats.levels
                 <<",\"ordinary_batches\":"<<tree_stats.ordinary_batches
                 <<",\"sparse_batches\":"<<tree_stats.power_batches
                 <<",\"nonunit_divisions\":"<<tree_stats.nonunit_divisions
                 <<",\"prime_bases\":[3,5,239],\"maximum_r\":257}\n";
    } catch(const std::exception& error) {
        std::cerr<<"FAIL: "<<error.what()<<'\n'; return 1;
    }
}
