// Independent exact checks for the bounded-word optimized implementation.
// Build from the repository root with:
// g++ -O2 -std=c++17 validation/weighted_optimized_arithmetic_checks.cpp -lgmpxx -lgmp -o tmp/weighted_optimized_arithmetic_checks
// No benchmark timings are collected by this executable.
#define main weighted_optimized_program_main
#include "../decimal_extract_weighted_optimized.cpp"
#undef main

static void require(bool condition,const char* message) {
    if(!condition) throw std::runtime_error(message);
}
static u64 reference_power(u64 base,u64 exponent,u64 modulus) {
    u64 result=1%modulus;
    base%=modulus;
    while(exponent) {
        if(exponent&1) result=static_cast<u128>(result)*base%modulus;
        base=static_cast<u128>(base)*base%modulus;
        exponent>>=1;
    }
    return result;
}
static mpz_class gmp_inverse(const mpz_class& denominator,const mpz_class& modulus) {
    mpz_class result;
    require(mpz_invert(result.get_mpz_t(),denominator.get_mpz_t(),modulus.get_mpz_t())!=0,
            "GMP reference inverse does not exist");
    return result;
}
static mpz_class reference_words(u64 N,i64 K,i64 limit,u64 guard,u64& count) {
    mpz_class result=0;
    count=0;
    for(const auto f:families) for(i64 k=0;k<=K;++k) {
        const i64 exponent=static_cast<i64>(N-1)+f.c-4*k;
        const u64 modulus=f.a*k+f.b;
        if(modulus==1) continue;
        u64 residue=0;
        if(exponent>=0)
            residue=static_cast<u128>(reference_power(2,exponent,modulus))*reference_power(5,N-1,modulus)%modulus;
        else if(-exponent<limit)
            residue=static_cast<u128>(reference_power(5,N-1,modulus))*reference_power((modulus+1)/2,-exponent,modulus)%modulus;
        if(residue) {
            const mpz_class scaled=integer(residue)<<bitcount(guard);
            const mpz_class term=floor_div(scaled,integer(modulus));
            if(f.sign>0) result+=term; else result-=term;
            ++count;
        }
    }
    reduce(result,guard);
    return result;
}

int main() {
    try {
        std::mt19937_64 rng(20261024);
        unsigned inverses=0,block_cases=0,inverse2_cases=0,word_cases=0,core_cases=0;
        constexpr unsigned arithmetic_cases=20000,joint_cases=2000,radix_cases=20;
        for(unsigned i=0;i<arithmetic_cases;++i) {
            const u64 modulus=i%2 ? UINT32_MAX-2*(rng()%65536) : (3+(rng()%(UINT32_MAX-2)))|1;
            const u64 x=rng()&UINT32_MAX,y=rng()&UINT32_MAX;
            const WordModulus reducer(modulus);
            require(reducer.multiply(x,y)==static_cast<u128>(x)*y%modulus,"reciprocal multiplication mismatch");
            require(mulmod(x,y,modulus)==static_cast<u128>(x)*y%modulus,"bounded multiplication mismatch");
            if(std::gcd(x,modulus)==1) {
                const u64 inverse=invmod(x,modulus);
                require(inverse<modulus && static_cast<u128>(x)*inverse%modulus==1,"uint32-boundary inverse mismatch");
                ++inverses;
            }
        }
        // Explicit endpoints and long Euclid chains complement random samples.
        for(const auto modulus:{3ULL,65537ULL,2147483647ULL,4294967291ULL,4294967295ULL}) {
            for(const auto x:{1ULL,2ULL,123456789ULL,4294967294ULL,4294967295ULL}) {
                if(std::gcd(x,modulus)!=1) continue;
                require(integer(invmod(x,modulus))==gmp_inverse(integer(x),integer(modulus)),"GMP word inverse mismatch");
                ++inverses;
            }
        }
        for(unsigned i=0;i<joint_cases;++i) {
            const u64 modulus=(3+(rng()%(UINT32_MAX-2)))|1;
            // These deliberately exceed 32-bit exponents, with both signs.
            const u64 n=(1ULL<<40)|(rng()&((1ULL<<40)-1));
            const i64 magnitude=static_cast<i64>((1ULL<<36)|(rng()&((1ULL<<36)-1)));
            const i64 exponent=i&1 ? magnitude : -magnitude;
            const u64 two=reference_power(exponent<0 ? (modulus+1)/2 : 2,magnitude,modulus);
            const u64 expected=static_cast<u128>(reference_power(5,n,modulus))*two%modulus;
            require(weighted_words::joint_power(n,exponent,modulus)==expected,"joint signed exponent mismatch");
        }
        for(const u64 precision:{1ULL,2ULL,3ULL,7ULL,31ULL,32ULL,33ULL,63ULL,64ULL,65ULL,127ULL,128ULL,129ULL,191ULL,192ULL,193ULL,255ULL,257ULL,511ULL,1000ULL,1025ULL}) {
            const mpz_class modulus=mpz_class(1)<<bitcount(precision);
            for(unsigned sample=0;sample<8;++sample) {
                mpz_class denominator=integer(rng()|1);
                denominator+=(integer(rng())<<bitcount(precision/2+3));
                if(sample&1) denominator=-denominator;
                require(inverse_2exp(denominator,precision)==gmp_inverse(denominator,modulus),"Newton inverse mismatch");
                ++inverse2_cases;
            }
        }
        for(unsigned i=0;i<radix_cases;++i) {
            std::vector<u64> records;
            for(u64 j=0;j<257;++j) records.push_back(((rng()&UINT32_MAX)<<32)|j);
            auto expected=records;
            std::sort(expected.begin(),expected.end());
            sort_residual_primes(records,UINT32_MAX);
            require(records==expected,"radix ordering mismatch");
        }
        const auto primes=primes_to(isqrt(UINT32_MAX));
        for(unsigned sample=0;sample<32;++sample) {
            const u64 stride=1ULL<<(1+rng()%4);
            const u64 length=sample<16 ? 1+rng()%61 : 128+rng()%129;
            const u64 start=sample%2 ? (UINT32_MAX-stride*(length-1)-2*(rng()%65536)) : 2*(rng()%1000000)+1;
            Stats stats;
            const auto actual=weighted_block(start,stride,length,primes,stats);
            mpq_class expected=0;
            for(u64 j=0;j<length;++j)
                expected+=mpq_class(mpz_class(1)<<bitcount(4*(length-1-j)),integer(start+stride*j));
            require(actual.n*expected.get_den()==actual.d*expected.get_num(),"packed weighted rational mismatch");
            ++block_cases;
        }
        for(const u64 N:{1ULL,2ULL,3ULL,31ULL,32ULL,33ULL,127ULL,257ULL,997ULL}) {
            mpz_class X;
            mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(N-1));
            for(const u64 guard:{1ULL,31ULL,65ULL,192ULL}) {
                const i64 K=cutoff(N,guard),limit=bits(X)+guard+8;
                u64 expected_count=0,actual_count=0;
                const mpz_class expected=reference_words(N,K,limit,guard,expected_count);
                const mpz_class actual=weighted_words::accumulate(N,K,limit,guard,&actual_count);
                require(actual==expected && actual_count==expected_count,"shared word accumulator mismatch");
                ++word_cases;
            }
        }
        for(const u64 N:{1ULL,17ULL,129ULL,257ULL}) {
            mpz_class X;
            mpz_ui_pow_ui(X.get_mpz_t(),5,static_cast<unsigned long>(N-1));
            for(const u64 length:{0ULL,1ULL,7ULL,128ULL,257ULL}) {
                Stats stats;
                const mpz_class actual=block_core(N,192,bits(X),stats,length);
                const mpz_class modulus=mpz_class(1)<<bitcount(stats.Jmax);
                mpz_class expected=0;
                for(const auto f:families) for(i64 k=0;k<=cutoff(N,192);++k) {
                    const i64 exponent=static_cast<i64>(N-1)+f.c-4*k;
                    if(exponent>=0 || -exponent>=static_cast<i64>(bits(X)+200)) continue;
                    mpz_class term=gmp_inverse(integer(f.a*k+f.b),modulus);
                    term<<=bitcount(stats.Jmax+exponent);
                    if(f.sign>0) expected+=term; else expected-=term;
                }
                reduce(expected,stats.Jmax);
                require(actual==expected,"independent exact modular core mismatch");
                ++core_cases;
            }
        }
        std::cout<<"{\"status\":\"PASS\",\"modular_product_cases\":"<<arithmetic_cases
            <<",\"word_inverse_cases\":"<<inverses<<",\"joint_power_cases\":"<<joint_cases
            <<",\"power_of_two_inverse_cases\":"<<inverse2_cases<<",\"radix_cases\":"<<radix_cases
            <<",\"exact_block_cases\":"<<block_cases<<",\"word_accumulator_cases\":"<<word_cases
            <<",\"exact_core_cases\":"<<core_cases<<"}\n";
    } catch(const std::exception& error) {
        std::cerr<<"FAIL: "<<error.what()<<"\n";
        return 1;
    }
}
