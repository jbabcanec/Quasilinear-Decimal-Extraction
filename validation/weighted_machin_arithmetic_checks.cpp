// Independent checks for arithmetic shared by weighted and control Machin.
#define main weighted_machin_program_main
#include "../decimal_extract_weighted_machin.cpp"
#undef main

int main() {
    try {
        std::mt19937_64 rng(20261024);
        unsigned inverse_cases=0,block_cases=0,scale_cases=0;
        for(u64 p:{3ULL,5ULL,239ULL}) {
            for(u64 exponent:{1ULL,2ULL,3ULL,7ULL,31ULL,32ULL,33ULL,63ULL,64ULL,65ULL,127ULL,128ULL,129ULL,257ULL,511ULL,513ULL,1025ULL}) {
                const mpz_class modulus=machin::power(p,exponent);
                for(unsigned sample=0;sample<12;++sample) {
                    mpz_class denominator;
                    if(sample<2) denominator=1;
                    else if(sample<4) denominator=modulus+1;
                    else {
                        denominator=integer(rng());
                        denominator+=(integer(rng())<<bitcount(bits(modulus)+113));
                        if(mpz_divisible_ui_p(denominator.get_mpz_t(),static_cast<unsigned long>(p))) ++denominator;
                    }
                    if(sample&1) denominator=-denominator;
                    mpz_class expected;
                    if(!mpz_invert(expected.get_mpz_t(),denominator.get_mpz_t(),modulus.get_mpz_t()))
                        throw std::runtime_error("reference denominator was not a unit");
                    const mpz_class actual=machin::lift_inverse(denominator,p,exponent);
                    if(actual!=expected || actual<0 || actual>=modulus)
                        throw std::runtime_error("p-adic inverse differs from mpz_invert");
                    ++inverse_cases;
                }
            }
            for(unsigned sample=0;sample<24;++sample) {
                const u64 length=sample<12 ? 1+rng()%32 : 128+rng()%129;
                const u64 A=sample%3==0 ? 1 : sample%3==1 ? UINT32_MAX-2*(length-1)-2*(rng()%10000) : 2*(rng()%100000)+1;
                Stats stats;
                const auto actual=weighted_machin_block(p,A,2,length,primes_to(isqrt(A+2*(length-1))),stats);
                mpq_class expected=0;
                for(u64 j=0;j<length;++j) {
                    mpz_class numerator=prime_power(p,2*(length-1-j));
                    if((length-1-j)&1) numerator=-numerator;
                    expected+=mpq_class(numerator,integer(A+2*j));
                }
                expected*=prime_power(p,actual.valuation);
                if(actual.value.n*expected.get_den()!=actual.value.d*expected.get_num())
                    throw std::runtime_error("signed weighted Machin block differs from exact rational sum");
                ++block_cases;
            }
        }
        // Exercise the scale optimization independently against exact signed
        // multiplication, including zero and multiples of the odd modulus.
        std::mt19937_64 scale_rng(20261025);
        for(u64 p:{3ULL,5ULL,239ULL}) for(u64 exponent:{1ULL,2ULL,31ULL,64ULL,129ULL,257ULL}) {
            const mpz_class modulus=machin::power(p,exponent);
            for(unsigned sample=0;sample<30;++sample) {
                mpz_class value=integer(scale_rng())<<(sample*7);
                value+=integer(scale_rng());
                if(sample&1) value=-value;
                mpz_class Y=integer(scale_rng())<<(sample*31);
                if(sample%5==0) Y=0;
                if(sample%5==1) Y=modulus*integer(scale_rng());
                if(sample&2) Y=-Y;
                const mpz_class expected=machin::modulo(mpz_class(value*Y),modulus);
                const mpz_class actual=machin::scale_residue(value,Y,modulus);
                if(actual!=expected || actual<0 || actual>=modulus)
                    throw std::runtime_error("stripped scale differs from exact signed multiplication");
                ++scale_cases;
            }
        }
        std::cout<<"{\"status\":\"PASS\",\"prime_power_inverse_cases\":"<<inverse_cases
            <<",\"exact_signed_block_cases\":"<<block_cases
            <<",\"signed_scale_cases\":"<<scale_cases
            <<",\"prime_bases\":[3,5,239],\"maximum_inverse_exponent\":1025,\"maximum_denominator\":4294967295}\n";
    } catch(const std::exception& error) {
        std::cerr<<"FAIL: "<<error.what()<<'\n';
        return 1;
    }
}

