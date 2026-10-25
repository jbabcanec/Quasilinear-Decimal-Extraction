// Independent exact checks for the prospective fixed-prime real correction.
// The oracle obtains H_r from the original direct PQT and projects its reduced
// denominator onto the selected prime power. No timing assertions are made.
#define main weighted_machin_lucas_program_main
#include "../decimal_extract_weighted_machin_lucas.cpp"
#undef main
#include "machin_lucas_adic.hpp"
#include <set>

static mpz_class adic_check_power(u64 prime,u64 exponent) {
    mpz_class result;
    mpz_ui_pow_ui(result.get_mpz_t(),static_cast<unsigned long>(prime),
                  static_cast<unsigned long>(exponent));
    return result;
}

static u64 adic_check_valuation(const mpz_class& value,u64 prime) {
    if(value==0) throw std::runtime_error("valuation of zero in checker");
    mpz_class unit;
    const mpz_class factor=integer(prime);
    return static_cast<u64>(mpz_remove(unit.get_mpz_t(),value.get_mpz_t(),factor.get_mpz_t()));
}

static Rational adic_check_reduce(Rational value) {
    mpz_class divisor;
    mpz_gcd(divisor.get_mpz_t(),value.n.get_mpz_t(),value.d.get_mpz_t());
    mpz_divexact(value.n.get_mpz_t(),value.n.get_mpz_t(),divisor.get_mpz_t());
    mpz_divexact(value.d.get_mpz_t(),value.d.get_mpz_t(),divisor.get_mpz_t());
    return value;
}

static Rational adic_check_oracle(const Rational& direct,u64 prime) {
    mpz_class unit;
    const mpz_class factor=integer(prime);
    const u64 exponent=static_cast<u64>(
        mpz_remove(unit.get_mpz_t(),direct.d.get_mpz_t(),factor.get_mpz_t()));
    const mpz_class denominator=adic_check_power(prime,exponent);
    if(denominator==1) return {0,1};
    mpz_class inverse,numerator;
    if(!mpz_invert(inverse.get_mpz_t(),unit.get_mpz_t(),denominator.get_mpz_t()))
        throw std::runtime_error("oracle unit was not invertible");
    numerator=direct.n*inverse;
    mpz_fdiv_r(numerator.get_mpz_t(),numerator.get_mpz_t(),denominator.get_mpz_t());
    return {std::move(numerator),denominator};
}

int main() {
    try {
        constexpr u64 selected_prime=13;
        u64 cases=0,integral_cases=0,cancellation_cases=0,high_valuation_cases=0;
        u64 maximum_denominator_valuation=0,maximum_input_denominator_valuation=0;
        const std::vector<u64> counts={1,2,3,5,13,14,27,50,127,256,500};
        for(const u64 machin_prime:{5ULL,239ULL}) {
            const u64 C=machin_prime*machin_prime+1;
            for(const u64 r:counts) {
                std::set<u64> origins={1,3,9,11,13,15,25,101,167,169,171,2000001};
                u64 power=1;
                for(u64 exponent=1;exponent<=6;++exponent) {
                    power*=selected_prime;
                    origins.insert(power);
                    // Force a selected denominator A+2k to equal 13^exponent
                    // at both endpoints and inside the active truncation.
                    for(const u64 k:{1ULL,r/2,r-1})
                        if(k<r && power>2*k) origins.insert(power-2*k);
                }
                for(const u64 A:origins) {
                    const auto raw=machin_real_bs(A,2,C,0,r);
                    const Rational exact=adic_check_reduce({raw.t,raw.q});
                    const auto expected=adic_check_oracle(exact,selected_prime);
                    const auto actual=machin_lucas::fixed_prime_fraction(A,2,C,r,selected_prime);
                    const auto fail=[&](const char* message) {
                        throw std::runtime_error(std::string(message)+" (p0="+
                            std::to_string(machin_prime)+", A="+std::to_string(A)+
                            ", r="+std::to_string(r)+")");
                    };
                    if(actual.d<=0 || actual.n<0 || actual.n>=actual.d)
                        fail("fixed-prime output is outside [0,1)");
                    mpz_class denominator_unit;
                    const mpz_class factor=integer(selected_prime);
                    mpz_remove(denominator_unit.get_mpz_t(),actual.d.get_mpz_t(),factor.get_mpz_t());
                    if(denominator_unit!=1) fail("output denominator is not a pure selected-prime power");
                    if(actual.n*expected.d!=expected.n*actual.d)
                        fail("primary component differs from exact direct-PQT oracle");
                    const mpz_class difference_n=exact.n*actual.d-actual.n*exact.d;
                    const mpz_class difference_d=exact.d*actual.d;
                    if(difference_n!=0 && adic_check_valuation(difference_n,selected_prime)<
                                              adic_check_valuation(difference_d,selected_prime))
                        fail("difference from H_r is not selected-prime integral");
                    const u64 raw_valuation=adic_check_valuation(raw.q,selected_prime);
                    const u64 reduced_valuation=adic_check_valuation(exact.d,selected_prime);
                    if(raw_valuation>reduced_valuation) ++cancellation_cases;
                    if(expected.n==0) ++integral_cases;
                    maximum_denominator_valuation=std::max(maximum_denominator_valuation,reduced_valuation);
                    u64 largest_input=0;
                    for(u64 k=1;k<r;++k) {
                        u64 denominator=A+2*k,valuation=0;
                        while(denominator%selected_prime==0) { denominator/=selected_prime; ++valuation; }
                        largest_input=std::max(largest_input,valuation);
                    }
                    maximum_input_denominator_valuation=std::max(maximum_input_denominator_valuation,largest_input);
                    if(largest_input>=4) ++high_valuation_cases;
                    ++cases;
                }
            }
        }
        if(!integral_cases || !cancellation_cases || !high_valuation_cases ||
           maximum_input_denominator_valuation<6)
            throw std::runtime_error("checker missed a required fixed-prime route");
        std::cout<<"{\"status\":\"PASS\",\"exact_primary_component_cases\":"<<cases
                 <<",\"integral_cases\":"<<integral_cases
                 <<",\"raw_denominator_cancellation_cases\":"<<cancellation_cases
                 <<",\"high_input_valuation_cases\":"<<high_valuation_cases
                 <<",\"maximum_reduced_denominator_valuation\":"<<maximum_denominator_valuation
                 <<",\"maximum_input_denominator_valuation\":"<<maximum_input_denominator_valuation
                 <<",\"fixed_prime\":13,\"machin_primes\":[5,239],\"a\":2,\"maximum_r\":500}\n";
    } catch(const std::exception& error) {
        std::cerr<<"FAIL: "<<error.what()<<'\n';
        return 1;
    }
}
