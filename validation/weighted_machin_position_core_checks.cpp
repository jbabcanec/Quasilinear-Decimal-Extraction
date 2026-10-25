// Compare both complete prime-power residues, not only their decimal windows.
// Build with -O3 -march=native -std=c++17 -lgmpxx -lgmp.
#define main weighted_machin_program_main
#include "../decimal_extract_weighted_machin.cpp"
#undef main

int main() {
    try {
        unsigned checks=0;
        u64 maximum_bits=0;
        for(u64 N:{1ULL,2ULL,97ULL,997ULL,10000ULL,100000ULL,1000000ULL}) {
            const u64 n=N-1;
            const auto ten=machin::power(10,n);
            for(u64 p:{5ULL,239ULL}) {
                const u64 coefficient=p==5 ? 16 : 4;
                const u64 cancelled=p==5 ? n : 0;
                const int sign=p==5 ? 1 : -1;
                const auto K=machin::cutoff(p,coefficient,n,192,ten);
                const mpz_class Y=p==5 ? mpz_class(mpz_class(1)<<(n+4)) : mpz_class(4*ten);
                machin::Stats profile;
                const auto got=machin::prime_core(p,K,cancelled,Y,sign,&profile);
                const auto want=machin::reference_prime_core(p,K,cancelled,Y,sign);
                if(got.exponent!=want.exponent || got.modulus!=want.modulus ||
                   got.residue!=want.residue || got.terms!=want.terms)
                    throw std::runtime_error("complete position core mismatch");
                maximum_bits=std::max(maximum_bits,bits(got.modulus));
                ++checks;
            }
        }
        std::cout<<"{\"status\":\"PASS\",\"complete_position_cores\":"<<checks
                 <<",\"maximum_position\":1000000,\"largest_modulus_bits\":"<<maximum_bits<<"}\n";
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n'; return 1;
    }
}
