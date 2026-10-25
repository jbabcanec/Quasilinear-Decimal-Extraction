// Full-prefix comparison control: compute pi to N+8 decimal places,
// form the scaled integer prefix, and retain its final nine digits.
// This is a fixed-guard numerical calculation, not an interval certificate.
#include <gmp.h>
#include <mpfr.h>
#include <cstdio>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: mpfr_prefix POSITION");
        unsigned long n = std::stoul(argv[1]);
        if (n < 1 || n > 1000000) throw std::runtime_error("position outside 1..1000000");
        mpz_t scale, prefix;
        mpz_inits(scale, prefix, nullptr);
        mpz_ui_pow_ui(scale, 10, n + 8);
        mpfr_t pi;
        mpfr_init2(pi, mpz_sizeinbase(scale, 2) + 192);
        mpfr_const_pi(pi, MPFR_RNDN);
        mpfr_mul_z(pi, pi, scale, MPFR_RNDN);
        mpfr_get_z(prefix, pi, MPFR_RNDD);
        std::printf("%09lu\n", mpz_fdiv_ui(prefix, 1000000000UL));
        mpfr_clear(pi);
        mpz_clears(scale, prefix, nullptr);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
}
