/**
 * High-Performance Decimal Digit Extraction for π via Machin's Formula
 * =====================================================================
 *
 * This implementation uses:
 * - GMP/MPFR for arbitrary precision arithmetic
 * - OpenMP for parallel computation of P₁ and Q₁ sums
 * - Modular exponentiation for O(log N) per-term computation
 *
 * Compile (macOS with Homebrew):
 *   clang++ -std=c++17 -O3 -Xpreprocessor -fopenmp \
 *     -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/include \
 *     -L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/lib \
 *     -o machin_digit_extract machin_digit_extract.cpp -lmpfr -lgmp -lomp
 *
 * Compile (Linux):
 *   g++ -std=c++17 -O3 -fopenmp -o machin_digit_extract machin_digit_extract.cpp -lmpfr -lgmp
 *
 * Author: Joseph Babcanec (Benedict College)
 * Date: 2026
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>
#include <string>
#include <algorithm>
#include <mpfr.h>
#include <gmp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// Modular exponentiation: compute base^exp mod mod
inline unsigned long mod_pow(unsigned long base, unsigned long exp, unsigned long mod) {
    if (mod == 1) return 0;
    if (exp == 0) return 1;

    unsigned long result = 1;
    base = base % mod;

    while (exp > 0) {
        if (exp & 1) {
            result = (__uint128_t)result * base % mod;
        }
        exp >>= 1;
        base = (__uint128_t)base * base % mod;
    }
    return result;
}

class MachineDigitExtractor {
private:
    mpfr_prec_t precision;
    int guard_digits;

    // Compute P₁ = Σ 1/((4j+1)·625^j) directly
    void compute_P1_direct(mpfr_t result) {
        mpfr_set_zero(result, 1);
        mpfr_t term, power;
        mpfr_init2(term, precision + 100);
        mpfr_init2(power, precision + 100);

        mpfr_set_ui(power, 1, MPFR_RNDN);  // 625^0 = 1

        for (unsigned long j = 0; j < 200; j++) {
            unsigned long denom = 4 * j + 1;
            mpfr_ui_div(term, 1, power, MPFR_RNDN);
            mpfr_div_ui(term, term, denom, MPFR_RNDN);
            mpfr_add(result, result, term, MPFR_RNDN);

            mpfr_mul_ui(power, power, 625, MPFR_RNDN);

            if (mpfr_cmp_d(term, pow(10.0, -(double)(precision + 10))) < 0) break;
        }

        mpfr_clear(term);
        mpfr_clear(power);
    }

    // Compute Q₁ = Σ 1/((4j+3)·625^j) directly
    void compute_Q1_direct(mpfr_t result) {
        mpfr_set_zero(result, 1);
        mpfr_t term, power;
        mpfr_init2(term, precision + 100);
        mpfr_init2(power, precision + 100);

        mpfr_set_ui(power, 1, MPFR_RNDN);

        for (unsigned long j = 0; j < 200; j++) {
            unsigned long denom = 4 * j + 3;
            mpfr_ui_div(term, 1, power, MPFR_RNDN);
            mpfr_div_ui(term, term, denom, MPFR_RNDN);
            mpfr_add(result, result, term, MPFR_RNDN);

            mpfr_mul_ui(power, power, 625, MPFR_RNDN);

            if (mpfr_cmp_d(term, pow(10.0, -(double)(precision + 10))) < 0) break;
        }

        mpfr_clear(term);
        mpfr_clear(power);
    }

    // Compute {5^N · P₁} where P₁ = Σ 1/((4j+1)·625^j)
    void compute_frac_P1(mpfr_t result, long N) {
        mpfr_set_zero(result, 1);

        unsigned long max_j = (N > 0 ? N / 4 : 0) + precision / 2 + 50;

        int num_threads = 1;
        #ifdef _OPENMP
        num_threads = omp_get_max_threads();
        #endif

        mpfr_t* partial_sums = new mpfr_t[num_threads];
        for (int i = 0; i < num_threads; i++) {
            mpfr_init2(partial_sums[i], precision + 100);
            mpfr_set_zero(partial_sums[i], 1);
        }

        #pragma omp parallel
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif

            mpfr_t term;
            mpfr_init2(term, precision + 100);

            #pragma omp for schedule(dynamic, 100)
            for (unsigned long j = 0; j <= max_j; j++) {
                unsigned long denom = 4 * j + 1;
                long exp = N - 4 * (long)j;

                if (exp >= 0) {
                    unsigned long numerator = mod_pow(5, (unsigned long)exp, denom);
                    mpfr_set_ui(term, numerator, MPFR_RNDN);
                    mpfr_div_ui(term, term, denom, MPFR_RNDN);
                } else {
                    mpfr_set_ui(term, 5, MPFR_RNDN);
                    mpfr_pow_ui(term, term, (unsigned long)(-exp), MPFR_RNDN);
                    mpfr_mul_ui(term, term, denom, MPFR_RNDN);
                    mpfr_ui_div(term, 1, term, MPFR_RNDN);
                }

                mpfr_add(partial_sums[tid], partial_sums[tid], term, MPFR_RNDN);

                if (mpfr_cmp_ui(partial_sums[tid], 1000) > 0) {
                    mpfr_frac(partial_sums[tid], partial_sums[tid], MPFR_RNDN);
                }
            }

            mpfr_clear(term);
        }

        for (int i = 0; i < num_threads; i++) {
            mpfr_add(result, result, partial_sums[i], MPFR_RNDN);
            mpfr_clear(partial_sums[i]);
        }
        delete[] partial_sums;

        mpfr_frac(result, result, MPFR_RNDN);
        if (mpfr_sgn(result) < 0) {
            mpfr_add_ui(result, result, 1, MPFR_RNDN);
        }
    }

    // Compute {5^N · Q₁} where Q₁ = Σ 1/((4j+3)·625^j)
    void compute_frac_Q1(mpfr_t result, long N) {
        mpfr_set_zero(result, 1);

        unsigned long max_j = (N > 0 ? N / 4 : 0) + precision / 2 + 50;

        int num_threads = 1;
        #ifdef _OPENMP
        num_threads = omp_get_max_threads();
        #endif

        mpfr_t* partial_sums = new mpfr_t[num_threads];
        for (int i = 0; i < num_threads; i++) {
            mpfr_init2(partial_sums[i], precision + 100);
            mpfr_set_zero(partial_sums[i], 1);
        }

        #pragma omp parallel
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif

            mpfr_t term;
            mpfr_init2(term, precision + 100);

            #pragma omp for schedule(dynamic, 100)
            for (unsigned long j = 0; j <= max_j; j++) {
                unsigned long denom = 4 * j + 3;
                long exp = N - 4 * (long)j;

                if (exp >= 0) {
                    unsigned long numerator = mod_pow(5, (unsigned long)exp, denom);
                    mpfr_set_ui(term, numerator, MPFR_RNDN);
                    mpfr_div_ui(term, term, denom, MPFR_RNDN);
                } else {
                    mpfr_set_ui(term, 5, MPFR_RNDN);
                    mpfr_pow_ui(term, term, (unsigned long)(-exp), MPFR_RNDN);
                    mpfr_mul_ui(term, term, denom, MPFR_RNDN);
                    mpfr_ui_div(term, 1, term, MPFR_RNDN);
                }

                mpfr_add(partial_sums[tid], partial_sums[tid], term, MPFR_RNDN);

                if (mpfr_cmp_ui(partial_sums[tid], 1000) > 0) {
                    mpfr_frac(partial_sums[tid], partial_sums[tid], MPFR_RNDN);
                }
            }

            mpfr_clear(term);
        }

        for (int i = 0; i < num_threads; i++) {
            mpfr_add(result, result, partial_sums[i], MPFR_RNDN);
            mpfr_clear(partial_sums[i]);
        }
        delete[] partial_sums;

        mpfr_frac(result, result, MPFR_RNDN);
        if (mpfr_sgn(result) < 0) {
            mpfr_add_ui(result, result, 1, MPFR_RNDN);
        }
    }

    // Compute {5^N · arctan(1/239)}
    void compute_frac_arctan239(mpfr_t result, unsigned long N) {
        mpfr_t five_N, power, term, inv_239_sq;
        mpfr_init2(five_N, precision + 100);
        mpfr_init2(power, precision + 100);
        mpfr_init2(term, precision + 100);
        mpfr_init2(inv_239_sq, precision + 100);

        mpfr_set_ui(five_N, 5, MPFR_RNDN);
        mpfr_pow_ui(five_N, five_N, N, MPFR_RNDN);

        mpfr_set_ui(power, 1, MPFR_RNDN);
        mpfr_div_ui(power, power, 239, MPFR_RNDN);

        mpfr_set_ui(inv_239_sq, 1, MPFR_RNDN);
        mpfr_div_ui(inv_239_sq, inv_239_sq, 239*239, MPFR_RNDN);

        mpfr_set_zero(result, 1);

        unsigned long max_terms = (unsigned long)(0.42 * N) + 100;

        for (unsigned long k = 0; k <= max_terms; k++) {
            mpfr_mul(term, five_N, power, MPFR_RNDN);
            mpfr_div_ui(term, term, 2 * k + 1, MPFR_RNDN);

            if (k % 2 == 0) {
                mpfr_add(result, result, term, MPFR_RNDN);
            } else {
                mpfr_sub(result, result, term, MPFR_RNDN);
            }

            mpfr_abs(term, term, MPFR_RNDN);
            if (mpfr_cmp_d(term, pow(10.0, -(double)(precision + 10))) < 0) {
                break;
            }

            mpfr_mul(power, power, inv_239_sq, MPFR_RNDN);

            if (mpfr_cmpabs_ui(result, 1000) > 0) {
                mpfr_frac(result, result, MPFR_RNDN);
            }
        }

        mpfr_frac(result, result, MPFR_RNDN);
        if (mpfr_sgn(result) < 0) {
            mpfr_add_ui(result, result, 1, MPFR_RNDN);
        }

        mpfr_clear(five_N);
        mpfr_clear(power);
        mpfr_clear(term);
        mpfr_clear(inv_239_sq);
    }

    // Helper: compute {x} ensuring result in [0,1)
    void frac_positive(mpfr_t x) {
        mpfr_frac(x, x, MPFR_RNDN);
        if (mpfr_sgn(x) < 0) {
            mpfr_add_ui(x, x, 1, MPFR_RNDN);
        }
    }

public:
    MachineDigitExtractor(int guard = 100) : guard_digits(guard) {}

    int extract_digit(unsigned long N) {
        // Precision in bits: need ~N*log2(10) + guard for the 2^(N-1) multiplication
        precision = (mpfr_prec_t)(N * 3.5) + guard_digits + 200;

        mpfr_t frac_A, frac_B, frac_C, frac_5pi, frac_10pi, temp, P1, Q1, five_power;
        mpfr_init2(frac_A, precision + 100);
        mpfr_init2(frac_B, precision + 100);
        mpfr_init2(frac_C, precision + 100);
        mpfr_init2(frac_5pi, precision + 100);
        mpfr_init2(frac_10pi, precision + 100);
        mpfr_init2(temp, precision + 100);
        mpfr_init2(P1, precision + 100);
        mpfr_init2(Q1, precision + 100);
        mpfr_init2(five_power, precision + 100);

        // Compute {16 * 5^(N-2) * P₁}
        if (N >= 2) {
            compute_frac_P1(temp, (long)N - 2);
        } else {
            // N=1: compute P1 directly, then multiply by 5^(N-2) = 5^(-1) = 0.2
            compute_P1_direct(P1);
            mpfr_set_ui(five_power, 5, MPFR_RNDN);
            mpfr_pow_si(five_power, five_power, (long)N - 2, MPFR_RNDN);
            mpfr_mul(temp, P1, five_power, MPFR_RNDN);
            frac_positive(temp);
        }
        mpfr_mul_ui(frac_A, temp, 16, MPFR_RNDN);
        frac_positive(frac_A);

        // Compute {16 * 5^(N-4) * Q₁}
        if (N >= 4) {
            compute_frac_Q1(temp, (long)N - 4);
        } else {
            // N<4: compute Q1 directly, then multiply by 5^(N-4)
            compute_Q1_direct(Q1);
            mpfr_set_ui(five_power, 5, MPFR_RNDN);
            mpfr_pow_si(five_power, five_power, (long)N - 4, MPFR_RNDN);
            mpfr_mul(temp, Q1, five_power, MPFR_RNDN);
            frac_positive(temp);
        }
        mpfr_mul_ui(frac_B, temp, 16, MPFR_RNDN);
        frac_positive(frac_B);

        // Compute {4 * 5^(N-1) * arctan(1/239)}
        compute_frac_arctan239(temp, N - 1);
        mpfr_mul_ui(frac_C, temp, 4, MPFR_RNDN);
        frac_positive(frac_C);

        // Borrow tracking: {A - B - C}
        int eps1 = (mpfr_cmp(frac_A, frac_B) < 0) ? 1 : 0;
        mpfr_sub(frac_5pi, frac_A, frac_B, MPFR_RNDN);
        mpfr_add_ui(frac_5pi, frac_5pi, eps1, MPFR_RNDN);

        int eps2 = (mpfr_cmp(frac_5pi, frac_C) < 0) ? 1 : 0;
        mpfr_sub(frac_5pi, frac_5pi, frac_C, MPFR_RNDN);
        mpfr_add_ui(frac_5pi, frac_5pi, eps2, MPFR_RNDN);

        frac_positive(frac_5pi);

        // Base conversion: {10^(N-1) · π} = {2^(N-1) · {5^(N-1) · π}}
        mpfr_set_ui(temp, 2, MPFR_RNDN);
        mpfr_pow_ui(temp, temp, N - 1, MPFR_RNDN);
        mpfr_mul(frac_10pi, frac_5pi, temp, MPFR_RNDN);
        frac_positive(frac_10pi);

        // Extract digit
        mpfr_mul_ui(frac_10pi, frac_10pi, 10, MPFR_RNDN);
        int digit = (int)mpfr_get_ui(frac_10pi, MPFR_RNDZ);

        mpfr_clear(frac_A);
        mpfr_clear(frac_B);
        mpfr_clear(frac_C);
        mpfr_clear(frac_5pi);
        mpfr_clear(frac_10pi);
        mpfr_clear(temp);
        mpfr_clear(P1);
        mpfr_clear(Q1);
        mpfr_clear(five_power);

        return digit;
    }
};

struct KnownDigit {
    unsigned long pos;
    int digit;
};

KnownDigit KNOWN_DIGITS[] = {
    {1, 1}, {2, 4}, {3, 1}, {4, 5}, {5, 9}, {6, 2}, {7, 6}, {8, 5}, {9, 3}, {10, 5},
    {100, 9}, {500, 2}, {1000, 9}
};
const int NUM_KNOWN = sizeof(KNOWN_DIGITS) / sizeof(KnownDigit);

void print_header() {
    std::cout << std::string(70, '=') << "\n";
    std::cout << "  MACHIN-TYPE DECIMAL DIGIT EXTRACTION FOR π\n";
    std::cout << "  High-Performance C++ Implementation with OpenMP\n";
    std::cout << std::string(70, '=') << "\n\n";

    #ifdef _OPENMP
    std::cout << "OpenMP enabled with " << omp_get_max_threads() << " threads\n";
    #else
    std::cout << "OpenMP not available (single-threaded)\n";
    #endif

    std::cout << "MPFR version: " << mpfr_get_version() << "\n";
    std::cout << "GMP version: " << gmp_version << "\n\n";
}

int main(int argc, char* argv[]) {
    print_header();

    MachineDigitExtractor extractor;

    // Verify first 20 digits
    std::cout << "Verifying first 20 digits of π:\n";
    std::cout << "π = 3.";
    for (int i = 1; i <= 20; i++) {
        std::cout << extractor.extract_digit(i);
    }
    std::cout << "...\n";
    std::cout << "Expected: 3.14159265358979323846...\n\n";

    // Verify all 1000 digits
    std::cout << "Verifying first 1000 digits...\n";

    // Compute reference pi
    mpfr_t pi;
    mpfr_init2(pi, 4000);
    mpfr_const_pi(pi, MPFR_RNDN);
    char* pi_str = new char[1100];
    mpfr_sprintf(pi_str, "%.1010Rf", pi);

    // Find decimal point and get digits after it
    std::string pi_digits;
    bool found_dot = false;
    for (int i = 0; pi_str[i] && pi_digits.size() < 1000; i++) {
        if (found_dot && pi_str[i] >= '0' && pi_str[i] <= '9') {
            pi_digits += pi_str[i];
        }
        if (pi_str[i] == '.') found_dot = true;
    }

    auto start = std::chrono::high_resolution_clock::now();

    int errors = 0;
    for (int n = 1; n <= 1000; n++) {
        int computed = extractor.extract_digit(n);
        int expected = pi_digits[n-1] - '0';
        if (computed != expected) {
            errors++;
            if (errors <= 10) {
                std::cout << "  ERROR at position " << n << ": got " << computed
                          << ", expected " << expected << "\n";
            }
        }
        if (n % 100 == 0) {
            std::cout << "  Verified " << n << "/1000...\n";
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "\n";
    std::cout << std::string(50, '-') << "\n";
    std::cout << "Total time: " << std::fixed << std::setprecision(2) << elapsed.count() << " seconds\n";

    if (errors == 0) {
        std::cout << "All 1000 digits verified correctly!\n";
    } else {
        std::cout << errors << " errors found!\n";
    }

    delete[] pi_str;
    mpfr_clear(pi);

    return errors == 0 ? 0 : 1;
}
