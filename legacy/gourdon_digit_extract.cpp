/**
 * Gourdon (2003) n-th decimal digit of pi, low memory
 * ===================================================
 *
 * Implements Algorithm 1 + Algorithm 2 from:
 *   Xavier Gourdon, "Computation of the n-th decimal digit of π with low memory", 2003.
 *   https://plouffe.fr/simon/articles/nthdecimaldigit.pdf
 *
 * This is a *verification-grade* implementation:
 * - Uses only modular arithmetic + floating accumulation mod 1.
 * - Includes optional MPFR reference check for small/medium positions.
 *
 * Build (macOS Homebrew):
 *   clang++ -std=c++17 -O3 -march=native -o gourdon_digit_extract gourdon_digit_extract.cpp -lmpfr -lgmp
 *
 * Build (Linux):
 *   g++ -std=c++17 -O3 -march=native -o gourdon_digit_extract gourdon_digit_extract.cpp -lmpfr -lgmp
 *
 * Usage:
 *   ./gourdon_digit_extract <pos> [guard_digits]
 *
 * Where:
 *   pos = 1-indexed decimal digit position after the decimal point.
 *   guard_digits = n0 in Gourdon's paper (default: 12).
 *
 * Notes:
 * - This implementation is intended for research/benchmarking and correctness testing,
 *   not yet heavily optimized (e.g., factoring/modular-power caching is not implemented).
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// MPFR for optional verification (optional dependency).
#if defined(__has_include)
#  if __has_include(<mpfr.h>)
#    include <mpfr.h>
#    define HAVE_MPFR 1
#  endif
#endif
#ifndef HAVE_MPFR
#  define HAVE_MPFR 0
#endif

static inline uint64_t mod_mul_u64(uint64_t a, uint64_t b, uint64_t mod) {
    return (uint64_t)((__uint128_t)a * b % mod);
}

static uint64_t mod_pow_u64(uint64_t base, uint64_t exp, uint64_t mod) {
    if (mod == 1) return 0;
    uint64_t res = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) res = mod_mul_u64(res, base, mod);
        exp >>= 1;
        if (exp) base = mod_mul_u64(base, base, mod);
    }
    return res;
}

static int64_t egcd(int64_t a, int64_t b, int64_t& x, int64_t& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int64_t x1 = 0, y1 = 0;
    int64_t g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

static uint64_t mod_inv_u64(uint64_t a, uint64_t mod) {
    int64_t x = 0, y = 0;
    int64_t g = egcd((int64_t)a, (int64_t)mod, x, y);
    if (g != 1 && g != -1) {
        throw std::runtime_error("mod_inv: non-invertible");
    }
    int64_t r = x % (int64_t)mod;
    if (r < 0) r += (int64_t)mod;
    return (uint64_t)r;
}

static std::vector<uint32_t> sieve_primes_upto(uint32_t n) {
    if (n < 2) return {};
    std::vector<bool> is_prime(n + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (uint32_t p = 2; (uint64_t)p * p <= n; p++) {
        if (!is_prime[p]) continue;
        for (uint64_t q = (uint64_t)p * p; q <= n; q += p) is_prime[(size_t)q] = false;
    }
    std::vector<uint32_t> primes;
    primes.reserve((size_t)(n / std::log((double)std::max<uint32_t>(n, 3))));
    for (uint32_t i = 2; i <= n; i++) if (is_prime[i]) primes.push_back(i);
    return primes;
}

static inline long double frac_pos(long double x) {
    long double f = x - std::floor(x);
    if (f < 0) f += 1.0L;
    // normalize tiny negatives due to fp error
    if (f >= 1.0L) f -= 1.0L;
    return f;
}

// Borrow-tracking for fractional parts: {A - B} from {A}, {B}
static inline long double frac_diff(long double fa, long double fb) {
    // assumes fa, fb are already in [0,1)
    long double x = fa - fb;
    if (x < 0) x += 1.0L;
    // normalize
    if (x >= 1.0L) x -= 1.0L;
    return x;
}

// Kahan compensated summation (for long double).
struct KahanSum {
    long double s = 0.0L;
    long double c = 0.0L;
    void add(long double x) {
        long double y = x - c;
        long double t = s + y;
        c = (t - s) - y;
        s = t;
    }
};

// Fixed-point accumulator modulo 1 with scale 2^64.
// Stores v in [0, 2^64) representing v / 2^64 in [0,1).
struct FixedFrac64 {
    uint64_t v = 0;

    static uint64_t term(uint64_t numer, uint64_t denom) {
        // Return round(numer/denom * 2^64) modulo 2^64.
        // Precondition: denom > 0 and numer < denom (true for our mod residues).
        __uint128_t x = ((__uint128_t)numer << 64);
        x += (denom / 2); // round-to-nearest
        return (uint64_t)(x / denom);
    }

    void add_term(uint64_t numer, uint64_t denom) { v += term(numer, denom); }
    void sub_term(uint64_t numer, uint64_t denom) { v -= term(numer, denom); }

    long double to_ld() const { return (long double)v / std::ldexp(1.0L, 64); }
};

static inline FixedFrac64 fixed_diff(FixedFrac64 a, FixedFrac64 b) {
    FixedFrac64 r;
    r.v = a.v - b.v; // wraps mod 2^64
    return r;
}

static std::tuple<int, bool> digit_from_frac_with_error(long double frac, long double err) {
    // frac in [0,1), err >= 0 absolute error bound on frac.
    frac = frac_pos(frac);
    long double y = 10.0L * frac;
    int d = (int)std::floor(y);
    if (d < 0) d = 0;
    if (d > 9) d = 9;
    long double lo = (long double)d / 10.0L;
    long double hi = (long double)(d + 1) / 10.0L;
    // Sufficient condition (no wraparound because err is tiny in practice).
    bool ok = (frac - err >= lo) && (frac + err < hi);
    return {d, ok};
}

static uint64_t choose_u64(const uint64_t n, const uint64_t k) {
    // Exact binomial for small n (used only in self-tests).
    if (k > n) return 0;
    uint64_t kk = k;
    if (kk > n - kk) kk = n - kk;
    __uint128_t num = 1;
    __uint128_t den = 1;
    for (uint64_t i = 1; i <= kk; i++) {
        num *= (n - kk + i);
        den *= i;
        // exact division at the end is safe for small n; keep values bounded
    }
    return (uint64_t)(num / den);
}

static uint64_t binom_prefix_sum_mod_bruteforce(uint64_t N, uint64_t k, uint64_t m) {
    if (m == 1) return 0;
    if (k > N) k = N;
    uint64_t s = 0;
    for (uint64_t j = 0; j <= k; j++) {
        uint64_t c = choose_u64(N, j);
        s += (c % m);
        s %= m;
    }
    return s;
}

// Algorithm 2 in Gourdon (2003): compute S = sum_{j=0}^k binom(N,j) mod m
static uint64_t binom_prefix_sum_mod(uint64_t N, uint64_t k, uint64_t m, const std::vector<uint32_t>& primes) {
    if (m == 1) return 0;
    if (k >= N) {
        // sum_{j=0}^N binom(N,j) = 2^N
        return mod_pow_u64(2, N, m);
    }
    // Optimization noted by Gourdon: if k > N/2, compute 2^N - sum_{j=0}^{N-k-1} binom(N,j)
    if (k > N / 2) {
        uint64_t other = binom_prefix_sum_mod(N, (N - k - 1), m, primes);
        uint64_t twoN = mod_pow_u64(2, N, m);
        return (twoN + m - other) % m;
    }

    // Step 1: factor m, but only keep distinct prime factors p <= k
    std::vector<uint32_t> pf;
    uint64_t tmp = m;
    for (uint32_t p : primes) {
        if (p > k) break;
        if ((uint64_t)p * p > tmp) break;
        if (tmp % p == 0) {
            pf.push_back(p);
            while (tmp % p == 0) tmp /= p;
        }
    }
    if (tmp > 1 && tmp <= k) {
        // remaining prime factor <= k
        pf.push_back((uint32_t)tmp);
    }

    const size_t L = pf.size();
    std::vector<uint64_t> R(L, 1); // exact prime powers (as integers)

    uint64_t A = 1 % m;
    uint64_t B = 1 % m;
    uint64_t C = 1 % m; // corresponds to sum up to j=0, since binom(N,0)=1

    for (uint64_t j = 1; j <= k; j++) {
        uint64_t a = N - j + 1;
        uint64_t b = j;

        // Step 3b/c: remove prime factors <=k from a and b, updating R_i
        for (size_t idx = 0; idx < L; idx++) {
            uint32_t p = pf[idx];
            uint64_t alpha = 0;
            while (a % p == 0) {
                a /= p;
                alpha++;
            }
            uint64_t beta = 0;
            while (b % p == 0) {
                b /= p;
                beta++;
            }
            if (alpha == beta) continue;
            if (alpha > beta) {
                uint64_t d = alpha - beta;
                // multiply R[idx] by p^d
                for (uint64_t t = 0; t < d; t++) R[idx] *= (uint64_t)p;
            } else {
                uint64_t d = beta - alpha;
                // divide R[idx] by p^d (exact)
                for (uint64_t t = 0; t < d; t++) R[idx] /= (uint64_t)p;
            }
        }

        uint64_t a_star_mod = a % m;
        uint64_t b_star_mod = b % m;

        A = mod_mul_u64(A, a_star_mod, m);
        B = mod_mul_u64(B, b_star_mod, m);

        // Compute prod_i R_i (mod m)
        uint64_t Rprod_mod = 1 % m;
        for (size_t idx = 0; idx < L; idx++) {
            Rprod_mod = mod_mul_u64(Rprod_mod, (R[idx] % m), m);
        }

        // C := C*b* + A*Rprod  (mod m)
        uint64_t term1 = mod_mul_u64(C, b_star_mod, m);
        uint64_t term2 = mod_mul_u64(A, Rprod_mod, m);
        C = term1 + term2;
        if (C >= m) C %= m;
    }

    // Final step: S = C / B (mod m) i.e. C * inv(B) mod m
    uint64_t invB = mod_inv_u64(B % m, m);
    return mod_mul_u64(C % m, invB, m);
}

struct GourdonParams {
    uint64_t n;   // exponent for 10^n
    uint64_t n0;  // guard digits
    uint64_t M;
    uint64_t N;
};

static GourdonParams choose_params(uint64_t n, uint64_t n0) {
    if (n < 10) {
        throw std::runtime_error("n too small for Gourdon parameter heuristics; use n>=10");
    }
    long double ln_n = std::log((long double)n);
    long double ln_n3 = ln_n * ln_n * ln_n;
    long double M_real = 2.0L * std::ceil((long double)n / ln_n3);
    uint64_t M = (uint64_t)M_real;
    if (M < 4) M = 4;
    if (M & 1) M++; // ensure even

    // N = ceil((n+n0+1)*ln(10) / ln(2 e M))
    long double denom = std::log(2.0L * std::exp(1.0L) * (long double)M);
    long double N_real = std::ceil(((long double)(n + n0 + 1) * std::log(10.0L)) / denom);
    uint64_t Ncap = (uint64_t)std::max<long double>(1.0L, N_real);

    return GourdonParams{n, n0, M, Ncap};
}

struct CVZParamsAB {
    uint64_t n;   // exponent for 10^n
    uint64_t n0;  // guard digits
    uint64_t a;   // x^a
    uint64_t b;   // (1-x)^b
    uint64_t N;   // outer power
};

static long double cvz_logG(uint64_t a, uint64_t b) {
    // logG = ln( 2^b * (a+b)^(a+b) / (a^a b^b) )
    // must be > 0 for acceleration.
    if (a == 0 || b == 0) return 0.0L;
    long double aa = (long double)a;
    long double bb = (long double)b;
    long double apb = aa + bb;
    return (bb * std::log(2.0L))
        + (apb * std::log(apb))
        - (aa * std::log(aa))
        - (bb * std::log(bb));
}

static uint64_t cvz_choose_N(uint64_t n, uint64_t n0, uint64_t a, uint64_t b) {
    long double logG = cvz_logG(a, b);
    if (!(logG > 0)) throw std::runtime_error("invalid (a,b): acceleration factor not > 1");
    long double N_real = std::ceil(((long double)(n + n0 + 1) * std::log(10.0L)) / logG);
    return (uint64_t)std::max<long double>(1.0L, N_real);
}

static CVZParamsAB choose_params_ab(uint64_t n, uint64_t n0, uint64_t a_in, uint64_t b_in) {
    if (n < 10) throw std::runtime_error("n too small for parameter heuristics; use n>=10");
    if (b_in < 1) throw std::runtime_error("b must be >= 1");
    uint64_t a = a_in;
    uint64_t b = b_in;
    if (a < 1) a = 1;
    if (a & 1) a++; // keep even like Gourdon (not required, but tends to help symmetry)

    uint64_t N = cvz_choose_N(n, n0, a, b);
    // Need exponent 10^(n-bN+2) to be nonnegative (Gourdon-like integer numerator condition).
    // We allow +2 slack like in Gourdon's proof.
    for (int it = 0; it < 64 && b * N > n + 2; it++) {
        a += 2;
        N = cvz_choose_N(n, n0, a, b);
    }
    if (b * N > n + 2) {
        throw std::runtime_error("failed to find parameters with b*N <= n+2; try smaller b or larger n");
    }

    return CVZParamsAB{n, n0, a, b, N};
}

static CVZParamsAB choose_best_params_ab(uint64_t n, uint64_t n0, uint64_t b,
                                        uint64_t a_hint = 0,
                                        uint64_t window_mul = 4,
                                        uint64_t window_add = 2000) {
    if (b < 1) throw std::runtime_error("b must be >= 1");
    if (n < 10) throw std::runtime_error("n too small for parameter heuristics; use n>=10");

    // Default hint: Gourdon's M ~ 2 * ceil(n / log^3 n)
    if (a_hint == 0) {
        long double ln_n = std::log((long double)n);
        long double ln_n3 = ln_n * ln_n * ln_n;
        long double M_real = 2.0L * std::ceil((long double)n / ln_n3);
        a_hint = (uint64_t)std::max<long double>(4.0L, M_real);
    }
    if (a_hint & 1) a_hint++;

    uint64_t a_lo = std::max<uint64_t>(2, a_hint / window_mul);
    uint64_t a_hi = std::max<uint64_t>(a_lo + 2, a_hint * window_mul + window_add);
    if (a_lo & 1) a_lo++;
    if (a_hi & 1) a_hi--;

    bool found = false;
    CVZParamsAB best{n, n0, 0, b, 0};
    __uint128_t best_cost = 0;

    auto try_range = [&](uint64_t lo, uint64_t hi) {
        for (uint64_t a = lo; a <= hi; a += 2) {
            long double logG = cvz_logG(a, b);
            if (!(logG > 0)) continue;
            uint64_t N = cvz_choose_N(n, n0, a, b);
            if ((__uint128_t)b * N > (__uint128_t)n + 2) continue;
            uint64_t m = (a + b) * N;
            __uint128_t cost = (__uint128_t)m + (__uint128_t)b * N;
            if (!found || cost < best_cost) {
                found = true;
                best = CVZParamsAB{n, n0, a, b, N};
                best_cost = cost;
            }
        }
    };

    try_range(a_lo, a_hi);
    if (!found) {
        // Expand search if needed
        try_range(2, a_hi + 8000);
    }
    if (!found) {
        throw std::runtime_error("no feasible (a,N) found for given b and n; try smaller b or larger n");
    }
    return best;
}

static int extract_digit_gourdon(uint64_t pos, uint64_t n0) {
    if (pos < 1) throw std::runtime_error("pos must be >= 1");
    // We need {10^(pos-1) * pi}. In Gourdon paper, n denotes exponent on 10^n.
    uint64_t n = pos - 1;
    if (n < 10) {
        // For tiny n, avoid parameter heuristics. long double precision is enough here.
        long double pi = acosl(-1.0L);
        long double x = pi * std::powl(10.0L, (long double)n);
        long double f = frac_pos(x);
        int d = (int)std::floor(10.0L * f + 1e-18L);
        if (d < 0) d = 0;
        if (d > 9) d = 9;
        return d;
    }

    GourdonParams p = choose_params(n, n0);
    if (p.N > p.n + 2) {
        // The paper's key condition for integer numerators is N <= n+2 (for n >= 4n0).
        // If this fails, we increase M slightly to increase ln(2eM) and shrink N.
        // (Simple fallback; a tighter implementation would tune M until the condition holds.)
        for (int t = 0; t < 32 && p.N > p.n + 2; t++) {
            p.M += 2;
            long double denom = std::log(2.0L * std::exp(1.0L) * (long double)p.M);
            long double N_real = std::ceil(((long double)(p.n + p.n0 + 1) * std::log(10.0L)) / denom);
            p.N = (uint64_t)std::max<long double>(1.0L, N_real);
        }
    }

    const uint64_t M = p.M;
    const uint64_t Ncap = p.N;

    // Pre-sieve primes up to Ncap (sufficient for Algorithm 2 factoring cutoff p<=k<=Ncap)
    if (Ncap > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("N too large for this implementation (prime sieve limit).");
    }
    std::vector<uint32_t> primes = sieve_primes_upto((uint32_t)Ncap);

    // Compute B
    const uint64_t Kmax = (M + 1) * Ncap;
    // Split alternating sum into two positive sums (even/odd indices).
    // This is exactly your "series splitting" trick, now applied to Gourdon's alternating sums,
    // which makes the accumulation numerically safer.
    KahanSum b_even, b_odd;
    KahanSum c_even, c_odd;

    for (uint64_t k = 0; k < Kmax; k++) {
        uint64_t denom = 2 * k + 1;
        uint64_t pow10 = mod_pow_u64(10, p.n, denom);
        uint64_t x = mod_mul_u64(4 % denom, pow10, denom);
        long double term = (long double)x / (long double)denom;
        if (k & 1) {
            b_odd.add(term);
        } else {
            b_even.add(term);
        }

        // keep bounded (reduce both occasionally)
        if ((k & 8191ULL) == 8191ULL) {
            b_even.s = frac_pos(b_even.s);
            b_odd.s = frac_pos(b_odd.s);
        }
    }
    long double b_frac = frac_diff(frac_pos(b_even.s), frac_pos(b_odd.s));

    // Compute C
    // Numerator is: 5^(N-2) * 10^(n-N+2) * s_k  (mod denom)
    // where denom = 2MN + 2k + 1, s_k = sum_{j=0}^k binom(N,j) mod denom
    int64_t exp10 = (int64_t)p.n - (int64_t)Ncap + 2;
    if (exp10 < 0) {
        // For our intended use (n >= 4n0), this shouldn't happen, but guard anyway.
        throw std::runtime_error("Parameter choice produced negative exponent for 10^(n-N+2).");
    }
    for (uint64_t k = 0; k < Ncap; k++) {
        uint64_t denom = 2 * M * Ncap + 2 * k + 1;
        uint64_t s_k = binom_prefix_sum_mod(Ncap, k, denom, primes);

        uint64_t a = mod_pow_u64(5, (Ncap >= 2 ? (Ncap - 2) : 0), denom);
        uint64_t b10 = mod_pow_u64(10, (uint64_t)exp10, denom);
        uint64_t ymod = mod_mul_u64(a, b10, denom);
        ymod = mod_mul_u64(ymod, s_k, denom);

        long double term = (long double)ymod / (long double)denom;
        if (k & 1) {
            c_odd.add(term);
        } else {
            c_even.add(term);
        }
        if ((k & 8191ULL) == 8191ULL) {
            c_even.s = frac_pos(c_even.s);
            c_odd.s = frac_pos(c_odd.s);
        }
    }
    long double c_frac = frac_diff(frac_pos(c_even.s), frac_pos(c_odd.s));

    // Final subtraction (again safe via frac_diff)
    long double xfrac = frac_diff(b_frac, c_frac);
    long double tenx = 10.0L * xfrac;
    int digit = (int)std::floor(tenx + 1e-18L); // tiny bias to counter fp underflow at integer boundaries
    if (digit < 0) digit = 0;
    if (digit > 9) digit = 9;
    return digit;
}

// Generalized CVZ/Gourdon family: P(x) = (x^a (1-x)^b)^N
// Derivation yields:
//   π ≈ 4 * Σ_{k=0}^{m-1} (-1)^k/(2k+1)  -  (4/2^{bN}) * Σ_{k=0}^{bN-1} (-1)^{aN+k} s_k /(2aN+2k+1)
// where m = (a+b)N and s_k = Σ_{j=0}^{k} C(bN, j).
// Multiply by 10^n, take fractional part using modular arithmetic.
static int extract_digit_cvz_ab(uint64_t pos, uint64_t n0, uint64_t a, uint64_t b) {
    if (pos < 1) throw std::runtime_error("pos must be >= 1");
    uint64_t n = pos - 1;
    if (n < 10) {
        long double pi = acosl(-1.0L);
        long double x = pi * std::powl(10.0L, (long double)n);
        int d = (int)std::floor(10.0L * frac_pos(x) + 1e-18L);
        if (d < 0) d = 0;
        if (d > 9) d = 9;
        return d;
    }

    CVZParamsAB p = choose_params_ab(n, n0, a, b);
    const uint64_t aN = p.a * p.N;
    const uint64_t bN = p.b * p.N;
    const uint64_t m = (p.a + p.b) * p.N;

    if (m > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("m too large for this implementation (prime sieve limit).");
    }
    // Need primes up to max(bN, m) for factoring in binomial prefix sum algorithm; we cap at m.
    std::vector<uint32_t> primes = sieve_primes_upto((uint32_t)std::max<uint64_t>(bN, (uint64_t)100));

    // B = Σ_{k=0}^{m-1} (-1)^k * (4*10^n mod (2k+1))/(2k+1)
    KahanSum b_even, b_odd;
    for (uint64_t k = 0; k < m; k++) {
        uint64_t denom = 2 * k + 1;
        uint64_t pow10 = mod_pow_u64(10, p.n, denom);
        uint64_t xmod = mod_mul_u64(4 % denom, pow10, denom);
        long double term = (long double)xmod / (long double)denom;
        if (k & 1) b_odd.add(term);
        else b_even.add(term);
        if ((k & 8191ULL) == 8191ULL) {
            b_even.s = frac_pos(b_even.s);
            b_odd.s = frac_pos(b_odd.s);
        }
    }
    long double B = frac_diff(frac_pos(b_even.s), frac_pos(b_odd.s));

    // C = Σ_{k=0}^{bN-1} (-1)^{aN+k} * (4 * 10^n / 2^{bN}) * s_k /(2aN+2k+1)
    // Scale: 10^n / 2^{bN} = 5^{bN} * 10^{n-bN}
    int64_t exp10 = (int64_t)p.n - (int64_t)bN;
    if (exp10 < 0) throw std::runtime_error("internal error: exp10 < 0 despite parameter constraint");
    KahanSum c_even, c_odd;
    for (uint64_t k = 0; k < bN; k++) {
        uint64_t denom = 2 * aN + 2 * k + 1;
        uint64_t s_k = binom_prefix_sum_mod(bN, k, denom, primes);
        uint64_t a5 = mod_pow_u64(5, bN, denom);
        uint64_t t10 = mod_pow_u64(10, (uint64_t)exp10, denom);
        uint64_t y = mod_mul_u64(4 % denom, a5, denom);
        y = mod_mul_u64(y, t10, denom);
        y = mod_mul_u64(y, s_k, denom);
        long double term = (long double)y / (long double)denom;
        // We will apply the (-1)^{aN+k} sign by even/odd splitting afterwards.
        if (k & 1) c_odd.add(term);
        else c_even.add(term);
        if ((k & 8191ULL) == 8191ULL) {
            c_even.s = frac_pos(c_even.s);
            c_odd.s = frac_pos(c_odd.s);
        }
    }
    // Apply (-1)^{aN+k}:
    // if aN even: sum_{even} term - sum_{odd} term
    // if aN odd : -(sum_{even} term - sum_{odd} term) = sum_{odd} term - sum_{even} term
    long double Ce = frac_pos(c_even.s);
    long double Co = frac_pos(c_odd.s);
    long double C = ((aN & 1ULL) == 0) ? frac_diff(Ce, Co) : frac_diff(Co, Ce);

    long double xfrac = frac_diff(B, C);
    int digit = (int)std::floor(10.0L * xfrac + 1e-18L);
    if (digit < 0) digit = 0;
    if (digit > 9) digit = 9;
    return digit;
}

static int extract_digit_cvz_ab_params(uint64_t pos, const CVZParamsAB& p) {
    // Same as extract_digit_cvz_ab, but uses a fixed (a,b,N) chosen externally.
    if (pos < 1) throw std::runtime_error("pos must be >= 1");
    uint64_t n = pos - 1;
    if (n != p.n) throw std::runtime_error("internal: pos does not match params n");
    const uint64_t aN = p.a * p.N;
    const uint64_t bN = p.b * p.N;
    const uint64_t m = (p.a + p.b) * p.N;

    if (n < 10) {
        long double pi = acosl(-1.0L);
        long double x = pi * std::powl(10.0L, (long double)n);
        int d = (int)std::floor(10.0L * frac_pos(x) + 1e-18L);
        if (d < 0) d = 0;
        if (d > 9) d = 9;
        return d;
    }
    if ((__uint128_t)p.b * p.N > (__uint128_t)p.n + 2) {
        throw std::runtime_error("params violate b*N <= n+2 condition");
    }

    std::vector<uint32_t> primes = sieve_primes_upto((uint32_t)std::max<uint64_t>(bN, (uint64_t)100));

    // B
    KahanSum b_even, b_odd;
    for (uint64_t k = 0; k < m; k++) {
        uint64_t denom = 2 * k + 1;
        uint64_t pow10 = mod_pow_u64(10, p.n, denom);
        uint64_t xmod = mod_mul_u64(4 % denom, pow10, denom);
        long double term = (long double)xmod / (long double)denom;
        if (k & 1) b_odd.add(term);
        else b_even.add(term);
        if ((k & 8191ULL) == 8191ULL) {
            b_even.s = frac_pos(b_even.s);
            b_odd.s = frac_pos(b_odd.s);
        }
    }
    long double B = frac_diff(frac_pos(b_even.s), frac_pos(b_odd.s));

    // C
    int64_t exp10 = (int64_t)p.n - (int64_t)bN;
    if (exp10 < 0) throw std::runtime_error("internal error: exp10 < 0");
    KahanSum c_even, c_odd;
    for (uint64_t k = 0; k < bN; k++) {
        uint64_t denom = 2 * aN + 2 * k + 1;
        uint64_t s_k = binom_prefix_sum_mod(bN, k, denom, primes);
        uint64_t a5 = mod_pow_u64(5, bN, denom);
        uint64_t t10 = mod_pow_u64(10, (uint64_t)exp10, denom);
        uint64_t y = mod_mul_u64(4 % denom, a5, denom);
        y = mod_mul_u64(y, t10, denom);
        y = mod_mul_u64(y, s_k, denom);
        long double term = (long double)y / (long double)denom;
        if (k & 1) c_odd.add(term);
        else c_even.add(term);
        if ((k & 8191ULL) == 8191ULL) {
            c_even.s = frac_pos(c_even.s);
            c_odd.s = frac_pos(c_odd.s);
        }
    }
    long double Ce = frac_pos(c_even.s);
    long double Co = frac_pos(c_odd.s);
    long double C = ((aN & 1ULL) == 0) ? frac_diff(Ce, Co) : frac_diff(Co, Ce);

    long double xfrac = frac_diff(B, C);
    int digit = (int)std::floor(10.0L * xfrac + 1e-18L);
    if (digit < 0) digit = 0;
    if (digit > 9) digit = 9;
    return digit;
}

static std::tuple<int, bool> extract_digit_cvz_ab_fixed64(uint64_t pos, uint64_t n0, uint64_t a, uint64_t b) {
    if (pos < 1) throw std::runtime_error("pos must be >= 1");
    uint64_t n = pos - 1;
    if (n < 10) {
        long double pi = acosl(-1.0L);
        return digit_from_frac_with_error(frac_pos(pi * std::powl(10.0L, (long double)n)), 0.0L);
    }
    CVZParamsAB p = choose_params_ab(n, n0, a, b);
    const uint64_t aN = p.a * p.N;
    const uint64_t bN = p.b * p.N;
    const uint64_t m = (p.a + p.b) * p.N;

    std::vector<uint32_t> primes = sieve_primes_upto((uint32_t)std::max<uint64_t>(bN, (uint64_t)100));

    FixedFrac64 b_even, b_odd;
    for (uint64_t k = 0; k < m; k++) {
        uint64_t denom = 2 * k + 1;
        uint64_t pow10 = mod_pow_u64(10, p.n, denom);
        uint64_t xmod = mod_mul_u64(4 % denom, pow10, denom);
        if (k & 1) b_odd.add_term(xmod, denom);
        else b_even.add_term(xmod, denom);
    }
    FixedFrac64 Bfx = fixed_diff(b_even, b_odd);

    int64_t exp10 = (int64_t)p.n - (int64_t)bN;
    if (exp10 < 0) throw std::runtime_error("internal error: exp10 < 0");

    FixedFrac64 c_even, c_odd;
    for (uint64_t k = 0; k < bN; k++) {
        uint64_t denom = 2 * aN + 2 * k + 1;
        uint64_t s_k = binom_prefix_sum_mod(bN, k, denom, primes);
        uint64_t a5 = mod_pow_u64(5, bN, denom);
        uint64_t t10 = mod_pow_u64(10, (uint64_t)exp10, denom);
        uint64_t y = mod_mul_u64(4 % denom, a5, denom);
        y = mod_mul_u64(y, t10, denom);
        y = mod_mul_u64(y, s_k, denom);
        if (k & 1) c_odd.add_term(y, denom);
        else c_even.add_term(y, denom);
    }
    FixedFrac64 Cfx = ((aN & 1ULL) == 0) ? fixed_diff(c_even, c_odd) : fixed_diff(c_odd, c_even);

    FixedFrac64 Xfx = fixed_diff(Bfx, Cfx);
    long double frac = Xfx.to_ld();
    long double T = (long double)(m + bN);
    long double err = T * std::ldexp(1.0L, -65);
    return digit_from_frac_with_error(frac, err);
}

// Same algorithm, but using FixedFrac64 accumulation with a provable error bound
// coming only from per-term rounding to 2^-64 (round-to-nearest => <= 2^-65 each).
static std::tuple<int, bool> extract_digit_gourdon_fixed64(uint64_t pos, uint64_t n0) {
    if (pos < 1) throw std::runtime_error("pos must be >= 1");
    uint64_t n = pos - 1;
    if (n < 10) {
        long double pi = acosl(-1.0L);
        long double x = pi * std::powl(10.0L, (long double)n);
        return digit_from_frac_with_error(frac_pos(x), 0.0L);
    }

    GourdonParams p = choose_params(n, n0);
    if (p.N > p.n + 2) {
        for (int t = 0; t < 32 && p.N > p.n + 2; t++) {
            p.M += 2;
            long double denom = std::log(2.0L * std::exp(1.0L) * (long double)p.M);
            long double N_real = std::ceil(((long double)(p.n + p.n0 + 1) * std::log(10.0L)) / denom);
            p.N = (uint64_t)std::max<long double>(1.0L, N_real);
        }
    }

    const uint64_t M = p.M;
    const uint64_t Ncap = p.N;

    if (Ncap > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("N too large for this implementation (prime sieve limit).");
    }
    std::vector<uint32_t> primes = sieve_primes_upto((uint32_t)Ncap);

    // B = even - odd
    const uint64_t Kmax = (M + 1) * Ncap;
    FixedFrac64 b_even, b_odd;
    for (uint64_t k = 0; k < Kmax; k++) {
        uint64_t denom = 2 * k + 1;
        uint64_t pow10 = mod_pow_u64(10, p.n, denom);
        uint64_t xmod = mod_mul_u64(4 % denom, pow10, denom);
        if (k & 1) b_odd.add_term(xmod, denom);
        else b_even.add_term(xmod, denom);
    }
    FixedFrac64 b_frac_fx = fixed_diff(b_even, b_odd);

    // C = even - odd
    int64_t exp10 = (int64_t)p.n - (int64_t)Ncap + 2;
    if (exp10 < 0) throw std::runtime_error("Parameter choice produced negative exponent for 10^(n-N+2).");
    FixedFrac64 c_even, c_odd;
    for (uint64_t k = 0; k < Ncap; k++) {
        uint64_t denom = 2 * M * Ncap + 2 * k + 1;
        uint64_t s_k = binom_prefix_sum_mod(Ncap, k, denom, primes);
        uint64_t a = mod_pow_u64(5, (Ncap >= 2 ? (Ncap - 2) : 0), denom);
        uint64_t b10 = mod_pow_u64(10, (uint64_t)exp10, denom);
        uint64_t ymod = mod_mul_u64(a, b10, denom);
        ymod = mod_mul_u64(ymod, s_k, denom);
        if (k & 1) c_odd.add_term(ymod, denom);
        else c_even.add_term(ymod, denom);
    }
    FixedFrac64 c_frac_fx = fixed_diff(c_even, c_odd);

    FixedFrac64 x_fx = fixed_diff(b_frac_fx, c_frac_fx);
    long double frac = x_fx.to_ld();

    long double T = (long double)(Kmax + Ncap);
    long double err = T * std::ldexp(1.0L, -65);
    return digit_from_frac_with_error(frac, err);
}

#if HAVE_MPFR
static int reference_digit_mpfr(uint64_t pos) {
    if (pos < 1) throw std::runtime_error("pos must be >=1");
    uint64_t n = pos - 1;
    // bits ~ (n+5) * log2(10) + safety
    mpfr_prec_t prec = (mpfr_prec_t)((long double)(n + 25) * 3.32192809488736234787L) + 256;
    mpfr_t pi, scale, x;
    mpfr_init2(pi, prec);
    mpfr_init2(scale, prec);
    mpfr_init2(x, prec);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_set_ui(scale, 10, MPFR_RNDN);
    mpfr_pow_ui(scale, scale, n, MPFR_RNDN);
    mpfr_mul(x, pi, scale, MPFR_RNDN);
    mpfr_frac(x, x, MPFR_RNDN);
    mpfr_mul_ui(x, x, 10, MPFR_RNDN);
    uint64_t d = mpfr_get_ui(x, MPFR_RNDZ);
    mpfr_clear(pi);
    mpfr_clear(scale);
    mpfr_clear(x);
    return (int)d;
}
#endif

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " <pos> [guard_digits]\n"
        << "  " << prog << " --ab <a> <b> <pos> [guard_digits]\n"
        << "  " << prog << " --ab-auto <b> <pos> [guard_digits]\n"
        << "  " << prog << " --ab-auto-fixed <b> <pos> [guard_digits]\n"
        << "  " << prog << " --ab-fixed <a> <b> <pos> [guard_digits]\n"
        << "  " << prog << " --ab-scan <bmin> <bmax> <pos> [guard_digits]\n"
        << "  " << prog << " --ab-scan-top <bmin> <bmax> <pos> [guard_digits] [topk]\n"
        << "  " << prog << " --ab-scan-save <bmin> <bmax> <pos> [guard_digits] <outfile>\n"
        << "  " << prog << " --fixed <pos> [guard_digits]\n"
        << "  " << prog << " --selftest\n"
        << "\n"
        << "pos: 1-indexed decimal digit position after decimal point.\n"
        << "guard_digits: n0 in Gourdon's paper (default 12).\n";
}

struct ScanRow {
    uint64_t b = 0;
    uint64_t a = 0;
    uint64_t N = 0;
    uint64_t m = 0;
    uint64_t bN = 0;
    long double logG = 0.0L;
    __uint128_t cost_terms = 0;
    int digit = -1;
    bool certified = false;
    bool error = false;
    std::string err;
};

static std::string u128_to_string(__uint128_t x) {
    if (x == 0) return "0";
    std::string s;
    while (x > 0) {
        int d = (int)(x % 10);
        s.push_back((char)('0' + d));
        x /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

static void print_scan_csv(std::ostream& os,
                           const std::vector<ScanRow>& rows,
                           uint64_t bmin, uint64_t bmax,
                           uint64_t pos, uint64_t n0,
                           bool mpfr_available, int mpfr_digit,
                           size_t limit_rows = 0) {
    os << "scan pos=" << pos << " n0=" << n0 << " b=[" << bmin << "," << bmax << "]\n";
    if (mpfr_available) {
        if (mpfr_digit >= 0) os << "mpfr_digit=" << mpfr_digit << "\n";
        else os << "mpfr_digit=skipped\n";
    } else {
        os << "mpfr_digit=unavailable\n";
    }
    os << "b,a,N,m,bN,logG,cost_terms,digit,certified,status\n";
    size_t nrows = rows.size();
    size_t take = (limit_rows == 0) ? nrows : std::min(nrows, limit_rows);
    for (size_t i = 0; i < take; i++) {
        const auto& r = rows[i];
        os << r.b << "," << r.a << "," << r.N << "," << r.m << "," << r.bN << ","
           << std::setprecision(8) << (double)r.logG << ","
           << u128_to_string(r.cost_terms) << ","
           << r.digit << "," << (r.certified ? "yes" : "no") << ","
           << (r.error ? r.err : "OK")
           << "\n";
    }
}

static std::vector<ScanRow> run_ab_scan(uint64_t bmin, uint64_t bmax, uint64_t pos, uint64_t n0) {
    if (bmin < 1 || bmax < 1 || bmin > bmax) {
        throw std::runtime_error("require 1 <= bmin <= bmax");
    }
    if (pos < 1) throw std::runtime_error("pos must be >= 1");
    uint64_t n = pos - 1;
    std::vector<ScanRow> rows;
    rows.reserve((size_t)(bmax - bmin + 1));
    for (uint64_t b = bmin; b <= bmax; b++) {
        if (((b - bmin) % 10) == 0) {
            std::cerr << "scan progress: b=" << b << "/" << bmax << "\n";
        }
        ScanRow r;
        r.b = b;
        try {
            CVZParamsAB best = choose_best_params_ab(n, n0, b);
            r.a = best.a;
            r.N = best.N;
            r.m = (best.a + best.b) * best.N;
            r.bN = best.b * best.N;
            r.logG = cvz_logG(best.a, best.b);
            r.cost_terms = (__uint128_t)r.m + (__uint128_t)r.bN;
            auto [d, ok] = extract_digit_cvz_ab_fixed64(pos, n0, best.a, best.b);
            r.digit = d;
            r.certified = ok;
        } catch (const std::exception& e) {
            r.error = true;
            r.err = e.what();
        }
        rows.push_back(r);
    }
    std::stable_sort(rows.begin(), rows.end(), [](const ScanRow& x, const ScanRow& y) {
        if (x.error != y.error) return !x.error;
        if (x.cost_terms != y.cost_terms) return x.cost_terms < y.cost_terms;
        return x.b < y.b;
    });
    return rows;
}

static bool selftest_core_math() {
    // Checkpoint 1: fractional arithmetic invariants for frac_pos / frac_diff
    {
        const std::array<long double, 7> xs = {
            0.0L, 0.1L, 0.9999999999999999L, -0.1L, -1.9L, 10.25L, -10.25L
        };
        for (long double x : xs) {
            long double f = frac_pos(x);
            if (!(f >= 0.0L && f < 1.0L)) return false;
        }
        // frac_diff should behave like (a-b) mod 1 for fractional inputs.
        const std::array<long double, 5> fs = {0.0L, 0.2L, 0.5L, 0.9L, 0.999999999999L};
        for (long double a : fs) {
            for (long double b : fs) {
                long double d = frac_diff(a, b);
                long double ref = frac_pos(a - b);
                if (fabsl(d - ref) > 1e-18L) return false;
            }
        }
    }

    // Checkpoint 2: modular exponentiation sanity
    {
        for (uint64_t mod : {3ULL, 5ULL, 17ULL, 101ULL}) {
            for (uint64_t e = 0; e < 200; e++) {
                uint64_t p = mod_pow_u64(10, e, mod);
                // Fermat sanity when mod is prime and gcd(10,mod)=1:
                // 10^(mod-1) ≡ 1 (mod mod). Not always applicable; just check p in range.
                if (p >= mod) return false;
            }
        }
    }

    // Checkpoint 3: Algorithm 2 (binomial prefix sum mod) vs brute force on small cases
    {
        auto primes = sieve_primes_upto(500);
        for (uint64_t N : {5ULL, 10ULL, 20ULL, 40ULL}) {
            for (uint64_t k = 0; k <= N; k++) {
                for (uint64_t m : {3ULL, 5ULL, 7ULL, 9ULL, 11ULL, 25ULL, 49ULL, 97ULL, 101ULL}) {
                    uint64_t got = binom_prefix_sum_mod(N, k, m, primes);
                    uint64_t ref = binom_prefix_sum_mod_bruteforce(N, k, m);
                    if (got != ref) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

static bool selftest_digits(uint64_t n0) {
    // Checkpoint 4: known small digits (cheap structural sanity)
    struct Pair { uint64_t pos; int digit; };
    const std::array<Pair, 6> known = {{
        {1, 1}, {2, 4}, {3, 1}, {10, 5}, {100, 9}, {1000, 9}
    }};
    for (auto [pos, dig] : known) {
        int got = extract_digit_gourdon(pos, n0);
        if (got != dig) return false;
    }

    // Checkpoint 5: guard stability on a deterministic set of positions
    const std::array<uint64_t, 8> posset = {1, 2, 3, 10, 100, 1000, 5000, 10000};
    for (uint64_t pos : posset) {
        int d1 = extract_digit_gourdon(pos, n0);
        int d2 = extract_digit_gourdon(pos, n0 + 4);
        int d3 = extract_digit_gourdon(pos, n0 + 8);
        if (!(d1 == d2 && d2 == d3)) return false;
    }

    // Checkpoint 5b: fixed-point backend matches and (usually) certifies on the same ladder
    for (uint64_t pos : posset) {
        auto [df, ok] = extract_digit_gourdon_fixed64(pos, n0);
        int d = extract_digit_gourdon(pos, n0);
        if (df != d) return false;
        if (pos >= 10 && !ok) return false;
    }

#if HAVE_MPFR
    // Checkpoint 6: MPFR spot-check on a deterministic ladder (moderate cost)
    const std::array<uint64_t, 7> refset = {1, 10, 100, 1000, 5000, 10000, 20000};
    for (uint64_t pos : refset) {
        int got = extract_digit_gourdon(pos, n0);
        int ref = reference_digit_mpfr(pos);
        if (got != ref) return false;
    }
#endif

    // Checkpoint 7: generalized CVZ family sanity + fixed-point certification
    {
        // A couple of (a,b) pairs that should work and match MPFR at moderate positions.
        const std::array<std::tuple<uint64_t, uint64_t, uint64_t>, 2> cases = {{
            {6, 1, 10000},
            {200, 2, 1000},
        }};
        for (auto [a, b, pos] : cases) {
            int d = extract_digit_cvz_ab(pos, n0, a, b);
            auto [df, ok] = extract_digit_cvz_ab_fixed64(pos, n0, a, b);
            if (d != df) return false;
            if (!ok) return false;
#if HAVE_MPFR
            int ref = reference_digit_mpfr(pos);
            if (d != ref) return false;
#endif
        }
    }

    return true;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv[0]);
        return 2;
    }
    std::string a1(argv[1]);
    if (a1 == "--help" || a1 == "-h") {
        usage(argv[0]);
        return 0;
    }
    if (a1 == "--selftest") {
        uint64_t n0 = 12;
        bool ok = selftest_core_math() && selftest_digits(n0);
        std::cout << "selftest=" << (ok ? "PASS" : "FAIL") << "\n";
#if HAVE_MPFR
        std::cout << "mpfr=enabled\n";
#else
        std::cout << "mpfr=disabled\n";
#endif
        return ok ? 0 : 1;
    }
    if (a1 == "--fixed") {
        if (argc < 3) {
            usage(argv[0]);
            return 2;
        }
        uint64_t pos = 0;
        uint64_t n0 = 12;
        try {
            pos = std::stoull(std::string(argv[2]));
            if (argc >= 4) n0 = std::stoull(std::string(argv[3]));
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            auto [d, ok] = extract_digit_gourdon_fixed64(pos, n0);
            std::cout << "pos=" << pos << " digit=" << d << " certified=" << (ok ? "yes" : "no") << "\n";
#if HAVE_MPFR
            if (pos <= 20000) {
                int ref = reference_digit_mpfr(pos);
                std::cout << "reference=" << ref << (ref == d ? "  OK" : "  MISMATCH") << "\n";
                if (ref != d) return 1;
            } else {
                std::cout << "reference=skipped (pos too large)\n";
            }
#endif
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (a1 == "--ab") {
        if (argc < 5) {
            usage(argv[0]);
            return 2;
        }
        uint64_t a = 0, b = 0, pos = 0;
        uint64_t n0 = 12;
        try {
            a = std::stoull(std::string(argv[2]));
            b = std::stoull(std::string(argv[3]));
            pos = std::stoull(std::string(argv[4]));
            if (argc >= 6) n0 = std::stoull(std::string(argv[5]));
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            int d = extract_digit_cvz_ab(pos, n0, a, b);
            std::cout << "pos=" << pos << " digit=" << d << " (ab=" << a << "," << b << ")\n";
#if HAVE_MPFR
            if (pos <= 20000) {
                int ref = reference_digit_mpfr(pos);
                std::cout << "reference=" << ref << (ref == d ? "  OK" : "  MISMATCH") << "\n";
                if (ref != d) return 1;
            } else {
                std::cout << "reference=skipped (pos too large)\n";
            }
#endif
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (a1 == "--ab-auto") {
        if (argc < 4) {
            usage(argv[0]);
            return 2;
        }
        uint64_t b = 0, pos = 0, n0 = 12;
        try {
            b = std::stoull(std::string(argv[2]));
            pos = std::stoull(std::string(argv[3]));
            if (argc >= 5) n0 = std::stoull(std::string(argv[4]));
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            uint64_t n = pos - 1;
            CVZParamsAB best = choose_best_params_ab(n, n0, b);
            int d = extract_digit_cvz_ab_params(pos, best);
            long double logG = cvz_logG(best.a, best.b);
            uint64_t m = (best.a + best.b) * best.N;
            __uint128_t cost = (__uint128_t)m + (__uint128_t)best.b * best.N;
            std::cout << "pos=" << pos << " digit=" << d
                      << " (auto ab=" << best.a << "," << best.b
                      << " N=" << best.N
                      << " m=" << m
                      << " bN=" << (best.b * best.N)
                      << " logG=" << std::setprecision(6) << (double)logG
                      << " cost_terms=" << (unsigned long long)cost
                      << ")\n";
#if HAVE_MPFR
            if (pos <= 20000) {
                int ref = reference_digit_mpfr(pos);
                std::cout << "reference=" << ref << (ref == d ? "  OK" : "  MISMATCH") << "\n";
                if (ref != d) return 1;
            } else {
                std::cout << "reference=skipped (pos too large)\n";
            }
#endif
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (a1 == "--ab-auto-fixed") {
        if (argc < 4) {
            usage(argv[0]);
            return 2;
        }
        uint64_t b = 0, pos = 0, n0 = 12;
        try {
            b = std::stoull(std::string(argv[2]));
            pos = std::stoull(std::string(argv[3]));
            if (argc >= 5) n0 = std::stoull(std::string(argv[4]));
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            uint64_t n = pos - 1;
            CVZParamsAB best = choose_best_params_ab(n, n0, b);
            auto [d, ok] = extract_digit_cvz_ab_fixed64(pos, n0, best.a, best.b);
            long double logG = cvz_logG(best.a, best.b);
            uint64_t m = (best.a + best.b) * best.N;
            __uint128_t cost = (__uint128_t)m + (__uint128_t)best.b * best.N;
            std::cout << "pos=" << pos << " digit=" << d << " certified=" << (ok ? "yes" : "no")
                      << " (auto ab=" << best.a << "," << best.b
                      << " N=" << best.N
                      << " m=" << m
                      << " bN=" << (best.b * best.N)
                      << " logG=" << std::setprecision(6) << (double)logG
                      << " cost_terms=" << (unsigned long long)cost
                      << ")\n";
#if HAVE_MPFR
            if (pos <= 20000) {
                int ref = reference_digit_mpfr(pos);
                std::cout << "reference=" << ref << (ref == d ? "  OK" : "  MISMATCH") << "\n";
                if (ref != d) return 1;
            } else {
                std::cout << "reference=skipped (pos too large)\n";
            }
#endif
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (a1 == "--ab-fixed") {
        if (argc < 5) {
            usage(argv[0]);
            return 2;
        }
        uint64_t a = 0, b = 0, pos = 0;
        uint64_t n0 = 12;
        try {
            a = std::stoull(std::string(argv[2]));
            b = std::stoull(std::string(argv[3]));
            pos = std::stoull(std::string(argv[4]));
            if (argc >= 6) n0 = std::stoull(std::string(argv[5]));
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            auto [d, ok] = extract_digit_cvz_ab_fixed64(pos, n0, a, b);
            std::cout << "pos=" << pos << " digit=" << d << " certified=" << (ok ? "yes" : "no")
                      << " (ab=" << a << "," << b << ")\n";
#if HAVE_MPFR
            if (pos <= 20000) {
                int ref = reference_digit_mpfr(pos);
                std::cout << "reference=" << ref << (ref == d ? "  OK" : "  MISMATCH") << "\n";
                if (ref != d) return 1;
            } else {
                std::cout << "reference=skipped (pos too large)\n";
            }
#endif
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (a1 == "--ab-scan") {
        if (argc < 5) {
            usage(argv[0]);
            return 2;
        }
        uint64_t bmin = 0, bmax = 0, pos = 0, n0 = 12;
        try {
            bmin = std::stoull(std::string(argv[2]));
            bmax = std::stoull(std::string(argv[3]));
            pos = std::stoull(std::string(argv[4]));
            if (argc >= 6) n0 = std::stoull(std::string(argv[5]));
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            auto rows = run_ab_scan(bmin, bmax, pos, n0);
#if HAVE_MPFR
            int ref = (pos <= 20000) ? reference_digit_mpfr(pos) : -1;
            if (ref >= 0) {
                for (auto& r : rows) {
                    if (!r.error && r.digit != ref) {
                        r.error = true;
                        r.err = "MPFR mismatch";
                    }
                }
                std::stable_sort(rows.begin(), rows.end(), [](const ScanRow& x, const ScanRow& y) {
                    if (x.error != y.error) return !x.error;
                    if (x.cost_terms != y.cost_terms) return x.cost_terms < y.cost_terms;
                    return x.b < y.b;
                });
            }
            print_scan_csv(std::cout, rows, bmin, bmax, pos, n0, true, ref, 0);
#else
            print_scan_csv(std::cout, rows, bmin, bmax, pos, n0, false, -1, 0);
#endif
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (a1 == "--ab-scan-top") {
        if (argc < 5) {
            usage(argv[0]);
            return 2;
        }
        uint64_t bmin = 0, bmax = 0, pos = 0, n0 = 12;
        size_t topk = 10;
        try {
            bmin = std::stoull(std::string(argv[2]));
            bmax = std::stoull(std::string(argv[3]));
            pos = std::stoull(std::string(argv[4]));
            if (argc >= 6) n0 = std::stoull(std::string(argv[5]));
            if (argc >= 7) topk = (size_t)std::stoull(std::string(argv[6]));
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            auto rows = run_ab_scan(bmin, bmax, pos, n0);
#if HAVE_MPFR
            int ref = (pos <= 20000) ? reference_digit_mpfr(pos) : -1;
            if (ref >= 0) {
                for (auto& r : rows) {
                    if (!r.error && r.digit != ref) {
                        r.error = true;
                        r.err = "MPFR mismatch";
                    }
                }
                std::stable_sort(rows.begin(), rows.end(), [](const ScanRow& x, const ScanRow& y) {
                    if (x.error != y.error) return !x.error;
                    if (x.cost_terms != y.cost_terms) return x.cost_terms < y.cost_terms;
                    return x.b < y.b;
                });
            }
            print_scan_csv(std::cout, rows, bmin, bmax, pos, n0, true, ref, topk);
#else
            print_scan_csv(std::cout, rows, bmin, bmax, pos, n0, false, -1, topk);
#endif
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (a1 == "--ab-scan-save") {
        if (argc < 6) {
            usage(argv[0]);
            return 2;
        }
        uint64_t bmin = 0, bmax = 0, pos = 0, n0 = 12;
        std::string outpath;
        try {
            bmin = std::stoull(std::string(argv[2]));
            bmax = std::stoull(std::string(argv[3]));
            pos = std::stoull(std::string(argv[4]));
            // If argv[5] is numeric, treat it as n0 and argv[6] as outfile; else argv[5] is outfile.
            try {
                n0 = std::stoull(std::string(argv[5]));
                if (argc < 7) {
                    usage(argv[0]);
                    return 2;
                }
                outpath = std::string(argv[6]);
            } catch (const std::exception&) {
                n0 = 12;
                outpath = std::string(argv[5]);
            }
        } catch (const std::exception&) {
            usage(argv[0]);
            return 2;
        }
        try {
            auto rows = run_ab_scan(bmin, bmax, pos, n0);
#if HAVE_MPFR
            int ref = (pos <= 20000) ? reference_digit_mpfr(pos) : -1;
            if (ref >= 0) {
                for (auto& r : rows) {
                    if (!r.error && r.digit != ref) {
                        r.error = true;
                        r.err = "MPFR mismatch";
                    }
                }
                std::stable_sort(rows.begin(), rows.end(), [](const ScanRow& x, const ScanRow& y) {
                    if (x.error != y.error) return !x.error;
                    if (x.cost_terms != y.cost_terms) return x.cost_terms < y.cost_terms;
                    return x.b < y.b;
                });
            }
#endif
            std::ofstream out(outpath);
            if (!out) throw std::runtime_error("failed to open outfile");
#if HAVE_MPFR
            print_scan_csv(out, rows, bmin, bmax, pos, n0, true, (pos <= 20000 ? reference_digit_mpfr(pos) : -1), 0);
#else
            print_scan_csv(out, rows, bmin, bmax, pos, n0, false, -1, 0);
#endif
            std::cout << "wrote " << outpath << "\n";
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (!a1.empty() && a1[0] == '-') {
        // Unknown flag
        usage(argv[0]);
        return 2;
    }

    uint64_t n0 = 12;
    uint64_t pos = 0;
    try {
        pos = std::stoull(a1);
        if (argc >= 3) n0 = std::stoull(argv[2]);
    } catch (const std::exception&) {
        usage(argv[0]);
        return 2;
    }

    try {
        // Two-pass stability check: run with n0 and n0+4; if different, warn.
        int d1 = extract_digit_gourdon(pos, n0);
        int d2 = extract_digit_gourdon(pos, n0 + 4);

        std::cout << "pos=" << pos << " digit=" << d1;
        if (d1 != d2) {
            std::cout << "  (WARNING: guard instability: " << d1 << " vs " << d2 << ")";
        }
        std::cout << "\n";

        // Reference check for moderate positions (MPFR cost grows with pos)
        #if HAVE_MPFR
        if (pos <= 20000) {
            int ref = reference_digit_mpfr(pos);
            std::cout << "reference=" << ref << (ref == d1 ? "  OK" : "  MISMATCH") << "\n";
            if (ref != d1) return 1;
        } else {
            std::cout << "reference=skipped (pos too large)\n";
        }
        #else
        std::cout << "reference=skipped (mpfr not available at build time)\n";
        #endif
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

