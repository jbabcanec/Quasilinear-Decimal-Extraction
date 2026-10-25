/*
 * Nine-digit benchmark adapter for Gourdon (2003), Algorithms 1 and 2.
 * Adapted from this repository's legacy/gourdon_digit_extract.cpp.
 * This is the manuscript author's reimplementation, not Gourdon's pidec.
 * Algorithms/parameters unchanged; one fixed-point extraction, no MPFR.
 * The prime sieve is cached; this executable does not enforce O(log^2 N)
 * total memory. See README.md for provenance and comparison scope.
 */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
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

static uint64_t extract_window_gourdon(uint64_t pos, uint64_t n0) {
    if (pos < 1) throw std::runtime_error("pos must be >= 1");
    uint64_t n = pos - 1;
    if (n < 10) throw std::runtime_error("position must be at least 11");

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
    // Nine digits from the SAME accumulated fractional part.
    // No repeated extraction and no claim of a full output certificate.
    return (uint64_t)(((__uint128_t)x_fx.v * 1000000000ULL) >> 64);
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: gourdon_window POSITION");
        uint64_t pos = std::stoull(argv[1]);
        // Restrict to the audited experiment, avoiding wide-input overflow.
        if (pos < 11 || pos > 1000000) throw std::runtime_error("position outside 11..1000000");
        auto digits = extract_window_gourdon(pos, 12);
        std::printf("%09llu\n", (unsigned long long)digits);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
}
