/**
 * three_phase_extract.cpp
 * =======================
 * N-th decimal digit of pi via the three-phase split of the hex BBP formula
 * (research_log.md Part 16; Python prototype: decimal_bbp_extract.py).
 *
 * Every term of 10^(N-1)*pi is +/- 2^E * 5^(N-1) / mu with mu odd:
 *   E >= 0 ("easy"): (2^E * 5^(N-1) mod mu)/mu via word-sized modexps.
 *   E <  0 ("band"): X/(mu*2^J), X = 5^(N-1), J = -E, and exactly
 *        frac(X/(mu*2^J)) = (h + u)/mu,
 *        h = floor(X/2^J) mod mu,   u = frac(X/2^J) (bit-window of X at J).
 *   All band terms share the single precomputed byte string of X.
 *   h comes from a Horner reduction over the shorter of the prefix [0,J) or
 *   suffix [J,LX) of X's bits (v1: naive per-term reduction, O(N^2/w) total;
 *   the batched prefix-remainder tree of Part 16.4 is future work).
 *
 * Build (MSYS2 g++):
 *   g++ -O2 -march=native -std=c++17 -fopenmp -o three_phase_extract three_phase_extract.cpp
 *
 * Usage:
 *   ./three_phase_extract <N>                     # digits of pi at position N
 *   ./three_phase_extract <N> --x5 <file>         # load 5^(N-1) bytes (little endian)
 *   ./three_phase_extract --check <digitsfile>    # verify N=1..60,100,1000 vs reference digits
 *   Optional: --threads <T>
 *
 * The digits file for --check is the decimal expansion of pi after the point
 * (one long line), e.g. produced by:
 *   python decimal_bbp_extract.py  (pi_reference)  -- see run notes.
 */
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>
#if defined(_OPENMP)
#  include <omp.h>
#endif

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef int64_t i64;

static const int G = 96;   // accumulator fraction bits
static const int GP = 80;  // window read bits
static const u128 MASKG = (((u128)1) << G) - 1;

static inline u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }

static u64 powmod(u64 b, u64 e, u64 m) {
    if (m == 1) return 0;
    u64 r = 1;
    b %= m;
    while (e) {
        if (e & 1) r = mulmod(r, b, m);
        b = mulmod(b, b, m);
        e >>= 1;
    }
    return r;
}

// floor((h*2^96 + U)/mu), h < mu < 2^63, U < 2^96. Exact, base-2^48 long division.
static inline u128 q96(u64 h, u128 U, u64 mu) {
    const u64 M48 = ((u64)1 << 48) - 1;
    u64 U_hi = (u64)(U >> 48);          // < 2^48
    u64 U_lo = (u64)U & M48;
    u128 t1 = (((u128)h) << 48) + U_hi;
    u128 q1 = t1 / mu;
    u64 r1 = (u64)(t1 % mu);
    u128 t2 = (((u128)r1) << 48) + U_lo;
    u128 q2 = t2 / mu;
    return (q1 << 48) + q2;
}

// ---------------------------------------------------------------------------
// X = 5^(N-1) as a little-endian byte string
// ---------------------------------------------------------------------------
struct XString {
    std::vector<uint8_t> b;
    i64 bits = 0;

    // read `len` (<= 96) bits at bit position `pos` (zero-fill outside [0,bits);
    // negative pos shifts the value left, i.e. bits below 0 are zeros)
    inline u128 get(i64 pos, int len) const {
        if (len <= 0) return 0;
        if (pos < 0) {
            int sh = (int)(-pos);
            if (sh >= len) return 0;
            return get(0, len - sh) << sh;
        }
        if (pos >= bits) return 0;
        size_t lo = (size_t)(pos >> 3);
        int off = (int)(pos & 7);
        u128 v = 0;
        int nb = (len + off + 7) / 8;
        for (int k = 0; k < nb; k++) {
            size_t idx = lo + (size_t)k;
            if (idx >= b.size()) break;
            v |= ((u128)b[idx]) << (8 * k);
        }
        v >>= off;
        if (len < 128) v &= ((((u128)1) << len) - 1);
        return v;
    }
};

// compute 5^(n) by repeated multiply with 5^27 (< 2^63) word passes
static XString pow5_bytes(i64 n) {
    std::vector<u64> w(1, 1);
    i64 left = n;
    const u64 CH = 7450580596923828125ull; // 5^27
    auto mul_small = [&](u64 m) {
        u128 carry = 0;
        for (size_t i = 0; i < w.size(); i++) {
            u128 t = (u128)w[i] * m + carry;
            w[i] = (u64)t;
            carry = t >> 64;
        }
        while (carry) { w.push_back((u64)carry); carry >>= 64; }
    };
    while (left >= 27) { mul_small(CH); left -= 27; }
    u64 rest = 1;
    for (i64 i = 0; i < left; i++) rest *= 5;
    if (rest > 1) mul_small(rest);

    XString X;
    X.b.resize(w.size() * 8);
    memcpy(X.b.data(), w.data(), w.size() * 8);
    while (X.b.size() > 1 && X.b.back() == 0) X.b.pop_back();
    i64 nb = (i64)X.b.size();
    uint8_t top = X.b.back();
    int tb = 0;
    while (top) { tb++; top >>= 1; }
    X.bits = (nb - 1) * 8 + tb;
    return X;
}

static bool load_x5(const char *path, XString &X) {
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    X.b.resize((size_t)sz);
    if (fread(X.b.data(), 1, (size_t)sz, f) != (size_t)sz) { fclose(f); return false; }
    fclose(f);
    while (X.b.size() > 1 && X.b.back() == 0) X.b.pop_back();
    uint8_t top = X.b.back();
    int tb = 0;
    while (top) { tb++; top >>= 1; }
    X.bits = ((i64)X.b.size() - 1) * 8 + tb;
    return true;
}

// ---------------------------------------------------------------------------
// Horner reduction of X's bit range [lo, hi) (MSB first) modulo mu, Barrett
// ---------------------------------------------------------------------------
struct Barrett {
    u64 mu, inv;
    explicit Barrett(u64 m) : mu(m), inv(m > 1 ? ~(u64)0 / m : ~(u64)0) {}
    inline u64 reduce(u64 t) const {           // t < 2^64 -> t mod mu
        if (mu == 1) return 0;
        u64 q = (u64)(((u128)t * inv) >> 64);
        u64 r = t - q * mu;
        while (r >= mu) r -= mu;
        return r;
    }
};

static u64 horner_range(const XString &X, i64 lo, i64 hi, const Barrett &bar) {
    if (hi > X.bits) hi = X.bits;
    if (hi <= lo) return 0;
    u64 r = 0;
    i64 pos = hi;
    while (pos > lo) {
        int c = (int)((pos - lo) < 32 ? (pos - lo) : 32);
        pos -= c;
        u64 chunk = (u64)X.get(pos, c);
        r = bar.reduce((r << c) | chunk);
    }
    return r;
}

// ---------------------------------------------------------------------------
// The extractor
// ---------------------------------------------------------------------------
struct Result {
    std::string digits;
    double t_pre = 0, t_loop = 0;
    i64 mem_bytes = 0;
};

static Result extract_pi(i64 N, const XString *Xpre, int nthreads, int window = 12) {
    Result res;
    auto tp0 = std::chrono::steady_clock::now();
    XString Xlocal;
    const XString *X = Xpre;
    if (!X) { Xlocal = pow5_bytes(N - 1); X = &Xlocal; }
    res.mem_bytes = (i64)X->b.size();
    auto tp1 = std::chrono::steady_clock::now();
    res.t_pre = std::chrono::duration<double>(tp1 - tp0).count();

    const i64 n5 = N - 1;
    const i64 LX = X->bits;
    const i64 K = (i64)(((double)(N - 1) * 3.321928094887362 + G + 24) / 4) + 2;

    if (8 * K + 6 >= ((i64)1 << 62)) { fprintf(stderr, "N too large\n"); exit(1); }

#if defined(_OPENMP)
    if (nthreads > 0) omp_set_num_threads(nthreads);
    int T = nthreads > 0 ? nthreads : omp_get_max_threads();
#else
    int T = 1; (void)nthreads;
#endif

    std::vector<u128> partial((size_t)(T > 0 ? T : 1), 0);

#if defined(_OPENMP)
#pragma omp parallel
#endif
    {
#if defined(_OPENMP)
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        u128 acc = 0;
#if defined(_OPENMP)
#pragma omp for schedule(dynamic, 64)
#endif
        for (i64 k = 0; k <= K; k++) {
            i64 A = N - 1 - 4 * k;
            const i64 Es[4] = { A + 2, A - 1, A, A - 1 };
            const u64 mus[4] = { (u64)(8 * k + 1), (u64)(2 * k + 1),
                                 (u64)(8 * k + 5), (u64)(4 * k + 3) };
            const int sgn[4] = { +1, -1, -1, -1 };
            for (int t = 0; t < 4; t++) {
                i64 E = Es[t];
                u64 mu = mus[t];
                u128 f;
                if (E >= 0) {
                    if (mu == 1) continue;
                    u64 r = mulmod(powmod(2, (u64)E, mu), powmod(5, (u64)n5, mu), mu);
                    if (!r) continue;
                    f = q96(r, 0, mu);
                } else {
                    i64 J = -E;
                    if (J >= LX + G + 8) continue;
                    u64 h;
                    Barrett bar(mu);
                    if (LX - J <= J) {
                        h = horner_range(*X, J, LX, bar);      // floor(X/2^J) mod mu
                    } else {
                        u64 Yt = horner_range(*X, 0, J, bar);  // (X mod 2^J) mod mu
                        u64 R = powmod(5, (u64)n5, mu);
                        u64 inv2J = powmod((mu + 1) >> 1, (u64)J, mu);
                        h = mulmod((R + mu - Yt) % mu, inv2J, mu);
                    }
                    u128 uG = X->get(J - GP, GP);
                    f = q96(h, uG << (G - GP), mu);
                }
                if (sgn[t] > 0) acc = (acc + f) & MASKG;
                else            acc = (acc + ((MASKG + 1) - (f & MASKG))) & MASKG;
            }
        }
        partial[(size_t)tid] = acc;
    }

    u128 acc = 0;
    for (u128 p : partial) acc = (acc + p) & MASKG;

    auto tp2 = std::chrono::steady_clock::now();
    res.t_loop = std::chrono::duration<double>(tp2 - tp1).count();

    std::string d;
    for (int i = 0; i < window; i++) {
        acc *= 10;
        d.push_back((char)('0' + (int)(acc >> G)));
        acc &= MASKG;
    }
    res.digits = d;
    return res;
}

// ---------------------------------------------------------------------------
static int run_check(const char *digitsfile, int nthreads) {
    FILE *f = fopen(digitsfile, "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", digitsfile); return 1; }
    std::string ref;
    int c;
    while ((c = fgetc(f)) != EOF) if (c >= '0' && c <= '9') ref.push_back((char)c);
    fclose(f);
    printf("reference digits: %zu\n", ref.size());

    bool ok = true;
    for (i64 N = 1; N <= 60; N++) {
        Result r = extract_pi(N, nullptr, nthreads, 8);
        std::string want = ref.substr((size_t)(N - 1), 8);
        if (r.digits.substr(0, 8) != want) {
            printf("FAIL N=%lld got=%s want=%s\n", (long long)N, r.digits.c_str(), want.c_str());
            ok = false;
        }
    }
    printf("N=1..60: %s\n", ok ? "all OK" : "FAILURES");
    for (i64 N : { (i64)100, (i64)1000 }) {
        if ((size_t)(N + 8) > ref.size()) break;
        Result r = extract_pi(N, nullptr, nthreads, 8);
        std::string want = ref.substr((size_t)(N - 1), 8);
        bool g = r.digits.substr(0, 8) == want;
        ok = ok && g;
        printf("N=%-6lld got=%s want=%s  %s\n", (long long)N, r.digits.c_str(),
               want.c_str(), g ? "OK" : "FAIL");
    }
    printf("CHECK %s\n", ok ? "PASSED" : "FAILED");
    return ok ? 0 : 1;
}

int main(int argc, char **argv) {
    i64 N = 0;
    const char *x5file = nullptr;
    const char *checkfile = nullptr;
    int nthreads = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--x5") && i + 1 < argc) x5file = argv[++i];
        else if (!strcmp(argv[i], "--check") && i + 1 < argc) checkfile = argv[++i];
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else N = atoll(argv[i]);
    }
    if (checkfile) return run_check(checkfile, nthreads);
    if (N <= 0) {
        fprintf(stderr, "usage: %s <N> [--x5 file] [--threads T] | --check <digitsfile>\n", argv[0]);
        return 1;
    }
    XString X;
    XString *Xp = nullptr;
    if (x5file) {
        if (!load_x5(x5file, X)) { fprintf(stderr, "cannot load %s\n", x5file); return 1; }
        Xp = &X;
    }
    Result r = extract_pi(N, Xp, nthreads);
    printf("pi digits at %lld..%lld: %s\n", (long long)N, (long long)(N + 11), r.digits.c_str());
    printf("precompute %.3fs, loop %.3fs, total %.3fs, aux mem %.2f MB\n",
           r.t_pre, r.t_loop, r.t_pre + r.t_loop, (double)r.mem_bytes / 1e6);
    return 0;
}
