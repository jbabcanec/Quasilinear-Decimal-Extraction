// decimal_extract_v3.cpp
// Quasi-linear decimal digit extraction of pi via the division-free 2-adic
// pipeline (engine "v3"), a C++/GMP port of decimal_bbp_extract.py::extract_v3,
// which is verified against public digit corpora through decimal position 1e7.
//
// Build (needs GMP + its C++ wrapper):
//     g++ -O3 -march=native decimal_extract_v3.cpp -lgmpxx -lgmp -o extract_v3
//
// Run:
//     ./extract_v3 1000000       # 12 decimal digits of pi at positions N..N+11
//
// Validated against the verified Python reference and public digit corpora:
//     N = 1000000  -> 130927562832   (~3.9 s)
//     N = 10000000 -> 725915133612   (~54 s, vs 111.6 s for the Python reference)
// single-threaded on one desktop core; both phases scale quasi-linearly.
// Set the environment variable V3_TIMING=1 to print per-phase timing on stderr.
//
// Limits: uses a 32-bit exponent for 5^(N-1), so N < ~4.29e9 (well past the
// ~1e9 memory ceiling of a single machine). Word moduli mu stay < 2^64.

#include <gmpxx.h>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
using namespace std;
static std::chrono::steady_clock::time_point _T0;
static const bool _timing = std::getenv("V3_TIMING") != nullptr;   // set V3_TIMING=1 for per-phase timing
static void lap(const char* s){ if(!_timing) return; auto n=std::chrono::steady_clock::now();
    fprintf(stderr,"  [%-6s %6.2fs]\n", s, std::chrono::duration<double>(n-_T0).count()); _T0=n; }

typedef uint64_t u64;
typedef unsigned __int128 u128;

static inline u64 mulmod(u64 a, u64 b, u64 m) { return (u128)a * b % m; }
static u64 powmod(u64 b, u64 e, u64 m) {
    u64 r = 1 % m; b %= m;
    while (e) { if (e & 1) r = mulmod(r, b, m); b = mulmod(b, b, m); e >>= 1; }
    return r;
}
// mpz from a full 64-bit value (gmpxx has no unsigned-long-long ctor on LLP64)
static inline mpz_class from_u64(u64 v) {
    mpz_class r; mpz_import(r.get_mpz_t(), 1, 1, sizeof(u64), 0, 0, &v); return r;
}

// a node of the division-free merge: value = 2^shift * num/den, kept mod 2^(Jmax-shift)
struct Node { mpz_class num, den; long long shift; };

int main(int argc, char** argv) {
    const long long N  = (argc > 1) ? atoll(argv[1]) : 1000000;
    const int  G = 192, window = 12;           // accumulator bits, output digits
    const long long n5 = N - 1;
    const long long K  = (long long)(((double)n5 * log2(10.0) + G + 24) / 4) + 2;

    // shared string X = 5^(N-1)
    _T0 = std::chrono::steady_clock::now();
    mpz_class X;
    mpz_ui_pow_ui(X.get_mpz_t(), 5, (unsigned long)n5);
    const long long LX  = (long long)mpz_sizeinbase(X.get_mpz_t(), 2);
    const long long Jcut = LX + G + 8;
    lap("X");

    mpz_class acc = 0;
    mpz_class MOD = mpz_class(1) << G;          // 2^G

    // pi = sum 16^-k [ 4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6) ], shifted by 10^(N-1):
    // each family gives mu = a*k + b (odd), exponent E = e0 - 4k, sign +/-.
    struct Fam { int sign; u64 a, b; long long e0; };
    const Fam fams[4] = {
        {  1, 8, 1, n5 + 2 },   //  4/(8k+1)         = 2^(A+2) X/(8k+1)
        { -1, 2, 1, n5 - 1 },   // -2/(8k+4)         = -2^(A-1) X/(2k+1)
        { -1, 8, 5, n5     },   // -1/(8k+5)         = -2^A     X/(8k+5)
        { -1, 4, 3, n5 - 1 },   // -1/(8k+6)         = -2^(A-1) X/(4k+3)
    };

    struct Band { int sign; long long J; u64 mu; };
    vector<Band> band;

    mpz_class t;
    for (const Fam& f : fams) {
        for (long long k = 0; k <= K; k++) {
            const long long E  = f.e0 - 4 * k;
            const u64       mu = f.a * (u64)k + f.b;
            if (E >= 0) {                                   // easy term
                if (mu == 1) continue;
                const u64 r = mulmod(powmod(2, (u64)E, mu), powmod(5, (u64)n5, mu), mu);
                if (r) {
                    mpz_class mz = from_u64(mu);
                    t = from_u64(r) << G;
                    mpz_fdiv_q(t.get_mpz_t(), t.get_mpz_t(), mz.get_mpz_t());   // (r<<G)/mu
                    if (f.sign > 0) acc += t; else acc -= t;
                }
            } else {                                        // band term
                const long long J = -E;
                if (J < Jcut) {
                    band.push_back({ f.sign, J, mu });
                    if (mu > 1) {                            // Bezout odd half -> word phase
                        const u64 R  = powmod(5, (u64)n5, mu);
                        const u64 iv = powmod((mu + 1) >> 1, (u64)J, mu);   // 2^-J mod mu
                        const u64 ro = mulmod(R, iv, mu);
                        if (ro) {
                            mpz_class mz = from_u64(mu);
                            t = from_u64(ro) << G;
                            mpz_fdiv_q(t.get_mpz_t(), t.get_mpz_t(), mz.get_mpz_t());
                            if (f.sign > 0) acc += t; else acc -= t;
                        }
                    }
                }
            }
        }
    }
    mpz_mod(acc.get_mpz_t(), acc.get_mpz_t(), MOD.get_mpz_t());
    lap("word");

    // ---- V phase: consolidate the dyadic halves into one 2-adic integer ----
    if (!band.empty()) {
        long long Jmax = 0;
        for (const Band& b : band) Jmax = max(Jmax, b.J);
        sort(band.begin(), band.end(), [](const Band& x, const Band& y) { return x.J < y.J; });

        auto merge = [&](Node a, Node b) -> Node {
            const long long s   = min(a.shift, b.shift);
            const long long cap = Jmax - s;
            mpz_class n_ = ((a.num * b.den) << (unsigned long)(a.shift - s))
                         + ((b.num * a.den) << (unsigned long)(b.shift - s));
            mpz_class d_ = a.den * b.den;
            // mask only when the magnitude exceeds cap (as in the Python reference):
            // masking a small negative would inflate it to 2^cap - |n|, a cap-bit giant,
            // and those propagate through the tree, turning the merge quadratic.
            if ((long long)mpz_sizeinbase(n_.get_mpz_t(), 2) > cap)
                mpz_fdiv_r_2exp(n_.get_mpz_t(), n_.get_mpz_t(), (mp_bitcnt_t)cap);
            if ((long long)mpz_sizeinbase(d_.get_mpz_t(), 2) > cap)
                mpz_fdiv_r_2exp(d_.get_mpz_t(), d_.get_mpz_t(), (mp_bitcnt_t)cap);
            return Node{ std::move(n_), std::move(d_), s };
        };
        auto reduce_nodes = [&](vector<Node> nodes) -> Node {
            while (nodes.size() > 1) {
                vector<Node> nxt; nxt.reserve((nodes.size() + 1) / 2);
                size_t i = 0;
                for (; i + 1 < nodes.size(); i += 2)
                    nxt.push_back(merge(std::move(nodes[i]), std::move(nodes[i + 1])));
                if (nodes.size() & 1) nxt.push_back(std::move(nodes.back()));
                nodes = std::move(nxt);
            }
            return std::move(nodes[0]);
        };

        const size_t BLK = 1u << 15;
        vector<Node> roots;
        for (size_t lo = 0; lo < band.size(); lo += BLK) {
            const size_t hi = min(lo + BLK, band.size());
            vector<Node> leaves; leaves.reserve(hi - lo);
            for (size_t i = lo; i < hi; i++) {
                Node nd;
                nd.num   = band[i].sign;                    // +/-1 (residue taken in merge)
                nd.den   = from_u64(band[i].mu);
                nd.shift = Jmax - band[i].J;
                leaves.push_back(std::move(nd));
            }
            roots.push_back(reduce_nodes(std::move(leaves)));
        }
        Node root = reduce_nodes(std::move(roots));
        const long long cap = Jmax - root.shift;

        // Newton-Hensel inverse of the odd denominator, mod 2^cap
        mpz_class inv = 1;
        for (long long prec = 1; prec < cap; ) {
            prec = min(2 * prec, cap);
            inv = inv * (mpz_class(2) - root.den * inv);
            mpz_fdiv_r_2exp(inv.get_mpz_t(), inv.get_mpz_t(), (mp_bitcnt_t)prec);
        }
        mpz_class V = root.num * inv;
        mpz_fdiv_r_2exp(V.get_mpz_t(), V.get_mpz_t(), (mp_bitcnt_t)cap);
        V <<= (unsigned long)root.shift;
        mpz_fdiv_r_2exp(V.get_mpz_t(), V.get_mpz_t(), (mp_bitcnt_t)Jmax);

        mpz_class W = X * V;
        mpz_fdiv_r_2exp(W.get_mpz_t(), W.get_mpz_t(), (mp_bitcnt_t)Jmax);   // X*V mod 2^Jmax

        mpz_class dy;
        if (Jmax >= G) dy = W >> (unsigned long)(Jmax - G);                 // top G bits
        else           dy = W << (unsigned long)(G - Jmax);
        acc += dy;
        mpz_mod(acc.get_mpz_t(), acc.get_mpz_t(), MOD.get_mpz_t());
    }
    lap("V");

    // digits = floor(acc * 10^window / 2^G), zero-padded to `window` places
    mpz_class out = acc, tp;
    mpz_ui_pow_ui(tp.get_mpz_t(), 10, window);
    out *= tp;
    out >>= G;
    char fmt[16]; snprintf(fmt, sizeof(fmt), "%%0%dZd", window);
    gmp_printf(fmt, out.get_mpz_t());
    gmp_printf("\n");
    return 0;
}
