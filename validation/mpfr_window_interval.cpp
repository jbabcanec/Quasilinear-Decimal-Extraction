// Independent verification only: directed-rounding MPFR prefix enclosure.
// This program is never called from the no-prefix extraction algorithms.
#include <gmpxx.h>
#include <mpfr.h>
#include <cstdio>
#include <stdexcept>
#include <string>

int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::invalid_argument("usage: mpfr_window_interval POSITION");
        std::size_t used=0;
        const std::string arg=argv[1];
        const unsigned long N=std::stoul(arg,&used);
        if(used!=arg.size() || N<1 || N>10000000 || arg[0]=='-')
            throw std::invalid_argument("position must be in 1..10000000");
        mpz_class scale,low_integer,high_integer,modulus;
        mpz_ui_pow_ui(scale.get_mpz_t(),10,N+11);
        mpz_ui_pow_ui(modulus.get_mpz_t(),10,12);
        const mp_bitcnt_t scale_bits=mpz_sizeinbase(scale.get_mpz_t(),2);
        for(unsigned guard=128;guard<=4096;guard*=2) {
            mpfr_t lo,hi;
            mpfr_inits2(scale_bits+guard,lo,hi,(mpfr_ptr)0);
            mpfr_const_pi(lo,MPFR_RNDD);
            mpfr_const_pi(hi,MPFR_RNDU);
            mpfr_mul_z(lo,lo,scale.get_mpz_t(),MPFR_RNDD);
            mpfr_mul_z(hi,hi,scale.get_mpz_t(),MPFR_RNDU);
            mpfr_get_z(low_integer.get_mpz_t(),lo,MPFR_RNDD);
            mpfr_get_z(high_integer.get_mpz_t(),hi,MPFR_RNDD);
            mpfr_clears(lo,hi,(mpfr_ptr)0);
            if(low_integer==high_integer) {
                low_integer%=modulus;
                gmp_printf("{\"N\":%lu,\"digits\":\"%012Zd\",\"guard_bits\":%u,"
                           "\"method\":\"MPFR directed-rounding enclosure\"}\n",
                           N,low_integer.get_mpz_t(),guard);
                return 0;
            }
        }
        throw std::runtime_error("reference interval crosses an integer boundary");
    } catch(const std::exception& error) {
        std::fprintf(stderr,"error: %s\n",error.what());
        return 1;
    }
}
