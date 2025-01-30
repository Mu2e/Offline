// Ed Callaghan
// A collection of routines using gmp functions and types
// August 2024

#include "Offline/Blinding/inc/GMPRoutines.hh"

namespace mu2e{
  namespace gmp{
    // test if x is a quadratic resiue of N = p * q
    bool is_quadratic_residue(mpz_t x,
                              const mpz_t& p,
                              const mpz_t& q,
                              const mpz_t& N){
        mpz_t gcd;
        mpz_init(gcd);

        mpz_gcd(gcd, x, N);
        int lp = mpz_legendre(x, p);
        int lq = mpz_legendre(x, q);

        bool rv = true;
        rv &= (mpz_cmp_si(gcd, 1) == 0);
        rv &= (lp == 1);
        rv &= (lq == 1);

        mpz_clear(gcd);
        return rv;
    }

    // identify a non-quadratic residue of N = p * q
    void sample_nonresidue(mpz_t x,
                          const mpz_t& p,
                          const mpz_t& q,
                          const mpz_t& N,
                          gmp_randstate_t state){
        mpz_set(x, N);

        do{
            mpz_urandomm(x, state, N);
        } while (is_quadratic_residue(x, p, q, N));
    }

    // identify x s.t. x == 1 (mod N)
    void sample_unit_modulo(mpz_t x,
                            const mpz_t& N,
                            gmp_randstate_t state){
        mpz_set(x, N);

        mpz_t gcd;
        mpz_init_set_ui(gcd, 0);
        while (mpz_cmp_ui(gcd, 1) != 0){
            mpz_urandomm(x, state, N);
            mpz_gcd(gcd, x, N);
        }
        mpz_clear(gcd);
    }

    // apply Goldwasser-Micali encryption to a single bit
    void gm_encrypt(mpz_t encrypted,
                    unsigned int set,
                    const mpz_t& x,
                    const mpz_t& N,
                    gmp_randstate_t state){
        mpz_t y;
        mpz_init(y);
        sample_unit_modulo(y, N, state);

        mpz_t rest;
        mpz_init(rest);
        mpz_pow_ui(rest, x, set);
        mpz_mul(encrypted, y, y);
        mpz_mul(encrypted, encrypted, rest);
        mpz_mod(encrypted, encrypted, N);

        mpz_clears(y, rest, NULL);
    }

    // apply Goldwasser-Micali decryption for a single bit; N = p * q
    unsigned int gm_decrypt(mpz_t encrypted,
                            const mpz_t& p,
                            const mpz_t& q,
                            const mpz_t& N){
        unsigned int rv = 1;
        if (is_quadratic_residue(encrypted, p, q, N)){
            rv = 0;
        }
        return rv;
    }
  } // namespace gmp
} // namespace mu2e
