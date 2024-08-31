// Ed Callaghan
// A collection of routines using gmp functions and types
// August 2024

#ifndef GMPRoutines_hh
#define GMPRoutines_hh

#include <gmp.h>

namespace mu2e{
  namespace gmp{
    // test if x is a quadratic resiue of N = p * q
    bool is_quadratic_residue(mpz_t x,
                              const mpz_t& p,
                              const mpz_t& q,
                              const mpz_t& N);

    // identify a non-quadratic residue of N = p * q
    void sample_nonresidue(mpz_t x,
                          const mpz_t& p,
                          const mpz_t& q,
                          const mpz_t& N,
                          gmp_randstate_t state);

    // identify x s.t. x == 1 (mod N)
    void sample_unit_modulo(mpz_t x,
                            const mpz_t& N,
                            gmp_randstate_t state);

    // apply Goldwasser-Micali encryption to a single bit
    void gm_encrypt(mpz_t encrypted,
                    unsigned int set,
                    const mpz_t& x,
                    const mpz_t& N,
                    gmp_randstate_t state);

    // apply Goldwasser-Micali decryption for a single bit; N = p * q
    unsigned int gm_decrypt(mpz_t encrypted,
                            const mpz_t& p,
                            const mpz_t& q,
                            const mpz_t& N);
  } // namespace gmp
} // namespace mu2e

#endif
