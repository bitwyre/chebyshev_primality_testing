/*
 * By Olexa Bilaniuk
 */
/* Include Guards */
#ifndef INTEGERFACTORING_H
#define INTEGERFACTORING_H


/* Includes */
#include <stdio.h>
#include <stdint.h>

/* Defines */


/* C++ Extern "C" Guard */
#ifdef __cplusplus
extern "C" {
#endif

/* Functions */

/**
 * @brief Checks whether an integer is prime.
 *
 * @param [in] n   The integer whose primality is to be checked.
 * @return 1 if prime; 0 if not prime.
 *
 * NB: This is *not* a probabilistic primality checker. For all integers it can
 *     be given as input, it will correctly report "prime" or "composite".
 * NB: Internally, this function uses the Miller-Rabin test, which *is*
 *     probabilistic, and may falsely report a number as prime when in fact it
 *     is composite. However, this function uses a deterministic set of
 *     Miller-Rabin "witnesses", which ensures that there are no strong
 *     probable primes equal to or below 2^64-1 (the size of the input
 *     argument). This set of witnesses is
 *
 *         $$a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37$$
 *
 *     See https://oeis.org/A014233
 */

int      gaIIsPrime(uint64_t n);
/* End C++ Extern "C" Guard */
#ifdef __cplusplus
}
#endif


/* End Include Guards */
#endif
