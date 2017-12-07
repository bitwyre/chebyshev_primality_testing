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

/**
 * @brief Count trailing zeros of a 64-bit integer.
 *
 * @param [in] n  The integer whose trailing zero count is to be computed.
 * @return     If n != 0, returns trailing zero count; Else returns 64.
 */

int      gaICtz(uint64_t n);

/**
 * @brief Count leading zeros of a 64-bit integer.
 *
 * @param [in] n  The integer whose leading zero count is to be computed.
 * @return     If n != 0, returns leading zero count; Else returns 64.
 */

int      gaIClz(uint64_t n);

/**
 * @brief Integer Modular Addition.
 *
 * Computes
 *
 *     $$a+b \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

uint64_t gaIAddMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Subtraction.
 *
 * Computes
 *
 *     $$a-b \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

uint64_t gaISubMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Average.
 *
 * Computes
 *
 *     $$\frac{a+b}{2} \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

uint64_t gaIAvgMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Multiplication.
 *
 * Computes
 *
 *     $$a*b \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

uint64_t gaIMulMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Exponentiation.
 *
 * Computes
 *
 *     $$x^a \pmod m$$
 *
 * efficiently for 64-bit unsigned integers x, a, m.
 */

uint64_t gaIPowMod    (uint64_t x, uint64_t a, uint64_t m);

/**
 * @brief Jacobi Symbol
 *
 * Computes the Jacobi symbol, notated
 *
 *     $$(a/n)$$
 *
 * efficiently for 64-bit unsigned integers a, n.
 */

int      gaIJacobiSymbol(uint64_t a, uint64_t n);

/**
 * @brief Strong Fermat base-a probable prime test.
 *
 * @param [in] n  An odd integer >= 3.
 * @param [in] a  A witness integer > 0.
 * @return Non-zero if n is a strong probable prime to base a and zero if n is
 *         composite.
 */

int      gaIIsPrimeStrongFermat(uint64_t n, uint64_t a);

/**
 * @brief Strong Lucas probable prime test.
 *
 * The function uses Selfridge's Method A for selecting D,P,Q.
 *
 * @param [in] n  An odd integer >= 3.
 * @return Non-zero if n is a strong probable prime and zero if n is composite.
 */
/* End C++ Extern "C" Guard */
#ifdef __cplusplus
}
#endif


/* End Include Guards */
#endif
