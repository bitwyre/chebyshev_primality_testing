/* 
 * By Olexa Bilaniuk
 */

/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "integer-factoring-baseline.h"


/* Detect when to avoid VLAs. */
#if defined(_MSC_VER) || defined(__STDC_NO_VLA__)
#define GA_USING_MALLOC_FOR_VLA 1
#endif


/* Defines */
#define GA_IS_COMPOSITE      0
#define GA_IS_PRIME          1
#define GA_IS_PROBABLY_PRIME 2


/**
 * Static Function Prototypes
 */

/**
 * @brief Count trailing zeros of a 64-bit integer.
 *
 * @param [in] n  The integer whose trailing zero count is to be computed.
 * @return     If n != 0, returns trailing zero count; Else returns 64.
 */

static int      gaICtz(uint64_t n);

/**
 * @brief Count leading zeros of a 64-bit integer.
 *
 * @param [in] n  The integer whose leading zero count is to be computed.
 * @return     If n != 0, returns leading zero count; Else returns 64.
 */

static int      gaIClz(uint64_t n);

/**
 * @brief Integer Modular Addition.
 *
 * Computes
 *
 *     $$a+b \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

static uint64_t gaIAddMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Subtraction.
 *
 * Computes
 *
 *     $$a-b \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

static uint64_t gaISubMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Average.
 *
 * Computes
 *
 *     $$\frac{a+b}{2} \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

static uint64_t gaIAvgMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Multiplication.
 *
 * Computes
 *
 *     $$a*b \pmod m$$
 *
 * efficiently for 64-bit unsigned integers a, b, m.
 */

static uint64_t gaIMulMod    (uint64_t a, uint64_t b, uint64_t m);

/**
 * @brief Integer Modular Exponentiation.
 *
 * Computes
 *
 *     $$x^a \pmod m$$
 *
 * efficiently for 64-bit unsigned integers x, a, m.
 */

static uint64_t gaIPowMod    (uint64_t x, uint64_t a, uint64_t m);

/**
 * @brief Jacobi Symbol
 *
 * Computes the Jacobi symbol, notated
 *
 *     $$(a/n)$$
 *
 * efficiently for 64-bit unsigned integers a, n.
 */

static int      gaIJacobiSymbol(uint64_t a, uint64_t n);

/**
 * @brief Strong Fermat base-a probable prime test.
 *
 * @param [in] n  An odd integer >= 3.
 * @param [in] a  A witness integer > 0.
 * @return Non-zero if n is a strong probable prime to base a and zero if n is
 *         composite.
 */

static int      gaIIsPrimeStrongFermat(uint64_t n, uint64_t a);

/**
 * @brief Strong Lucas probable prime test.
 *
 * The function uses Selfridge's Method A for selecting D,P,Q.
 *
 * @param [in] n  An odd integer >= 3.
 * @return Non-zero if n is a strong probable prime and zero if n is composite.
 */

/**
 * Function Definitions
 */

static int      gaICtz       (uint64_t n){
#if __GNUC__ >= 4
	return n ? __builtin_ctzll(n) : 64;
#else
	int z;

	for(z=0;z<64;z++){
		if((n>>z) & 1){break;}
	}

	return z;
#endif
}

static int      gaIClz       (uint64_t n){
#if __GNUC__ >= 4
	return n ? __builtin_clzll(n) : 64;
#else
	int z;

	for(z=63;z>=0;z--){
		if((n>>z) & 1){break;}
	}

	return 63-z;
#endif
}

static uint64_t gaIAddMod    (uint64_t a, uint64_t b, uint64_t m){
	a %= m;
	b %= m;

	if(m-a > b){
		return a+b;
	}else{
		return a+b-m;
	}
}

static uint64_t gaISubMod    (uint64_t a, uint64_t b, uint64_t m){
	a %= m;
	b %= m;

	if(a >= b){
		return a-b;
	}else{
		return a-b+m;
	}
}

static uint64_t gaIAvgMod    (uint64_t a, uint64_t b, uint64_t m){
	uint64_t s = gaIAddMod(a,b,m);

	if(s&1){
		return (s>>1)+(m>>1)+(s&m&1);
	}else{
		return s>>1;
	}
}

static uint64_t gaIMulMod    (uint64_t a, uint64_t b, uint64_t m){
#if (__GNUC__ >= 4) && defined(__x86_64__) && !defined(__STRICT_ANSI__)
	uint64_t r;

	asm(
	    "mul %2\n\t"
	    "div %3\n\t"
	    : "=&d"(r), "+a"(a)   /* Outputs */
	    : "r"(b),  "r"(m)     /* Inputs */
	    : "cc"
	);

	return r;
#elif ((__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)) && defined(__SIZEOF_INT128__) && __SIZEOF_INT128__ >= 16
	/* Hardcore GCC 4.6+ optimization jazz */
	return ((unsigned __int128)a * (unsigned __int128)b) % m;
#else
	const uint64_t TWOPOW32 = (uint64_t)1<<32;
	int i;

	a %= m;
	b %= m;

	if(m <= TWOPOW32){
		/**
		 * Fast path: When performing modulo arithmetic on values <= 2^32,
		 * (a*b) % m gives the correct answer.
		 */

		return (a*b) % m;
	}else{
		/**
		 * Slow path: Have to simulate 128-bit arithmetic long division.
		 */

		uint64_t ah   = a>>32;
		uint64_t al   = (uint32_t)a;
		uint64_t bh   = b>>32;
		uint64_t bl   = (uint32_t)b;

		uint64_t ahbh = ah*bh;
		uint64_t ahbl = ah*bl;
		uint64_t albh = al*bh;
		uint64_t albl = al*bl;

		uint64_t md   = ahbl+albh;

		uint64_t lo   = albl + (md<<32);
		uint64_t hi   = ahbh + (md>>32);

		/* Propagate carry-outs from `md` and `lo` into `hi` */
		if(lo < albl){hi++;}
		if(md < ahbl){hi+=TWOPOW32;}

		/**
		 * Begin 128-bit-by-64-bit remainder.
		 *
		 * 1) Cut down `hi` mod `m`. This implements the first few iterations
		 *    of a shift-and-subtract loop, leaving only 64 iterations to go.
		 * 2) Iterate 64 times:
		 *     2.1) Shift left [hi:lo] by 1 bit, into [newHi:newLo].
		 *     2.2) If:
		 *         2.2.1) newHi < hi, then there was an overflow into bit 128.
		 *                The value [1:newHi:newLo] is definitely larger than
		 *                m, so we subtract. This situation can only occur if
		 *                m > 2^63.
		 *         2.2.2) newHi > m, then we must subtract m out of newHi in
		 *                order to bring back newHi within the range [0, m).
		 * 3) The modulo is in hi.
		 */

		hi %= m;
		for(i=0;i<64;i++){
			uint64_t newLo = (lo<<1);
			uint64_t newHi = (hi<<1) + (newLo<lo);

			if(newHi < hi || newHi > m){newHi -= m;}

			hi = newHi;
			lo = newLo;
		}

		return hi;
	}
#endif
}

static uint64_t gaIPowMod    (uint64_t x, uint64_t a, uint64_t m){
	uint64_t r;

	/**
	 * Special cases (order matters!):
	 * - A modulo of 0 makes no sense and a modulo of 1 implies a return value
	 *   of 0, since the result must be integer.
	 * - An exponent of 0 requires a return value of 1.
	 * - A base of 0 or 1 requires a return value of 0 or 1.
	 * - An exponent of 1 requires a return value of x.
	 * - An exponent of 2 can be handled by the modulo multiplication directly.
	 */

	if(m<=1){
		return 0;
	}

	x %= m;

	if(a==0){
		return 1;
	}else if(x<=1){
		return x;
	}else if(a==1){
		return x;
	}else if(a==2){
		return gaIMulMod(x,x,m);
	}

	/**
	 * Otherwise, perform modular exponentiation by squaring.
	 */

	r = 1;
	while(a){
		if(a&1){
			r = gaIMulMod(r, x, m);
		}

		x = gaIMulMod(x, x, m);
		a >>= 1;
	}

	return r;
}

static int      gaIJacobiSymbol(uint64_t a, uint64_t n){
	int      s=0;
	uint64_t e, a1, n1;

	a %= n;

	if(a == 1 || n == 1){
		return 1;
	}

	if(a == 0){
		return 0;
	}

	e  = gaICtz(a);
	a1 = a >> e;

	if(e%2 == 0){
		s =  1;
	}else if(n%8 == 1 || n%8 == 7){
		s =  1;
	}else if(n%8 == 3 || n%8 == 5){
		s = -1;
	}

	if(n%4 == 3 && a1%4 == 3){
		s = -s;
	}

	n1 = n%a1;
	return s*gaIJacobiSymbol(n1,a1);
}

static int      gaIIsPrimeStrongFermat(uint64_t n, uint64_t a){
	/**
	 * The Fermat strong probable prime test the Miller-Rabin test relies upon
	 * uses integer "witnesses" in an attempt at proving the number composite.
	 * Should it fail to prove an integer composite, it reports the number as
	 * "probably prime". However, if the witnesses are chosen carefully, the
	 * Miller-Rabin test can be made deterministic below a chosen threshold.
	 *
	 * One can use the primes 2 to 37 in order to ensure the correctness of the
	 * identifications for integers under 2^64.
	 *
	 * Jim Sinclair has found that the seven witnesses
	 *     2, 325, 9375, 28178, 450775, 9780504, 1795265022
	 * also deterministically classify all integers <2^64.
	 *
	 *
	 * The Fermat strong probable prime test states that, for integers
	 *             n = d*2^s+1,  d odd, s integer >= 0
	 *             a             integer (chosen witness)
	 * n is a Fermat strong probable prime if
	 *     a^(d    ) =  1 mod n       or
	 *     a^(d*2^r) = -1 mod n       for any integer r, 0 <= r < s.
	 *
	 *
	 * The justification for this comes from Fermat's Little Theorem: If n is
	 * prime and a is any integer, then the following always holds:
	 *           a^n =  a mod n
	 * If n is prime and a is coprime to n, then the following always holds:
	 *       a^(n-1) =  1 mod n
	 *
	 *
	 * In effect, the logic goes
	 *
	 *   A:   The number  n  is prime.                               (Statement)
	 *   B:   The number  n  does not divide a.                      (Statement)
	 *   C:   a^(  n-1)       =  1 mod n                             (Statement)
	 *   D:   The commutative ring Z/nZ is a finite field.           (Statement)
	 *   E:   Finite fields are unique factorization domains.        (Statement)
	 *   F:   x^2 = 1 mod n factorizes as (x+1)(x-1) = 0 mod n.      (Statement)
	 *   G:   x^2 mod n only has the trivial square roots 1 and -1   (Statement)
	 *   H:   The number  n  is odd and >= 3.                        (Statement)
	 *   I:   The number n-1 equals d*2^s, with d,s int > 0, d odd.  (Statement)
	 *   J:   a^(    d)       =   1 mod n                            (Statement)
	 *   K:   a^(d*2^r)       =  -1 mod n   for some 0 <= r < s.     (Statement)
	 *   L:   a^(d*2^(r+1))   =   1 mod n   for some 0 <= r < s.     (Statement)
	 *   M:   a^(d*2^r)      != +-1 mod n   AND                      (Statement)
	 *        a^(d*2^(r+1))   =   1 mod n   for some 0 <= r < s.
	 *
	 *   A&B           -->  C                 (Proposition:     Fermat's Little Theorem)
	 *   !C            -->  !(A&B) = !A|!B    (Contrapositive:  Fermat's Little Theorem)
	 *   A             <->  D                 (Proposition)
	 *   E                                    (Proposition:     By definition)
	 *   F                                    (Proposition:     x^2-x+x-1 = x^2-1 mod n)
	 *   D&E&F         -->  G                 (Proposition:     (x+1)(x-1) is the only
	 *                                                           factorization)
	 *   !G            -->  !D|!E|!F          (Contrapositive:  See above)
	 *   H&I&J         -->  C                 (Proposition:     Squaring  1 gives 1)
	 *   H&I&K         -->  L                 (Proposition:     Squaring -1 gives 1)
	 *   H&I&L         -->  C                 (Proposition:     1, squared or not, gives 1)
	 *   H&I&K         -->  C                 (Hypothetical Syllogism)
	 *   H&I&(J|K)     -->  C                 (Union)
	 *   H&I&!(J|K)    -->  M|!C              (Proposition:     Either squaring
	 *                                                            a^(d*2^(s-1)) != +-1 mod n
	 *                                                          gives a 1, in which case
	 *                                                          M holds, or it does not
	 *                                                          give 1 and therefore
	 *                                                            a^(n-1) != 1 mod n)
	 *                                                          and thus !C holds.
	 *   H&I&!(J|K)    -->  H&I&M | !A | !B   (Absorbtion, Hypothetical Syllogism)
	 *   H&I&M         -->  !G                (Proposition:     x^2 = 1 mod n but x!=+1,
	 *                                                          so x^2 - 1 has roots
	 *                                                          other than +-1)
	 *   H&I&M         -->  !D|!E|!F          (Modus Tollens)
	 *   H&I&M         -->  !D                (Disjunctive Syllogism)
	 *   H&I&M         -->  !A                (Biconditional)
	 *   H&I&!(J|K)    -->  !A | !A | !B      (Hypothethical Syllogism)
	 *   H&I&!(J|K)&B  -->  !A | !A           (Absorbtion)
	 *   H&I&!(J|K)&B  -->  !A | !A           (Disjunctive Syllogism)
	 *   H&I&!(J|K)&B  -->  !A                (Disjunctive Simplification)
	 *                           ***** Conclusions: *****
	 *                            H&I&M         -->  !A
	 *                            H&I&!(J|K)&B  -->  !A
	 *
	 * Broadly speaking, what the above tells us is:
	 *   - We can't prove n prime (A), but we can prove it composite (!A).
	 *   - Either H&I&M or H&I&!(J|K)&B prove compositeness.
	 *   - If H&I&(J|K) for any r, then we've proven C true. If we prove C true,
	 *     we can't use the contrapositive of Fermat's Little Theorem, so no
	 *     conclusions about the truth-value of A can be made. The test is
	 *     inconclusive. Thus this function returns "probably prime".
	 */

	uint64_t d, x;
	int64_t  s, r;

	a %= n;
	if(a==0){
		return GA_IS_PROBABLY_PRIME;
	}

	s  = gaICtz(n-1);
	d  = (n-1) >> s;
	x  = gaIPowMod(a,d,n);

	if(x==1 || x==n-1){
		return GA_IS_PROBABLY_PRIME;
	}

	for(r=0;r<s-1;r++){
		x = gaIMulMod(x,x,n);
		if(x==1){
			return GA_IS_COMPOSITE;
		}else if(x == n-1){
			return GA_IS_PROBABLY_PRIME;
		}
	}

	return GA_IS_COMPOSITE;
}

static int      gaIIsPrimeStrongLucas(uint64_t n){
	uint64_t Dp, Dm, D, K, U, Ut, V, Vt;
	int      J, r, i;

	/**
	 * FIPS 186-4 C.3.3 (General) Lucas Probabilistic Primality Test
	 *
	 * 1. Test if n is perfect square. If so, return "composite".
	 *
	 *     NOTE: The only strong base-2 Fermat pseudoprime squares are
	 *           1194649 and 12327121;
	 */

	if(n==1194649 || n==12327121){
		return GA_IS_COMPOSITE;
	}

	/**
	 * 2. Find first D in sequence 5,-7,9,-11,... s.t. Jacobi symbol (D/n) < 1.
	 *     Iff Jacobi symbol is 0, return "composite".
	 */

	Dp = gaIAddMod(0, 5, n);
	Dm = gaISubMod(0, 7, n);
	while(1){
		J = gaIJacobiSymbol(Dp, n);
		if     (J ==  0){return GA_IS_COMPOSITE;}
		else if(J == -1){D = Dp;break;}

		J = gaIJacobiSymbol(Dm, n);
		if     (J ==  0){return GA_IS_COMPOSITE;}
		else if(J == -1){D = Dm;break;}

		Dp = gaIAddMod(Dp, 4, n);
		Dm = gaISubMod(Dm, 4, n);
	}

	/**
	 * 3. K = n+1
	 *
	 *     NOTE: Cannot overflow, since 2^64-1 is eliminated by strong Fermat
	 *           base-2 test.
	 */

	K = n+1;

	/**
	 * 4. Let Kr, Kr–1, ..., K0 be the binary expansion of K, with Kr = 1.
	 */

	r = 63-gaIClz(K);

	/**
	 * 5. Set Ur = 1 and Vr = 1.
	 */

	U = V = 1;

	/**
	 * 6. For i=r–1 to 0, do
	 */

	for(i=r-1;i>=0;i--){
		Ut = gaIMulMod(U,V,n);
		Vt = gaIAvgMod(gaIMulMod(V,V,n), gaIMulMod(D,gaIMulMod(U,U,n),n), n);
		if((K>>i)&1){
			U = gaIAvgMod(Ut,Vt,n);
			V = gaIAvgMod(Vt,gaIMulMod(D,Ut,n),n);
		}else{
			U = Ut;
			V = Vt;
		}
	}

	/**
	 * 7. If U0==0, then return "probably prime". Otherwise, return "composite".
	 */

	return U==0 ? GA_IS_PROBABLY_PRIME : GA_IS_COMPOSITE;
}

int      gaIIsPrime   (uint64_t n){
	int            hasNoSmallFactors, hasSmallFactors;

	/**
	 * Check if it is 2, the oddest prime.
	 */

	if(n==2){return GA_IS_PRIME;}

	/**
	 * Check if it is an even integer.
	 */

	if((n&1) == 0){return GA_IS_COMPOSITE;}

	/**
	 * For small integers, read directly the answer in a table.
	 */

	if(n<256){
        /*
         * This is a lookup table for the first 256 integer composites and primes
         * n means it's a composite and y means it's a prime
         */

		return "nnyynynynnnynynnnynynnnynnnnnyny"
		       "nnnnnynnnynynnnynnnnnynnnnnynynn"
		       "nnnynnnynynnnnnynnnynnnnnynnnnnn"
		       "nynnnynynnnynynnnynnnnnnnnnnnnny"
		       "nnnynnnnnynynnnnnnnnnynynnnnnynn"
		       "nnnynnnynnnnnynnnnnynynnnnnnnnny"
		       "nynnnynynnnnnnnnnnnynnnnnnnnnnny"
		       "nnnynynnnynnnnnynynnnnnnnnnynnnn"[n] == 'y';
	}

	/**
	 * Test small prime factors.
	 */

	hasNoSmallFactors = n% 3 && n% 5 && n% 7 && n%11 && n%13 && n%17 && n%19 &&
	                    n%23 && n%29 && n%31 && n%37 && n%41 && n%43 && n%47 &&
	                    n%53 && n%59 && n%61 && n%67 && n%71 && n%73 && n%79;
	hasSmallFactors   = !hasNoSmallFactors;
	if(hasSmallFactors){
		return GA_IS_COMPOSITE;
	}

	/**
	 * We implement the Baillie-Pomerance-Selfridge-Wagstaff primality checker.
	 *   1) A Fermat base-2 strong probable prime that is also
	 *   2) A Lucas strong probable prime is
	 *   3) Prime.
	 * The BPSW test has no known failure cases and is proven to have no failures
	 * for all numbers under 2^64. It is expected to have failures (composites
	 * classified as "probably prime") but they are expected to be enormous.
	 *
	 * We begin with the Fermat base-2 strong primality test
	 * (Miller-Rabin test with one witness only, a=2).
	 */

	return gaIIsPrimeStrongFermat(n,          2) &&

	/**
	 * Assuming this is one of the base-2 Fermat strong probable primes, we run
	 * the Lucas primality test with Selfridge's Method A for selecting D.
	 */

	       gaIIsPrimeStrongLucas (n            );
}

