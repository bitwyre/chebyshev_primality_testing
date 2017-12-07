# Chebyshev polynomials of the first kind and primality testing

Claim:

If r is a prime number that does not divide n and if

(X − 1)^n = X^n − 1 (mod X^r − 1, n),

then either n is prime or n^2 = 1 (mod r).

If this conjecture is true, we can modify the algorithm slightly to first search for an r which does
not divide n^2 − 1. Such an r can assuredly be found in the range [2, 4 log n]. This is because the
product of prime numbers less than x is at least e^x (see [Apo97]). Thereafter we can test whether
the congruence (6) holds or not. Verifying the congruence takes time O∼(r log^2 n). This gives a time
complexity of O∼(log^3 n).

Recently, Hendrik Lenstra and Carl Pomerance [LP03b] have given a heuristic argument that suggests
that the above conjecture is false. However, some variant of the conjecture may still be true (for example,
if we force r > log n).

See [Chebyshev polynomials of the first kind and primality testing](https://mathoverflow.net/questions/286304/chebyshev-polynomials-of-the-first-kind-and-primality-testing) for details.
