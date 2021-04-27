#include <iostream>
#include <vector>
#include "../include/benchmark.h"
#include "../include/primality-test-baseline.h"

using namespace std;

typedef struct polynomial
{
    polynomial(uint64_t r, uint64_t n) : n(n), p(r) {}
    uint64_t          n;
    std::vector <uint64_t> p;

    // implementation using Galois field
    // Galoid field (x^r)^2x2
    // Finite field arithmetic for lookup
    polynomial operator+ (const polynomial& other) {
        uint64_t r = p.size();
        polynomial ret(r, n);
        for (int i=0; i< this->p.size(); i++) {
            ret.p[i] = gaIAddMod(this->p[i], other.p[i], this->n);
        }
        return ret;
    }
    
    // this is only there 
    // because the finite field is the polynomial
    // is mod x^r -1
    polynomial operator* (const polynomial& other) {
        uint64_t r = this->p.size();
        polynomial ret(r, n);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < r; j++) {
                 ret.p[(i+j) % r] = gaIAddMod(ret.p[(i+j) % r], gaIMulMod(this->p[i], other.p[j], this->n), this->n);
            }
        }
        return ret;
    }
} polynomial;


typedef struct matrix
{
    matrix(uint64_t r, uint64_t n) : n(n), p00(r,n), p01(r,n), p10(r,n), p11(r,n) {}
    
    uint64_t   n;
    polynomial p00;
    polynomial p01;
    polynomial p10;
    polynomial p11;

    // | p00 p01 | * | q00 q01 | = | p00*q00+p01*q10 p00*q01+p01q11 |
    // | p10 p11 |   | q10 q11 |   | p10*q00+p11*q10 p10*q01+p11q11 |

    // operator for fast exponentiation
    matrix operator* (const matrix& other){
        uint64_t r = p00.p.size();
        matrix ret(r, n);
        ret.p00 = this->p00*other.p00 + this->p01*other.p10;
        ret.p01 = this->p00*other.p01 + this->p01*other.p11;
        ret.p10 = this->p10*other.p00 + this->p11*other.p10;
        ret.p11 = this->p10*other.p01 + this->p11*other.p11;
        return ret;
    }

} matrix;

bool isprime_chebyshev(uint64_t n)
{   
    // we asssume that the prime is false
    // at first
    bool result = false;
    
    if( n<2){return false;}
    if( n<4){return true;}
    if(~n&1){return false;}


    /*
     * If this conjecture is true, we can modify the algorithm 
     * slightly to first search for an r which does not divide n cong 2 − 1. 
     * Such an r can assuredly be found in the range [2, 4 log n]. 
     * This is because the product of prime numbers less than x is at least e x
     * (see [Apo97]). Thereafter we can test whether
     * the congruence (6) holds or not. Verifying the congruence takes time O∼(r log2n). 
     * This gives a time complexity of O∼(log3 n). return result;
     */

    const int PRIMES[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 
                          37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 
                          79, 83, 89, 97, 101, 103, 107, 109, 113,
                          127, 131, 137, 139, 149, 151, 157, 163, 
                          167, 173};
    uint64_t s, x;
    for (int i=0; i < sizeof(PRIMES)/sizeof(*PRIMES); i++){
        s = PRIMES[i];
        if(  n   == s){return true;}
        if(n % s == 0){return false;}
        x = n % s;
        if(x*x % s != 1){break;}
    }
    const uint64_t r = s;

    /* 
     * We have selected the r that satisfies the conditions above. 
     * Let n be a natural number greater than two . 
     * Let r be the smallest odd prime number such that r∤nr∤n and n2 \neq 1(modr)n2≢1(modr). 
     * Let Tn(x)Tn(x) be Chebyshev polynomial of the first kind, 
     * then nn is a prime number 
     * if and only if Tn(x) \eq x^n(mod x^r−1,n)
     */

    matrix poly(r, n);
    
    poly.p00.p[1] =  2 ;
    poly.p01.p[0] = n-1;
    poly.p10.p[0] =  1 ;

    /* now since we already have the exponent which is n -1
     * now we could just do fast exponentiation
     * square the matrix, then check if the n-1 is 
     * shifting integer to the right
     * if it's 
     */
    
    matrix powered(r, n);

    /**
     * Initialize powered to an 
     * Identity matrix 
     */

    powered.p00.p[0] = 1;
    powered.p11.p[0] = 1;
   
    // fast exponentiation starts here 
    // see gaIMod in Olexa's code
    
    x = n-1;

    while(x){
        if(x & 1){
            // 
            powered = powered*poly;
        }
        poly = poly*poly;
        x >>= 1;
    }

    // Powered is poly**(n-1);
    //
    
    polynomial v0(r,n), v1(r,n);
    v0.p[1] = 1;// x
    v1.p[0] = 1;// 1

    polynomial Tn = powered.p00*v0 + powered.p01*v1;

    // Is Tn === x^n (mod x^r - 1)
    // This means
    //   1) Tn.p[n % r ] == 1
    //   2) Tn.p[others] == 0
    //

    for(int i=0; i<r; i++){
        if(i == n%r){
            if(Tn.p[i] != 1){
                return false;
            }
        }else{
            if(Tn.p[i] != 0){
                return false;
            }
        }

        //if(Tn.p[i] == (i == n%r)){
        //    return false;
        //}
    }

    return true;
}

static void BM_chebyshev(benchmark::State& state) {

  for (auto _ : state) {
    bool prime = isprime_chebyshev(state.range(0));
    bool sanity_test = gaIIsPrime(state.range(0));
    state.counters["IS PRIME"] = prime; 
    if (sanity_test != prime) {
        std::cout << "Sanity check failed for " << state.range(0) << "\n";
        break;
    }
  }
  state.SetComplexityN(state.range(0));
}

BENCHMARK(BM_chebyshev)->DenseRange(1, std::stol(std::getenv("MAX_INT_CHEBYSHEV") ) )->Complexity();
BENCHMARK_MAIN();
