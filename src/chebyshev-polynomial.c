/* 
 * By Dendi Suhubdy, 2017
 */

#ifndef __CHEBYSHEV_H__
#define __CHEBYSHEV_H__

double inline T0(const &x)
{
    return 1.0;
}

double inline T1(const &x)
{
    return x;
}

double inline T2(const &x)
{
    return (2.0 * x * x) - 1.0;
}

double inline Tn(unsigned long n, const &x)
{
    if (n == 0)
    {
        return T0(x);
    }
    else if (n == 1)
    {
        return T1(x);
    }
    else if (n==2)
    {
        return T2(x);
    }
    /* we ignore the use of recurrence here for speed
     * return (2.0 * x * Tn(n - 1, x)) - Tn(n - 2, x)
     * We use external memory to speed up the process
     */
    
    double inline tnm1(T2(x));
    double inline tnm2(T1(x));
    double inline tn(tnm1);

    for (unsigned int l=3 ; l <= n; l++)
    {
        tn = (2 * x * tnm1) - tnm2;
        tnm2 = tnm1;
        tnm1 = tn;
    }

    return tn;
}

#endif
