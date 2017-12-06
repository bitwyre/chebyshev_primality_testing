/* 
 * By Dendi Suhubdy, 2017
 */

#ifndef __CHEBYSHEV_H__
#define __CHEBYSHEV_H__

double T0(const long x)
{
    return 1.0;
}

double T1(const long x)
{
    return x;
}

double T2(const long x)
{
    return (2.0 * x * x) - 1.0;
}

double Tn(unsigned long n, const long x)
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
    
    double tnm1 = T2(x);
    double tnm2 = T1(x);
    double tn = tnm1;

    for (unsigned int l=3 ; l <= n; l++)
    {
        tn = (2 * x * tnm1) - tnm2;
        tnm2 = tnm1;
        tnm1 = tn;
    }

    return tn;
}

#endif
