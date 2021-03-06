#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdbool.h>

#ifndef NTP_OFFSET
#define NTP_OFFSET 2208988800
#endif

/*
 * double polyval(const double *p, size_t N, double x)
 * --------------------------------------------------------------------------------
 * Computes the polynomial evaluation 
 *  p[0] * x^(N-1) + p[1] * x^(N-2) + ... + p[N-1]
 * p - a list of coefficients
 * N - the length of the array
 * x - the value to compute
 */
double polyval(const double *p, size_t N, double x);

static inline bool nearly_equal(double a, double b, double epsilon)
{
    double diff = a - b;
    return (diff < epsilon && diff > -epsilon);
}

short int ntp_month(double t);
#endif /* __UTILS_H___*/
