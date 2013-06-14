#include <math.h>
#include <stddef.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include "utils.h"


double polyval(const double *p, size_t N, double x)
{
    double total=0;
    size_t i=0;
    for(i=0;i<N;i++) {
        total += p[i] * pow(x, N-(i+1));
    }
    return total;
}

short int ntp_month(double t)
{
    time_t tval = lround(t - NTP_OFFSET);
    struct tm cal;

    gmtime_r(&tval, &cal);
    return cal.tm_mon;
}

