#include <math.h>
#include <stddef.h>
#include <stdbool.h>
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

