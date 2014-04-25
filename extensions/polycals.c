#include "polycals.h"
#include <math.h>

static inline size_t min(size_t a, size_t b) 
{
    return ((a < b) ? a : b);
}

static double polyval(const coeff_vector *cal, const double v)
{
    double retval = 0.0;
    size_t i=0;
    for(i=0; i < cal->N; i++) {
        retval += pow(v, cal->N - (i+1)) * cal->coeff[i];
    }
    return retval;
}

size_t search_sorted(size_t *out, double *a, size_t a_len, double *v, size_t v_len)
{
       size_t i=0;
       size_t j=0;
       for(i=0;i<v_len;i++) 
       {
           for(; j < a_len && a[j] < v[i]; j++);
           out[i] = j;
       }
       return i;
}



size_t polycal(double *out, 
               coeff_vector *cals, 
               double *cal_t,
               size_t cal_len, 
               double *t,
               double *x,
               size_t x_len) 
{
    size_t a[cal_len];

    search_sorted(a, t, x_len, cal_t, cal_len);

    size_t i=0;
    size_t j=0;
    double x0;
    double x1;
    double w;

    for(i=0; i<x_len; i++)
    {
        if((j+1) >= cal_len) 
        {
            out[i] = x[i];
            continue;
        } 
        else if(i > a[j+1]) 
        {
            j = min(j+1, cal_len-1);

            if( (j+1 < cal_len) && (a[j] <= i && i <= a[j+1]) )
            {
                x0 = polyval(&cals[j], x[i]);
                x1 = polyval(&cals[j], x[i]);
                w = (t[i] - t[a[j]] * 1.0) / ( t[a[j+1]] - t[a[j]] * 1.0);
                out[i] = x0 * (1 - w) + x1 * w;
            }
            else
            {
                out[i] = x[i];
            }
        }
        else if(a[j] <= i && i <= a[j+1])
        {
            x0 = polyval(&cals[j], x[i]);
            x1 = polyval(&cals[j], x[i]);
            w = (t[i] - t[a[j]] * 1.0) / ( t[a[j+1]] - t[a[j]] * 1.0);
            out[i] = x0 * (1 - w) + x1 * w;
        }
        else
        {
            out[i] = x[i];
        }
    }
}




