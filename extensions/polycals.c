#include <stdio.h>
#include <math.h>
#include "polycals.h"

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



/*
 * polycal
 *
 * The algorithm combines a polynominal evaluation and a linear interpolation
 * to produce a calibrated vector of data for the input vector x. At each point
 * of x, if the point lies between two calibration values, defined by the cal_t
 * vector and t, then the output is the linear interpolation of a polynominal
 * evaluation of x using the calibration coefficients between the two
 * calibration points.
 *         { x_i, if i is not between calibrations
 * out_i = { (1-w) * u_j + w * u_(j+1), if i is between calibrations
 *
 * u is the polynomial evaluation of x using cals
 *
 *     N 
 *     __             i
 * u = >   c_(i-N) * x
 *     --
 *     i
 */
size_t polycal(double *out,         /* output vector */
               coeff_vector *cals,  /* calibration vector (ragged array) */
               double *cal_t,       /* timestamps for calibrations */
               size_t cal_len,      /* length of cals and cal_t */
               double *x,           /* vector of data, x */
               double *t,           /* timestamps for x */
               size_t x_len)        /* length of x and t */
{
    size_t a[cal_len];

    search_sorted(a, t, x_len, cal_t, cal_len);

    size_t i=0; // step in x
    size_t j=0; // step in the calibrations
    double x0;
    double x1;
    double w;

    /* Iterate through every x */
    for(i=0; i<x_len; i++)
    {
        /* If there are no more calibrations, set out_i = x_i */
        if((j+1) >= cal_len) 
        {
            out[i] = x[i];
            continue;
        } 
        /* x is no longer bounded between this segment of calibrations */
        else if(i > a[j+1]) 
        {
            /* set the calibration step to the next segment boundary */
            j = min(j+1, cal_len-1);

            /* There is a case where after advacing j that this could be the last
             * segment. If this is not the last segment, interpolate. */
            if( (j+1 < cal_len) && (a[j] <= i && i <= a[j+1]) )
            {
                x0 = polyval(&cals[j], x[i]);
                x1 = polyval(&cals[j+1], x[i]);
                w = (t[i] - t[a[j]] * 1.0) / ( t[a[j+1]] - t[a[j]] * 1.0);
                out[i] = x0 * (1 - w) + x1 * w;
            }
            /* This is now one step after the last segment, so we don't calibrate */
            else
            {
                out[i] = x[i];
            }
        }
        /* x lies between two calibration boundaries, so we calibrate based on the calibrations */
        else if(a[j] <= i && i <= a[j+1])
        {
            x0 = polyval(&cals[j], x[i]);
            x1 = polyval(&cals[j+1], x[i]);
            w = (t[i] - t[a[j]] * 1.0) / ( t[a[j+1]] - t[a[j]] * 1.0);
            out[i] = x0 * (1 - w) + x1 * w;
        }
        /* x lies outside the boundaries, we don't calibrate */
        else
        {
            out[i] = x[i];
        }
    }
    return i;
}




