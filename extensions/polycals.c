/*
 * polycals.c -- Implementation of the polycals algorithm
 *
 * Description:
 *
 *   Secondary Calibrations of data streams.
 *
 * Implemented by:
 *   
 *   2014-04-30: Luke Campbell, initial implementation.
 */
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

/*
 * search_sorted
 *
 * The function finds indices into a sorted array `a` such that, if the
 * corresponding elements in `v` were inserted before the indices, the order
 * would be preserved.
 *
 * Arguments:
 *   size_t *out  - The array of indices that gets set.
 *   double *a    - The haystack
 *   size_t a_len - Length of haystack
 *   double *v    - The needle
 *   size_t v_len - Length of needle
 */
size_t search_sorted(size_t *out, double *a, size_t a_len, double *v, size_t v_len)
{
       size_t i=0;
       size_t j=0;
       for(i=0;i<v_len;i++) {
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
 * Arguments:
 *
 * double *out        - Output vector
 * coeff_vector *cals - Calibration vector (ragged array)
 * double *cal_t      - Timestamps for calibrations
 * size_t cal_len     - Length of cals and cal_t
 * double *x          - Vector of data
 * double *t          - Timestamps for x
 * size_t x_len       - Lenght of x and t
 *
 */
size_t polycal(double *out,
               coeff_vector *cals,
               double *cal_t,
               size_t cal_len,
               double *x,
               double *t,
               size_t x_len)
{
    size_t a[cal_len];

    size_t i=0;     /* step in x */
    size_t lower=0; /* lower bound in calibrations */
    size_t upper=1; /* the upper bound in calibrations  */
    size_t a_lower;
    size_t a_upper;

    double x0;
    double x1;
    double w;

    /*
     * Fill a with the indexes where cal_t fits in t
     */
    search_sorted(a, t, x_len, cal_t, cal_len);
    /* 
     * Iterate through every x 
     */
    for(i=0; i < x_len; i++) {
        if ( upper >= cal_len ) {
            /* 
             * If there are no more calibrations, set out_i = x_i 
             */
            out[i] = x[i];
            continue;
        } else if (i > a[upper]) {
            /* 
             * x is no longer bounded between this segment of calibrations 
             * set the calibration step to the next segment boundary 
             */
            lower = min(upper, cal_len - 1);
            upper = lower + 1;
        }

        if ( (upper < cal_len) && (a[lower] <= i && i <= a[upper])) {
            /* 
             * x lies between two calibration boundaries, so we calibrate based on
             * the calibrations 
             */
            x0 = polyval(&cals[lower], x[i]);
            x1 = polyval(&cals[upper], x[i]);
            a_lower = a[lower];
            a_upper = a[upper];
            w = (t[i] - t[a_lower]) / (t[a_upper] - t[a_lower]);
            out[i] = x0 * (1 - w) + x1 * w;
        } else {
            out[i] = x[i];
        }
    }
    return i;
}




