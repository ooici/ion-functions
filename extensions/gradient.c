/*
 * gradient.c -- Implementation for the gradient algorithm
 * 
 *
 * Description:
 *  
 *   A quality control algorithm test that identifies if changes between
 *   successive data points fall within a certain range.
 *
 * Implemented by:
 *
 *   2012-07-17: DPS authored by Mathias Lankhorst. Example code provided
 *               for Matlab.
 *   2013-04-06: Christopher Wingard. Initial python implementation.
 *   2014-04-01: Luke Campbell optimized python implementation by porting to C.
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "gradient.h"


static inline int tolerance(double a, double b, double tolerance)
{
    return (fabs(b - a) <= tolerance);
}


/*
 * gradient
 * 
 * The algorithm scans the data points in dat successively, one-by-one.
 * Starting from a data point assumed to be good, if the change from this point
 * to the next, divided by the respective advance in x, is within the range
 * defined by ddatdx, the next point will be assumed good, else bad. Once a
 * data point is identified as bad, successive data points will be considered
 * bad until a data point falls within toldat of the last known good data
 * point.
 *
 * As a default starting point, the first data point is considered good. This
 * default can be overridden by defining startdat, in which case the first data
 * point that falls within toldat of startdat will be considered good (and
 * potentially earlier ones as bad).
 *
 * x must be strictly increasing. In order to avoid small x steps, which can
 * lead to exaggerated gradients because the x difference is used in the
 * denominator, setting mindx to a value greater than zero is an option to
 * remove all data points dat (and x) for which x is separated by mindx or
 * less.
 *
 * Note: the out array should be initialized to 1s by the client
 *
 * Arguments:
 * signed char *out           - An array of values that represent the
 *                              evaluation of this QC.  Only bad values and
 *                              skipped values are set, all others are 
 *                              ignored. The implication here is that an array
 *                              of ones, should be passed in.
 * const double *dat          - An array of data values to apply QC against.
 * const double *x            - The axis for the data values, must be in
 *                              increasing order.
 * size_t len                 - Length of out, dat and x. Each of these arrays
 *                              MUST be the same size and match len.
 * double grad_min            - Minimum change in dat over the change in x.
 * double grad_max            - Maximum change in dat over the change in x.
 * double mindx               - If dx is less than mindx, the value is skipped.
 * double startdat            - Value to apply as a starting value (known good).
 * double toldat              - The tolerance of d(dat)/dx to reach a good 
 *                              value.
 * const double skipped_value - Caller specified value to apply when qc
 *                              evaluation is "skipped".
 */

int gradient(
        signed char *out, 
        const double *dat, 
        const double *x, 
        size_t len, 
        double grad_min, 
        double grad_max, 
        double mindx, 
        double startdat, 
        double toldat, 
        const signed char skipped_value) 
{ 
    double ddatdx; 
    size_t i=0; 
    int j=0; 
    int skipped=0; 
    bool bad=false; 

    /*
     * If startdat is not set, set it to dat[0]. Otherwise, use it to evaluate 
     * the tolerance of dat[0].
     */
    if (startdat == 0) { 
        startdat = dat[0];
    } else { 
        if( !tolerance(dat[0], startdat, toldat) ) {
            bad=true;
            out[0] = 0;
        }
    }
    for(i = 1; i < len; i++) {

        /* 
         * Check if dx < mindx and skip if it's not.
         */
        if ( tolerance(x[i], x[i - (1 + skipped)], mindx) ) {
            skipped++;
            out[i] = skipped_value;
            continue;
        }

        /*
         * If the last value was bad
         */

        if(bad) { /* Only start again if dat[i] is within toldat of startdat */
            if (tolerance(dat[i], startdat, toldat)) { /* It's good again */
                bad = false;
            } else {
                out[i] = 0; /* still bad mark it and move on */
            }
            continue; 
        }
        ddatdx = (dat[i] - dat[i-(1+skipped)])/(x[i] - x[i-(1+skipped)]);
        /* Calculate the rate of change */
        if(ddatdx < grad_min || ddatdx > grad_max) {
            /* If the differential is outside of the min/max */
            for(j=1;j<=skipped;j++) {
                /* Set all the skipped to 0 as well */
                out[i-j] = 0;
            }
            skipped = 0; /* Reset the skipped */
            out[i] = 0;  /* Set the output to false */
            bad = true;  /* Mark a bad */
        }
        else {
            /* Continue on our way and update startdat */
            startdat = dat[i];
            /* Reset the skipped count, we're done skiping for the moment */
            skipped = 0;
        }
    }
    return 0;
}



