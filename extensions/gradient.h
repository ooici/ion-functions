#ifndef __GRADIENT_H__
#define __GRADIENT_H__

#include <stddef.h>

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
int gradient(signed char *out, 
        const double *dat, 
        const double *x, 
        size_t len, 
        double grad_min, 
        double grad_max, 
        double mindx, 
        double startdat, 
        double toldat, 
        const signed char skipped_value);

#endif /* __GRADIENT_H__ */
