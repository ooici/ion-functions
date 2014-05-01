#ifndef __POLYCALS_H__
#define __POLYCALS_H__

#include <stddef.h>

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
size_t search_sorted(size_t *out, double *a, size_t a_len, double *v, size_t v_len);

typedef struct _coeff_vector 
{
    size_t N;
    double *coeff; /* Coefficients */
} coeff_vector;

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
               size_t x_len);
#endif /* __POLYCALS_H__ */
