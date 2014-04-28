#ifndef __POLYCALS_H__
#define __POLYCALS_H__

#include <stddef.h>
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
 */
size_t polycal(double *out,         /* output vector */
               coeff_vector *cals,  /* calibration vector (ragged array) */
               double *cal_t,       /* timestamps for calibrations */
               size_t cal_len,      /* length of cals and cal_t */
               double *x,           /* vector of data, x */
               double *t,           /* timestamps for x */
               size_t x_len);       /* length of x and t */
#endif /* __POLYCALS_H__ */
