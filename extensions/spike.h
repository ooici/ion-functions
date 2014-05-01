#ifndef __SPIKE_H__
#define __SPIKE_H__

/*
 * spike
 *
 * This algorithm generates flags for data values according to whether a single
 * data value deviates significantly from surrounding data values. The purpose
 * of this document is to serve as a reference in order to document which
 * processing steps have been applied to a data product.
 *
 * Arguments:
 * signed char *out  - The output array of flags
 * const double *dat - The data vector
 * size_t len        - Length of the data vector
 * int L             - Window Length
 * double N          - Range multipier
 * double ACC        - accuracy
 */
int spike(signed char *out, const double *dat, size_t len, int L, double N, double acc);

#endif /* __SPIKE_H__ */
