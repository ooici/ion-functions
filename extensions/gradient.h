#ifndef __GRADIENT_H__
#define __GRADIENT_H__

#include <stddef.h>

int gradient(signed char *out, const double *dat, const double *x, size_t len, double grad_min, double grad_max, double mindx, double startdat, double toldat, const double skipped_value);

#endif /* __GRADIENT_H__ */
