#ifndef __POLYCALS_H__
#define __POLYCALS_H__

#include <stddef.h>
size_t search_sorted(size_t *out, double *a, size_t a_len, double *v, size_t v_len);

typedef struct _coeff_vector 
{
    size_t N;
    double *coeff; /* Coefficients */
} coeff_vector;


#endif /* __POLYCALS_H__ */
