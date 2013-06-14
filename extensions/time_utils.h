#ifndef __TIME_UTILS_H__
#define __TIME_UTILS_H__

int ntp_month_vector(short int *out, const double *in, size_t len);

#define MONTH_S (3600*24*28)

#endif /* __TIME_UTILS_H__ */
