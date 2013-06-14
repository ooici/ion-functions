#include <stddef.h>
#include "utils.h"
#include "time_utils.h"


int ntp_month_vector(short int *out, const double *in, size_t len)
{
    size_t i;
    /* We're going to chop it up baby! */
    short int i0, i1;
    double dt;
    int ret1;
    int ret2;
    i0 = ntp_month(in[0]);
    i1 = ntp_month(in[len-1]);
    dt = in[len-1] - in[0];
    if((i0==i1) && (dt < MONTH_S)) {/* The first and last value are in the same month */
        for(i=0;i<len;i++) {
            out[i] = i0;
        }
    }
    else {
        /* Binary traversal */
        if(!(len%2)) { /* len is even */
            ret1 = ntp_month_vector(out, in, len/2);
            ret2 = ntp_month_vector(out + (len/2), in + (len/2), len/2);
            return ret1 + ret2;
        }
        else {
            ret1 = ntp_month_vector(out, in, len/2);
            ret2 = ntp_month_vector(out + (len/2), in + (len/2), len/2 + 1);
        }
    }
    return len;
}

