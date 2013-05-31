#include <math.h>
#include <stdlib.h>



static int all(const int *arr, int len);

double 
d_abs(double a) {
    if(a>0)
        return a;
    return -a;
}

int 
i_min(int a, int b) {
    if(a<=b) 
        return a;
    return b;
}


int
comp_res(double a, double b, double res) {
    return (d_abs(b-a) < res);
}


int 
stuck(signed char *out, int out_len, const double *dat, int dat_len, double reso, int num)
{
    int i,j,k;
    int i_max;
    int *tmp = malloc(sizeof(int) * (num-1));
    for(i=0;i<(dat_len - num + 1);i++) {
        i_max = i_min(i+num, dat_len)-1;
        
        if(comp_res(dat[i], dat[i_max], reso)) {
            for(j=i;j<i_max;j++) {
                tmp[(i_max)-1-j] = comp_res(dat[j], dat[i_max], reso);
            }
            if(all(tmp, num-1)) {
                for(k=i;k<=i_max;k++) {
                    out[k] = 0;
                }
            }
        }
    }
    free(tmp);
    return 0;
}
                

static int 
all(const int *arr, int len) {
    int i=0;
    for(i=0;i<len;i++) {
        if(!arr[i])
            return 0;
    }
    return 1;
}

