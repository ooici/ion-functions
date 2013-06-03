#include <math.h>
#include <stdlib.h>

static int all(const int *arr, size_t len);
static double double_abs(double a);
static int int_min(int a, int b);
static int comp_res(double a, double b, double res);

/*
 * stuck
 * sets out[i] to 0 where i in dat is a stuck value
 */
int stuck(signed char *out, const double *dat, size_t len, double reso, int num)
{
    int i,j,k;
    int i_max;
    int *tmp = malloc(sizeof(int) * (num-1));
    for(i=0;i<(len - num + 1);i++) {
        i_max = int_min(i+num, len)-1;
        
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
                
/*
 * Performs a boolean evaluation on an entire array
 */
static int all(const int *arr, size_t len) 
{
    int i=0;
    for(i=0;i<len;i++) {
        if(!arr[i])
            return 0;
    }
    return 1;
}

/* double_abs
 * Returns the absolute value of a
 */
static double double_abs(double a) 
{
    if(a>0)
        return a;
    return -a;
}
/*
 * int_min
 * Returns the minimum value between a and b
 */
static int int_min(int a, int b) 
{
    if(a<=b) 
        return a;
    return b;
}

/*
 * comp_res
 * Returns true if the difference between a and b is less than res
 */
static int comp_res(double a, double b, double res) 
{
    return (double_abs(b-a) < res);
}

