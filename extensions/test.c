#include <stdio.h>
#include <strings.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include "spike.h"
#include "utils.h"
#include "gradient.h"

void arange(double *arr, size_t len);
signed char all(signed char *, size_t);
void print_array(signed char *, size_t);
void print_double_array(double *, size_t);
char test_spike_simple(void);
char test_spike_l(void);
char test_spike_long(void);
char test_polyval(void);
char test_gradient(void);
char test_gradient2(void);
char test_gradient3(void);
void test(char (*func)(void));

extern bool nearly_equal(double, double, double);

static const char *message=NULL;

int main(int argc, char *argv[])
{
    test(&test_spike_simple);
    test(&test_spike_l);
    test(&test_spike_long);
    test(&test_polyval);
    test(&test_gradient);
    test(&test_gradient2);
    test(&test_gradient3);
    return 0;
}

void test(char (*func)(void))
{
    message=NULL;
    if(func())
        printf("ok\n");
    else {
        printf("FAIL\n");
        if(message)
            printf("%s\n", message);
    }
}

char test_gradient()
{
    double dat[] = {3., 5., 98., 99., 4.};
    double x[] = {1., 2., 3., 4., 5. };
    double grad_min = -50;
    double grad_max = 50;
    double mindx = 0;
    double startdat = 0;
    double toldat = 5;
    size_t len = 5;
    signed char expected[] = {1, 1, 0, 0, 1};
    signed char out[] = {1, 1, 1, 1, 1};
    printf("test_gradient... ");
    gradient(out, dat, x, len, grad_min, grad_max, mindx, startdat, toldat);
    for(int i=0;i<len;i++) {
        if(!(expected[i]==out[i])) {
            message = "Expected doesn't match received";
            printf("\n");
            print_double_array(dat,5);
            print_array(expected,5);
            print_array(out,5);
            return false;
        }
    }
    return true;
}

char test_gradient2()
{
    double dat[] = {3., 5., 98., 99., 4.};
    double x[] = {1., 2., 3., 4., 5. };
    double grad_min = -50;
    double grad_max = 50;
    double mindx = 0;
    double startdat = 100;
    double toldat = 5;
    size_t len = 5;
    signed char expected[] = {0, 0, 1, 1, 0};
    signed char out[] = {1, 1, 1, 1, 1};
    printf("test_gradient2... ");
    gradient(out, dat, x, len, grad_min, grad_max, mindx, startdat, toldat);
    for(int i=0;i<len;i++) {
        if(!(expected[i]==out[i])) {
            message = "Expected doesn't match received";
            printf("\n");
            print_double_array(dat,5);
            print_array(expected,5);
            print_array(out,5);
            return false;
        }
    }
    return true;
}
char test_gradient3()
{
    double dat[] = {3., 5., 98., 99., 4.};
    double x[] = {1., 2., 3., 3.1, 4. };
    double grad_min = -50;
    double grad_max = 50;
    double mindx = 0.2;
    double startdat = 0;
    double toldat = 5;
    size_t len = 5;
    signed char expected[] = {1, 1, 0, 0, 1};
    signed char out[] = {1, 1, 1, 1, 1};
    printf("test_gradient3... ");
    gradient(out, dat, x, len, grad_min, grad_max, mindx, startdat, toldat);
    for(int i=0;i<len;i++) {
        if(!(expected[i]==out[i])) {
            message = "Expected doesn't match received";
            printf("\n");
            print_double_array(dat,5);
            print_array(expected,5);
            print_array(out,5);
            return false;
        }
    }
    return true;
}



char test_polyval()
{
    double p[] = {1., 2., -3.};
    double inputs[] = {0., 1., 2., 3.};
    double expected[] = {-3., 0., 5., 12.};
    double e = 0.00001;
    int i=0;
    printf("test_polyval... ");

    for(i=0;i<4;i++)
        if(!nearly_equal(expected[i], polyval(p,3,inputs[i]), e)) 
            return false;
    return true;
}

char test_spike_simple()
{
    double example[20];
    size_t i=0;
    signed char output[20];
    printf("test_spike_simple... ");
    arange(example,20);
    for(i=0;i<20;i++)
        output[i] = 1;
    example[0] = 100.0; /* Make a spike at the beginning */
    example[19] = 100.0; /* Make a spike at the end */
    example[10] = 100.0; /* Make a spike in the middle */
    spike(output, example, 20, 5, 5, 0.1);
    if(all(output,20))
        return 0;
    for(i=0;i<20;i++) {
        if(i==0 || i==19 || i==10) {
            if(output[i]!=0)
                return 0;
        }
        else {
            if(output[i]!=1)
                return 0;
        }
    }
    return 1;
}

char test_spike_l()
{
    double example[20];
    size_t i=0;
    signed char output[20];
    printf("test_spike_l simple... ");
    arange(example,20);
    for(i=0;i<20;i++)
        output[i] = 1;
    example[0] = 100.0; /* Make a spike at the beginning */
    example[19] = 100.0; /* Make a spike at the end */
    example[10] = 100.0; /* Make a spike in the middle */
    spike(output, example, 20, 3, 5, 0.1);
    if(all(output,20))
        return 0;
    for(i=0;i<20;i++) {
        if(i==0 || i==19 || i==10) {
            if(output[i]!=0)
                return 0;
        }
        else {
            if(output[i]!=1)
                return 0;
        }
    }
    return 1;
}
char test_spike_long()
{
    double dat[] = { -1 , 3 , 40 , -1 , 1 , -6 , -6 , 1 , 2 , 4 , 3 , 1 , -1 , 40 , 1 , 1 , 4 , 2 , 2 , 2 , 1 , 2 , 100 };
    signed char expected[] = { 1 , 1 , 1  , 1  , 1 , 1  , 1  , 1 , 1 , 1 , 1 , 1 , 1  , 0  , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 };
    size_t len = sizeof(dat)/sizeof(double);
    size_t i=0;
    signed char output[len];
    printf("test_spike_long... ");

    spike(output, dat, len, 7, 5, 0.1);
    for(i=0;i<len;i++) {
        if(expected[i] != output[i]) {
            message = "Expected does not match received.";
            printf("\n");
            print_array(expected, len);
            print_array(output, len);
            return 0;
        }
    }
    return 1;
}

signed char all(signed char *input, size_t len) 
{
    size_t i=0;
    for(i=0;i<len;i++) {
        if(!input[i]) 
            return 0;
    }
    return 1;
}

void arange(double *arr, size_t len)
{
    size_t i=0;
    for(i=0;i<len;i++) {
        arr[i] = (double) i;
    }
}

void print_array(signed char *array, size_t len)
{
    size_t i=0;
    printf("[");
    for(i=0;i<len-1;i++) {
        printf("%d ", (int)array[i]);
    }
    printf("%d]\n", (int)array[i]);
}

void print_double_array(double *array, size_t len)
{
    size_t i=0;
    printf("[");
    for(i=0;i<len-1;i++) {
        printf("%.2f ", (double)array[i]);
    }
    printf("%.2f]\n", (double)array[i]);
}

