#include <stdio.h>
#include <strings.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include "spike.h"

void arange(double *arr, size_t len);
signed char all(signed char *, size_t);
void print_array(signed char *, size_t);
void print_double_array(double *, size_t);

int main(int argc, char *argv[])
{
    double example[20];
    size_t i=0;
    signed char output[20];
    arange(example,20);
    for(i=0;i<20;i++)
        output[i] = 1;
    printf("Dat\n");
    example[0] = 100.0; /* Make a spike at the beginning */
    example[19] = 100.0; /* Make a spike at the end */
    example[10] = 100.0; /* Make a spike in the middle */
    print_double_array(example, 20);
    spike(output, example, 20, 5, 5, 0.1);
    print_array(output, 20);
    assert(!all(output,20));
    for(i=0;i<20;i++) {
        if(i==0 || i==19 || i==10) {
            assert(output[i]==0);
        }
        else {
            assert(output[i]==1);
        }
    }
    printf("Tests pass...\n");
    return 0;
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

