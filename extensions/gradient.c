#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "gradient.h"

int gradient(signed char *out, const double *dat, const double *x, size_t len, double grad_min, double grad_max, double mindx, double startdat, double toldat)
{
    double ddatdx;
    size_t i=0;
    int j=0;
    int skipped=0;
    bool bad=false;
    if(startdat==0) { 
        /* If there was no startdat then set it to dat[0] */
        startdat = dat[0];
    }
    else {
        /* startdat was specified */
        if(fabs(dat[0] - startdat) > toldat) {
            /* If dat[0] is not within toldat of startdat then it's a bad value */
            bad=true;
            out[0] = 0;
        }
    }
    for(i=1;i<len;i++) {
        if(bad) { /* Only start again if dat[i] is within toldat of startdat */
            if(fabs(dat[i] - startdat) <= toldat) { /* It's good again */
                bad = false;
            }
            else {
                out[i] = 0; /* still bad mark it and move on */
            }
            /* The test for this element is simple, if it's within the tolerance
             * set it to startdat and continue on
             */
            continue; 
        }

        if (fabs((x[i]-x[i-(1+skipped)])) > mindx) {
            /* If the change in x is greater than mindx */
            ddatdx = (dat[i] - dat[i-(1+skipped)])/(x[i] - x[i-(1+skipped)]);
            /* Calculate the rate of change */
            if(ddatdx < grad_min || ddatdx > grad_max) {
                /* If the differential is outside of the min/max */
                for(j=1;j<=skipped;j++) {
                    /* Set all the skipped to 0 as well */
                    out[i-j] = 0;
                }
                skipped = 0; /* Reset the skipped */
                out[i] = 0;  /* Set the output to false */
                bad=true;    /* Mark a bad */
            }
            else {
                /* Continue on our way and update startdat */
                startdat = dat[i];
            }
        }
        else {
            /* dx was <= mindx */
            skipped++;
        }
    }
    return 0;
}



