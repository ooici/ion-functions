#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <libgen.h>
#include "wmm.h"
#include "EGM9615.h"

// Add support for internal caching of initialized models.  Alloc storate here,
// and ignore external calls to free.
#define WMM_MAX_NCACHED 20
WMM_ModelCache wmm_mcache[WMM_MAX_NCACHED];


static bool fexists(char *filename)
{
    FILE *f;
    if( (f = fopen(filename,"r")) ) {
        fclose(f);
        return true;
    }
    return false;
}



int wmm_initialize(char *filename, WMM_Model **model)
{
    int nMax = 0;
    int NumTerms;
    char *bname;
    static int firstPass = 1;
    int cIdx, tIdx;       // cache indexs
    struct timeval tv;
    struct timezone tz;
    time_t oldest;

    // if we get here, we need to build a model, with storage pointed to by wwm_this
    if(!fexists(filename)) {
        wmm_errmsg = "File does not exist";
        return 1;
    }


    // index our cache by the basename of the config file, so that FOO1234.COF results in the
    // same model, regardless of the directory it is loaded from
    bname = basename(filename);

    // on first access, init the cache fields
    if ( firstPass == 1 ) {
        for ( cIdx = 0; cIdx < WMM_MAX_NCACHED; cIdx++ ) {
            wmm_mcache[cIdx].inited = 0;
        }
        firstPass = 0;
    }
   
    // check to see if we have already built a model based upon the supplied .COF file
    for ( cIdx = 0; cIdx < WMM_MAX_NCACHED; cIdx++ ) {
        if (wmm_mcache[cIdx].inited == 1) {
            if (strcmp(wmm_mcache[cIdx].filename, bname) == 0)  {
                // found an entry associated with the same config file,
                // return that model
                *model = &wmm_mcache[cIdx].model;
                // printf("Found existing cache entry %d\n",cIdx);
                // update the last touched time for this cache site
                gettimeofday(&tv,&tz);
                wmm_mcache[cIdx].lastTouch = tv.tv_sec;
                return 0;
            }
        }
    }


    // not already built, pick a slot
    tIdx = -1;  // not found flag
    for ( cIdx = 0; cIdx < WMM_MAX_NCACHED; cIdx++ ) {
        if (wmm_mcache[cIdx].inited == 0) {
            tIdx = cIdx;
            break;
        }
    }

    if ( tIdx == -1 ) {
        // didn't find an empty slot, pick the oldest
        oldest = wmm_mcache[0].lastTouch;
        tIdx = 0;
        for ( cIdx = 1; cIdx < WMM_MAX_NCACHED; cIdx++ ) {
            if ( wmm_mcache[tIdx].lastTouch < oldest ) {
                oldest = wmm_mcache[tIdx].lastTouch;
                tIdx = cIdx;
            }
        }
    }

    
    // printf("Building into cache entry %d\n",tIdx);

    if (wmm_mcache[tIdx].inited == 1) {
        // this slot was already used, free prealloc'd model memory to prevent memory leak
        // when MAG_robustReadMagModels internally allocs memory
        MAG_FreeMagneticModelMemory(wmm_mcache[tIdx].model.TimedMagneticModel);
        MAG_FreeMagneticModelMemory(wmm_mcache[tIdx].model.MagneticModel);
        wmm_mcache[tIdx].model.initialized = 0;
    }
   
                 
    if(!MAG_robustReadMagModels(filename, (MAGtype_MagneticModel *(*)[]) &wmm_mcache[tIdx].model.MagneticModel, 1)) {
        wmm_errmsg = "Failed to read and initialize WMM";
        return 1;
    }

    if(nMax < wmm_mcache[tIdx].model.MagneticModel->nMax)
        nMax = wmm_mcache[tIdx].model.MagneticModel->nMax;

    NumTerms = ((nMax + 1) * (nMax + 2)/2);
    wmm_mcache[tIdx].model.TimedMagneticModel = MAG_AllocateModelMemory(NumTerms);
    if(wmm_mcache[tIdx].model.MagneticModel == NULL || wmm_mcache[tIdx].model.TimedMagneticModel == NULL)
    {
        wmm_errmsg = "Failed to allocate WMM models";
        return 1;
    }

    // mark internal model as initialized
    wmm_mcache[tIdx].model.initialized = 1;

    // mark cache site as initialized
    wmm_mcache[tIdx].inited = 1;

    // remember the filename associated with this config
    strcpy(wmm_mcache[tIdx].filename, bname);
            
    // update the last touched time for this cache site
    gettimeofday(&tv,&tz);
    wmm_mcache[tIdx].lastTouch = tv.tv_sec;

    // and set the return pointer to out current entry
    *model = &wmm_mcache[tIdx].model;

    return 0;
}

int wmm_free(WMM_Model *model)
{
    // There is now a persistant cache of built models.  There will be at most WMM_MAX_NCACHE
    // internal mallocs ( from the Geomagnatism MAG_robustReadMagModels and MAG_AllocateModelMemory
    // calls ).  When replacing a cached model, these mallocs are freed above prior to realloc.
    // The WMM_MAX_NCACHE mallocs held in the cache will be GC'd on program exit.
    return 0;
}


size_t wmm_velocity_correction(const velocity_profile *in, WMM_Model *model, velocity_profile *out)
{
    size_t i=0;
    double M[2];
    double theta, theta_deg;
    int y, m, d;
    struct tm time_info;
    if(in == NULL || out == NULL)
        return 0;
    if(in->len <= 0 || out->len <= 0)
        return 0;
    if(in->uu == NULL || in->vv == NULL || in->lat == NULL 
            || in->lon == NULL || in->z == NULL 
            || in->timestamp == NULL || out->uu == NULL 
            || out->vv == NULL)
        return 0;
    for(i=0;i<out->len;i++) {
        gmtime_r((time_t*) &in->timestamp[i], &time_info);
        y = time_info.tm_year + 1900;
        m = time_info.tm_mon + 1;
        d = time_info.tm_mday;

        theta_deg = wmm_declination(model, in->lat[i], in->lon[i], in->z[i], y, m, d);
        theta = theta_deg * M_PI / 180.;
        M[0] = cos(theta);
        M[1] = sin(theta);
        out->uu[i] =  in->uu[i] * M[0] + in->vv[i] * M[1];
        out->vv[i] = -in->uu[i] * M[1] + in->vv[i] * M[0];
    }
    return i;
}

double wmm_declination(WMM_Model *model, double lat, double lon, double z, int year, int month, int day)
{
    MAGtype_Ellipsoid Ellip;
    MAGtype_CoordSpherical CoordSpherical;
    MAGtype_CoordGeodetic CoordGeodetic;
    MAGtype_Date MagneticDate;
    MAGtype_GeoMagneticElements GeoMagneticElements;
    MAGtype_Geoid Geoid;
    if(model->initialized != 1)
    {
        wmm_errmsg = "Uninitialized model";
        return NaN;
    }

    MAG_SetDefaults(&Ellip, &Geoid);

    Geoid.GeoidHeightBuffer = GeoidHeightBuffer;
    Geoid.Geoid_Initialized = 1;

    /* In lieu of User Input */

    CoordGeodetic.phi = lat;
    CoordGeodetic.lambda = lon;
    /* Set the Height MSL */
    CoordGeodetic.HeightAboveGeoid = z;
    Geoid.UseGeoid = 1;
    MAG_ConvertGeoidToEllipsoidHeight(&CoordGeodetic, &Geoid);
    
    /* Set the date */
    MagneticDate.Year = year;
    MagneticDate.Month = month;
    MagneticDate.Day = day;
    MAG_DateToYear(&MagneticDate, NULL);

    /* Convert from geodetic to spherical equations */
    MAG_GeodeticToSpherical(Ellip, CoordGeodetic, &CoordSpherical);
    /* Time adjust the coefficients */
    MAG_TimelyModifyMagneticModel(MagneticDate, model->MagneticModel, model->TimedMagneticModel);
    MAG_Geomag(Ellip, CoordSpherical, CoordGeodetic, model->TimedMagneticModel, &GeoMagneticElements);

    return GeoMagneticElements.Decl;

}

