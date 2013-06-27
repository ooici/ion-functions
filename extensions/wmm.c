#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "wmm.h"
#include "EGM9615.h"

static bool fexists(char *filename)
{
    FILE *f;
    if(f = fopen(filename,"r")) {
        fclose(f);
        return true;
    }
    return false;
}



int wmm_initialize(char *filename, WMM_Model *model)
{
    int nMax = 0;
    int NumTerms;

    if(!fexists(filename)) {
        wmm_errmsg = "File does not exist";
        return 1;
    }

    if(!MAG_robustReadMagModels(filename, (MAGtype_MagneticModel *(*)[]) &model->MagneticModel, 1)) {
        wmm_errmsg = "Failed to read and initialize WMM";
        return 1;
    }

    if(nMax < model->MagneticModel->nMax)
        nMax = model->MagneticModel->nMax;

    NumTerms = ((nMax + 1) * (nMax + 2)/2);
    model->TimedMagneticModel = MAG_AllocateModelMemory(NumTerms);
    if(model->MagneticModel == NULL || model->TimedMagneticModel == NULL)
    {
        wmm_errmsg = "Failed to allocate WMM models";
        return 1;
    }
    model->initialized = 1;
    return 0;
}

int wmm_free(WMM_Model *model)
{
    MAG_FreeMagneticModelMemory(model->TimedMagneticModel);
    MAG_FreeMagneticModelMemory(model->MagneticModel);
    model->initialized = 0;
    return 0;
}




double wmm_declination(WMM_Model *model, double lat, double lon, double z, int year, int month, int day)
{
    int nMax = 0;
    int epochs = 1;
    int NumTerms;

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

