#ifndef __WMM_H__
#define __WMM_H__


#include "GeomagnetismHeader.h"
#include <inttypes.h>


static const char *wmm_errmsg;
typedef struct velocity_profile_ {
    size_t len;
    double *uu;
    double *vv;
    double *lat;
    double *lon;
    double *z;
    int64_t *timestamp;
} velocity_profile;

typedef struct wmm_model_ {
    MAGtype_MagneticModel *MagneticModel;
    MAGtype_MagneticModel *TimedMagneticModel;
    int initialized;
} WMM_Model;


// adding a associated array mapping configuration file with WMM_model,
// to prevent reinitialization when presented with a previously seen
// configuration. Arrange as a simple linked list
typedef struct WMM_PrevConf WMM_PrevConf;

typedef struct WMM_ModelCache {
    char filename[256];             // base filename of .COF file used to init mode
    WMM_Model model;                // model built from the above COF
    int inited;                     // flag to indicate a completed entry
    time_t lastTouch;               // last time the model was used, to prioritize replacements if cache is full
} WMM_ModelCache;


size_t wmm_velocity_correction(const velocity_profile *in, WMM_Model *model, velocity_profile *out);

int wmm_initialize(char *filename, WMM_Model **model);
int wmm_free(WMM_Model *model);
double wmm_declination(WMM_Model *model, double lat, double lon, double z, int year, int month, int day);

#ifndef NaN
#define NaN 0./0.
#endif 

#endif /* __WMM_H__ */
