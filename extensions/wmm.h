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

size_t wmm_velocity_correction(const velocity_profile *in, WMM_Model *model, velocity_profile *out);

int wmm_initialize(char *filename, WMM_Model *model);
int wmm_free(WMM_Model *model);
double wmm_declination(WMM_Model *model, double lat, double lon, double z, int year, int month, int day);

#ifndef NaN
#define NaN 0./0.
#endif 

#endif /* __WMM_H__ */
