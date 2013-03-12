#!/usr/bin/env python

"""
@package ion_functions.data.data_functions
@file ion_functions/data/data_functions.py
@author Christopher Mueller
@brief Module containing data-calculation functions.  Primarily used for calculating values in Parameter Functions
"""


def data_density(SP, p, t, lat, lon):
    from pygsw import vectors as gsw

    abs_sal = gsw.sa_from_sp(SP, p, lon, lat)
    cons_temp = gsw.ct_from_t(abs_sal, t, p)

    return gsw.rho(abs_sal, cons_temp, p)

def data_l2_density(conductivity, temp,pressure, lat, lon):
    '''
    Based on calculations done here: https://github.com/ooici/coi-services/blob/master/ion/processes/data/transforms/ctd/ctd_L2_density.py#L55
    '''
    from pygsw import vectors as gsw
    sp = gsw.sp_from_c(conductivity/42.914, temp, pressure)
    sa = gsw.sa_from_sp(sp, pressure, lon, lat)
    rho = gsw.rho(sa, temp, pressure)
    return rho

def data_l2_salinity(conductivity, temp, pressure):
    '''
    Based on calculations done here: https://github.com/ooici/coi-services/blob/master/ion/processes/data/transforms/ctd/ctd_L2_salinity.py#L58
    '''
    import pygsw.vectors as gsw
    sal_value = gsw.sp_from_c(conductivity/42.914, temp, pressure)
    return sal_value
