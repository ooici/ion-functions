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

