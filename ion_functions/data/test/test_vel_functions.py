#!/usr/bin/env python
"""
@package ion_functions.test.vel_functions
@file ion_functions/test/vel_functions.py
@author Stuart Pearce
@brief Unit tests for vel_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data.vel_functions import nobska_mag_corr_east, nobska_mag_corr_north
from ion_functions.data.generic_functions import error


@attr('UNIT', group='func')
class TestGenericFunctionsUnit(BaseUnitTestCase):

    ## No VELPTTU Nobska test data as of yet.
    def test_nobska_mag_correction(self):
    #    """
    #    Test the nobska_mag_correction function.
    #    
    #    Description: test values from DPS
    #    
    #    Implemented by Stuart Pearce, April 2013
    #    """
    #    uu = np.array([])
    #    vv = np.array([])
    #    ww = np.array([])
    #    lat = np.array([])
    #    lon = np.array([])
    #    z = np.array([])
    #    timestamp = np.array([])
    #    
    #    output = []
    #    for ii in range(len(lat)):
    #        out_ = velfunc.nobska_mag_correction(uu,vv,
    #                                              lat[ii],lon[ii],z[ii],
    #                                              timestamp[ii],zflag=-1)
    #        output.append(out_)
    #    output = np.array(output)
    #    
    #    check_values = np.array([])
    #    self.assertTrue(np.allclose(output, check_values,
    #                                rtol=0, atol=1e-2))
        pass

    def test_vel3d_nobska(self):
        lat = 14.6846
        lon = -51.044
        ts = np.array([3193419600, 3193423200, 3193426800, 3193430400,
            3193434000, 3193437600, 3193441200, 3193444800, 3193448400,
            3193452000],dtype=np.float)

        ve = np.array([ -3.2,  0.1,  0. ,  2.3, -0.1,  5.6,  5.1,  5.8,
            8.8, 10.3])

        vn = np.array([ 18.2,  9.9, 12. ,  6.6, 7.4,  3.4, -2.6,  0.2,
            -1.5,  4.1])
        vu = np.array([-1.1, -0.6, -1.4, -2, -1.7, -2, 1.3, -1.6, -1.1, -4.5])

        vn_cor = nobska_mag_corr_north(ve, vn, lat, lon, ts, vu)
        vn_expected = np.array([ 0.164012,  0.094738,  0.114471,  0.06986,  0.07029,
                    0.049237, -0.009499,  0.019311,  0.012096,  0.070017])
        self.assertTrue((np.abs(vn_cor - vn_expected) < 0.001).all())

        ve_cor = nobska_mag_corr_east(ve, vn, lat, lon, ts, vu)
        ve_expected = np.array([-0.085136, -0.028752, -0.036007, 0.002136,
            -0.023158, 0.043218, 0.056451, 0.054727, 0.088446, 0.085952])

        self.assertTrue((np.abs(ve_cor - ve_expected) < 0.001).all())



