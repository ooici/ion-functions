#!/usr/bin/env python
"""
@package ion_functions.test.velocity_functions
@file ion_functions/test/velocity_functions.py
@author Stuart Pearce
@brief Unit tests for velocity_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import velocity_functions as velfunc


@attr('UNIT', group='func')
class TestGenericFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_vel3d_b_mag_correction(self):
            """
            Test the vel3d_b_mag_correction function.
    
            Description: test values from DPS
            
            Implemented by Stuart Pearce, April 2013
            """
    
            lat = np.array([45.0, 45.0, 80.0, 0.0, -80.0, 80.0, 0.0, -80.0])
            lon = np.array([-128.0, -128.0, 0.0,120.0,
                            240.0, 0.0, 120.0, 240.0])
            z = np.array([0.0, 1000.0, 0.0, 0.0,
                          0.0, 100000.0, 100000.0, 100000.0])
            timestamp = np.array([3575053740.7382507,  # 2013-04-15 22:29:00
                                  3575053740.7382507,  # UTC
                                  3471292800.0,        # 2010-01-01 UTC
                                  3471292800.0,
                                  3471292800.0,
                                  3471292800.0,
                                  3471292800.0,
                                  3471292800.0])
            zflag = np.array([-1, -1, -1, -1, -1, 1, 1, 1])
            
            
            output = []
            for ii in range(len(lat)):
                out_ = velfunc.vel3d_b_mag_correction(uu,vv,ww,
                                                      lat[ii],lon[ii],z[ii],
                                                      timestamp[ii],zflag[ii])
                output.append(out_)
            output = np.array(output)
            
            check_values = np.array([16.46093044096720,
                                     16.46376239313584,
                                     -6.13,
                                     0.97,
                                     70.21,
                                     -6.57,
                                     0.94,
                                     69.62])
            self.assertTrue(np.allclose(output, check_values,
                                        rtol=0, atol=1e-2))