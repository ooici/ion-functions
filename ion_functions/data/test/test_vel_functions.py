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
from ion_functions.data import vel_functions as velfunc


@attr('UNIT', group='func')
class TestGenericFunctionsUnit(BaseUnitTestCase):
    pass

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
