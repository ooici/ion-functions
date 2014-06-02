#!/usr/bin/env python
'''
@file ion_functions/data/test/test_interpolation
'''
import numpy as np

from ion_functions.test.base_test import BaseUnitTestCase
from ion_functions.data.interpolation import polyval_calibration, secondary_interpolation

class TestInterpolation(BaseUnitTestCase):
    def test_basic_polyval(self):
        '''
        Ensures that the polyval evaluation for secondary calibrations work
        Note: It currently uses record arrays, that is not intended for release.
        '''
        # Note to readers: DONT USE RECORD ARRAYS!!!!!
        # This is a stop-gap measure and is only temporary!
        # WARNING! PLEASE READ THE ABOVE STATEMENT!
        calibrations = np.array([(0.0, 0.0, 0.0, 1.02, 2.0)] * 10, dtype=('<f4,<f4,<f4,<f4,<f4'))
        values = np.arange(10, dtype='<f4')

        output = polyval_calibration(calibrations, values)
        np.testing.assert_allclose(output, values * 1.02 + 2.0)

    def test_interpolation(self):
        '''
        Tests the basic functionality and interface for secondary_calibration
        '''
        signal0 = np.ones(20, dtype='<f4') * 30
        signal1 = np.ones(20, dtype='<f4') * 50

        time = np.arange(2000, 2020)
        starts = np.ones(20) * 2010
        ends = np.ones(20) * 2015

        interpolated = secondary_interpolation(time, signal0, signal1, starts, ends)
        expected = np.array([ np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  np.nan,  
                              np.nan,  np.nan,  np.nan,  np.nan,  30.,     34.,  
                              38.,     42.,     46.,     50.,     np.nan,  np.nan,  
                              np.nan,  np.nan], dtype='<f4')
        np.testing.assert_allclose(interpolated,expected)




