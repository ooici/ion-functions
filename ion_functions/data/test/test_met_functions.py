#!/usr/bin/env python
"""
@package ion_functions.test.met_functions
@file ion_functions/test/met_functions.py
@author Stuart Pearce
@brief Unit tests for met_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np

from ion_functions.data.met_functions import  # TODO: fill this out with the functions to be tested


@attr('UNIT', group='func')
class TestMetFunctionsUnit(BaseUnitTestCase):  # TODO: Name the class appropriately

    # Nobska instrument outputs velocities in cm/s
    def test_vel3d_nobska(self):

        lat = 14.6846
        lon = -51.044

        # timestamp in seconds since 1900-01-01
        ts = np.array([
            3319563600, 3319567200, 3319570800, 3319574400,
            3319578000, 3319581600, 3319585200, 3319588800, 3319592400,
            3319596000], dtype=np.float)

        # input velocities in cm/s
        ve = np.array([
            -3.2,  0.1,  0., 2.3, -0.1,  5.6,  5.1,  5.8,
            8.8, 10.3])
        vn = np.array([
            18.2,  9.9, 12., 6.6, 7.4,  3.4, -2.6,  0.2,
            -1.5,  4.1])
        vu = np.array([-1.1, -0.6, -1.4, -2, -1.7, -2, 1.3, -1.6, -1.1, -4.5])

        # convert from cm/s to m/s for metbk
        ve /= 100.0
        vn /= 100.0
        vu /= 100.0

        # expected output velocities in m/s
        ve_expected = np.array([
            -0.085136, -0.028752, -0.036007, 0.002136, -0.023158,
            0.043218, 0.056451, 0.054727, 0.088446, 0.085952])
        vn_expected = np.array([
            0.164012,  0.094738,  0.114471,  0.06986,  0.07029,
            0.049237, -0.009499,  0.019311,  0.012096,  0.070017])
        vu_expected = np.array([
            -0.011, -0.006, -0.014, -0.02, -0.017, -0.02,
            0.013, -0.016, -0.011, -0.045])
        ve_cor = metbk_mag_corr_east(ve, vn, lat, lon, ts, 6)
        vn_cor = metbk_mag_corr_north(ve, vn, lat, lon, ts, 6)
        vu_cor = vu

        np.testing.assert_array_almost_equal(ve_cor, ve_expected)
        np.testing.assert_array_almost_equal(vn_cor, vn_expected)
        np.testing.assert_array_almost_equal(vu_cor, vu_expected)
