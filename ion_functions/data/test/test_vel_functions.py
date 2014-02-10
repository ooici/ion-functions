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
from ion_functions.data.vel_functions import nortek_mag_corr_east, nortek_mag_corr_north, valid_lat, valid_lon


# these test data definitions are used in the test methods
# for both the Nobska and Nortek VEL3D instruments
LAT = 14.6846
LON = -51.044

# timestamp in seconds since 1900-01-01
TS = np.array([
    3319563600, 3319567200, 3319570800, 3319574400, 3319578000,
    3319581600, 3319585200, 3319588800, 3319592400, 3319596000],
    dtype=np.float)

# input velocities in cm/s (Nobsak outputs cm/s, Nortek mm/s)
VE = np.array([-3.2, 0.1, 0., 2.3, -0.1, 5.6, 5.1, 5.8, 8.8, 10.3])
VN = np.array([18.2, 9.9, 12., 6.6, 7.4, 3.4, -2.6, 0.2, -1.5, 4.1])
VU = np.array([-1.1, -0.6, -1.4, -2, -1.7, -2, 1.3, -1.6, -1.1, -4.5])

# expected output velocities in m/s
VE_EXPECTED = np.array([
    -0.085136, -0.028752, -0.036007, 0.002136, -0.023158,
    0.043218, 0.056451, 0.054727, 0.088446, 0.085952])
VN_EXPECTED = np.array([
    0.164012,  0.094738,  0.114471,  0.06986,  0.07029,
    0.049237, -0.009499,  0.019311,  0.012096,  0.070017])
VU_EXPECTED = np.array([
    -0.011, -0.006, -0.014, -0.02, -0.017, -0.02,
    0.013, -0.016, -0.011, -0.045])


@attr('UNIT', group='func')
class TestVelFunctionsUnit(BaseUnitTestCase):

    # Nobska instrument outputs velocities in cm/s
    def test_vel3d_nobska(self):
        ve_cor = nobska_mag_corr_east(
            VE, VN, LAT, LON, TS, 6)
        vn_cor = nobska_mag_corr_north(
            VE, VN, LAT, LON, TS, 6)
        vu_cor = VU / 100.0

        np.testing.assert_array_almost_equal(ve_cor, VE_EXPECTED)
        np.testing.assert_array_almost_equal(vn_cor, VN_EXPECTED)
        np.testing.assert_array_almost_equal(vu_cor, VU_EXPECTED)

    # Nortek instrument outputs velocities in mm/s
    def test_vel3d_nortek(self):
        # convert velocities from cm/s to mm/s
        ve = VE * 10.
        vn = VN * 10.
        vu = VU * 10.

        ve_cor = nortek_mag_corr_east(
            ve, vn, LAT, LON, TS, 6)
        vn_cor = nortek_mag_corr_north(
            ve, vn, LAT, LON, TS, 6)
        vu_cor = vu / 1000.0

        np.testing.assert_array_almost_equal(ve_cor, VE_EXPECTED)
        np.testing.assert_array_almost_equal(vn_cor, VN_EXPECTED)
        np.testing.assert_array_almost_equal(vu_cor, VU_EXPECTED)

    def test_latlon_valid(self):
        lat = np.array([-90, -80, 0, 80, 90])
        self.assertTrue(valid_lat(lat))

        lat = np.array([-9999, -9999])
        self.assertFalse(valid_lat(lat))

        self.assertTrue(valid_lat(30) and valid_lat(-30))

        self.assertFalse(valid_lat(-9999) or valid_lat(-99) or valid_lat(92))

        lon = np.array([-180, -80, 0, 80, 180])
        self.assertTrue(valid_lon(lon))

        lon = np.array([-9999, -9999])
        self.assertFalse(valid_lon(lon))

        self.assertTrue(valid_lon(30) and valid_lon(-30))

        self.assertFalse(valid_lon(-9999) or valid_lon(-181) or valid_lon(181))

    def test_invalids_correct(self):
        ve_cor = nobska_mag_corr_east(
            VE, VN, -9999, -9999, TS, 6)
        vn_cor = nobska_mag_corr_north(
            VE, VN, -9999, -9999, TS, 6)

        np.testing.assert_array_almost_equal(ve_cor, np.ones(VU.shape) * -9999)
        np.testing.assert_array_almost_equal(vn_cor, np.ones(VU.shape) * -9999)
