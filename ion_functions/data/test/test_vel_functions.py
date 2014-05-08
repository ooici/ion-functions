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
import pdb
from ion_functions.data.vel_functions import nobska_mag_corr_east, nobska_mag_corr_north
from ion_functions.data.vel_functions import nortek_mag_corr_east, nortek_mag_corr_north
from ion_functions.data.vel_functions import vel3dk_east, vel3dk_north
from ion_functions.data.vel_functions import vel3dk_up
from ion_functions.data.vel_functions import valid_lat, valid_lon
from ion_functions.data.vel_functions import vel_mag_correction
from exceptions import ValueError

# these test data definitions are used in the test methods
# for both the Nobska and Nortek instruments
LAT = 14.6846 * np.ones(10)
LON = -51.044 * np.ones(10)
DEPTH = 6 * np.ones(10)  # meters

# timestamp in seconds since 1900-01-01
TS = np.array([
    3319563600, 3319567200, 3319570800, 3319574400, 3319578000,
    3319581600, 3319585200, 3319588800, 3319592400, 3319596000],
    dtype=np.float)

# input velocities (Nobska outputs cm/s, Nortek m/s)
VE_NOBSKA = np.array([-3.2, 0.1, 0., 2.3, -0.1, 5.6, 5.1, 5.8, 8.8, 10.3])
VN_NOBSKA = np.array([18.2, 9.9, 12., 6.6, 7.4, 3.4, -2.6, 0.2, -1.5, 4.1])
VU_NOBSKA = np.array([-1.1, -0.6, -1.4, -2, -1.7, -2, 1.3, -1.6, -1.1, -4.5])
VE_NORTEK = VE_NOBSKA / 100.
VN_NORTEK = VN_NOBSKA / 100.
VU_NORTEK = VU_NOBSKA / 100.

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

    # Nobska instrument unit test
    def test_vel3d_nobska(self):
        """
        Tests functions nobska_mag_corr_east and nobska_mag_corr_north
        from the ion_functions.data.vel_functions module using test data
        found in the VELPTMN DPS.

        References:

            OOI (2012). Data Product Specification for Turbulent Point Water
                Velocity. Document Control Number 1341-00781.
                https://alfresco.oceanobservatories.org/ (See: Company Home
                >> OOI >> Controlled >> 1000 System Level >>
                1341-00781_Data_Product_SPEC_VELPTTU_Nobska_OOI.pdf

            OOI (2012). Data Product Specification for Mean Point Water
                Velocity. Document Control Number 1341-00790.
                https://alfresco.oceanobservatories.org/ (See: Company Home
                >> OOI >> Controlled >> 1000 System Level >>
                1341-00790_Data_Product_SPEC_VELPTMN_OOI.pdf)
        """
        ve_cor = nobska_mag_corr_east(
            VE_NOBSKA, VN_NOBSKA, LAT, LON, TS, DEPTH)
        vn_cor = nobska_mag_corr_north(
            VE_NOBSKA, VN_NOBSKA, LAT, LON, TS, DEPTH)

        np.testing.assert_array_almost_equal(ve_cor, VE_EXPECTED)
        np.testing.assert_array_almost_equal(vn_cor, VN_EXPECTED)

    def test_vel3d_nortek(self):
        """
        Tests functions nortek_mag_corr_east and nortek_mag_corr_north
        from the ion_functions.data.vel_functions module using test data
        from the VELPTMN DPS.

        Implemented by:

        2013-04-17: Stuart Pearce. Initial code.
        2013-04-24: Stuart Pearce. Changed to be general for all velocity
                    instruments.
        2014-02-05: Christopher Wingard. Edited to use magnetic corrections in
                    the generic_functions module.

        References:

            OOI (2012). Data Product Specification for Turbulent Point Water
                Velocity. Document Control Number 1341-00780.
                https://alfresco.oceanobservatories.org/ (See: Company Home
                >> OOI >> Controlled >> 1000 System Level >>
                1341-00780_Data_Product_SPEC_VELPTTU_Nortek_OOI.pdf

            OOI (2012). Data Product Specification for Mean Point Water
                Velocity. Document Control Number 1341-00790.
                https://alfresco.oceanobservatories.org/ (See: Company Home
                >> OOI >> Controlled >> 1000 System Level >>
                1341-00790_Data_Product_SPEC_VELPTMN_OOI.pdf)
        """

        ve_cor = nortek_mag_corr_east(
            VE_NORTEK, VN_NORTEK, LAT, LON, TS, DEPTH)
        vn_cor = nortek_mag_corr_north(
            VE_NORTEK, VN_NORTEK, LAT, LON, TS, DEPTH)

        np.testing.assert_array_almost_equal(ve_cor, VE_EXPECTED)
        np.testing.assert_array_almost_equal(vn_cor, VN_EXPECTED)

    def test_vel3dk(self):
        """
        Tests functions vel3dk_mag_corr_east and vel3dk_mag_corr_north
        from the ion_functions.data.vel_functions module using test data
        from the VELPTMN DPS.

        Implemented by:

        2013-04-17: Stuart Pearce. Initial code.
        2013-04-24: Stuart Pearce. Changed to be general for all velocity
                    instruments.
        2014-02-05: Christopher Wingard. Edited to use magnetic corrections in
                    the generic_functions module.

        References:

            OOI (2012). Data Product Specification for Mean Point Water
                Velocity. Document Control Number 1341-00790.
                https://alfresco.oceanobservatories.org/ (See: Company Home
                >> OOI >> Controlled >> 1000 System Level >>
                1341-00790_Data_Product_SPEC_VELPTMN_OOI.pdf)
        """

        beam1 = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2])
        beam2 = np.array([2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
        beam3 = np.array([4, 4, 4, 4, 4, 4, 4, 4, 4, 4])
        beam4 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        hdg = np.array([0., 36., 72., 108., 144., 180., 216., 252., 288., 324.])

        ptch = np.array([
            -20.10, -5.82, 0.0, -10.26, -8.09, -1.36, 5.33, 9.96, 13.68, 26.77])

        rll = np.array([
            -15.0, -11.67, -8.33, -5.0, 0.0, 1.67, 5.0, 8.33, 11.67, 15.0])

        vel0 = np.array([
            0, 6305, 6000, 7810, 10048, 11440, -10341, -10597, 9123, -12657])

        vel1 = np.array([
            0, 1050, 856, -1672, 3593, -2487, -5731, -3085, -1987, 2382])

        vel2 = np.array([
            7628, 0, 4974, -4204, 4896, 5937, 6480, -6376, -7271, -7576])

        Vscale = -4

        ve_cor = vel3dk_east(
            vel0, vel1, vel2, hdg, ptch, rll, beam1, beam2, beam3,
            beam4, LAT, LON, TS, DEPTH, Vscale)
        vn_cor = vel3dk_north(
            vel0, vel1, vel2, hdg, ptch, rll, beam1, beam2, beam3,
            beam4, LAT, LON, TS, DEPTH, Vscale)

        VE_expected = np.array([
            0.34404501, -0.01039404,  0.64049184, -0.17489265, -0.0739631,
            -1.09305797, -0.47947474,  0.11710443,  1.97369869, -1.6466505])

        VN_expected = np.array([
            0.91742235,  0.04629215,  0.06132321,  0.56597656, -0.35874325,
            0.37553716, -1.87672302, -1.12589293,  0.0720366, -0.6617893])

        np.testing.assert_array_almost_equal(ve_cor, VE_expected, decimal=7)
        np.testing.assert_array_almost_equal(vn_cor, VN_expected, decimal=7)

    def test_vel3dk_up(self):
        """
        Tests functions vel3dk_scale_up_vel from the
        ion_functions.data.vel_functions module using test data from the
        VELPTMN DPS.

        Implemented by:

        2014-03-06: Stuart Pearce. Initial code.

        References:

            VEL3D-K IDD (2014) (No DPS as of 2014-03-03)
                https://confluence.oceanobservatories.org/display/
                instruments/VEL3D-K__stc_imodem+-+Telemetered
        """
        beam1 = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2])
        beam2 = np.array([2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
        beam3 = np.array([4, 4, 4, 4, 4, 4, 4, 4, 4, 4])
        beam4 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        hdg = np.array([0., 36., 72., 108., 144., 180., 216., 252., 288., 324.])

        ptch = np.array([
            -20.10, -5.82, 0.0, -10.26, -8.09, -1.36, 5.33, 9.96, 13.68, 26.77])

        rll = np.array([
            -15.0, -11.67, -8.33, -5.0, 0.0, 1.67, 5.0, 8.33, 11.67, 15.0])

        vel0 = np.array([
            0, 6305, 6000, 7810, 10048, 11440, -10341, -10597, 9123, -12657])

        vel1 = np.array([
            0, 1050, 856, -1672, 3593, -2487, -5731, -3085, -1987, 2382])

        vel2 = np.array([
            7628, 0, 4974, -4204, 4896, 5937, 6480, -6376, -7271, -7576])

        Vscale = -4

        vu_calcd = vel3dk_up(
            vel0, vel1, vel2, hdg, ptch, rll,
            beam1, beam2, beam3, beam4, Vscale)

        vu_expected = np.array([[
            -0.42452079, 0.8123358, 0.443022, 1.29753572, 0.99032304,
            1.21870746, 0.76989652, -0.19519593, -0.05982637, -0.51904823]])

        #np.testing.assert_array_almost_equal(vu_calcd, vu_expected)

    def test_zero_case(self, ):
        """
        Tests the case of all zero inputs to the main function
        vel_mag_correction.

        Implemented by:
            2014-02-26: Stuart Pearce
        """
        v0 = vel_mag_correction(0., 0., 0., 0., 0., 0., zflag=-1)
        np.testing.assert_allclose(v0, 0.0)

    def test_latlon_valid(self):
        """
        Tests for sub-functions valid_lat and valid_lon in the
        ion_functions.data.vel_functions module

        Implemented by:
            Luke Campbell sometime in late 2013
        """
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
        """
        Tests that if the sub-functions valid_lat and valid_lon return
        False within the main functions, that a ValueError exception is
        raised.

        Implemented by:
            2014-02-26: Stuart Pearce.
        """
        # try with Latitudes out of bounds
        self.assertRaises(ValueError, nobska_mag_corr_east, VE_NOBSKA[0], VN_NOBSKA[0], 91.0, -125.0, TS[0], DEPTH)
        self.assertRaises(ValueError, nobska_mag_corr_north, VE_NOBSKA[0], VN_NOBSKA[0], -100, -125.0, TS[0], DEPTH)
        self.assertRaises(ValueError, nortek_mag_corr_east, VE_NOBSKA[0], VN_NOBSKA[0], -90.01, -125.0, TS[0], DEPTH)
        self.assertRaises(ValueError, nortek_mag_corr_north, VE_NOBSKA[0], VN_NOBSKA[0], 120, -125.0, TS[0], DEPTH)
        # try with Longitudes out of bounds
        self.assertRaises(ValueError, nobska_mag_corr_east, VE_NOBSKA[0], VN_NOBSKA[0], 45.0, 180.01, TS[0], DEPTH)
        self.assertRaises(ValueError, nobska_mag_corr_north, VE_NOBSKA[0], VN_NOBSKA[0], 45.0, -180.1, TS[0], DEPTH)
        self.assertRaises(ValueError, nortek_mag_corr_east, VE_NOBSKA[0], VN_NOBSKA[0], 45.0, -210.0, TS[0], DEPTH)
        self.assertRaises(ValueError, nortek_mag_corr_north, VE_NOBSKA[0], VN_NOBSKA[0], 45.0, 200.0, TS[0], DEPTH)
