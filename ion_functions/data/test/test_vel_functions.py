#!/usr/bin/env python
"""
@package ion_functions.test.vel_functions
@file ion_functions/test/vel_functions.py
@author Stuart Pearce, Russell Desiderio
@brief Unit tests for vel_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
import pdb
from ion_functions.data.vel_functions import nobska_mag_corr_east, nobska_mag_corr_north
from ion_functions.data.vel_functions import nortek_mag_corr_east, nortek_mag_corr_north
from ion_functions.data.vel_functions import velpt_mag_corr_east, velpt_mag_corr_north
from ion_functions.data.vel_functions import velpt_up_vel
from ion_functions.data.vel_functions import vel3dk_east, vel3dk_north
from ion_functions.data.vel_functions import vel3dk_up
from ion_functions.data.vel_functions import valid_lat, valid_lon
from ion_functions.data.vel_functions import vel_mag_correction
from ion_functions.data.vel_functions import fsi_acm_rsn_east, fsi_acm_rsn_north
from ion_functions.data.vel_functions import fsi_acm_sio_east, fsi_acm_sio_north
from ion_functions.data.vel_functions import fsi_acm_up_profiler_ascending
from ion_functions.data.vel_functions import fsi_acm_up_profiler_descending
from ion_functions.data.vel_functions import fsi_acm_nautical_heading, fsi_acm_compass_cal
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
VE_VELPT = VE_NOBSKA * 10.
VN_VELPT = VN_NOBSKA * 10.
VU_VELPT = VU_NOBSKA * 10.

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

    def test_velpt(self):
        """
        Tests functions velpt_mag_corr_east and velpt_mag_corr_north
        from the ion_functions.data.vel_functions module using test data
        from the VELPTMN DPS.

        Implemented by:

        2013-04-17: Stuart Pearce. Initial code.
        2013-04-24: Stuart Pearce. Changed to be general for all velocity
                    instruments.
        2014-02-05: Christopher Wingard. Edited to use magnetic corrections in
                    the generic_functions module.
        2014-10-29: Stuart Pearce. Adds the velpt functions.

        References:

            OOI (2012). Data Product Specification for Mean Point Water
                Velocity. Document Control Number 1341-00790.
                https://alfresco.oceanobservatories.org/ (See: Company Home
                >> OOI >> Controlled >> 1000 System Level >>
                1341-00790_Data_Product_SPEC_VELPTMN_OOI.pdf)
        """

        ve_cor = velpt_mag_corr_east(
            VE_VELPT, VN_VELPT, LAT, LON, TS, DEPTH)
        vn_cor = velpt_mag_corr_north(
            VE_VELPT, VN_VELPT, LAT, LON, TS, DEPTH)

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

        beams = np.array([
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0]])

        hdg = 10. * np.array(
            [0., 36., 72., 108., 144., 180., 216., 252., 288., 324.])

        ptch = 10. * np.array([
            -20.10, -5.82, 0.0, -10.26, -8.09, -1.36, 5.33, 9.96, 13.68, 26.77])

        rll = 10. * np.array([
            -15.0, -11.67, -8.33, -5.0, 0.0, 1.67, 5.0, 8.33, 11.67, 15.0])

        vel0 = np.array([
            0, 6305, 6000, 7810, 10048, 11440, -10341, -10597, 9123, -12657])

        vel1 = np.array([
            0, 1050, 856, -1672, 3593, -2487, -5731, -3085, -1987, 2382])

        vel2 = np.array([
            7628, 0, 4974, -4204, 4896, 5937, 6480, -6376, -7271, -7576])

        Vscale = np.array([-4])

        ve_cor = vel3dk_east(
            vel0, vel1, vel2, hdg, ptch, rll, beams, LAT, LON, TS, Vscale, DEPTH)
        vn_cor = vel3dk_north(
            vel0, vel1, vel2, hdg, ptch, rll, beams, LAT, LON, TS, Vscale, DEPTH)

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
        beams = np.array([
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [1, 2, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0],
            [2, 3, 4, 0, 0]])

        hdg = 10. * np.array([0., 36., 72., 108., 144., 180., 216., 252., 288., 324.])

        ptch = 10. * np.array([
            -20.10, -5.82, 0.0, -10.26, -8.09, -1.36, 5.33, 9.96, 13.68, 26.77])/0.1

        rll = 10. * np.array([
            -15.0, -11.67, -8.33, -5.0, 0.0, 1.67, 5.0, 8.33, 11.67, 15.0])/0.1

        vel0 = np.array([
            0, 6305, 6000, 7810, 10048, 11440, -10341, -10597, 9123, -12657])

        vel1 = np.array([
            0, 1050, 856, -1672, 3593, -2487, -5731, -3085, -1987, 2382])

        vel2 = np.array([
            7628, 0, 4974, -4204, 4896, 5937, 6480, -6376, -7271, -7576])

        Vscale = -4

        vu_calcd = vel3dk_up(
            vel0, vel1, vel2, hdg, ptch, rll,
            beams, Vscale)

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

    # vel3d-a,l unit tests
    def test_fsi_acm(self):
        """
        Description:

            Tests fsi_acm functions which calculate velptmn data products from vel3d series A
            and L instruments. No velocity data for unit tests were provided. I constructed
            these unit tests with reference to the definitions and conventions regarding beam
            and instrument coordinates and heading designations so that the test results could
            be visualized and predicted geometrically.

        Implemented by:

            2015-02-20: Russell Desiderio. Initial code.
            2015-05-29: Russell Desiderio.
                        Time vectorized cal coeffs hdg_cal, hx_cal, hy_cal (fsi_acm_rsn functions):
                            Single time point case, shapes were changed from [8,] to [1,8].
                            Multiple time point case, these arrays were tiled in the time dimension.
                        Added two time point, different cals unit test.

        References:

            OOI (2015). Data Product Specification for Mean Point Water Velocity
                Data from FSI Acoustic Current Meters. Document Control Number
                1341-00792. https://alfresco.oceanobservatories.org/  (See:
                Company Home >> OOI >> Controlled >> 1000 System Level >>
                1341-00792_Data_Product_SPEC_VELPTMN_ACM_OOI.pdf)

            OOI (2015). 1341-00792_VELPTMN Artifact: McLane Moored Profiler User Manual.
                https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
                >> REFERENCE >> Data Product Specification Artifacts >> 1341-00792_VELPTMN >>
                MMP-User Manual-Rev-E-WEB.pdf)
        """
        ###
        ### test upwards velocity functions
        ###
        # for these tests the ascending and descending velocities will not always agree,
        # because w_asc uses vp4 and not vp2, w_dsc vp2 and not vp4.
        # the expected values have units of [m/s]
        xpctd_asc = np.array([-2.83, 0.00, -1.41, 1.41, -2.83, 0.00, -1.41, 1.41,
                              -1.41, 1.41, -0.00, 2.83, -1.41, 1.41, -0.00, 2.83])
        xpctd_dsc = np.array([2.83,  2.83,  1.41,  1.41, -0.00, -0.00, -1.41, -1.41,
                              1.41,  1.41,  0.00,  0.00, -1.41, -1.41, -2.83, -2.83])
        # the raw beam values have units of [cm/s]
        vp1 = 100.0 * np.array([+1, +1, +1, +1, +1, +1, +1, +1, -1, -1, -1, -1, -1, -1, -1, -1])
        vp2 = 100.0 * np.array([+1, +1, +1, +1, -1, -1, -1, -1, +1, +1, +1, +1, -1, -1, -1, -1])
        vp3 = 100.0 * np.array([+1, +1, -1, -1, +1, +1, -1, -1, +1, +1, -1, -1, +1, +1, -1, -1])
        vp4 = 100.0 * np.array([+1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1])
        calc = fsi_acm_up_profiler_ascending(vp1, vp3, vp4)
        np.testing.assert_allclose(calc, xpctd_asc, rtol=0.0, atol=1.e-2)
        calc = fsi_acm_up_profiler_descending(vp1, vp2, vp3)
        np.testing.assert_allclose(calc, xpctd_dsc, rtol=0.0, atol=1.e-2)

        #############################
        ### fsi_acm_compass_cal tests
        #############################
        #
        ### these are single time point cases:
        ### calcoeffs will be supplied by CI as row vectors, not 1D arrays.
        #
        # (a)  trivial perfect cal case
        # offsets are zero, scale factors are 1, bias is 0
        xpctd = (np.array([0.]), np.array([0.]), np.array([1.]), np.array([1.]), np.array([0.]))
        hdg_percal = np.array([[0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]])
        hx_percal = np.array([[1.0, 1.0, 0.0, -1.0, -1.0, -1.0, +0.0, +1.0]])
        hy_percal = np.array([[0.0, 1.0, 1.0, +1.0, +0.0, -1.0, -1.0, -1.0]])
        calc = fsi_acm_compass_cal(hdg_percal, hx_percal, hy_percal)
        np.testing.assert_allclose(calc, xpctd, rtol=0.0, atol=1.e-12)
        # (b)  using MMP manual example (p 8-19)
        #
        # the expected values are hx and hy offsets, hx and hy scaling factors, and compass bias.
        # as in case (a), the xpctd arguments are actually each 1-element 1D arrays; this is
        # "up-dimensioned" by the action of the assertion testing to a column vector, so,
        # equivalently:
        xpctd_MMP_cal = np.array([[0.01763750], [-0.04673750], [0.34686250],
                                  [0.33813750], [7.46424566]])
        # convert the manual's cartesian reference directions to nautical
        hdg_cal = np.array([[0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]])
        hx_cal = np.array([[0.3645, 0.2383, -0.0246, -0.2508, -0.3212, -0.1979, 0.0511, 0.2817]])
        hy_cal = np.array([[0.0008, 0.2321, 0.2914, 0.1515, -0.1064, -0.3218, -0.3652, -0.2563]])
        calc = fsi_acm_compass_cal(hdg_cal, hx_cal, hy_cal)
        np.testing.assert_allclose(calc, xpctd_MMP_cal, rtol=0.0, atol=1.e-6)

        #############################################################################
        # for all but the last of the remaining test sets, 4 time points will be used.
        # therefore, the calcoeff arrays will be time-vectorized by tiling so that
        # the lead dimension of the resulting 2D array will be 4.
        npoints = 4
        hdg_percal_tv = np.tile(hdg_percal, (npoints, 1))
        hx_percal_tv = np.tile(hx_percal, (npoints, 1))
        hy_percal_tv = np.tile(hy_percal, (npoints, 1))
        hdg_cal_tv = np.tile(hdg_cal, (npoints, 1))
        hx_cal_tv = np.tile(hx_cal, (npoints, 1))
        hy_cal_tv = np.tile(hy_cal, (npoints, 1))
        #############################################################################

        #################################
        ### fsi_acm_nautical_heading test
        #################################
        rt3 = np.sqrt(3.0)
        # pick a point in each quadrant.
        #
        # heading values given are of the direction of the x-axis in the instrument coordinate
        # system, before cal correction, ccw relative to magnetic east (cartesian) or cw relative
        # to magnetic north (nautical).
        #
        # note that hx and hy define the direction of magnetic north so that for (hx, hy) =
        # (rt3, 1), magnetic north lies in the 1st quadrant of the instrument coordinate system
        # 30 degrees ccw from the positive x-axis; therefore the positive x-axis lies 30 degrees
        # cw from magnetic north, which is the nautical heading, and 60 degrees ccw from magnetic
        # east, which is the cartesian heading.
        #
        # the 4 (hx,hy) values track in 90 degree steps ccw in the instrument coordinate system,
        # which means that magnetic north is also tracking ccw relative to the x-axis. the result is
        # that the nautical heading of the x-axis tracks cw in +90 steps, cartesian ccw in -90 steps.
        # cartesian      60   -30  -120   150
        # nautical       30   120   210   300
        hx = np.array([+rt3, -1.0, -rt3, +1.0])
        hy = np.array([+1.0, +rt3, -1.0, -rt3])
        # note: for this test these values are not normalized to lie on the unit circle.
        #
        #       the perfect cal unit test results are independent of scaling because the
        #       relevant calculation involves the atan2 of the ratio of the direction cosines
        #       so that any scaling of the raw (hx, hy) values divides out to 1.
        #
        #       real calibrations result in corrections to hx and hy which will change the
        #       value of the direction cosine ratio resulting in answers dependent on the
        #       scaling of the raw test values.
        #
        # (a)  trivial perfect cal case
        xpctd = np.array([30.0, 120.0, 210.0, 300.0])
        calc = fsi_acm_nautical_heading(hx, hy, hdg_percal_tv, hx_percal_tv, hy_percal_tv)
        np.testing.assert_allclose(calc, xpctd, rtol=0.0, atol=1.e-12)
        # (b)  using calibration data from the example in the MMP manual
        xpctd_MMP_hdg = np.array([24.59490, 111.68440, 201.73555, 292.14246])
        calc = fsi_acm_nautical_heading(hx, hy, hdg_cal_tv, hx_cal_tv, hy_cal_tv)
        np.testing.assert_allclose(calc, xpctd_MMP_hdg, rtol=0.0, atol=1.e-4)

        ####################################
        ### fsi_acm_sio east and north tests
        ####################################
        #
        # the sio instruments output heading directly, and do not correct the field
        # data using compass cals. therefore, there are no calcoeffs to time-vectorize
        # in the tests in this section.
        #
        # figure out OOI timestamp for Feb 1, 2011 at midnight [sec since 1900_01_01]
        ts_unix_time = 1296518400.0  # midnight, Feb 1, 2011, sec since 1970_01_01
        delta_epoch = 2208988800.0  # OOI time for midnight 1970_01_01 [sec since 1900_01_01]
        ref_ts = ts_unix_time + delta_epoch
        # pick a test lat and lon which with this timestamp gives a magnetic variation of
        # close to 15 degrees.
        ref_lat = 37.768812
        ref_lon = -135.0
        # this will give a magnetic variation value of magvar = 15.000000171 (WMM).

        # the standard convention is that a positive magnetic variation delta indicates that
        # magnetic north lies delta degrees to the east (clockwise) of true north. this means
        # that a magnetic (nautical) heading of 0 degrees is at a true geographic (nautical)
        # heading of +delta degrees. therefore, to correct nautical headings, the magnetic
        # declination value magvar is added.
        #
        # because cartesian headings are defined as positive in the ccw direction, correcting
        # cartesian headings for magnetic declination involves subtracting the magvar value.

        # the MMP manual uses the opposite sign convention in their definition of magnetic variation.

        # set up 4 test raw beam velocity pairs (vp1,vp3) with water speeds of
        # sqrt(vp1*vp1 + vp3*vp3) = 100 cm/s = 1 m/s.
        # these correlate with the following cartesian angles from positive x-axis ccw to water
        # velocity vector in instrument coordinates:
        # ccw from x:           165   255   345    75
        vp1 = 100.0 * np.array([rt3, -1.0, -rt3, +1.0]) / 2.0
        vp3 = 100.0 * np.array([1.0, +rt3, -1.0, -rt3]) / 2.0
        # note that the v1 and vp3 beam paths are oriented at 45 degrees to the x-axis, so that
        # when a water velocity vector has 30 and/or 60 degree components in the (vp1,vp3)
        # beam coordinate system (so, sqrt(3)/2 and 1/2) this will result in 15 degree offsets
        # in the instrument coordinate system above [165 255 345 75] which also sets the choice
        # of value of magnetic variation for these tests (15 degrees) so that the corrected
        # heading values will be integral multiples of 30 degrees.

        # construct lat, lon, and ts vectors
        lat = np.zeros(4) + ref_lat
        lon = np.zeros(4) + ref_lon
        ts = np.zeros(4) + ref_ts

        # (a) test when x-axis is always at magnetic north
        hdg = np.array([0., 0., 0., 0.])
        # because magvar = +15 degrees, the x-axis is 15 degrees cw from true (geographic)
        # north and 75 degrees ccw from true east.
        # therefore the xpctd (u,v) ENU cartesian directions are 75 + [165, 255, 345, 75]
        xpctd_angles = np.array([240.0, 330.0, 60.0, 150.0])
        # knowing these angles, the velocity component magnitudes can be calculated because
        # the raw beam velocities were set up to have water velocity speeds of 1 m/s.
        # the units of the data product velocities are m/s, so:
        u_xpctd_hdg_0 = np.cos(np.deg2rad(xpctd_angles))
        v_xpctd_hdg_0 = np.sin(np.deg2rad(xpctd_angles))
        # which should equal
        u_xpctd = np.array([-1.0, +rt3, 1.0, -rt3]) / 2.0
        v_xpctd = np.array([-rt3, -1.0, rt3, +1.0]) / 2.0
        # make sure
        np.testing.assert_allclose(u_xpctd_hdg_0, u_xpctd, rtol=0.0, atol=1.e-8)
        np.testing.assert_allclose(v_xpctd_hdg_0, v_xpctd, rtol=0.0, atol=1.e-8)
        # unit test:
        u_calc = fsi_acm_sio_east(vp1, vp3, hdg, lat, lon, ts)
        v_calc = fsi_acm_sio_north(vp1, vp3, hdg, lat, lon, ts)
        calc_angles = np.mod(360.0 + np.rad2deg(np.arctan2(v_calc, u_calc)), 360.0)
        np.testing.assert_allclose(calc_angles, xpctd_angles, rtol=0.0, atol=1.e-4)
        np.testing.assert_allclose(u_calc, u_xpctd, rtol=0.0, atol=1.e-6)
        np.testing.assert_allclose(v_calc, v_xpctd, rtol=0.0, atol=1.e-6)
        # (b) non-zero hdg test
        # pick the heading which should rotate the instrument coordinate system so that the
        # water velocity vectors will lie along the cardinal compass points N, W, S, E.
        hdg = np.array([150., 150., 150., 150])
        # this should result in ENU angles of 90, 180, 270, 360. so:
        u_xpctd = np.array([0.0, -1.0, +0.0, 1.0])  # cos(angle)
        v_xpctd = np.array([1.0, +0.0, -1.0, 0.0])  # sin(angle)
        u_calc = fsi_acm_sio_east(vp1, vp3, hdg, lat, lon, ts)
        v_calc = fsi_acm_sio_north(vp1, vp3, hdg, lat, lon, ts)
        np.testing.assert_allclose(u_calc, u_xpctd, rtol=0.0, atol=1.e-6)
        np.testing.assert_allclose(v_calc, v_xpctd, rtol=0.0, atol=1.e-6)

        ####################################
        ### fsi_acm_rsn east and north tests
        ####################################
        # (a) use fsi_acm_compass_cal and fsi_acm_nautical_heading test set-up along with sio
        #     vp1 and vp3 values and a perfect cal so that the expected values can be predicted
        #     on geometrical grounds.
        #
        #     the 4 (vp1, vp3) pairs have cartesian angles from the x-axis of [165, 255,  345,  75]
        #     the 4 (hx, hy) pairs result in cartesian headings of            [ 60, -30, -120, 150]
        #     the magnetic variation is 15 degrees, so subtract (cartesian)   [-15, -15,  -15, -15]
        #     add columns; so, all 4 angles should be 210 degrees             [210, 210,  210, 210]
        u_xpctd = np.zeros(4) - rt3/2.0  # cos(210)
        v_xpctd = np.zeros(4) - 0.5      # sin(210)
        u_calc = fsi_acm_rsn_east(vp1, vp3, hx, hy, hdg_percal_tv, hx_percal_tv, hy_percal_tv, lat, lon, ts)
        v_calc = fsi_acm_rsn_north(vp1, vp3, hx, hy, hdg_percal_tv, hx_percal_tv, hy_percal_tv, lat, lon, ts)
        np.testing.assert_allclose(u_calc, u_xpctd, rtol=0.0, atol=1.e-6)
        np.testing.assert_allclose(v_calc, v_xpctd, rtol=0.0, atol=1.e-6)
        # (b) fix (vp1, vp3) and vary (hx, hy), perfect cal case
        # ccw from x:           165  165  165  165
        vp1 = 100.0 * np.array([rt3, rt3, rt3, rt3]) / 2.0
        vp3 = 100.0 * np.array([1.0, 1.0, 1.0, 1.0]) / 2.0
        # cartesian     90    0,  -90,  180
        # nautical       0,  90,  180,  270
        hx = np.array([1.0, 0.0, -1.0, +0.0])
        hy = np.array([0.0, 1.0, +0.0, -1.0])
        # predicted angles:
        #     (vp1, vp3) cartesian angles: [ 165, 165, 165, 165]
        #     ( hx,  hy) cartesian angles: [  90,   0, -90, 180]
        #     subtract magnetic variation: [ -15, -15, -15, -15]
        #     predicted:                   [-120, 150,  60, -30]
        u_xpctd = np.array([-1.0, -rt3, 1.0, +rt3]) / 2.0
        v_xpctd = np.array([-rt3, +1.0, rt3, -1.0]) / 2.0
        u_calc = fsi_acm_rsn_east(vp1, vp3, hx, hy, hdg_percal_tv, hx_percal_tv, hy_percal_tv, lat, lon, ts)
        v_calc = fsi_acm_rsn_north(vp1, vp3, hx, hy, hdg_percal_tv, hx_percal_tv, hy_percal_tv, lat, lon, ts)
        np.testing.assert_allclose(u_calc, u_xpctd, rtol=0.0, atol=1.e-6)
        np.testing.assert_allclose(v_calc, v_xpctd, rtol=0.0, atol=1.e-6)
        # (c) fix (hx, hy) and vary (vp1, vp3), perfect cal
        # ccw from x:           165   255   345    75
        vp1 = 100.0 * np.array([rt3, -1.0, -rt3, +1.0]) / 2.0
        vp3 = 100.0 * np.array([1.0, +rt3, -1.0, -rt3]) / 2.0
        # cartesian   -120 -120, -120, -120
        # nautical     210, 210,  210,  210
        hx = np.zeros(4) - rt3
        hy = np.zeros(4) - 1.0
        # predicted angles:
        #     (vp1, vp3) cartesian angles: [ 165,  255,  345,   75]
        #     ( hx,  hy) cartesian angles: [-120, -120, -120, -120]
        #     subtract magnetic variation: [ -15,  -15,  -15,  -15]
        #     predicted:                   [  30,  120,  210,  300]
        u_xpctd = np.array([rt3, -1.0, -rt3, +1.0]) / 2.0
        v_xpctd = np.array([1.0, +rt3, -1.0, -rt3]) / 2.0
        u_calc = fsi_acm_rsn_east(vp1, vp3, hx, hy, hdg_percal_tv, hx_percal_tv, hy_percal_tv, lat, lon, ts)
        v_calc = fsi_acm_rsn_north(vp1, vp3, hx, hy, hdg_percal_tv, hx_percal_tv, hy_percal_tv, lat, lon, ts)
        np.testing.assert_allclose(u_calc, u_xpctd, rtol=0.0, atol=1.e-6)
        np.testing.assert_allclose(v_calc, v_xpctd, rtol=0.0, atol=1.e-6)
        # (d) non-perfect cal case, work through by hand
        #     single time point.
        vp1 = 100.0 * rt3 / 2.0
        vp3 = 100.0 * 1.0 / 2.0
        hx = -rt3
        hy = -1.0
        #
        # make sure the shapes of the lat, lon, and ts variables are (1,) as CI will supply them.
        u_calc = fsi_acm_rsn_east(vp1, vp3, hx, hy, hdg_cal, hx_cal, hy_cal, lat[[0]], lon[[0]], ts[[0]])
        v_calc = fsi_acm_rsn_north(vp1, vp3, hx, hy, hdg_cal, hx_cal, hy_cal, lat[[0]], lon[[0]], ts[[0]])
        # step through by stages
        u_inst_coord = -(vp1 + vp3) / np.sqrt(2.0)
        v_inst_coord = (vp1 - vp3) / np.sqrt(2.0)
        angle_inst_coord = np.arctan2(v_inst_coord, u_inst_coord)
        hdg_naut = xpctd_MMP_hdg[2]   # calculated in a previous test; already contains compass bias
        magvar = 15.0                 # reference magnetic variation by design
        # using the cartesian angle convention, convert nautical heading and flip magvar sign.
        hdg_cart = np.deg2rad((90.0 - hdg_naut) - magvar)
        angle_cart = angle_inst_coord + hdg_cart
        #
        speed = np.sqrt(u_inst_coord * u_inst_coord + v_inst_coord * v_inst_coord) / 100.0
        u_xpctd = speed * np.cos(angle_cart)
        v_xpctd = speed * np.sin(angle_cart)
        np.testing.assert_allclose(u_calc, u_xpctd, rtol=0.0, atol=1.e-6)
        np.testing.assert_allclose(v_calc, v_xpctd, rtol=0.0, atol=1.e-6)

        # (e) two time point case, different cals.
        #         1st point: case (d) immediately above
        #         2nd point: 1st point of case (c)
        #     these cases have identical vp1, vp3, hx, hy, lat, lon, and ts, but different cals.
        vp1 = np.array([vp1, vp1])
        vp3 = np.array([vp3, vp3])
        hx = np.array([hx, hx])
        hy = np.array([hy, hy])
        hdg_2cals = np.vstack([hdg_cal, hdg_percal])
        hx_2cals = np.vstack([hx_cal, hx_percal])
        hy_2cals = np.vstack([hy_cal, hy_percal])
        u_xpctd = np.array([u_xpctd, rt3 / 2.0])
        v_xpctd = np.array([v_xpctd, 0.5])
        u_calc = fsi_acm_rsn_east(vp1, vp3, hx, hy, hdg_2cals, hx_2cals, hy_2cals, lat[0:2], lon[0:2], ts[0:2])
        v_calc = fsi_acm_rsn_north(vp1, vp3, hx, hy, hdg_2cals, hx_2cals, hy_2cals, lat[0:2], lon[0:2], ts[0:2])

        np.testing.assert_allclose(u_calc, u_xpctd, rtol=0.0, atol=1.e-6)
        np.testing.assert_allclose(v_calc, v_xpctd, rtol=0.0, atol=1.e-6)

