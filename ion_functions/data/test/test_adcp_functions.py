"""
@package ion_functions.test.adcp_functions
@file ion_functions/test/test_adcp_functions.py
@author Christopher Wingard, Russell Desiderio, Craig Risien
@brief Unit tests for adcp_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np

from ion_functions.data import adcp_functions as af
from ion_functions.data.adcp_functions import ADCP_FILLVALUE
from ion_functions.data.generic_functions import SYSTEM_FILLVALUE


@attr('UNIT', group='func')
class TestADCPFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        """
        Implemented by:
            2014-02-06: Christopher Wingard. Initial Code.
            2015-06-12: Russell Desiderio. Changed raw beam data to type int. This
                        change did not affect any previously written unit tests.

        """
        # set test inputs -- values from DPS
        self.b1 = np.array([[-0.0300, -0.2950, -0.5140, -0.2340, -0.1880,
                            0.2030, -0.3250,  0.3050, -0.2040, -0.2940]]) * 1000
        self.b2 = np.array([[0.1800, -0.1320,  0.2130,  0.3090,  0.2910,
                            0.0490,  0.1880,  0.3730, -0.0020, 0.1720]]) * 1000
        self.b3 = np.array([[-0.3980, -0.4360, -0.1310, -0.4730, -0.4430,
                            0.1880, -0.1680,  0.2910, -0.1790, 0.0080]]) * 1000
        self.b4 = np.array([[-0.2160, -0.6050, -0.0920, -0.0580,  0.4840,
                            -0.0050,  0.3380,  0.1750, -0.0800, -0.5490]]) * 1000
        # the data type of the raw beam velocities is int;
        # set b1-b4 to int so that fill replacement can be tested.
        self.b1 = self.b1.astype(int)
        self.b2 = self.b2.astype(int)
        self.b3 = self.b3.astype(int)
        self.b4 = self.b4.astype(int)
        #
        self.echo = np.array([[0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250]])
        self.sfactor = 0.45
        # units of compass data are in centidegrees.
        self.heading = 9841
        self.pitch = 69
        self.roll = -254
        self.orient = 1
        self.lat = 50.0000
        self.lon = -145.0000
        self.depth = 0.0
        self.ntp = 3545769600.0    # May 12, 2012

        # set expected results -- velocity profiles in earth coordinates
        # (values in DPS)
        self.uu = np.array([[0.2175, -0.2814, -0.1002, 0.4831, 1.2380,
                            -0.2455, 0.6218, -0.1807, 0.0992, -0.9063]])
        self.vv = np.array([[-0.3367, -0.1815, -1.0522, -0.8676, -0.8919,
                            0.2585, -0.8497, -0.0873, -0.3073, -0.5461]])
        self.ww = np.array([[0.1401,  0.3977,  0.1870,  0.1637,  0.0091,
                            -0.1290,  0.0334, -0.3017, 0.1384, 0.1966]])

        # set expected results -- magnetic variation correction applied
        # (computed in Matlab using above values and mag_var.m)
        self.uu_cor = np.array([[0.1099, -0.3221, -0.4025, 0.2092, 0.9243,
                                -0.1595, 0.3471, -0.1983, 0.0053, -1.0261]])
        self.vv_cor = np.array([[-0.3855, -0.0916, -0.9773, -0.9707, -1.2140,
                                0.3188, -0.9940, -0.0308, -0.3229, -0.2582]])

        # set the expected results -- error velocity
        self.ee = np.array([[0.789762, 0.634704, -0.080630, 0.626434, 0.064090,
                             0.071326, -0.317352, 0.219148, 0.054787, 0.433129]])
        # set the expected results -- echo intensity conversion from counts to dB
        self.dB = np.array([[0.00, 11.25, 22.50, 33.75, 45.00, 56.25, 67.50,
                            78.75, 90.00, 101.25, 112.50]])

    def test_adcp_beam(self):
        """
        Directly tests DPA functions adcp_beam_eastward, adcp_beam_northward,
        adcp_beam_vertical, and adcp_beam_error.

        Tests adcp_beam2ins, adcp_ins2earth and magnetic_correction functions
        for ADCPs that output data in beam coordinates. All three functions
        must return the correct output for final tests cases to work.

        Values based on those defined in DPS:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by:

            2013-04-10: Christopher Wingard. Initial code.
            2014-02-06: Christopher Wingard. Added tests to confirm arrays of
                        arrays can be processed (in other words, vectorized the
                        code).
            2015-06-23: Russell Desiderio. Revised documentation. Added unit test
                        for the function adcp_beam_error.

        Notes:

            The original suite of tests within this function did not provide a
            test for adcp_beam_error. However, adcp_beam_error and vadcp_beam_error
            are identical functions, and vadcp_beam_error is implicitly tested in the
            test_vadcp_beam function when the 4th output argument of adcp_beam2inst
            is tested. Therefore values to directly test adcp_beam_error were
            then derived from the function itself and included as part of the unit
            test within this code (test_adcp_beam).
        """
        # single record case
        got_uu_cor = af.adcp_beam_eastward(self.b1, self.b2, self.b3, self.b4,
                                           self.heading, self.pitch, self.roll, self.orient,
                                           self.lat, self.lon, self.depth, self.ntp)
        got_vv_cor = af.adcp_beam_northward(self.b1, self.b2, self.b3, self.b4,
                                            self.heading, self.pitch, self.roll, self.orient,
                                            self.lat, self.lon, self.depth, self.ntp)
        got_ww = af.adcp_beam_vertical(self.b1, self.b2, self.b3, self.b4,
                                       self.heading, self.pitch, self.roll, self.orient)
        got_ee = af.adcp_beam_error(self.b1, self.b2, self.b3, self.b4)

        # test results
        np.testing.assert_array_almost_equal(got_uu_cor, self.uu_cor, 4)
        np.testing.assert_array_almost_equal(got_vv_cor, self.vv_cor, 4)
        np.testing.assert_array_almost_equal(got_ww, self.ww, 4)
        np.testing.assert_array_almost_equal(got_ee, self.ee, 4)

        # reset the test inputs for multiple records
        b1 = np.tile(self.b1, (24, 1))
        b2 = np.tile(self.b2, (24, 1))
        b3 = np.tile(self.b3, (24, 1))
        b4 = np.tile(self.b4, (24, 1))
        heading = np.ones(24, dtype=np.int) * self.heading
        pitch = np.ones(24, dtype=np.int) * self.pitch
        roll = np.ones(24, dtype=np.int) * self.roll
        orient = np.ones(24, dtype=np.int) * self.orient
        lat = np.ones(24) * self.lat
        lon = np.ones(24) * self.lon
        depth = np.ones(24) * self.depth
        ntp = np.ones(24) * self.ntp

        # reset outputs for multiple records
        uu_cor = np.tile(self.uu_cor, (24, 1))
        vv_cor = np.tile(self.vv_cor, (24, 1))
        ww = np.tile(self.ww, (24, 1))
        ee = np.tile(self.ee, (24, 1))

        # multiple record case
        got_uu_cor = af.adcp_beam_eastward(b1, b2, b3, b4,
                                           heading, pitch, roll, orient,
                                           lat, lon, depth, ntp)
        got_vv_cor = af.adcp_beam_northward(b1, b2, b3, b4,
                                            heading, pitch, roll, orient,
                                            lat, lon, depth, ntp)
        got_ww = af.adcp_beam_vertical(b1, b2, b3, b4,
                                       heading, pitch, roll, orient)
        got_ee = af.adcp_beam_error(b1, b2, b3, b4)

        # test results
        np.testing.assert_array_almost_equal(got_uu_cor, uu_cor, 4)
        np.testing.assert_array_almost_equal(got_vv_cor, vv_cor, 4)
        np.testing.assert_array_almost_equal(got_ww, ww, 4)
        np.testing.assert_array_almost_equal(got_ee, ee, 4)

    def test_adcp_beam_with_fill(self):
        """
        Directly tests DPA functions adcp_beam_eastward, adcp_beam_northward,
        adcp_beam_vertical, and adcp_beam_error when system fill values and
        ADCP fill values (bad value sentinels) are present in the data stream.

        Non-fill values are based on those used in test_adcp_beam in this module.

        Implemented by:

            2013-06-24: Russell Desiderio. Initial code.

        Notes:

        """
        # for convenience
        sfill = SYSTEM_FILLVALUE
        afill = ADCP_FILLVALUE

        ### set input data
        # units of compass data are in centidegrees.
        heading = np.array([9841])
        pitch = np.array([69])
        roll = np.array([-254])

        missingroll = np.array([sfill])

        orient = np.array([1])
        lat = np.array([50.0000])
        lon = np.array([-145.0000])
        depth = np.array([0.0])
        ntp = np.array([3545769600.0])  # May 12, 2012

        ###
        # for positional clarity, input beam and expected velocities will be explicitly
        # enumerated for each single time record test case.
        ###

        ### single time record case; missing roll data
        ## the ADCP does not use its bad flag sentinel for compass data, only beam data.
        ## however, it is possible that CI could supply the system fillvalue for missing compass data.
        # input data
        # beam velocity units are mm/s
        b1_x1 = np.array([[-30, -295, -514, -234, -188, 203, -325, 305, -204, -294]])
        b2_x1 = np.array([[180, -132, 213, 309,  291, 49, 188, 373, -2, 172]])
        b3_x1 = np.array([[-398, -436, -131, -473, -443, 188, -168,  291, -179, 8]])
        b4_x1 = np.array([[-216, -605, -92, -58, 484, -5, 338, 175, -80, -549]])
        # expected results if all good beam and compass data
        # these will be used later in the multiple time record test
        uu_x0 = np.array([[0.1099, -0.3221, -0.4025, 0.2092, 0.9243,
                           -0.1595, 0.3471, -0.1983, 0.0053, -1.0261]])
        vv_x0 = np.array([[-0.3855, -0.0916, -0.9773, -0.9707, -1.2140,
                           0.3188, -0.9940, -0.0308, -0.3229, -0.2582]])
        ww_x0 = np.array([[0.1401,  0.3977,  0.1870,  0.1637,  0.0091,
                           -0.1290,  0.0334, -0.3017, 0.1384, 0.1966]])
        ee_x0 = np.array([[0.789762, 0.634704, -0.080630, 0.626434, 0.064090,
                           0.071326, -0.317352, 0.219148, 0.054787, 0.433129]])
        # expected results for all good beam data, missing roll data;
        # nans for all results except for the error velocity, which does not depend on the compass
        uu_x1 = uu_x0 * np.nan
        vv_x1 = vv_x0 * np.nan
        ww_x1 = ww_x0 * np.nan
        ee_x1 = np.copy(ee_x0)

        uu_calc = af.adcp_beam_eastward(b1_x1, b2_x1, b3_x1, b4_x1, heading, pitch,
                                        missingroll,
                                        orient, lat, lon, depth, ntp)
        vv_calc = af.adcp_beam_northward(b1_x1, b2_x1, b3_x1, b4_x1, heading, pitch,
                                         missingroll,
                                         orient, lat, lon, depth, ntp)
        ww_calc = af.adcp_beam_vertical(b1_x1, b2_x1, b3_x1, b4_x1, heading, pitch,
                                        missingroll,
                                        orient)
        ee_calc = af.adcp_beam_error(b1_x1, b2_x1, b3_x1, b4_x1)
        # test results
        np.testing.assert_array_almost_equal(uu_calc, uu_x1, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_x1, 4)
        np.testing.assert_array_almost_equal(ww_calc, ww_x1, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_x1, 4)

        ### single time record case; missing and bad-flagged beam data, good compass data
        # input data
        b1_x2 = np.array([[sfill, -295, -514, -234,  -188, 203, -325, afill, -204,  -294]])
        b2_x2 = np.array([[sfill, -132,  213,  309,   291,  49,  188, afill,   -2, sfill]])
        b3_x2 = np.array([[sfill, -436, -131, -473,  -443, 188, -168, afill, -179,     8]])
        b4_x2 = np.array([[sfill, -605,  -92,  -58, afill,  -5,  338, afill,  -80,  -549]])
        # expected
        uu_x2 = np.array([[np.nan, -0.3221, -0.4025,  0.2092, np.nan,
                           -0.1595,  0.3471,  np.nan,  0.0053, np.nan]])
        vv_x2 = np.array([[np.nan, -0.0916, -0.9773, -0.9707, np.nan,
                           0.3188, -0.9940,  np.nan, -0.3229, np.nan]])
        ww_x2 = np.array([[np.nan,  0.3977,  0.1870,  0.1637, np.nan,
                           -0.1290,  0.0334,  np.nan,  0.1384, np.nan]])
        ee_x2 = np.array([[np.nan, 0.634704, -0.080630, 0.626434, np.nan,
                           0.071326, -0.317352, np.nan, 0.054787, np.nan]])
        # calculated
        uu_calc = af.adcp_beam_eastward(b1_x2, b2_x2, b3_x2, b4_x2,
                                        heading, pitch, roll, orient,
                                        lat, lon, depth, ntp)
        vv_calc = af.adcp_beam_northward(b1_x2, b2_x2, b3_x2, b4_x2,
                                         heading, pitch, roll, orient,
                                         lat, lon, depth, ntp)
        ww_calc = af.adcp_beam_vertical(b1_x2, b2_x2, b3_x2, b4_x2,
                                        heading, pitch, roll, orient)
        ee_calc = af.adcp_beam_error(b1_x2, b2_x2, b3_x2, b4_x2)
        # test results
        np.testing.assert_array_almost_equal(uu_calc, uu_x2, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_x2, 4)
        np.testing.assert_array_almost_equal(ww_calc, ww_x2, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_x2, 4)

        ### multiple (5) record case
        ## reset the test inputs for 5 time records
        #     1st time record is the bad/missing beam data case above
        #     2nd time record is a missing heading data case
        #     3rd time record is all good data
        #     4th time record is bad/missing beam and missing pitch data.
        #     5th time record is missing orientation data
        b1 = np.vstack((b1_x2, b1_x1, b1_x1, b1_x2, b1_x1))
        b2 = np.vstack((b2_x2, b2_x1, b2_x1, b2_x2, b2_x1))
        b3 = np.vstack((b3_x2, b3_x1, b3_x1, b3_x2, b3_x1))
        b4 = np.vstack((b4_x2, b4_x1, b4_x1, b4_x2, b4_x1))

        heading = np.hstack((heading, sfill, heading, heading, heading))
        pitch = np.hstack((pitch, pitch, pitch, sfill, pitch))
        roll = np.tile(roll, 5)
        orient = np.hstack((orient, orient, orient, orient, sfill))
        lat = np.tile(lat, 5)
        lon = np.tile(lon, 5)
        depth = np.tile(depth, 5)
        ntp = np.tile(ntp, 5)

        # set expected outputs for these 5 records
        # notes:
        #    (1) heading is not used in the calculation of vertical velocity,
        #        therefore the second entry to ww_xpctd is good data out (ww_x0),
        #        not nans as resulted from the missingroll test.
        #    (2) pitch is not used in the calculation of error velocity, so that
        #        in the mixed case (time record 4) the error velocity should be
        #        the same as that for the pure bad/missing beam case (ee_x2, 1st
        #        and 4th entries in ee_xpctd).
        #    (3) the orientation argument affects the roll calculation, so that
        #        when its value is missing (5th time record) the expected result
        #        would be the same as if the roll value were missing. therefore
        #        the 5th column entries are all x1 results.
        uu_xpctd = np.vstack((uu_x2, uu_x1, uu_x0, uu_x1, uu_x1))
        vv_xpctd = np.vstack((vv_x2, vv_x1, vv_x0, vv_x1, vv_x1))
        ww_xpctd = np.vstack((ww_x2, ww_x0, ww_x0, ww_x1, ww_x1))
        ee_xpctd = np.vstack((ee_x2, ee_x1, ee_x0, ee_x2, ee_x1))

        # calculated
        uu_calc = af.adcp_beam_eastward(b1, b2, b3, b4,
                                        heading, pitch, roll, orient,
                                        lat, lon, depth, ntp)
        vv_calc = af.adcp_beam_northward(b1, b2, b3, b4,
                                         heading, pitch, roll, orient,
                                         lat, lon, depth, ntp)
        ww_calc = af.adcp_beam_vertical(b1, b2, b3, b4,
                                        heading, pitch, roll, orient)
        ee_calc = af.adcp_beam_error(b1, b2, b3, b4)

        # test results
        np.testing.assert_array_almost_equal(uu_calc, uu_xpctd, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_xpctd, 4)
        np.testing.assert_array_almost_equal(ww_calc, ww_xpctd, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_xpctd, 4)

    def test_adcp_earth(self):
        """
        Tests magnetic_correction function for ADCPs set to output data in the
        Earth Coordinate system.

        Values were not defined in DPS, were recreated using test values above:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by:

            2014-02-06: Christopher Wingard. Initial code.
            2015-06-10: Russell Desiderio.
                            Changed adcp_ins2earth to require the units of the compass
                            data to be in centidegrees.
        """

        # set the test data
        u, v, w, e = af.adcp_beam2ins(self.b1, self.b2, self.b3, self.b4)

        ### old adcp_ins2earth returned 3 variables (CEW)
        # adcp_ins2earth now requires compass data in units of centidegrees (RAD)
        uu, vv, ww = af.adcp_ins2earth(u, v, w, self.heading, self.pitch,
                                       self.roll, self.orient)

        # test the magnetic variation correction
        got_uu_cor = af.adcp_earth_eastward(uu, vv, self.depth, self.lat, self.lon, self.ntp)
        got_vv_cor = af.adcp_earth_northward(uu, vv, self.depth, self.lat, self.lon, self.ntp)

        np.testing.assert_array_almost_equal(got_uu_cor, self.uu_cor, 4)
        np.testing.assert_array_almost_equal(got_vv_cor, self.vv_cor, 4)

        # reset the test inputs for multiple records using the integer inputs.
        uu = np.tile(uu, (24, 1))
        vv = np.tile(vv, (24, 1))
        depth = np.ones(24) * self.depth
        lat = np.ones(24) * self.lat
        lon = np.ones(24) * self.lon
        ntp = np.ones(24) * self.ntp

        # reset expected results for multiple records
        uu_cor = np.tile(self.uu_cor, (24, 1))
        vv_cor = np.tile(self.vv_cor, (24, 1))

        # compute the results for multiple records
        got_uu_cor = af.adcp_earth_eastward(uu, vv, depth, lat, lon, ntp)
        got_vv_cor = af.adcp_earth_northward(uu, vv, depth, lat, lon, ntp)

        # test the magnetic variation correction
        np.testing.assert_array_almost_equal(got_uu_cor, uu_cor, 4)
        np.testing.assert_array_almost_equal(got_vv_cor, vv_cor, 4)

    def test_adcp_earth_int_input_velocity_data(self):
        """
        Tests adcp_earth_eastward and adcp_earth_northward using int type raw velocity data,
        as will be supplied by CI. Also tests the almost trivial functions adcp_earth_vertical
        and adcp_earth_error (unit change).

        Input raw velocity values were derived from the float unit test in test_adcp_earth
        by rounding the uu and vv float output from adcp_ins2earth. These int inputs failed
        the assert_array_almost_equal unit tests (decimals=4) in test_adcp_earth because of
        round-off error but passed when the agreement precision was relaxed to decimals=3.
        This is taken as justification to more precisely calculate the expected values for
        unit tests in the current module from adcp_earth_eastward and adcp_earth_northward
        themselves (the very modules being tested), using as input the type int raw velocity
        data. Because these DPA functions were used to derive their own check data, the
        original (float type input velocity data) unit tests are retained in the
        test_adcp_earth function.

        The tests in this module will be used to derive unit tests checking the replacement
        of ADCP int bad value sentinels (-32768) with Nans; these tests require that the
        raw velocity data be of type int.

        Implemented by:

            2014-06-16: Russell Desiderio. Initial code.
        """
        # set the input test data [mm/sec]
        uu = np.array([[218, -281, -100, 483, 1238, -245, 622, -181, 99, -906]])
        vv = np.array([[-337, -182, -1052, -868, -892, 258, -850, -87, -307, -546]])
        ww = np.array([[140, 398, 187, 164, 9, -129, 33, -302, 138, 197]])
        ee = np.array([[790, 635, 81, 626, 64, 71, -317, 219, 55, 433]])

        # expected values, calculated using adcp_earth_eastward and adcp_earth_northward
        uu_cor = np.array([[0.11031103, -0.32184604, -0.40227939,  0.20903718,  0.92426103,
                            -0.15916447,  0.34724837, -0.19849871,  0.00522179, -1.02580274]])

        vv_cor = np.array([[-0.38590734, -0.09219615, -0.97717720, -0.97109035, -1.21410442,
                            0.31820696, -0.99438552, -0.03046741, -0.32252555, -0.25822614]])
        # expected values, calculated by changing units from mm/s to m/s
        ww_vel = ww / 1000.0
        ee_vel = ee / 1000.0

        # test the magnetic variation correction using type integer inputs for the velocities.
        got_uu_cor = af.adcp_earth_eastward(uu, vv, self.depth, self.lat, self.lon, self.ntp)
        got_vv_cor = af.adcp_earth_northward(uu, vv, self.depth, self.lat, self.lon, self.ntp)
        # and the unit change functions
        got_ww_vel = af.adcp_earth_vertical(ww)
        got_ee_vel = af.adcp_earth_error(ee)

        np.testing.assert_array_almost_equal(got_uu_cor, uu_cor, 4)
        np.testing.assert_array_almost_equal(got_vv_cor, vv_cor, 4)
        np.testing.assert_array_almost_equal(got_ww_vel, ww_vel, 4)
        np.testing.assert_array_almost_equal(got_ee_vel, ee_vel, 4)

        # reset the test inputs for multiple records using the integer inputs.
        uu = np.tile(uu, (24, 1))
        vv = np.tile(vv, (24, 1))
        ww = np.tile(ww, (24, 1))
        ee = np.tile(ee, (24, 1))
        depth = np.ones(24) * self.depth
        lat = np.ones(24) * self.lat
        lon = np.ones(24) * self.lon
        ntp = np.ones(24) * self.ntp

        # reset expected results for multiple records
        uu_cor = np.tile(uu_cor, (24, 1))
        vv_cor = np.tile(vv_cor, (24, 1))
        ww_vel = np.tile(ww_vel, (24, 1))
        ee_vel = np.tile(ee_vel, (24, 1))

        # compute the results for multiple records
        got_uu_cor = af.adcp_earth_eastward(uu, vv, depth, lat, lon, ntp)
        got_vv_cor = af.adcp_earth_northward(uu, vv, depth, lat, lon, ntp)
        got_ww_vel = af.adcp_earth_vertical(ww)
        got_ee_vel = af.adcp_earth_error(ee)

        # test the magnetic variation correction
        np.testing.assert_array_almost_equal(got_uu_cor, uu_cor, 4)
        np.testing.assert_array_almost_equal(got_vv_cor, vv_cor, 4)
        # and the unit change functions
        np.testing.assert_array_almost_equal(got_ww_vel, ww_vel, 4)
        np.testing.assert_array_almost_equal(got_ee_vel, ee_vel, 4)

    def test_adcp_earth_with_fill(self):
        """
        Tests adcp_earth_eastward, adcp_earth_northward, adcp_earth_vertical and
        adcp_earth_error when system fill values and ADCP fill values (bad value
        sentinels) are present in the data stream.

        Non-fill test values come from the function test_adcp_earth_int_input_velocity_data
        in this module.

        Implemented by:

            2014-06-25: Russell Desiderio. Initial code.
        """
        # for convenience
        sfill = SYSTEM_FILLVALUE
        afill = ADCP_FILLVALUE

        ### scalar time case
        # set the input test data
        lat = np.array([50.0000])
        lon = np.array([-145.0000])
        depth = np.array([0.0])
        ntp = np.array([3545769600.0])    # May 12, 2012
        # input velocities [mm/sec]
        uu_in0 = np.array([[218, sfill, -100, 483, afill, -245]])
        vv_in0 = np.array([[sfill, -182, -1052, -868, -892, afill]])
        ww_in0 = np.array([[sfill, 398, afill, 164, 9, -129]])
        ee_in0 = np.array([[afill, 635, 81, 626, sfill, 71]])
        # expected values [m/sec]
        uu_x0 = np.array([[np.nan, np.nan, -0.40227,  0.20903, np.nan, np.nan]])
        vv_x0 = np.array([[np.nan, np.nan, -0.97717, -0.97109, np.nan, np.nan]])
        ww_x0 = np.array([[np.nan, 0.398, np.nan, 0.164, 0.009, -0.129]])
        ee_x0 = np.array([[np.nan, 0.635, 0.081, 0.626, np.nan, 0.071]])
        # calculated
        uu_calc = af.adcp_earth_eastward(uu_in0, vv_in0, depth, lat, lon, ntp)
        vv_calc = af.adcp_earth_northward(uu_in0, vv_in0, depth, lat, lon, ntp)
        ww_calc = af.adcp_earth_vertical(ww_in0)
        ee_calc = af.adcp_earth_error(ee_in0)
        # test
        np.testing.assert_array_almost_equal(uu_calc, uu_x0, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_x0, 4)
        np.testing.assert_array_almost_equal(ww_calc, ww_x0, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_x0, 4)

        ### multiple time record case
        # set the input test data
        lat = np.tile(lat, 5)
        lon = np.tile(lon, 5)
        depth = np.tile(depth, 5)
        ntp = np.tile(ntp, 5)
        uu_in0 = np.tile(uu_in0, (5, 1))
        vv_in0 = np.tile(vv_in0, (5, 1))
        ww_in0 = np.tile(ww_in0, (5, 1))
        ee_in0 = np.tile(ee_in0, (5, 1))
        # expected
        uu_x0 = np.tile(uu_x0, (5, 1))
        vv_x0 = np.tile(vv_x0, (5, 1))
        ww_x0 = np.tile(ww_x0, (5, 1))
        ee_x0 = np.tile(ee_x0, (5, 1))
        # calculated
        uu_calc = af.adcp_earth_eastward(uu_in0, vv_in0, depth, lat, lon, ntp)
        vv_calc = af.adcp_earth_northward(uu_in0, vv_in0, depth, lat, lon, ntp)
        ww_calc = af.adcp_earth_vertical(ww_in0)
        ee_calc = af.adcp_earth_error(ee_in0)
        # test
        np.testing.assert_array_almost_equal(uu_calc, uu_x0, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_x0, 4)
        np.testing.assert_array_almost_equal(ww_calc, ww_x0, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_x0, 4)

    def test_adcp_backscatter(self):
        """
        Tests echo intensity scaling function (adcp_backscatter) for ADCPs
        in order to convert from echo intensity in counts to dB.

        Values were not defined in DPS, were created using test values above:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by Christopher Wingard, 2014-02-06
                       Russell Desiderio, 2015-06-25. Added tests for fill values.
        """
        # the single record case
        got = af.adcp_backscatter(self.echo, self.sfactor)
        np.testing.assert_array_almost_equal(got, self.dB, 4)

        # the multi-record case -- inputs
        raw = np.tile(self.echo, (24, 1))
        sf = np.ones(24) * self.sfactor

        # the multi-record case -- outputs
        dB = np.tile(self.dB, (24, 1))
        got = af.adcp_backscatter(raw, sf)
        np.testing.assert_array_almost_equal(got, dB, 4)

        ### test fill value replacement with nan
        # for convenience
        sfill = SYSTEM_FILLVALUE
        # the adcp bad sentinel fillvalue (requires 2 bytes) is not used for echo
        # intensity, which is stored in 1 byte.

        # the single time record case
        echo_with_fill, xpctd = np.copy(self.echo), np.copy(self.dB)
        echo_with_fill[0, 3], xpctd[0, 3] = sfill, np.nan
        echo_with_fill[0, 6], xpctd[0, 6] = sfill, np.nan
        echo_with_fill[0, 7], xpctd[0, 7] = sfill, np.nan

        got = af.adcp_backscatter(echo_with_fill, self.sfactor)
        np.testing.assert_array_almost_equal(got, xpctd, 4)

        # the multiple time record case
        echo_with_fill = np.vstack((echo_with_fill, self.echo, echo_with_fill))
        xpctd = np.vstack((xpctd, self.dB, xpctd))
        sfactor = np.tile(self.sfactor, (3, 1))

        got = af.adcp_backscatter(echo_with_fill, sfactor)
        np.testing.assert_array_almost_equal(got, xpctd, 4)

    def test_vadcp_beam(self):
        """
        Indirectly tests vadcp_beam_eastward, vadcp_beam_northward,
        vadcp_beam_vertical_est, and vadcp_beam_vertical_true functions (which
        call adcp_beam2ins and adcp_ins2earth) and vadcp_beam_error (which only
        calls adcp_beam2ins) for the specialized 5-beam ADCP. Application of
        the magnetic correction and conversion from mm/s to m/s is not applied.

        Values based on those defined in DPS:

        OOI (2012). Data Product Specification for Turbulent Velocity Profile
            and Echo Intensity. Document Control Number 1341-00760.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00760_Data_Product_SPEC_VELTURB_OOI.pdf)

        Implemented by:

            2014-07-24: Christopher Wingard. Initial code.
            2015-06-10: Russell Desiderio.
                           adcp_ins2earth now requires the units of the compass
                           data to be in centidegrees.
        """
        # test inputs
        b1 = np.ones((10, 10)).astype(int) * -325
        b2 = np.ones((10, 10)).astype(int) * 188
        b3 = np.ones((10, 10)).astype(int) * 168
        b4 = np.ones((10, 10)).astype(int) * -338
        b5 = np.ones((10, 10)).astype(int) * -70

        # units of centidegrees
        heading = np.array([30, 30, 30, 30, 30,
                            32, 32, 32, 32, 32]) * 100
        pitch = np.array([0, 2, 3, 3, 1, 2, 2, 3, 3, 1]) * 100
        roll = np.array([0, 4, 3, 4, 3, 3, 4, 3, 4, 3]) * 100
        orient = np.ones(10, dtype=np.int)

        # expected outputs
        vle = np.array([279.6195, 282.6881, 281.8311, 282.7147,
                        282.1188, 246.2155, 246.9874, 246.1226,
                        247.0156, 246.4276]).reshape(-1, 1)
        vle = np.reshape(np.tile(vle, 10), (10, 10))

        vln = np.array([-1015.5964, -1018.0226, -1018.2595, -1017.9765,
                        -1017.7612, -1027.3264, -1027.2681, -1027.4749,
                        -1027.2230, -1026.9870]).reshape(-1, 1)
        vln = np.reshape(np.tile(vln, 10), (10, 10))

        vlu = np.array([81.6756, 3.3916, 3.5950, -9.4974,
                        29.4154, 16.5077, 3.3916, 3.5950,
                        -9.4974, 29.4154]).reshape(-1, 1)
        vlu = np.reshape(np.tile(vlu, 10), (10, 10))

        evl = np.array([34.1128, 34.1128, 34.1128, 34.1128,
                        34.1128, 34.1128, 34.1128, 34.1128,
                        34.1128, 34.1128]).reshape(-1, 1)
        evl = np.reshape(np.tile(evl, 10), (10, 10))

        w5 = np.array([70.0000, -8.2485, -8.0487, -21.1287,
                       17.7575, 4.8552, -8.2485, -8.0487,
                       -21.1287, 17.7575]).reshape(-1, 1)
        w5 = np.reshape(np.tile(w5, 10), (10, 10))

        # test the transformations
        u, v, w_est, e = af.adcp_beam2ins(b1, b2, b3, b4)
        uu, vv, ww_est = af.adcp_ins2earth(u, v, w_est, heading, pitch, roll, orient)
        _, _, ww_true = af.adcp_ins2earth(u, v, b5, heading, pitch, roll, orient)

        # compare the results
        np.testing.assert_array_almost_equal(uu, vle, 4)
        np.testing.assert_array_almost_equal(vv, vln, 4)
        np.testing.assert_array_almost_equal(ww_est, vlu, 4)
        np.testing.assert_array_almost_equal(e, evl, 4)
        np.testing.assert_array_almost_equal(ww_true, w5, 4)

        #### KEEP: RAD 2015-06-22:
        """
        ## Given that these unit tests have validated the VADCP DPA functions, use these
        ## vadcp functions to generate values for unit tests with (a) type integer inputs
        ## (b) that use the vadcp functions themselves, instead of their constituent sub-
        ## routines, so that unit tests checking the trapping of CI fill values (-999999999)
        ## and ADCP instrument bad value sentinels (-32768) can be constructed.
        #lat = np.ones(10) * self.lat
        #lon = np.ones(10) * self.lon
        #z = np.ones(10) * self.depth
        #dt = np.ones(10) * self.ntp
        #
        #vle = af.vadcp_beam_eastward(b1, b2, b3, b4, heading, pitch, roll, orient, lat, lon, z, dt)
        #vln = af.vadcp_beam_northward(b1, b2, b3, b4, heading, pitch, roll, orient, lat, lon, z, dt)
        #vlu_4bm = af.vadcp_beam_vertical_est(b1, b2, b3, b4, heading, pitch, roll, orient)
        #vlu_5bm = af.vadcp_beam_vertical_true(b1, b2, b3, b4, b5, heading, pitch, roll, orient)
        #err = af.vadcp_beam_error(b1, b2, b3, b4)
        #
        #print vle.T
        #print vln.T
        #print vlu_4bm.T
        #print vlu_5bm.T
        #print err.T
        """
        #### RAD 2015-06-22

    def test_vadcp_beam_int_input_velocity_data(self):
        """
        Tests vadcp_beam_eastward, vadcp_beam_northward, vadcp_beam_vertical_est,
        vadcp_beam_vertical_true and vadcp_beam_error functions for the specialized 5-beam ADCP
        using int type raw velocity data, as will be supplied by CI.

        Test values come from the function test_vadcp_beam, in this module.

        The tests in this module will be used to derive unit tests checking the replacement
        of ADCP int bad value sentinels (-32768) with Nans; these tests require that the
        raw velocity data be of type int.

        Implemented by:

            2014-06-22: Russell Desiderio. Initial code.
        """
        # inputs
        b1 = np.ones((10, 10), dtype=np.int) * -325
        b2 = np.ones((10, 10), dtype=np.int) * 188
        b3 = np.ones((10, 10), dtype=np.int) * 168
        b4 = np.ones((10, 10), dtype=np.int) * -338
        b5 = np.ones((10, 10), dtype=np.int) * -70

        # units of centidegrees
        heading = np.array([30, 30, 30, 30, 30,
                            32, 32, 32, 32, 32]) * 100
        pitch = np.array([0, 2, 3, 3, 1, 2, 2, 3, 3, 1]) * 100
        roll = np.array([0, 4, 3, 4, 3, 3, 4, 3, 4, 3]) * 100
        orient = np.ones(10, dtype=np.int)

        lat = np.ones(10) * self.lat
        lon = np.ones(10) * self.lon
        z = np.ones(10) * self.depth
        dt = np.ones(10) * self.ntp

        # expected outputs from test_vadcp_beam
        vle_xpctd = np.array([[-0.02853200, -0.02630381, -0.02719268, -0.02626496, -0.02677222,
                               -0.06390457, -0.06314916, -0.06403668, -0.06310908, -0.06360277]])
        vle_xpctd = np.tile(vle_xpctd.T, (1, 10))

        vln_xpctd = np.array([[-1.05300003, -1.05621525, -1.05619207, -1.05617896, -1.05579924,
                               -1.05448459, -1.05465384, -1.05459965, -1.05461893, -1.05422174]])
        vln_xpctd = np.tile(vln_xpctd.T, (1, 10))

        vlu_4bm_xpctd = np.array([[0.08167564, 0.0033916, 0.00359505, -0.0094974, 0.02941538,
                                   0.01650774, 0.0033916, 0.00359505, -0.0094974, 0.02941538]])
        vlu_4bm_xpctd = np.tile(vlu_4bm_xpctd.T, (1, 10))

        vlu_5bm_xpctd = np.array([[0.07000000, -0.00824854, -0.00804866, -0.02112871, 0.01775751,
                                   0.00485518, -0.00824854, -0.00804866, -0.02112871, 0.01775751]])
        vlu_5bm_xpctd = np.tile(vlu_5bm_xpctd.T, (1, 10))

        err_vel_xpctd = np.tile(0.03411279, (10, 10))

        vle_calc = af.vadcp_beam_eastward(
            b1, b2, b3, b4, heading, pitch, roll, orient, lat, lon, z, dt)
        vln_calc = af.vadcp_beam_northward(
            b1, b2, b3, b4, heading, pitch, roll, orient, lat, lon, z, dt)
        vlu_4bm_calc = af.vadcp_beam_vertical_est(b1, b2, b3, b4, heading, pitch, roll, orient)
        vlu_5bm_calc = af.vadcp_beam_vertical_true(b1, b2, b3, b4, b5, heading, pitch, roll, orient)
        err_vel_calc = af.vadcp_beam_error(b1, b2, b3, b4)

        np.testing.assert_array_almost_equal(vle_calc, vle_xpctd, 6)
        np.testing.assert_array_almost_equal(vln_calc, vln_xpctd, 6)
        np.testing.assert_array_almost_equal(vlu_4bm_calc, vlu_4bm_xpctd, 6)
        np.testing.assert_array_almost_equal(vlu_5bm_calc, vlu_5bm_xpctd, 6)
        np.testing.assert_array_almost_equal(err_vel_calc, err_vel_xpctd, 6)

    def test_vadcp_beam_with_fill(self):
        """
        Tests vadcp_beam_eastward, vadcp_beam_northward, vadcp_beam_vertical_est,
        vadcp_beam_vertical_true and vadcp_beam_error functions for the specialized
        5-beam ADCP when system fill values and ADCP fill values (bad value sentinels)
        are present in the data stream.

        Non-fill test values come from the function test_vadcp_beam_int_input_velocity_data
        in this module.

        Implemented by:

            2014-06-25: Russell Desiderio. Initial code.

        Notes:

            Before this time there have been no scalar time tests for the vadcp functions.
            Therefore, scalar time tests are included in this test function.
        """
        # for convenience
        sfill = SYSTEM_FILLVALUE
        afill = ADCP_FILLVALUE

        ### scalar tests with all good data
        # inputs
        b1_x0 = np.tile(-325, (1, 6))
        b2_x0 = np.tile(188, (1, 6))
        b3_x0 = np.tile(168, (1, 6))
        b4_x0 = np.tile(-338, (1, 6))
        b5_x0 = np.tile(-70, (1, 6))

        # units of centidegrees
        heading = np.array([3000])
        pitch = np.array([200])
        roll = np.array([400])
        orient = np.array([1])  # vertical orientation

        lat = np.array([50.0000])
        lon = np.array([-145.0000])
        z = np.array([0.0])
        dt = np.array([3545769600.0])  # May 12, 2012

        # expected outputs from test_vadcp_beam
        uu_x0 = np.tile(-0.02630, (1, 6))
        vv_x0 = np.tile(-1.05621, (1, 6))
        ww_4bm_x0 = np.tile(0.00330, (1, 6))
        ww_5bm_x0 = np.tile(-0.00824, (1, 6))
        ee_x0 = np.tile(0.03411, (1, 6))

        uu_calc = af.vadcp_beam_eastward(
            b1_x0, b2_x0, b3_x0, b4_x0, heading, pitch, roll, orient, lat, lon, z, dt)
        vv_calc = af.vadcp_beam_northward(
            b1_x0, b2_x0, b3_x0, b4_x0, heading, pitch, roll, orient, lat, lon, z, dt)
        ww_4bm_calc = af.vadcp_beam_vertical_est(
            b1_x0, b2_x0, b3_x0, b4_x0, heading, pitch, roll, orient)
        ww_5bm_calc = af.vadcp_beam_vertical_true(
            b1_x0, b2_x0, b3_x0, b4_x0, b5_x0, heading, pitch, roll, orient)
        ee_calc = af.vadcp_beam_error(b1_x0, b2_x0, b3_x0, b4_x0)

        np.testing.assert_array_almost_equal(uu_calc, uu_x0, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_x0, 4)
        np.testing.assert_array_almost_equal(ww_4bm_calc, ww_4bm_x0, 4)
        np.testing.assert_array_almost_equal(ww_5bm_calc, ww_5bm_x0, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_x0, 4)

        ### single time record case; missing roll data
        ## the ADCP does not use its bad flag sentinel for compass data, only beam data.
        ## however, it is possible that CI could supply the system fillvalue for missing compass data.
        # input data
        # beam inputs are the same as for the all good data scalar test above
        b1_x1 = np.array([[-325, -325, -325, -325, -325, -325]])
        b2_x1 = np.array([[188, 188, 188, 188, 188, 188]])
        b3_x1 = np.array([[168, 168, 168, 168, 168, 168]])
        b4_x1 = np.array([[-338, -338, -338, -338, -338, -338]])
        b5_x1 = np.array([[-70, -70, -70, -70, -70, -70]])
        # compass data as above, except the roll value from the instrument is missing:
        missingroll = sfill
        # expected results for all good beam data, missing roll data;
        # nans for all results except for the error velocity, which does not depend on the compass
        uu_x1 = uu_x0 * np.nan
        vv_x1 = vv_x0 * np.nan
        ww_4bm_x1 = ww_4bm_x0 * np.nan
        ww_5bm_x1 = ww_5bm_x0 * np.nan
        ee_x1 = np.copy(ee_x0)
        # calculated
        uu_calc = af.vadcp_beam_eastward(b1_x1, b2_x1, b3_x1, b4_x1, heading, pitch,
                                         missingroll,
                                         orient, lat, lon, z, dt)
        vv_calc = af.vadcp_beam_northward(b1_x1, b2_x1, b3_x1, b4_x1, heading, pitch,
                                          missingroll,
                                          orient, lat, lon, z, dt)
        ww_4bm_calc = af.vadcp_beam_vertical_est(b1_x1, b2_x1, b3_x1, b4_x1, heading, pitch,
                                                 missingroll,
                                                 orient)
        ww_5bm_calc = af.vadcp_beam_vertical_true(b1_x1, b2_x1, b3_x1, b4_x1, b5_x1, heading, pitch,
                                                  missingroll,
                                                  orient)
        ee_calc = af.vadcp_beam_error(b1_x1, b2_x1, b3_x1, b4_x1)
        # test results
        np.testing.assert_array_almost_equal(uu_calc, uu_x1, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_x1, 4)
        np.testing.assert_array_almost_equal(ww_4bm_calc, ww_4bm_x1, 4)
        np.testing.assert_array_almost_equal(ww_5bm_calc, ww_5bm_x1, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_x1, 4)

        ### single time record case; missing and bad-flagged beam data, good compass data
        # input data
        b1_x2 = np.array([[-325, -325, -325, sfill, -325, -325]])
        b2_x2 = np.array([[188, afill, 188, 188, 188, 188]])
        b3_x2 = np.array([[168, 168, 168, 168, 168, 168]])
        b4_x2 = np.array([[-338, -338, -338, -338, -338, -338]])
        b5_x2 = np.array([[sfill, sfill, -70, -70, afill, -70]])
        # expected
        uu_x2 = np.array([[-0.02630, np.nan, -0.02630, np.nan, -0.02630, -0.02630]])
        vv_x2 = np.array([[-1.05621, np.nan, -1.05621, np.nan, -1.05621, -1.05621]])
        ww_4bm_x2 = np.array([[0.00330, np.nan, 0.00330, np.nan, 0.00330, 0.00330]])
        ww_5bm_x2 = np.array([[np.nan, np.nan, -0.00824, np.nan, np.nan, -0.00824]])
        ee_x2 = np.array([[0.03411, np.nan, 0.03411, np.nan, 0.03411, 0.03411]])
        # calculated
        uu_calc = af.vadcp_beam_eastward(b1_x2, b2_x2, b3_x2, b4_x2,
                                         heading, pitch, roll, orient,
                                         lat, lon, z, dt)
        vv_calc = af.vadcp_beam_northward(b1_x2, b2_x2, b3_x2, b4_x2,
                                          heading, pitch, roll, orient,
                                          lat, lon, z, dt)
        ww_4bm_calc = af.vadcp_beam_vertical_est(b1_x2, b2_x2, b3_x2, b4_x2,
                                                 heading, pitch, roll, orient)
        ww_5bm_calc = af.vadcp_beam_vertical_true(b1_x2, b2_x2, b3_x2, b4_x2, b5_x2,
                                                  heading, pitch, roll, orient)
        ee_calc = af.vadcp_beam_error(b1_x2, b2_x2, b3_x2, b4_x2)
        # test results
        np.testing.assert_array_almost_equal(uu_calc, uu_x2, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_x2, 4)
        np.testing.assert_array_almost_equal(ww_4bm_calc, ww_4bm_x2, 4)
        np.testing.assert_array_almost_equal(ww_5bm_calc, ww_5bm_x2, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_x2, 4)

        ### multiple (5) record case
        ## reset the test inputs for 5 time records
        #     1st time record is the bad/missing beam data case above
        #     2nd time record is a missing heading data case
        #     3rd time record is all good data (note, b1_x1 = b1_x0, etc)
        #     4th time record is bad/missing beam and missing pitch data.
        #     5th time record is missing orientation data
        b1 = np.vstack((b1_x2, b1_x1, b1_x1, b1_x2, b1_x1))
        b2 = np.vstack((b2_x2, b2_x1, b2_x1, b2_x2, b2_x1))
        b3 = np.vstack((b3_x2, b3_x1, b3_x1, b3_x2, b3_x1))
        b4 = np.vstack((b4_x2, b4_x1, b4_x1, b4_x2, b4_x1))
        b5 = np.vstack((b5_x2, b5_x1, b5_x1, b5_x2, b5_x1))

        heading = np.hstack((heading, sfill, heading, heading, heading))
        pitch = np.hstack((pitch, pitch, pitch, sfill, pitch))
        roll = np.tile(roll, 5)
        orient = np.hstack((orient, orient, orient, orient, sfill))
        lat = np.tile(lat, 5)
        lon = np.tile(lon, 5)
        z = np.tile(z, 5)
        dt = np.tile(dt, 5)

        # set expected outputs for these 5 records
        # notes:
        #    (1) heading is not used in the calculation of vertical velocities,
        #        therefore the second entries to the ww_xpctd products are good
        #        data out (ww_x0), not nans as resulted from the missingroll test.
        #    (2) pitch is not used in the calculation of error velocity, so that
        #        in the mixed case (time record 4) the error velocity should be
        #        the same as that for the pure bad/missing beam case (ee_x2, 1st
        #        and 4th entries in ee_xpctd).
        #    (3) the orientation argument affects the roll calculation, so that
        #        when its value is missing (5th time record) the expected result
        #        would be the same as if the roll value were missing. therefore
        #        the 5th column entries are all x1 results.
        uu_xpctd = np.vstack((uu_x2, uu_x1, uu_x0, uu_x1, uu_x1))
        vv_xpctd = np.vstack((vv_x2, vv_x1, vv_x0, vv_x1, vv_x1))
        ww_4bm_xpctd = np.vstack((ww_4bm_x2, ww_4bm_x0, ww_4bm_x0, ww_4bm_x1, ww_4bm_x1))
        ww_5bm_xpctd = np.vstack((ww_5bm_x2, ww_5bm_x0, ww_5bm_x0, ww_5bm_x1, ww_5bm_x1))
        ee_xpctd = np.vstack((ee_x2, ee_x1, ee_x0, ee_x2, ee_x1))

        # calculated
        uu_calc = af.vadcp_beam_eastward(b1, b2, b3, b4,
                                         heading, pitch, roll, orient,
                                         lat, lon, z, dt)
        vv_calc = af.vadcp_beam_northward(b1, b2, b3, b4,
                                          heading, pitch, roll, orient,
                                          lat, lon, z, dt)
        ww_4bm_calc = af.vadcp_beam_vertical_est(b1, b2, b3, b4,
                                                 heading, pitch, roll, orient)
        ww_5bm_calc = af.vadcp_beam_vertical_true(b1, b2, b3, b4, b5,
                                                  heading, pitch, roll, orient)
        ee_calc = af.vadcp_beam_error(b1, b2, b3, b4)
        # test results
        np.testing.assert_array_almost_equal(uu_calc, uu_xpctd, 4)
        np.testing.assert_array_almost_equal(vv_calc, vv_xpctd, 4)
        np.testing.assert_array_almost_equal(ww_4bm_calc, ww_4bm_xpctd, 4)
        np.testing.assert_array_almost_equal(ww_5bm_calc, ww_5bm_xpctd, 4)
        np.testing.assert_array_almost_equal(ee_calc, ee_xpctd, 4)

    def test_adcp_ins2earth_orientation(self):
        """
        Test the adcp worker routine adcp_inst2earth when the vertical orientation
        toggle switch codes for downward-looking (orient = 0).

        The instrument to earth coordinate transformation was coded in matlab using
        the January 2010 version of the Teledyne RD Instruments "ADCP Coordinate
        Transformation" manual. The code was checked against the upwards looking
        unit tests using adcp_inst2earth in the test_vadcp_beam function in this module;
        the values agreed to beyond single precision. This matlab code was then used
        to generate the check values for the downward looking case.

        Implemented by:

            2014-07-02: Russell Desiderio. Initial code.
        """
        # input values: these are the output of adcp_beam2inst in test_vadcp_beam
        # (velocities are in instrument coordinates)
        u = np.array([[-749.95582864]])
        v = np.array([[-739.72251324]])
        w = np.array([[-81.67564404]])
        # units of centidegrees
        heading = np.array([3200])
        pitch = np.array([300])
        roll = np.array([400])

        ### check test: upwards looking case
        orient_1 = np.array([1])
        # expected outputs, earth coordinates, upwards case, from test_vadcp_beam which
        # agrees with the test matlab code values to (much) better than single precision.
        vle = np.array([[247.015599]])
        vln = np.array([[-1027.223026]])
        vlu = np.array([[-9.497397]])
        xpctd = np.hstack((vle, vln, vlu))
        # calculated upwards looking case
        uu, vv, ww = af.adcp_ins2earth(u, v, w, heading, pitch, roll, orient_1)
        calc = np.hstack((uu, vv, ww))
        # test results
        np.testing.assert_array_almost_equal(calc, xpctd, 6)

        ### primary test: downwards looking case.
        # change one input:
        orient_0 = np.array([0])
        # expected outputs, earth coordinates, downwards case, from matlab test code
        vle = np.array([[-1029.9328104]])
        vln = np.array([[-225.7064203]])
        vlu = np.array([[-67.7426771]])
        xpctd = np.hstack((vle, vln, vlu))
        # calculated upwards looking case
        uu, vv, ww = af.adcp_ins2earth(u, v, w, heading, pitch, roll, orient_0)
        calc = np.hstack((uu, vv, ww))
        # test results
        np.testing.assert_array_almost_equal(calc, xpctd, 6)

    def test_adcp_bin_depths_meters(self):
        """
        Test the adcp_bin_depths_meters function.

        Implemented by:
            Craig Risien, January 2015. Initial code.
            Russell Desiderio. 26-Jun-2015. Added time-vectorized unit test after modifying DPA.
                               30-Jun-2015. Added fill value unit test.
        """
        sfill = SYSTEM_FILLVALUE

        ### scalar time case (1) - adcp looking up
        # test inputs - note, CI will be sending these into the DPAs as ndarrays, not python scalars.
        adcp_orientation = 1
        bin_size = 400
        dist_first_bin = 900
        num_bins = 29
        sensor_depth = 450
        # expected outputs
        # note that the output should be a row vector, not a 1D array.
        xpctd_bins_up = np.array([[441., 437., 433., 429., 425., 421., 417., 413., 409., 405., 401., 397., 393., 389.,
                                  385., 381., 377., 373., 369., 365., 361., 357., 353., 349., 345., 341., 337., 333.,
                                  329.]])
         # calculate bin depths
        calc_bins_up = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)

        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins_up, xpctd_bins_up, rtol=0.000001, atol=0.000001)

        ### scalar time case (2) - adcp looking down
        # test inputs
        adcp_orientation = np.array([0])
        bin_size = np.array([400])
        dist_first_bin = np.array([900])
        num_bins = np.array([29])
        sensor_depth = np.array([7])
        # expected outputs
        xpctd_bins_down = np.array([[16., 20., 24., 28., 32., 36., 40., 44., 48., 52., 56., 60., 64., 68., 72., 76., 80.,
                                    84., 88., 92., 96., 100., 104., 108., 112., 116., 120., 124., 128.]])
        # calculate bin depths
        calc_bins_down = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)

        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins_down, xpctd_bins_down, rtol=0.000001, atol=0.000001)

        ### time-vectorized case; cat the above two scalar cases together.
        # inputs
        dist_first_bin = np.array([900, 900])
        bin_size = np.array([400, 400])
        num_bins = np.array([29, 29])
        sensor_depth = np.array([450, 7])
        adcp_orientation = np.array([1, 0])
        # expected
        xpctd_bins = np.vstack((xpctd_bins_up, xpctd_bins_down))
        # calculated
        calc_bins = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)
        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins, xpctd_bins, rtol=0.000001, atol=0.000001)

        ### time-vectorized fill cases - test the action on a fill value in each of the 5 input data streams,
        # plus one instance of all good data.
        num_bins = np.array([29, 29, 29, sfill, 29, 29])  # NOTE: DPA uses only first num_bins value
        dist_first_bin = np.array([900, sfill, 900, 900, 900, 900])
        bin_size = np.array([400, 400, sfill, 400, 400, 400])
        sensor_depth = np.array([450, 7, 450, 7, 450, sfill])
        adcp_orientation = np.array([1, 0, 1, 0, sfill, 0])
        # 1st and 4th rows will have non-Nan data.
        xpctd_bins = np.tile(np.nan, (6, 29))
        xpctd_bins[0, :] = xpctd_bins_up
        xpctd_bins[3, :] = xpctd_bins_down
        # calculated
        calc_bins = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)
        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins, xpctd_bins, rtol=0.000001, atol=0.000001)

    def test_adcp_bin_depths_dapa(self):
        """
        Test the adcp_bin_depths_dapa function.

        Values based on z_from_p check values.

        Implemented by:
            Craig Risien, January 2015. Initial code.
            Russell Desiderio. 26-Jun-2015. Corrected pressure type and units and z_from_p usage.
                               30-Jun-2015. Added fill value unit test.
        """
        sfill = SYSTEM_FILLVALUE

        ### scalar time case (1) - adcp looking up
        # test inputs - note, CI will be sending these into the DPAs as ndarrays, not python scalars.
        # test inputs
        adcp_orientation = np.array([1])
        bin_size = np.array([400])
        dist_first_bin = np.array([900])
        latitude = np.array([4.0])
        num_bins = np.array([10])
        # input adcp pressure has units of decaPascals
        pressure = np.array([600000])
        # according to the z_from_p check value at 600db, this gives a depth of 595.8262 m
        # expected outputs
        # note that the output should be a row vector, not a 1D array.
        expected_bins_up = 0.8253480 + np.array([[586, 582, 578, 574, 570, 566, 562, 558, 554, 550]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_up, rtol=0.0, atol=0.000001)

        ### scalar time case (2) - adcp looking down
        # test inputs - should also work with python scalars, but this is not necessary
        adcp_orientation = 0
        bin_size = 400
        dist_first_bin = 900
        latitude = 4.0
        num_bins = 10
        pressure = 10000
        # expected depth from a pressure of 10000 decapascals is 9.94460074 m
        # expected outputs
        expected_bins_down = 0.9445834 + np.array([[18, 22, 26, 30, 34, 38, 42, 46, 50, 54]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_down, rtol=0.0, atol=0.000001)

        ### time-vectorized test - combine the 2 above
        adcp_orientation = np.array([1, 0])
        bin_size = np.array([400, 400])
        dist_first_bin = np.array([900, 900])
        latitude = np.array([4.0, 4.0])
        num_bins = np.array([10, 10])
        pressure = np.array([600000, 10000])
        #
        expected_bins = np.vstack((expected_bins_up, expected_bins_down))
        #
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        #
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

        ### time-vectorized fill cases - test the action on a fill value in each of the 5 input data streams,
        # plus one instance of all good data.
        num_bins = np.array([10, 10, 10, sfill, 10, 10])  # NOTE: DPA uses only first num_bins value
        dist_first_bin = np.array([900, sfill, 900, 900, 900, 900])
        bin_size = np.array([400, 400, sfill, 400, 400, 400])
        pressure = np.array([600000, 10000, 600000, 10000, 600000, sfill])
        adcp_orientation = np.array([1, 0, 1, 0, sfill, 0])
        latitude = np.array([4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
        # 1st and 4th rows will have non-Nan data.
        expected_bins = np.tile(np.nan, (6, 10))
        expected_bins[0, :] = expected_bins_up
        expected_bins[3, :] = expected_bins_down
        # calculated
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

    def test_adcp_bin_depths_bar(self):
        """
        Test the adcp_bin_depths_bar function.

        Values based on z_from_p check values.

        Implemented by:
            Craig Risien, January 2015. Initial code.
            Russell Desiderio. 26-Jun-2015. Corrected pressure type and units and z_from_p usage.
                               30-Jun-2015. Added fill value unit test.
        """
        sfill = SYSTEM_FILLVALUE

        ### scalar time case (1) - adcp looking up
        # test inputs - note, CI will be sending these into the DPAs as ndarrays, not python scalars.
        # test inputs
        adcp_orientation = np.array([1])
        bin_size = np.array([400])
        dist_first_bin = np.array([900])
        latitude = np.array([4.0])
        num_bins = np.array([10])
        # water pressure of gliders has units of bar
        pressure = np.array([60])
        # according to the z_from_p check value at 600db, this gives a depth of 595.8262 m
        # expected outputs
        # note that the output should be a row vector, not a 1D array.
        expected_bins_up = 0.8253480 + np.array([[586, 582, 578, 574, 570, 566, 562, 558, 554, 550]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_up, rtol=0.0, atol=0.000001)

        ### scalar time case (2) - adcp looking down
        # test inputs - should also work with python scalars, but this is not necessary
        adcp_orientation = 0
        bin_size = 400
        dist_first_bin = 900
        latitude = 4.0
        num_bins = 10
        pressure = 1
        # expected depth from a pressure of 1 bar is 9.94460074 m
        # expected outputs
        expected_bins_down = 0.9445834 + np.array([[18, 22, 26, 30, 34, 38, 42, 46, 50, 54]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_down, rtol=0.0, atol=0.000001)

        ### time-vectorized test - combine the 2 above
        adcp_orientation = np.array([1, 0])
        bin_size = np.array([400, 400])
        dist_first_bin = np.array([900, 900])
        latitude = np.array([4.0, 4.0])
        num_bins = np.array([10, 10])
        pressure = np.array([60, 1])
        #
        expected_bins = np.vstack((expected_bins_up, expected_bins_down))
        #
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        #
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

        ### time-vectorized fill cases - test the action on a fill value in each of the 5 input data streams,
        # plus one instance of all good data.
        num_bins = np.array([10, 10, 10, sfill, 10, 10])  # NOTE: DPA uses only first num_bins value
        dist_first_bin = np.array([900, sfill, 900, 900, 900, 900])
        bin_size = np.array([400, 400, sfill, 400, 400, 400])
        pressure = np.array([60, 1, 60, 1, 60, sfill])
        adcp_orientation = np.array([1, 0, 1, 0, sfill, 0])
        latitude = np.array([4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
        # 1st and 4th rows will have non-Nan data.
        expected_bins = np.tile(np.nan, (6, 10))
        expected_bins[0, :] = expected_bins_up
        expected_bins[3, :] = expected_bins_down
        # calculated
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

    def test_z_from_p(self):
        """
        Test the z_from_p function, which calculates depth in meters from pressure in decibars
        assuming a TEOS-10 standard ocean salinity and a TEOS-10 conservative temperature of 0 deg_C.

        Check values are taken from the TEOS-10 version 3.05 matlab documentation (see the function
        z_from_p in the adcp_functions.py module).

        There are no time-vectorized unit tests in this test function; the DPA that calls this
        function does have time-vectorized unit tests, however.

        Implemented by:

            2014-06-26: Russell Desiderio. Initial code.
            2015-07-01: Russell Desiderio. Updated check values to TEOS-10 ver. 3.05.
        """
        # test inputs
        p = np.array([10.0, 50.0, 125.0, 250.0, 600.0, 1000.0])
        lat = np.ones(6) * 4.0
        # outputs
        xpctd = np.array([-9.9445834469453,  -49.7180897012550, -124.2726219409978,
                          -248.4700576548589, -595.8253480356214, -992.0919060719987])
        calc = af.z_from_p(p, lat)
        # test relative accuracy
        np.testing.assert_allclose(calc, xpctd, rtol=0.00000001, atol=0.0)
        # test absolute accuracy
        np.testing.assert_allclose(calc, xpctd, rtol=0.0, atol=0.00000001)
