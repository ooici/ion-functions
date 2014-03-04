"""
@package ion_functions.test.adcp_functions
@file ion_functions/test/test_adcp_functions.py
@author Christopher Wingard
@brief Unit tests for adcp_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np

from ion_functions.data import adcp_functions as af


@attr('UNIT', group='func')
class TestADCPFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        # set test inputs -- values from DPS
        self.b1 = np.array([-0.0300, -0.2950, -0.5140, -0.2340, -0.1880,
                            0.2030, -0.3250,  0.3050, -0.2040, -0.2940]) * 1000
        self.b2 = np.array([0.1800, -0.1320,  0.2130,  0.3090,  0.2910,
                            0.0490,  0.1880,  0.3730, -0.0020, 0.1720]) * 1000
        self.b3 = np.array([-0.3980, -0.4360, -0.1310, -0.4730, -0.4430,
                            0.1880, -0.1680,  0.2910, -0.1790, 0.0080]) * 1000
        self.b4 = np.array([-0.2160, -0.6050, -0.0920, -0.0580,  0.4840,
                            -0.0050,  0.3380,  0.1750, -0.0800, -0.5490]) * 1000
        self.echo = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250])
        self.sfactor = 0.45
        self.heading = 98.4100
        self.pitch = 0.6900
        self.roll = -2.5400
        self.orient = 1
        self.lat = 50.0000
        self.lon = -145.0000
        self.depth = 0.0
        self.ntp = 3545769600.0    # May 12, 2012

        # set expected results -- velocity profiles in earth coordinates
        # (values in DPS)
        self.uu = np.array([0.2175, -0.2814, -0.1002, 0.4831, 1.2380,
                            -0.2455, 0.6218, -0.1807, 0.0992, -0.9063])
        self.vv = np.array([-0.3367, -0.1815, -1.0522, -0.8676, -0.8919,
                            0.2585, -0.8497, -0.0873, -0.3073, -0.5461])
        self.ww = np.array([0.1401,  0.3977,  0.1870,  0.1637,  0.0091,
                            -0.1290,  0.0334, -0.3017, 0.1384, 0.1966])

        # set expected results -- magnetic variation correction applied
        # (computed in Matlab using above values and mag_var.m)
        self.uu_cor = np.array([0.1099, -0.3221, -0.4025, 0.2092, 0.9243,
                                -0.1595, 0.3471, -0.1983, 0.0053, -1.0261])
        self.vv_cor = np.array([-0.3855, -0.0916, -0.9773, -0.9707, -1.2140,
                                0.3188, -0.9940, -0.0308, -0.3229, -0.2582])

        # set the expecte results -- echo intensity conversion from counts to dB
        self.dB = np.array([0.00, 11.25, 22.50, 33.75, 45.00, 56.25, 67.50,
                            78.75, 90.00, 101.25, 112.50])

    def test_adcp_beam(self):
        """
        Tests adcp_beam2ins, adcp_ins2earth and magentic_correction functions
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
        """
        ins2earth = np.vectorize(af.adcp_ins2earth)

        # single record case
        u, v, w, e = af.adcp_beam2ins(self.b1, self.b2, self.b3, self.b4)
        got_uu, got_vv, got_ww = ins2earth(u, v, w, self.heading,
                                           self.pitch, self.roll, self.orient)

        # test the beam to earth coordinate transforms
        np.testing.assert_array_almost_equal(got_uu / 1000, self.uu, 4)
        np.testing.assert_array_almost_equal(got_vv / 1000, self.vv, 4)
        np.testing.assert_array_almost_equal(got_ww / 1000, self.ww, 4)

        # reset the test inputs for multiple records
        b1 = np.tile(self.b1, (24, 1))
        b2 = np.tile(self.b2, (24, 1))
        b3 = np.tile(self.b3, (24, 1))
        b4 = np.tile(self.b4, (24, 1))
        heading = np.ones(24) * self.heading
        pitch = np.ones(24) * self.pitch
        roll = np.ones(24) * self.roll
        orient = np.ones(24) * self.orient

        # reset expected results for multiple records
        uu = np.tile(self.uu, (24, 1))
        vv = np.tile(self.vv, (24, 1))
        ww = np.tile(self.ww, (24, 1))

        u, v, w, e = af.adcp_beam2ins(b1, b2, b3, b4)
        h = heading.reshape(heading.shape[0], 1)
        p = pitch.reshape(pitch.shape[0], 1)
        r = roll.reshape(roll.shape[0], 1)
        vf = orient.reshape(orient.shape[0], 1)

        got_uu, got_vv, got_ww = ins2earth(u, v, w, h, p, r, vf)

        # test the beam to earth coordinate transforms for multiple records
        np.testing.assert_array_almost_equal(got_uu / 1000, uu, 4)
        np.testing.assert_array_almost_equal(got_vv / 1000, vv, 4)
        np.testing.assert_array_almost_equal(got_ww / 1000, ww, 4)

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

        Implemented by Christopher Wingard, 2014-02-06
        """
        ins2earth = np.vectorize(af.adcp_ins2earth)

        # single record case
        u, v, w, e = af.adcp_beam2ins(self.b1, self.b2, self.b3, self.b4)
        uu, vv, ww = ins2earth(u, v, w, self.heading,
                               self.pitch, self.roll, self.orient)
        got_uu_cor = af.adcp_earth_eastward(uu, vv, self.depth, self.lat, self.lon, self.ntp)
        got_vv_cor = af.adcp_earth_northward(uu, vv, self.depth, self.lat, self.lon, self.ntp)

        # test the magnetic variation correction
        np.testing.assert_array_almost_equal(got_uu_cor, self.uu_cor.reshape(1, 10), 4)
        np.testing.assert_array_almost_equal(got_vv_cor, self.vv_cor.reshape(1, 10), 4)

        # reset the test inputs for multiple records
        b1 = np.tile(self.b1, (24, 1))
        b2 = np.tile(self.b2, (24, 1))
        b3 = np.tile(self.b3, (24, 1))
        b4 = np.tile(self.b4, (24, 1))
        heading = np.ones(24) * self.heading
        pitch = np.ones(24) * self.pitch
        roll = np.ones(24) * self.roll
        orient = np.ones(24) * self.orient
        depth = np.ones(24) * self.depth
        lat = np.ones(24) * self.lat
        lon = np.ones(24) * self.lon
        ntp = np.ones(24) * self.ntp

        # reset expected results for multiple records
        uu_cor = np.tile(self.uu_cor, (24, 1))
        vv_cor = np.tile(self.vv_cor, (24, 1))

        # compute the results for multiple records
        u, v, w, e = af.adcp_beam2ins(b1, b2, b3, b4)
        h = heading.reshape(heading.shape[0], 1)
        p = pitch.reshape(pitch.shape[0], 1)
        r = roll.reshape(roll.shape[0], 1)
        vf = orient.reshape(orient.shape[0], 1)
        uu, vv, ww = ins2earth(u, v, w, h, p, r, vf)

        got_uu_cor = af.adcp_earth_eastward(uu, vv, depth, lat, lon, ntp)
        got_vv_cor = af.adcp_earth_northward(uu, vv, depth, lat, lon, ntp)

        # test the magnetic variation correction
        np.testing.assert_array_almost_equal(got_uu_cor, uu_cor, 4)
        np.testing.assert_array_almost_equal(got_vv_cor, vv_cor, 4)

    def test_adcp_echo(self):
        """
        Tests echo intensity scaling function for ADCPs in order to convert
        from echo intensity in counts to dB.

        Values were not defined in DPS, were created using test values above:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by Christopher Wingard, 2014-02-06
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
