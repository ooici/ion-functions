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
from ion_functions.data.do2_functions import do2_dofst_frequency
from ion_functions.data.do2_functions import do2_dofst_volt
#import pdb


@attr('UNIT', group='func')
class TestDo2FunctionsUnit(BaseUnitTestCase):

    def test_dofst_frequency(self):
        """ DOFST frequency test
        Unit Test of the do2_dofst_frequency function in
        ion_functions/data/do2_functions.py for calculation of
        oxygen from the frequency output of an SBE 43 instrument (DOFST).

        NOTE:
        SBE43F oxygen sensors are connected to a SBE 52-MP profiling CTD
        which provides the inputs of salinity, temperature, and
        pressure.
        """

        # FUNCTION CALIBRATION COEFFICIENTS
        #   These will come from a lookup table
        A = -4.1168e-3
        B = 2.4818e-4
        C = -3.8820e-6
        E = 0.036
        Foffset = -839.55
        Soc = 2.9968e-4

        # FUNCTION INPUTS
        lat = 45.0
        lon = -125.0

        salt = np.array([  # practical salinity in psu
            34.1145, 34.2845, 33.2464, 33.5524, 33.5619,
            33.2512, 33.2609, 33.2716, 33.4191, 33.2710,
            33.2808, 33.5483, 33.5424, 33.3458, 0.0000,
            37.7843, 35.7594, 33.3313, 33.3132, 33.3132,
            33.3132, 33.3132, 33.3132, 33.3132, 33.3132])

        temp = np.array([  # temperature in deg C
            15.5257, 15.3317, 11.9239, 12.8940, 12.9011,
            11.9350, 11.9715, 12.0110, 12.4553, 11.9932,
            12.0196, 12.8647, 12.8448, 12.2084, 12.0996,
            -10.1230, 0.0000, 12.0996, 12.0996, 12.0996,
            12.0996, 12.0996, 12.0996, 12.0996, 12.0996])

        pres = np.array([  # pressure in dbars
            60.5200, 72.5800, 31.4200, 70.8200, 74.8700,
            29.3300, 30.9500, 43.5800, 65.3700, 29.4600,
            31.0300, 74.5700, 75.0700, 57.9200, 42.9800,
            42.9800, 42.9800, 0.00000, 42.9800, 42.9800,
            42.9800, 42.9800, 42.9800, 42.9800, 42.9800])

        freq = np.array([  # 43F instrument output, frequency in Hz
            4354, 4143, 4583, 4476, 4481,
            4591, 4575, 4574, 4545, 4578,
            4572, 4505, 4383, 4555, 4569,
            4023, 4569, 4569, 0, 841,
            1000, 2000, 4000, 5000, 6000])

        # EXPECTED OUTPUTS
        #   intermediate calculation dissolved oxygen in ml/l
        #do_int_expected = np.array([
        #    5.89891032167396, 5.56780727769487, 6.76187243958794,
        #    6.458117534861, 6.46897458201929, 6.77275996877815,
        #    6.73969994032525, 6.74263221709132, 6.64102027182035,
        #    6.74036293148305, 6.7267587280842, 6.51674462650798,
        #    6.30302843255881, 6.6898131217667, 8.2830386610128,
        #    10.7809859878398, 8.95549253591715, 6.68181215593754,
        #    -1.51252046989329, 0.00261229787546345, 0.289064271805584,
        #    2.09064901350446, 5.6938184969022, 7.49540323860107,
        #    9.29698798029994])

        #   final output dissolved oxygen in micro-moles/kg
        do_expected = np.array([
            256.97434863158, 242.509215041926, 294.548757813511,
            281.302611659343, 281.773833754618, 295.022573240991,
            293.582249046689, 293.709591753566, 289.274637555853,
            293.610068502014, 293.01669510765, 283.855524405738,
            274.546683447082, 291.40284695473, 370.109515168768,
            467.353919671818, 388.785846175276, 291.052144126035,
            -65.8845095713829, 0.113790171971304, 12.5914710984794,
            91.0674517683417, 248.019413108066, 326.495393777929,
            404.971374447791])

        # CALCULATION
        do = do2_dofst_frequency(
            freq, Foffset, Soc, A, B, C, E, pres, temp, salt, lat, lon)

        # ASSERT INTERMEDIATE AND FINAL OUTPUTS
        #np.testing.assert_allclose(do_int, do_int_expected, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(do, do_expected, rtol=1e-6, atol=1e-6)

    def test_dofst_voltage(self):
        """ DOFST voltage test
        Unit Test of the do2_dofst_volt function in
        ion_functions/data/do2_functions.py for calculation of
        oxygen from the voltage output of an SBE 43 instrument (DOFST).

        NOTE:
        SBE43 oxygen sensors are connected to a SBE 16+ V2 CTD which
        provides the inputs of salinity, temperature, and pressure.
        """

        # FUNCTION CALIBRATION COEFFICIENTS
        #   These will come from a lookup table
        A = -3.1867e-3
        B = 1.7749e-4
        C = -3.5718e-6
        E = 0.036
        Voffset = -0.5186
        Soc = 0.4396

        # FUNCTION INPUTS
        volt_counts = np.array([  # SBE43 output, voltage counts
            0, 6798, 16384, 32768, 65535,
            0, 6798, 16384, 32768, 65535,
            0, 6798, 16384, 32768, 65535,
            0, 6798, 16384, 32768, 65535,
            0, 6798, 16384, 32768, 65535])

        salt = np.array([  # salinity in psu
            0.0, 33.4, 31.2, 20.1, 35.2, 35.2, 0.0, 31.2, 20.1, 33.4,
            35.2, 0.0, 20.1, 33.4, 31.2, 33.4, 35.2, 31.2, 20.1, 0.0,
            33.4, 31.2, 35.2, 20.1, 0.0])

        temp = np.array([  # temperature in deg C
            0.0, -30.1, 30.3, 10.1, 20.2, 20.2, 30.3, -30.1, 0.0, 10.1,
            20.2, 30.3, 10.1, 0.0, -30.1, -30.1, 20.2, 0.0, 10.1, 30.3,
            30.3, 0.0, -30.1, 10.1, 20.2])

        pres = np.array([  # pressure in dbars
            0.0, 307.5, 201.2, 5.2, 112.1, 5.2, 307.5, 201.2, 0.0,
            112.1, 5.2, 0.0, 112.1, 307.5, 201.2, 112.1, 5.2, 0.0,
            201.2, 307.5, 5.2, 0.0, 112.1, 201.2, 307.5])

        lat = np.array([  # latitude in decimal degrees N
            50.0, -42.0, 39.0, 60.0, 45.0, 60.0, 39.0, 45.0, 50.0,
            -42.0, 50.0, 39.0, 60.0, 45.0, -42.0, 60.0, 45.0, -42.0,
            50.0, 39.0, 60.0, 45.0, 39.0, -42.0, 50.0])

        lon = np.array([  # longitude in decimal degrees E
            145.0, -42.0, -70.5, 39.0, -125.0, 39.0, -70.5, -125.0,
            145.0, -42.0, 145.0, -70.5, 39.0, -125.0, -42.0, 39.0,
            -125.0, -42.0, 145.0, -70.5, 39.0, -125.0, -70.5,
            -42.0, 145.0])

        # EXPECTED OUTPUTS
        #   intermediate calculation dissolved oxygen in ml/l
        #do_int_expected = np.array([
        #    -2.332525266, 0.000797115, 1.412078813, 5.934280027,
        #    10.06589881, -1.149671963, 0.000125639, 10.82518961,
        #    7.744491469, 12.49523919, -1.149671963, 0.000121139,
        #    2.2205184, 7.347649726, 66.32586768, -7.415682793,
        #    0.000120053, 2.645019264, 6.083964669, 10.39697952,
        #    -0.966429789, 0.000195837, 10.2787545, 6.083964669,
        #    12.68706213])

        #   final output dissolved oxygen in micro-moles/kg
        do_expected = np.array([
            -104.1869283, 0.03494869, 61.89990653, 261.0228351,
            438.6325206, -50.09861089, 0.005635974, 475.5984302,
            340.3897211, 544.0600381, -50.09857466, 0.005434191,
            97.67068802, 319.5738329, 2914.002444, -325.155281,
            0.005231489, 115.240647, 267.6054819, 466.3908327,
            -42.29682113, 0.008532408, 449.9501918, 267.6060633,
            567.6400574])

        # CALCULATION
        do = do2_dofst_volt(
            volt_counts, Voffset, Soc, A, B, C, E, pres, temp, salt, lat, lon)

        # ASSERT INTERMEDIATE AND FINAL OUTPUTS
        #np.testing.assert_allclose(do_int, do_int_expected, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(do, do_expected, rtol=1e-6, atol=1e-6)
