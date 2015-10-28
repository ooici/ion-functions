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
from ion_functions.data.do2_functions import dosta_Topt_volt_to_degC
from ion_functions.data.do2_functions import dosta_phase_volt_to_degree
from ion_functions.data.do2_functions import do2_SVU
from ion_functions.data.do2_functions import o2_counts_to_uM
from ion_functions.data.do2_functions import do2_salinity_correction
from ion_functions.data.do2_functions import do2_dofst_frequency
from ion_functions.data.do2_functions import do2_dofst_volt
#import pdb

SYSTEM_FILLVALUE = -999999999


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

        2015-04-10: Russell Desiderio. Added test for implementation of
                    calibration coefficients as time-vectorized arguments.
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

        ### R. Desiderio
        # 10-Apr-2015: test new CI implementation of time-vectorized cal coeffs
        tval = freq.shape[0]
        A = np.tile(A, tval)
        B = np.tile(B, tval)
        C = np.tile(C, tval)
        E = np.tile(E, tval)
        Foffset = np.tile(Foffset, tval)
        Soc = np.tile(Soc, tval)

        # 05-Aug-2015: test replacement of type int system fillvalues with nan.
        freq[7], freq[14], freq[23] = SYSTEM_FILLVALUE, SYSTEM_FILLVALUE, SYSTEM_FILLVALUE
        do_expected[7], do_expected[14], do_expected[23] = np.nan, np.nan, np.nan

        # CALCULATION
        do = do2_dofst_frequency(
            freq, Foffset, Soc, A, B, C, E, pres, temp, salt, lat, lon)

        # ASSERT FINAL OUTPUTS
        np.testing.assert_allclose(do, do_expected, rtol=1e-6, atol=1e-6)

    def test_dofst_voltage(self):
        """ DOFST voltage test
        Unit Test of the do2_dofst_volt function in
        ion_functions/data/do2_functions.py for calculation of
        oxygen from the voltage output of an SBE 43 instrument (DOFST).

        NOTE:
        SBE43 oxygen sensors are connected to a SBE 16+ V2 CTD which
        provides the inputs of salinity, temperature, and pressure.

        2015-04-10: Russell Desiderio. Added test for implementation of
                    calibration coefficients as time-vectorized arguments.
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

        ### R. Desiderio
        # 10-Apr-2015: test new CI implementation of time-vectorized cal coeffs
        tval = volt_counts.shape[0]
        A = np.tile(A, tval)
        B = np.tile(B, tval)
        C = np.tile(C, tval)
        E = np.tile(E, tval)
        Voffset = np.tile(Voffset, tval)
        Soc = np.tile(Soc, tval)

        # 05-Aug-2015: test replacement of type int system fillvalues with nan.
        volt_counts[4], do_expected[4] = SYSTEM_FILLVALUE, np.nan
        volt_counts[11], do_expected[11] = SYSTEM_FILLVALUE, np.nan
        volt_counts[17], do_expected[17] = SYSTEM_FILLVALUE, np.nan

        # CALCULATION
        do = do2_dofst_volt(
            volt_counts, Voffset, Soc, A, B, C, E, pres, temp, salt, lat, lon)

        # ASSERT FINAL OUTPUTS
        np.testing.assert_allclose(do, do_expected, rtol=1e-6, atol=1e-6)

    def test_do2_SVU(self):
        """ DOSTA Stern-Volmer-Uchida (SVU) test
        Unit Test of the do2_SVU function in
        ion_functions/data/do2_functions.py for calculation of
        oxygen with the SVU equation from phase, temperature, and
        multipoint calibration coefficients from a Aanderaa oxygen
        optode.

        2015-04-10: Russell Desiderio. Added test for implementation of
                    calibration coefficients as time-vectorized arguments.
        2015-08-10: Russell Desiderio. Added test of conc_coef implementation.
        """
        calphase = np.array([
            32., 32., 32., 32., 39.825, 39.825, 39.825,
            39.825, 47.65, 47.65, 47.65, 47.65, 55.475, 55.475,
            55.475, 55.475, 63.3, 63.3, 63.3, 63.3])

        temp = np.array([
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.])

        csv = np.array([
            0.002848, 0.000114, 1.51e-6, 70.42301, -0.10302,
            -12.9462, 1.265377])

        svu_expected = np.array([
            3.67038772e+02,   2.89131248e+02,   2.32138540e+02,
            1.89376503e+02,   2.06105882e+02,   1.61517794e+02,
            1.28983515e+02,   1.04634998e+02,   1.12481219e+02,
            8.72771165e+01,   6.89718258e+01,   5.53355979e+01,
            5.12416074e+01,   3.87165116e+01,   2.97183595e+01,
            2.30890079e+01,   8.06153583e+00,   4.47641240e+00,
            2.04072696e+00,   3.51926228e-01])

        # SVU CALCULATION
        do_svu = do2_SVU(calphase, temp, csv)
        # SVU ASSERT
        np.testing.assert_array_almost_equal(do_svu, svu_expected, decimal=6)

        # R. Desiderio, 10-Apr-2015:
        # test new CI implementation of time-vectorized cal coeffs
        tval = calphase.shape[0]
        csv = np.tile(csv, (tval, 1))
        # SVU CALCULATION
        do_svu = do2_SVU(calphase, temp, csv)
        # SVU ASSERT
        np.testing.assert_array_almost_equal(do_svu, svu_expected, decimal=6)

        # R. Desiderio, 10-Aug-2015.  test implementation of conc_coef
        #               28-Oct-2015.  updated to handle 1D conc_coef arguments (case 3)

        offset, slope = 10.0, 0.9
        svu_expected = offset + slope * svu_expected

        # case 1:
        conc_coef = np.array([[offset, slope]])

        # SVU CALCULATION
        do_svu = do2_SVU(calphase, temp, csv, conc_coef)
        # SVU ASSERT
        np.testing.assert_array_almost_equal(do_svu, svu_expected, decimal=6)

        # case 2:
        conc_coef = np.tile(np.array([offset, slope]), (calphase.size, 1))

        # SVU CALCULATION
        do_svu = do2_SVU(calphase, temp, csv, conc_coef)
        # SVU ASSERT
        np.testing.assert_array_almost_equal(do_svu, svu_expected, decimal=6)

        # case 3 (added so that 1D entries into Omaha cal sheets won't crash):
        conc_coef = np.array([offset, slope])

        # SVU CALCULATION
        do_svu = do2_SVU(calphase, temp, csv, conc_coef)
        # SVU ASSERT
        np.testing.assert_array_almost_equal(do_svu, svu_expected, decimal=6)

    def test_salinity_correction(self):
        """ DOSTA salinity correction test
        Unit Test of the do2_salinity_correction function in
        ion_functions/data/do2_functions.py for correction of
        oxygen for pressure and salinity.

        2015-04-10: Russell Desiderio. No change, no calcoeffs in this
                    function's argument list.
        """
        do_svu = np.array([
            3.67038772e+02,   2.89131248e+02,   2.32138540e+02,
            1.89376503e+02,   2.06105882e+02,   1.61517794e+02,
            1.28983515e+02,   1.04634998e+02,   1.12481219e+02,
            8.72771165e+01,   6.89718258e+01,   5.53355979e+01,
            5.12416074e+01,   3.87165116e+01,   2.97183595e+01,
            2.30890079e+01,   8.06153583e+00,   4.47641240e+00,
            2.04072696e+00,   3.51926228e-01])

        pres = np.array([  # pressure in dbars
            60.52, 72.58, 31.42, 70.82, 74.87,
            29.33, 30.95, 43.58, 65.37, 29.46,
            31.03, 74.57, 75.07, 57.92, 42.98,
            42.98, 42.98, 0.00, 42.98, 42.98])

        salt = np.array([  # practical salinity in psu
            34.1145, 34.2845, 33.2464, 33.5524, 33.5619,
            33.2512, 33.2609, 33.2716, 33.4191, 33.2710,
            33.2808, 33.5483, 33.5424, 33.3458, 0.0000,
            37.7843, 35.7594, 33.3313, 33.3132, 33.3132])

        temp = np.array([
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.,
            10., 16.67, 23.33, 30.])

        lat = np.array([  # latitude in decimal degrees N
            50.0, -42.0, 39.0, 60.0, 45.0, 60.0, 39.0, 45.0, 50.0,
            -42.0, 50.0, 39.0, 60.0, 45.0, -42.0, 60.0, 45.0, -42.0,
            50.0, 39.0])

        lon = np.array([  # longitude in decimal degrees E
            145.0, -42.0, -70.5, 39.0, -125.0, 39.0, -70.5, -125.0,
            145.0, -42.0, 145.0, -70.5, 39.0, -125.0, -42.0, 39.0,
            -125.0, -42.0, 145.0, -70.5])

        do = do2_salinity_correction(do_svu, pres, temp, salt, lat, lon)

        do_expected = np.array([
            2.88296666e+02,   2.29721597e+02,   1.87629963e+02,
            1.54646004e+02,   1.62604386e+02,   1.29058268e+02,
            1.04241732e+02,   8.55210885e+01,   8.88041121e+01,
            6.97282873e+01,   5.57344581e+01,   4.51939082e+01,
            4.04323196e+01,   3.09440694e+01,   2.98348603e+01,
            1.83468213e+01,   6.25456746e+00,   3.57150705e+00,
            1.64934214e+00,   2.87559085e-01])

        np.testing.assert_array_almost_equal(do, do_expected, decimal=6)

    def test_o2_counts_to_uM(self):
        """
            2015-08-05: Russell Desiderio. Initial code.
        """
        counts = np.array([100000, SYSTEM_FILLVALUE, 2000000, 4000000, SYSTEM_FILLVALUE])
        do_xpctd = np.array([0.0, np.nan, 190.0, 390.0, np.nan])

        do_calc = o2_counts_to_uM(counts)

        np.testing.assert_array_almost_equal(do_calc, do_xpctd, decimal=6)

    def test_dosta_Topt_volt_to_degC(self):
        """
            2015-08-05: Russell Desiderio. Initial code.
        """
        volts = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        xpctd = np.array([-5.0, 3.0, 11.0, 19.0, 27.0, 35.0])

        calc = dosta_Topt_volt_to_degC(volts)

        np.testing.assert_array_almost_equal(calc, xpctd, decimal=6)

    def test_dosta_phase_volt_to_degree(self):
        """
            2015-08-05: Russell Desiderio. Initial code.
        """
        volts = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        xpctd = np.array([10.0, 22.0, 34.0, 46.0, 58.0, 70.0])

        calc = dosta_phase_volt_to_degree(volts)

        np.testing.assert_array_almost_equal(calc, xpctd, decimal=6)
