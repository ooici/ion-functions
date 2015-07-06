#!/usr/bin/env python

"""
@package ion_functions.test.generic_functions
@file ion_functions/test/generic_functions.py
@author Christopher Wingard, Stuart Pearce, Russell Desiderio
@brief Unit tests for generic_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import generic_functions as gfunc
from ion_functions.data.generic_functions import SYSTEM_FILLVALUE


# for test_replace_fill_with_nan
INST_FILLVALUE = -32768
ZERO_FILLVALUE = 0


@attr('UNIT', group='func')
class TestGenericFunctionsUnit(BaseUnitTestCase):

    def test_replace_fill_with_nan_with_vectortime_variables(self):
        """
        Description:

            Tests replace_fill_with_nan for variables that are time-vectorized
            with more than one time point. Functionally this means that the lead
            dimension of these test variables will have a shape greater than 1.

            Cases (instances of the system fill value are always processed):
                no instrument fill value
                one instrument fill value of 0
                two instrument fill values

        Implemented by:

            2013-06-11: Russell Desiderio. Initial code.
        """
        # constants for convenience in filling in test variables
        sfill = SYSTEM_FILLVALUE
        zfill = ZERO_FILLVALUE
        ifill = INST_FILLVALUE

        # set up 6 integer arrays to be tested, v1-v6.
        # this will also designate their ordering on input and output
        v1_good_1D = np.array([1, 2, 3, 4, 5])
        v2_good_2D = np.tile(v1_good_1D, (3, 1))
        v3_good_3D = np.tile(v2_good_2D, (3, 1, 1))

        v4_fill_1D = np.array([11, ifill, 13, zfill, sfill])
        v5_fill_2D = np.vstack((v4_fill_1D, v1_good_1D, v4_fill_1D))

        v6_fill_3D = np.copy(v3_good_3D)
        v6_fill_3D[0, 0, 0] = zfill
        v6_fill_3D[1, 1, 1] = sfill
        v6_fill_3D[2, 2, 2:5] = ifill  # last 3 elements in row '2'

        # and a float array which should pass through unchanged
        v7_float_array = v3_good_3D.astype(float)
        v7_float_array[0, 1, 2] = np.nan

        #SET INPUTS
        v1, v2, v3 = np.copy(v1_good_1D), np.copy(v2_good_2D), np.copy(v3_good_3D)
        v4, v5, v6 = np.copy(v4_fill_1D), np.copy(v5_fill_2D), np.copy(v6_fill_3D)
        v7 = np.copy(v7_float_array)

        ### CASE 1: no instrument fill values, only system fill values
        # set expected outputs
        x1 = v1_good_1D.astype(float)
        x2 = v2_good_2D.astype(float)
        x3 = v3_good_3D.astype(float)
        x4 = np.array([11.0, float(ifill), 13.0, float(zfill), np.nan])
        x5 = np.vstack((x4, x1, x4))
        x6 = v6_fill_3D.astype(float)
        x6[1, 1, 1] = np.nan
        x7 = np.copy(v7_float_array)

        # (1a): instrument_fillvalue = None
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            None, v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

        # (1b) instrument_fillvalue = empty list
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            [], v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

        # (1c) instrument_fillvalue: empty ndarray
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            np.array([]), v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

        ### CASE (2): one instrument fill value = 0
        # reset expected outputs; x1-x3 and x7 are unchanged
        x4[3] = np.nan  # x4 is now reset
        x5 = np.vstack((x4, x1, x4))
        x6[0, 0, 0] = np.nan  # x6 is now reset

        # (2a) instrument_fillvalue is a scalar
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            zfill, v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

        # (2b) instrument_fillvalue is a one-element list
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            [zfill], v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

        # (2c) instrument_fillvalue is a shape (1,) ndarray
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            np.array([zfill]), v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

        ### CASE (3): two instrument fill values
        # reset expected outputs; x1-x3 and x7 are unchanged
        x4[1] = np.nan  # x4 is now reset
        x5 = np.vstack((x4, x1, x4))
        x6[2, 2, 2:5] = np.nan  # x6 is now reset

        # (3a) instrument_fillvalues are contained in a list
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            [ifill, zfill], v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

        # (3b) instrument_fillvalues are contained in an ndarray
        (c1, c2, c3, c4, c5, c6, c7) = gfunc.replace_fill_with_nan(
            np.array([ifill, zfill]), v1, v2, v3, v4, v5, v6, v7)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)

    def test_replace_fill_with_nan_with_scalartime_variables(self):
        """
        Description:

            Tests replace_fill_with_nan for variables that are time-vectorized
            with one time point. Functionally this means that the lead dimension
            of these test variables will have shape 1.

            Cases (instances of the system fill value are always processed):
                no instrument fill value
                one instrument fill value of 0
                two instrument fill values

        Implemented by:

            2013-06-12: Russell Desiderio. Initial code.
        """
        # constants for convenience in filling in test variables
        sfill = SYSTEM_FILLVALUE
        zfill = ZERO_FILLVALUE
        ifill = INST_FILLVALUE

        # set up 8 integer arrays to be tested, v1-v8.
        # this will also designate their ordering on input and output
        v1 = np.array([500])
        v2 = np.array([ifill])
        v3 = np.array([zfill])
        v4 = np.array([sfill])
        v5 = np.array([[1, 2, 3, 4, 5]])
        v6 = np.array([[11, ifill, 13, zfill, sfill]])
        v7 = np.tile(v5, (3, 1))[np.newaxis, :]
        v8 = np.vstack((v6, v5, v6))[np.newaxis, :]

        ### CASE 1: no instrument fill values, only system fill values
        # set expected outputs
        x1, x2, x3 = v1.astype(float), v2.astype(float), v3.astype(float)
        x4 = np.array([np.nan])
        x5 = v5.astype(float)
        x6 = np.array([[11.0, float(ifill), 13.0, float(zfill), np.nan]])
        x7 = v7.astype(float)
        x8 = np.copy(v8.astype(float))
        x8[0, 0, 4] = np.nan
        x8[0, 2, 4] = np.nan

        # (1a) instrument_fillvalue = None
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            None, v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

        # (1b) instrument_fillvalue: empty list
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            [], v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

        # (1c) instrument_fillvalue: empty ndarray
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            np.array([]), v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

        ### CASE (2): one instrument fill value = 0
        # reset expected outputs; x1, x2, x4, x5 and x7 are unchanged
        x3 = np.array([np.nan])
        x6[0, 3] = np.nan  # x6 is now reset
        x8[0, 0, 3] = np.nan
        x8[0, 2, 3] = np.nan  # x8 is now reset

        # (2a) instrument_fillvalue is a scalar
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            zfill, v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

        # (2b) instrument_fillvalue is a one element list
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            [zfill], v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

        # (2c) instrument_fillvalue is a one element ndarray
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            np.array([zfill]), v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

        ### CASE (3): two instrument fill values
        # reset expected outputs; x1,x3-x5 and x7 are unchanged
        x2 = np.array([np.nan])
        x6[0, 1] = np.nan  # x6 is now reset
        x8[0, 0, 1] = np.nan
        x8[0, 2, 1] = np.nan  # x8 is now reset

        # (3a) instrument_fillvalues are contained in a list
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            [ifill, zfill], v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

        # (3b) instrument_fillvalues are contained in an ndarray
        (c1, c2, c3, c4, c5, c6, c7, c8) = gfunc.replace_fill_with_nan(
            np.array([ifill, zfill]), v1, v2, v3, v4, v5, v6, v7, v8)
        np.testing.assert_array_almost_equal(c1, x1, decimal=8)
        np.testing.assert_array_almost_equal(c2, x2, decimal=8)
        np.testing.assert_array_almost_equal(c3, x3, decimal=8)
        np.testing.assert_array_almost_equal(c4, x4, decimal=8)
        np.testing.assert_array_almost_equal(c5, x5, decimal=8)
        np.testing.assert_array_almost_equal(c6, x6, decimal=8)
        np.testing.assert_array_almost_equal(c7, x7, decimal=8)
        np.testing.assert_array_almost_equal(c8, x8, decimal=8)

    def test_replace_fill_with_nan_with_one_data_argument(self):
        """
        Description:

            Tests replace_fill_with_nan when called with only one data argument.
            Coding replace_fill_with_nan for this was initially not straightforward
            due to the core python\numpy "interface".

            Cases (instances of the system fill value are always processed):
                no instrument fill value
                one instrument fill value of 0
                two instrument fill values

        Implemented by:

            2013-06-23: Russell Desiderio. Initial code.
        """
        # constants for convenience in filling in test variables
        sfill = SYSTEM_FILLVALUE
        zfill = ZERO_FILLVALUE
        ifill = INST_FILLVALUE

        # set up a list of 8 input integer arrays to be tested one at a time
        v4 = np.array([[1, 2, 3, 4, 5]])
        v5 = np.array([[11, ifill, 13, zfill, sfill]])
        v = [np.array([500]),
             np.array([ifill]),
             np.array([zfill]),
             np.array([sfill]),
             v4,
             v5,
             np.tile(v4, (3, 1))[np.newaxis, :],
             np.vstack((v5, v4, v5))[np.newaxis, :]
             ]

        ### CASE(1): no instrument fill values
        # expected values
        x = [np.array([500.0]),
             np.array([ifill]).astype('float'),
             np.array([zfill]).astype('float'),
             np.array([np.nan]),
             np.array([[1.0, 2.0, 3.0, 4.0, 5.0]]),
             np.array([[11, ifill, 13, zfill, np.nan]]).astype('float'),
             v[6].astype(float),
             np.copy(v[7].astype(float))
             ]
        x[7][0, 0, 4] = np.nan
        x[7][0, 2, 4] = np.nan

        for ii in range(len(v)):
            c = gfunc.replace_fill_with_nan(None, v[ii])
            np.testing.assert_array_almost_equal(c, x[ii], decimal=8)

        ### CASE (2): one instrument fill value = 0
        # reset expected outputs; x[0], x[1], x[3], x[4] and x[6] are unchanged
        x[2] = np.array([np.nan])
        x[5][0, 3] = np.nan  # x[5] is now reset
        x[7][0, 0, 3] = np.nan
        x[7][0, 2, 3] = np.nan  # x[7] is now reset

        for ii in range(len(v)):
            c = gfunc.replace_fill_with_nan(zfill, v[ii])
            np.testing.assert_array_almost_equal(c, x[ii], decimal=8)

        ### CASE (3): two instrument fill values
        # reset expected outputs; only list elements 5 and 7 will change
        x[1] = np.array([np.nan])
        x[5][0, 1] = np.nan  # x[5] is now reset
        x[7][0, 0, 1] = np.nan
        x[7][0, 2, 1] = np.nan  # x[7] is now reset

        for ii in range(len(v)):
            c = gfunc.replace_fill_with_nan([ifill, zfill], v[ii])
            np.testing.assert_array_almost_equal(c, x[ii], decimal=8)

    def test_magnetic_declination(self):
        """
        Test magnetic_declination function.

        Some values based on those defined in the WMM document,
        WMM2010testvalues.pdf which accompanies the software.  Others
        were created and checked using online calculators.

        Implemented by Stuart Pearce, April 2013
        """

        lat = np.array([45.0, 45.0, 80.0, 0.0, -80.0, 80.0, 0.0, -80.0])
        lon = np.array([-128.0, -128.0, 0.0, 120.0,
                        240.0, 0.0, 120.0, 240.0])
        z = np.array([0.0, 1000.0, 0.0, 0.0,
                      0.0, 100000.0, 100000.0, 100000.0])
        timestamp = np.array([3575053740.7382507,  # 2013-04-15 22:29:00
                              3575053740.7382507,  # UTC
                              3471292800.0,        # 2010-01-01 UTC
                              3471292800.0,
                              3471292800.0,
                              3471292800.0,
                              3471292800.0,
                              3471292800.0])

        decln = np.array([16.46093044096720, 16.46376239313584, -6.13, 0.97,
                          70.21, -6.57, 0.94, 69.62])

        out = gfunc.magnetic_declination(lat, lon, timestamp, z, -1)

        self.assertTrue(np.allclose(out, decln, rtol=0, atol=1e-2))

    def test_magnetic_correction(self):
        """
        Test magentic_correction function.

        Input values based on those defined in DPS; output values calculated to
            more significant figures using matlab code specified in the DPS.

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        Modified by Russell Desiderio, April 07, 2014. Changed the rtol values
            from 1e4 to 1e-4 to get a fair test. Changed the output values by
            adding more significant figures.
        """
        # apply the magnetic declination correction.
        uu_cor, vv_cor = gfunc.magnetic_correction(16.9604, np.array([0.4413]),
                                                   np.array([0.1719]))

        # test the transform
        self.assertTrue(np.allclose(uu_cor, 0.472251, rtol=1e-4, atol=0))
        self.assertTrue(np.allclose(vv_cor, 0.035692, rtol=1e-4, atol=0))

    def test_ntp_to_unix_time(self):
        """
        Test ntp_to_unix_time function.

        Timestamp Values gathered from various internet sources
        including the NTP FAQ and HOWTO.

        Implemented by Stuart Pearce, April 2013
        """
        ntp_timestamps = np.array([3176736750.7358608,
                                   3359763506.2082224,
                                   3575049755.4380851])

        output = gfunc.ntp_to_unix_time(ntp_timestamps)

        check_values = np.array([967747950.735861,
                                 1150774706.2082224,
                                 1366060955.438085])
        self.assertTrue(np.allclose(output, check_values,
                                    rtol=0, atol=1e-6))

    def test_extract_parameters(self):
        """
        Test extract_parameter function.

        Array values created by author.

        Implemented by Christopher Wingard, April 2013
        """
        in_array = np.array([34, 67, 12, 15, 89, 100, 54, 36])
        self.assertTrue(np.equal(34, gfunc.extract_parameter(in_array, 0)))
        self.assertTrue(np.equal(67, gfunc.extract_parameter(in_array, 1)))
        self.assertTrue(np.equal(12, gfunc.extract_parameter(in_array, 2)))
        self.assertTrue(np.equal(15, gfunc.extract_parameter(in_array, 3)))
        self.assertTrue(np.equal(89, gfunc.extract_parameter(in_array, 4)))
        self.assertTrue(np.equal(100, gfunc.extract_parameter(in_array, 5)))
        self.assertTrue(np.equal(54, gfunc.extract_parameter(in_array, 6)))
        self.assertTrue(np.equal(36, gfunc.extract_parameter(in_array, 7)))
