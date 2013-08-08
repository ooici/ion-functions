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
import datetime as dt

from ion_functions.data.met_functions import windavg_mag_corr_east, windavg_mag_corr_north


@attr('UNIT', group='func')
class TestMetFunctionsUnit(BaseUnitTestCase):

    def test_metbk_windavg(self):
        # Date, lat, and lon provided approximately equal magnetic
        # declinations of -17 and +17 degrees

        # have to convert these dates to ntp timestamp (seonds since
        # 1900-01-01)
        date_str = np.array([
            '5/30/2013', '5/30/2013', '5/30/2013', '5/30/2013', '5/30/2013',
            '5/30/2013', '5/30/2013', '5/30/2013'])
        convert_to_ntp = lambda x: (
            dt.datetime.strptime(x, '%m/%d/%Y') -
            dt.datetime(1900, 1, 1)).total_seconds()
        date_ts = map(convert_to_ntp, date_str)

        lat = np.array([
            43.34, 43.34, 43.34, 43.34,
            47.767, 47.767, 47.767, 47.767])
        lon = np.array([
            -66, -66, -66, -66,
            -126, -126, -126, -126])

        ve = np.array([
            2.47, -2.47, -2.47, 2.47,
            2.47, -2.47, -2.47, 2.47])
        vn = np.array([
            6.52, 6.52, -6.52, -6.52,
            6.52, 6.52, -6.52, -6.52])

        ve_expected = np.array([
            0.46, -4.27, -0.46, 4.27,
            4.27, -0.46, -4.27, 0.46])
        vn_expected = np.array([
            6.96, 5.51, -6.96, -5.51,
            5.51, 6.96, -5.51, -6.96])

        ve_cor = windavg_mag_corr_east(ve, vn, lat, lon, date_ts)
        vn_cor = windavg_mag_corr_north(ve, vn, lat, lon, date_ts)

        # test data was only given to 2 decimals (despite that the
        # function can calculate to better accuracy based on comparison
        # of this function to the DPS Matlab function). So here the test
        # data is only tested to 2 decimal places
        np.testing.assert_array_almost_equal(ve_cor, ve_expected, decimal=2)
        np.testing.assert_array_almost_equal(vn_cor, vn_expected, decimal=2)

