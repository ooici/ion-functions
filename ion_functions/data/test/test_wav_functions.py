"""
@package ion_functions.test.test_wav_functions
@file ion_functions/test/test_wav_functions.py
@author Russell Desiderio
@brief Unit tests for wav_functions module
"""

import numpy as np

from ion_functions.test.base_test import BaseUnitTestCase
from nose.plugins.attrib import attr
from ion_functions.data import wav_functions as wv
from ion_functions.utils import fill_value as vfill


@attr('UNIT', group='func')
class TestWAVFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        # for calculating dataproducts from # of values, initial value, and interval
        self.nvalue = 5
        self.value0 = 0.01
        self.deltav = 0.001
        self.ivalue = np.array([0.010, 0.011, 0.012, 0.013, 0.014], ndmin=2)
        self.nrep = 10

        # for testing the correction of buoy displacement values for magnetic declination
        # (some values taken from the unit test module test_adcp_functions.py)
        self.lat = np.array([50.0, 45.0])
        self.lon = np.array([-145.0, -128.0])
        self.ntp = np.array([3545769600.0, 3575053740.0])
        self.xx = np.array([[0.2175, -0.2814, -0.1002, 0.4831, 1.2380],
                            [0.2455, 0.6218, -0.1807, 0.0992, -0.9063]])
        self.yy = np.array([[-0.3367, -0.1815, -1.0522, -0.8676, -0.8919],
                            [0.2585, -0.8497, -0.0873, 0.3073, 0.5461]])
        # set expected results -- magnetic variation correction applied
        # (computed in Matlab using above values and mag_var.m)
        self.xx_cor = np.array([[0.1099, -0.3221, -0.4025, 0.2092, 0.9243],
                                [0.3087, 0.3555, -0.1980, 0.1822, -0.7144]])
        self.yy_cor = np.array([[-0.3855, -0.0916, -0.9773, -0.9707, -1.2140],
                                [0.1783, -0.9911, -0.0325, 0.2666, 0.7805]])

    def test_wav_triaxys_nondir_freq(self):
        """
        Tests calculation of non-directional wave frequency bin values for WAVSS instruments.

        Values were not defined in DPS, are created above.

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

        Implemented by Russell Desiderio, 2014-04-08.
        """
        # the single input record case
        desired = self.ivalue
        actual = wv.wav_triaxys_nondir_freq(self.nvalue, self.value0, self.deltav)
        # test
        np.testing.assert_allclose(actual, desired, rtol=1e-8, atol=0)

        # the multi-record case -- inputs
        nvalues = np.repeat(self.nvalue, self.nrep)
        value0s = np.repeat(self.value0, self.nrep)
        deltavs = np.repeat(self.deltav, self.nrep)
        # the multi-record case -- outputs
        desired = np.tile(self.ivalue, (self.nrep, 1))
        actual = wv.wav_triaxys_nondir_freq(nvalues, value0s, deltavs)
        # test
        np.testing.assert_allclose(actual, desired, rtol=1e-8, atol=0)

    def test_wav_triaxys_dir_freq(self):
        """
        Tests calculation of directional wave frequency bin values for WAVSS instruments.

        Values were not defined in DPS, are created here and above.

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

        Implemented by Russell Desiderio, 2014-04-08.
        """
        # the single input record case
        desired = self.ivalue
        actual = wv.wav_triaxys_dir_freq(self.nvalue, self.nvalue, self.value0, self.deltav)
        # test
        np.testing.assert_allclose(actual, desired, rtol=1e-8, atol=0)

        # the multi-record case -- all nvalues_dir are equal -- inputs
        nvalues = np.repeat(self.nvalue, self.nrep)
        value0s = np.repeat(self.value0, self.nrep)
        deltavs = np.repeat(self.deltav, self.nrep)
        # the multi-record case -- all nvalues_dir are equal -- outputs
        desired = np.tile(self.ivalue, (self.nrep, 1))
        actual = wv.wav_triaxys_dir_freq(nvalues, nvalues, value0s, deltavs)
        # test
        np.testing.assert_allclose(actual, desired, rtol=1e-8, atol=0)

        # the multi-record case -- all nvalues_dir are not the same -- inputs
        nvalues_nondir = np.repeat(self.nvalue, self.nrep)
        nvalues_dir = np.array([4, 5, 3, 5, 4, 5, 1, 2, 5, 3])
        value0s = np.repeat(self.value0, self.nrep)
        deltavs = np.repeat(self.deltav, self.nrep)
        # the multi-record case -- all nvalues_dir are not the same -- outputs
        desired = np.array([[0.010, 0.011, 0.012, 0.013, vfill],
                            [0.010, 0.011, 0.012, 0.013, 0.014],
                            [0.010, 0.011, 0.012, vfill, vfill],
                            [0.010, 0.011, 0.012, 0.013, 0.014],
                            [0.010, 0.011, 0.012, 0.013, vfill],
                            [0.010, 0.011, 0.012, 0.013, 0.014],
                            [0.010, vfill, vfill, vfill, vfill],
                            [0.010, 0.011, vfill, vfill, vfill],
                            [0.010, 0.011, 0.012, 0.013, 0.014],
                            [0.010, 0.011, 0.012, vfill, vfill]])

        actual = wv.wav_triaxys_dir_freq(nvalues_nondir, nvalues_dir, value0s, deltavs)
        # test
        np.testing.assert_allclose(actual, desired, rtol=1e-8, atol=0)

    def test_wav_triaxys_buoymotion_time(self):
        """
        Tests calculation of times corresponding to (x,y,z) buoy displacement measurements
            for WAVSS instruments.

        Values were not defined in DPS, are created above.

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

        Implemented by Russell Desiderio, 2014-04-08.
        """
        # the single input record case
        ntp_timestamp = 3176736750.736
        desired = self.ivalue + ntp_timestamp
        actual = wv.wav_triaxys_buoymotion_time(ntp_timestamp, self.nvalue, self.value0, self.deltav)
        # test to msec
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.001)

        # the multi-record case -- inputs
        nvalues = np.repeat(self.nvalue, self.nrep)
        value0s = np.repeat(self.value0, self.nrep)
        deltavs = np.repeat(self.deltav, self.nrep)
        # the multi-record case -- outputs
        desired = np.tile(self.ivalue + ntp_timestamp, (self.nrep, 1))
        actual = wv.wav_triaxys_buoymotion_time(ntp_timestamp, nvalues, value0s, deltavs)
        # test to msec
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.001)

    def test_wav_triaxys_correct_mean_wave_direction(self):
        """
        Tests magnetic declination correction of mean wave direction WAVSTAT-D_L0
        from WAVSS instruments.

        Values were not defined in DPS. The values for the magnetic declination
        test are calculated directly from the magnetic_declination function in
        ion_functions/data/generic_functions.py.

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

        Implemented by Russell Desiderio, 2014-04-08.
        """
        # the single input record case
        dir_raw = 50.0
        lat = 45.0
        lon = -128.0
        ntp_ts = 3575053740.0
        desired = 50.0 + 16.461005
        actual = wv.wav_triaxys_correct_mean_wave_direction(dir_raw, lat, lon, ntp_ts)
        # test to first decimal place, in case model changes
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.1)

        # the multi-record case -- inputs
        # test "going around the corner" in both directions.
        dir_raw = np.array([50.0, 350.0, 1.0])
        lat = np.array([45.0, 45.0, 80.0])
        lon = np.array([-128.0, -128.0, 0.0])
        ntp_ts = np.array([3575053740.0, 3575053740.0, 3471292800.0])
        # the multi-record case -- outputs
        desired = np.array([66.461005, 366.461005 - 360.0, 1.0 - 6.133664 + 360.0])
        actual = wv.wav_triaxys_correct_mean_wave_direction(dir_raw, lat, lon, ntp_ts)
        # test
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.1)

    def test_wav_triaxys_correct_directional_wave_direction(self):
        """
        Tests magnetic declination correction of directional wave directions WAVSTAT-DDS_L0
        from WAVSS instruments.

        Values were not defined in DPS. The values for the magnetic declination
        test are calculated directly from the magnetic_declination function in
        ion_functions/data/generic_functions.py.

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

        Implemented by Russell Desiderio, 2014-04-10.
        """
        # the single input record case - no fill values (nfreq_dir = nfreq_nondir)
        dir_raw = np.array([50.0, 1.0, 359.0, 180.0, 245.0])
        lat = 45.0
        lon = -128.0
        ntp_ts = 3575053740.0
        # outputs
        desired = np.array([66.461, 17.461, 15.461, 196.461, 261.461], ndmin=2)
        actual = wv.wav_triaxys_correct_directional_wave_direction(dir_raw, lat, lon, ntp_ts)
        # test to first decimal place, in case model changes
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.1)

        # the single input record case - with fill values (nfreq_dir < nfreq_nondir)
        dir_raw = np.array([50.0, 1.0, 359.0, vfill, vfill])
        lat = 45.0
        lon = -128.0
        ntp_ts = 3575053740.0
        # outputs
        desired = np.array([66.461, 17.461, 15.461, vfill, vfill], ndmin=2)
        actual = wv.wav_triaxys_correct_directional_wave_direction(dir_raw, lat, lon, ntp_ts)
        # test to first decimal place, in case model changes
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.1)

        # the multi-record case -- inputs
        # test "going around the corner" in both directions.
        dir_raw = np.array([[50.0, 350.0, 1.0, 170.0, 240.0, vfill],
                            [150.0, 250.0, 11.0, vfill, vfill, vfill],
                            [50.0, 350.0, 1.0, 170.0, vfill, vfill]])
        lat = np.array([45.0, 45.0, 80.0])
        lon = np.array([-128.0, -128.0, 0.0])
        ntp_ts = np.array([3575053740.0, 3575053740.0, 3471292800.0])
        # the multi-record case -- outputs
        desired = np.array([[66.461, 6.461, 17.461, 186.461, 256.461, vfill],
                            [166.461, 266.461, 27.461, vfill, vfill, vfill],
                            [43.866, 343.866, 354.866, 163.866, vfill, vfill]])
        actual = wv.wav_triaxys_correct_directional_wave_direction(dir_raw, lat, lon, ntp_ts)
        # test
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.1)

    def test_wav_triaxys_magcor_buoymotion_x(self):
        """
        Tests calculation of magnetic corrections to eastward buoy displacements for WAVSS instruments.

        Values were not defined in DPS, are created as documented in setup module above.

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

        Implemented by Russell Desiderio, 2014-04-10.
        """
        # the single input record case
        lat = self.lat[0]
        lon = self.lon[0]
        ntp = self.ntp[0]
        xx = self.xx[0, :]
        yy = self.yy[0, :]
        # outputs
        desired = np.array(self.xx_cor[0, :], ndmin=2)
        actual = wv.wav_triaxys_magcor_buoymotion_x(xx, yy, lat, lon, ntp)
        # test
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.0001)

        # multiple records
        desired = np.array(self.xx_cor)
        actual = wv.wav_triaxys_magcor_buoymotion_x(self.xx, self.yy, self.lat, self.lon, self.ntp)
        # test
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.0001)

    def test_wav_triaxys_magcor_buoymotion_y(self):
        """
        Tests calculation of magnetic corrections to northward buoy displacements for WAVSS instruments.

        Values were not defined in DPS, are created as documented in setup module above.

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

        Implemented by Russell Desiderio, 2014-04-10.
        """
        # the single input record case
        lat = self.lat[0]
        lon = self.lon[0]
        ntp = self.ntp[0]
        xx = self.xx[0, :]
        yy = self.yy[0, :]
        # outputs
        desired = np.array(self.yy_cor[0, :], ndmin=2)
        actual = wv.wav_triaxys_magcor_buoymotion_y(xx, yy, lat, lon, ntp)
        # test
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.0001)

        # multiple records
        desired = np.array(self.yy_cor)
        actual = wv.wav_triaxys_magcor_buoymotion_y(self.xx, self.yy, self.lat, self.lon, self.ntp)
        # test
        np.testing.assert_allclose(actual, desired, rtol=0, atol=0.0001)
