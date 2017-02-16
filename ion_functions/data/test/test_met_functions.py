#!/usr/bin/env python
"""
@package ion_functions.test.met_functions
@file ion_functions/test/met_functions.py
@author Russell Desiderio, Chris Wingard
@brief Unit tests for met_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
import datetime as dt
import os

import ion_functions.data.met_functions as mb

"""
    List of tests:

        L2 on an hour time base (for netsirr see L1)
            test_met_buoyfls
            test_met_buoyflx
            test_met_frshflx
            test_met_heatflx
            test_met_latnflx
            test_met_mommflx
            test_met_netlirr
            test_met_rainflx
            test_met_sensflx
            test_met_sphum2m
            test_met_stablty
            test_met_tempa2m
            test_met_tempskn
            test_met_wind10m

        L2 on a minute time base
            test_met_heatflx_minute
            test_met_latnflx_minute
            test_met_netlirr_minute
            test_met_sensflx_minute

        L1
            test_met_barpres
            test_met_netsirr
            test_met_netsirr_hourly
            test_met_rainrte
            test_met_salsurf
            test_met_spechum

        Current and Wind
            test_met_current_direction
            test_met_current_speed
            test_met_relwind_direction
            test_met_relwind_speed
            test_met_windavg

        Time
            test_make_hourly_data
            test_met_timeflx
            test_multiple_days
            test_time_calcs_with_actual_metbk_data
            test_warmlayer_time_keys

        Misc
            test_rain_heat_flux
            test_time_vectorized_heights_and_switches
"""


@attr('UNIT', group='func')
class TestMetFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        """
            "SCALAR" TESTSET of 3-element vector data
                First values are selected for the scalar tests.
                All 3 values are selected for L2 perminute data product tests.

            These test values are adapted from array indices 1356, 1357, and
            1358 from the test dataset Revelle10minutesLeg3_r3.mat, from
            ftp://ftp.etl.noaa.gov/users/cfairall/bulkalg/cor3_5/.

            They were chosen because 1357 had significant Solar, IR, and Rainrate
            sensor values. The first time record was changed so that a local time
            (which is a function of longitude) earlier than 6:00 AM would result
            so that the warmlayer code would run and give non-trivial results; time
            records were initially spaced 600 seconds apart, but were changed to
            be 60 seconds apart so that unit tests for the L2 minute data products
            could be easily constructed. The timestamps are not used when scalar
            data are processed.

            The vectorized warmlayer tests producing hourly data do NOT use the
            3-element scalar testset data as inputs.

            Sensor heights are changed so that:
               (1) they are not all the same.
               (2) so that wind sensor height is not the same as that specified
                   for the derived wind data product reference height (10m).
               (3) ctd temperature sensor depth changed from 0.050 m to what
                   would be a more realistic OOI mooring deployment depth.
               (4) 40deg subtracted from longitude (to compare to original matlab
                   warmlayer code results; that code had a time artifact in the
                   calculation of local time).
        """
        tC_sea = np.array([31.125, 31.033, 31.022])
        wnd = np.array([4.805, 3.396, 3.843])
        tC_air = np.array([28.220, 27.961, 27.670])
        relhum = np.array([79.47, 82.30, 82.72])
        timestamp = np.array([2200.0, 2260.0, 2320.0])
        lon = np.array([80.5, 80.6, 80.7]) - 40.0
        ztmpwat = 1.5
        zwindsp = 8.0
        ztmpair = 5.0
        zhumair = 4.0
        lat = np.array([0.1, 1.1, 2.1])
        pr_air = np.array([1005.1, 1006.2, 1007.3])
        Rshort_down = np.array([572.7, 659.0, 634.0])
        Rlong_down = np.array([443.3, 456.0, 441.9])
        rain_rate = np.array([0.0, 17.3, 1.5])    # for scalar test input
        cm_prcp = np.array([0.0, 17.3, 18.8])/60  # for L2 minute product tests
        zinvpbl = 600.0

        """
            Even though the OOI METBK data products will always be computed with
            jwarm = jcool = 1, the more basic test cases (jwarm=0 especially) need
            to be checked first.

            Because jwarm=1 requires at least 2 values to give non-trivial output
            (with respect to the warmlayer calculation), no scalar tests can be
            done to test the warmlayer routine.

            Package the input arguments into tuples for convenience in calling
            the routines to be tested.
        """
        self.kk = 0  # so that any one of the 3 records can be easily selected
        #   NOTE that the correct scalar results for data products involving rain will
        #   only be obtained for self.kk=0, because by definition if there is only one
        #   value for cumulative precipitation, then the rain rate is 0.

        # for scalar unit tests
        self.args_scalar_inputs = (tC_sea[self.kk],
                                   wnd[self.kk],
                                   tC_air[self.kk],
                                   relhum[self.kk],
                                   timestamp[self.kk],
                                   lon[self.kk],
                                   ztmpwat,
                                   zwindsp,
                                   ztmpair,
                                   zhumair,
                                   lat[self.kk],
                                   pr_air[self.kk],
                                   Rshort_down[self.kk],
                                   Rlong_down[self.kk],
                                   rain_rate[self.kk],
                                   zinvpbl)

        # for L2 perminute unit tests
        self.args_permin_inputs = (tC_sea,
                                   wnd,
                                   tC_air,
                                   relhum,
                                   timestamp,
                                   lon,
                                   ztmpwat,
                                   zwindsp,
                                   ztmpair,
                                   zhumair,
                                   lat,
                                   pr_air,
                                   Rshort_down,
                                   Rlong_down,
                                   cm_prcp,
                                   zinvpbl)

        """
            "VECTOR" TESTSET for L2 products producing hourly data:

            The vector testset was constructed so that when processed into hourly data,
            it will give the scalar testset above. The timestamps were changed to span
            3 hours of data, such that the first 12 stamps occurred within the first hour,
            the next 8 the second, the last 10 the third. (The make_hourly_data routine was
            written to take into account the sporadic nature of the METBK sensor output,
            which is roughly once a minute but is usually greater than 60 sec with a value
            around 30 sec thrown in once every 5-ish values). The rain_rate data was replaced
            by cumulative precipitation data because this was required on input by the DPAs.
            (If the RAINRTE data product were to be specified as an input, the current CI
            model would broadcast this per hour data product at the timing of the raw inputs
            (per minute-ish); better to calculate rain rate fresh within the DPA itself).

            The test for the make_hourly_data routine takes as input the following vector
            testset data and should produce the scalar testset data.
        """
        # roughly 6-minute input, to be made into hourly testset.
        # Note that the time points are not evenly spaced; and, there are
        # 12 points in the 1st hour, 8 in the second, and 10 in the last.
        self.tC_sea = np.array([33.162, 30.446, 30.899, 29.540, 32.257,
                                31.351, 33.615, 29.993, 28.635, 31.804,
                                29.088, 32.710, 27.930, 34.136, 28.816,
                                31.476, 29.703, 33.250, 32.363, 30.590,
                                33.194, 30.298, 29.333, 31.746, 32.711,
                                28.850, 32.228, 29.816, 31.263, 30.781])

        self.wnd = np.array([5.189, 4.910, 4.421, 4.630, 4.700, 4.840, 4.980, 5.120, 4.490, 4.560,
                             5.050, 4.770, 3.639, 3.056, 3.542, 3.250, 3.347, 3.736, 3.153, 3.445,
                             4.112, 4.052, 3.813, 3.933, 3.634, 3.753, 3.694, 3.873, 3.992, 3.574])

        self.tC_air = np.array([26.373, 29.246, 29.657, 28.425, 27.604,
                                30.067, 28.015, 25.962, 27.194, 30.478,
                                26.783, 28.836, 28.360, 27.562, 25.165,
                                29.958, 26.763, 25.964, 30.757, 29.159,
                                26.164, 25.733, 27.455, 27.885, 28.316,
                                29.607, 28.746, 26.594, 29.176, 27.024])

        self.relhum = np.array([78.892, 74.268, 75.424, 83.516, 85.828,
                                80.048, 77.736, 82.360, 73.112, 84.672,
                                81.204, 76.580, 83.476, 90.530, 88.179,
                                81.124, 76.421, 85.827, 74.070, 78.773,
                                79.503, 85.937, 80.790, 83.363, 88.510,
                                87.224, 82.077, 84.650, 78.216, 76.930])

        self.timestamp = np.array([864400.0, 864760.0, 865120.0, 865300.0, 865840.0,
                                   865984.0, 866200.0, 866560.0, 867172.0, 867568.0,
                                   867640.0, 867820.0, 868360.0, 868828.0, 869224.0,
                                   869728.0, 870160.0, 870520.0, 871204.0, 871312.0,
                                   871960.0, 872320.0, 872680.0, 873040.0, 873400.0,
                                   873760.0, 874120.0, 874480.0, 874840.0, 875164.0])

        self.lon = np.array([40.205, 43.740, 38.438, 43.151, 39.616,
                             41.384, 37.849, 40.795, 42.562, 37.260,
                             41.973, 39.027, 37.700, 42.340, 41.180,
                             44.660, 36.540, 43.500, 38.860, 40.020,
                             39.750, 38.484, 40.383, 41.650, 43.549,
                             42.283, 37.851, 39.117, 42.916, 41.017])

        ztmpwat = 1.5
        zwindsp = 8.0
        ztmpair = 5.0
        zhumair = 4.0

        self.lat = np.array([0.098, 0.099, 0.095, 0.101, 0.102, 0.104, 0.096, 0.105, 0.107, 0.092,
                             0.093, 0.108, 1.084, 0.990, 1.210, 1.179, 1.116, 1.021, 1.147, 1.053,
                             2.051, 2.116, 1.986, 2.182, 2.018, 2.084, 1.953, 2.214, 2.149, 2.247])

        self.pr_air = np.array([1027.029, 968.551, 983.171, 1070.888, 997.790,
                                939.312, 1041.649, 924.692, 1085.508, 1012.410,
                                953.931, 1056.269, 934.329, 1020.574, 1106.820,
                                905.580, 1049.323, 991.826, 1078.071, 963.077,
                                983.796, 1077.811, 1046.473, 1015.135, 999.465,
                                936.789, 1062.142, 952.458, 968.127, 1030.804])

        self.Rshort_down = np.array([526.884, 568.535, 618.516, 576.865, 585.195,
                                     593.525, 610.186, 535.214, 543.544, 560.205,
                                     601.856, 551.875, 668.414, 630.757, 611.929,
                                     706.071, 649.586, 593.100, 724.900, 687.243,
                                     609.344, 668.518, 619.207, 658.656, 629.069,
                                     648.793, 599.482, 589.620, 678.380, 638.931])

        self.Rlong_down = np.array([407.836, 465.868, 459.420, 472.316, 433.628,
                                    414.284, 446.524, 420.732, 452.972, 478.764,
                                    440.076, 427.180, 423.429, 462.514, 436.457,
                                    475.543, 449.486, 410.400, 488.571, 501.600,
                                    410.967, 445.337, 424.715, 438.463, 459.085,
                                    465.959, 431.589, 472.833, 417.841, 452.211])

        self.cumu_prcp = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                                   0.00, 0.00, 16.25, 16.55, 16.85, 17.15, 17.45, 17.75, 18.05, 18.35,
                                   18.35, 18.45, 18.55, 18.65, 18.75, 18.85, 18.95, 19.05, 19.15, 19.25])

        zinvpbl = 600.0

        self.args_vector_inputs = (self.tC_sea,
                                   self.wnd,
                                   self.tC_air,
                                   self.relhum,
                                   self.timestamp,
                                   self.lon,
                                   ztmpwat,
                                   zwindsp,
                                   ztmpair,
                                   zhumair,
                                   self.lat,
                                   self.pr_air,
                                   self.Rshort_down,
                                   self.Rlong_down,
                                   self.cumu_prcp,
                                   zinvpbl)

        # construct argument list to test time_vectorization of sensor heights.
        # The DPAs were originally written without time vectorization of sensor
        # heights, nor of the jcool and jwarm switches.
        self.npts = self.lat.shape[0]
        self.args_vector_ztmvec = (self.tC_sea,
                                   self.wnd,
                                   self.tC_air,
                                   self.relhum,
                                   self.timestamp,
                                   self.lon,
                                   np.tile(ztmpwat, self.npts),
                                   np.tile(zwindsp, self.npts),
                                   np.tile(ztmpair, self.npts),
                                   np.tile(zhumair, self.npts),
                                   self.lat,
                                   self.pr_air,
                                   self.Rshort_down,
                                   self.Rlong_down,
                                   self.cumu_prcp,
                                   np.tile(zinvpbl, self.npts))

        # construct sets of arguments to test warmlayer_time_keys.py (placement of nans
        # in output when a day's data does not start before sunrise, taken to be 6AM).
        #
        # all 3 days represented by ts_3days do presently start before 6AM;
        # also, do not need to append values for jwarm and jcool
        ts_3days = np.hstack((self.timestamp, self.timestamp+86400.0, self.timestamp+2*86400.0))
        self.args_multiple_days = (np.tile(self.tC_sea, 3),
                                   np.tile(self.wnd, 3),
                                   np.tile(self.tC_air, 3),
                                   np.tile(self.relhum, 3),
                                   ts_3days,
                                   np.tile(self.lon, 3),
                                   ztmpwat,
                                   zwindsp,
                                   ztmpair,
                                   zhumair,
                                   np.tile(self.lat, 3),
                                   np.tile(self.pr_air, 3),
                                   np.tile(self.Rshort_down, 3),
                                   np.tile(self.Rlong_down, 3),
                                   np.tile(self.cumu_prcp, 3),
                                   zinvpbl)

        # the next 3 sets of values are not currently used in the unit tests.
        self.relwinddir = np.array([122.5, 118.5, 120.5, 114.5, 119.5, 116.5, 123.5, 124.5, 115.5, 125.5,
                                    121.5, 117.5, 206.5, 212.5, 208.5, 213.5, 207.5, 210.5, 209.5, 211.5,
                                    296.5, 300.5, 295.5, 303.5, 304.5, 299.5, 298.5, 302.5, 301.5, 297.5])

        self.vle_water = np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                   0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
                                   0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4])

        self.vln_water = self.vle_water + 0.3

    """
    Tests for data products requiring coolskin and warmlayer corrections.

    Description

        Matlab programs were written to generate unit test data for these data products.
        The programs compared values calculated from reference coare version 3.5 matlab code
        and a refactored matlab code. The refactored code calculates two versions of the
        data products: (1) 'old', which are calculated the same as in the reference code;
        and (2) 'new', where the matlab code incorporated the changes to be made in the OOI
        DPA implementation in python (examples: the celsius to kelvin conversion constant was
        corrected from 273.16 to 273.15; more significantly, the rain heat flux calculation
        was corrected). For each of these data products, the reference values generated
        by the original unchanged code were checked against the 'old' values given by the
        refactored code to make sure they were identical. The 'new' values calculated by
        the refactored code were then incorporated as the target test values for the python
        code to calculate.

        The coolskin\warmlayer 'xpctd' arrays are dimensioned as 4 x 3:
        {4 permutations of the [jwarm,jcool] switches} x {3 hourly data points}
        The switch assignments as a function of row are:
                jwarm   jcool
        row 1     0       0
        row 2     0       1
        row 3     1       0
        row 4     1       1

    Implemented by:

        2014-09-13: Russell Desiderio. Initial Code.
        2014-09-20: Russell Desiderio. Added tests of DPAs on sub-hourly testdata.
        2014-10-28: Russell Desiderio. New derivation of rain heat flux coded.
        2014-10-29: Russell Desiderio. Incorporated new unit test values for all data products
                    using warmlayer algorithm (which uses rain heat flux).
        2014-12-29: Russell Desiderio. Incorporated tests on Irminger METBK data.
        2015-07-13: Russell Desiderio. Incorporated tests using time-vectorized input for sensor
                    heights and algorithm switches; tested only met_latnflx as a proxy for all
                    the functions using the coolskin\warmlayer algorithm.
        2015-10-26: Russell Desiderio. Expanded tests developed on 7-13 to test all data products
                    using the coolskin\warmlayer algorithm (fixes problem exposed in redmine ticket
                    #8592).

                    Also fixed test_met_mommflx within this test module, which had been named
                    met_mommflx, so that up to this time the met_mommflx unit test function had
                    never been run. The unit test values turned out to be old values (probably
                    before the new rain heat flux code was written). The procedure detailed in
                    the Description section above was followed to generate the present unit test
                    values from the matlab code, which do agree with the python DPA code results.
        2017-02-03: Russell Desiderio. Added [jwarm,jcool] switch documentation as a function of row
                    in the Description section above.
                    Added test_met_netsirr_hourly.
        2017-02-14: Russell Desiderio. Added tests for the 4 new L2 data products that process
                    perminute data. Re-organized test placement, added documentation.

    References

        OOI (2014). 1341-00370_BULKFLX Artifacts. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> REFERENCE >> Data Product Specification Artifacts
            >> 1341-00370_BULKFLX
    """

    def test_met_buoyfls(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[36.10919371, 28.32153986, 33.56844146],
                          [31.81003830, 25.38513247, 29.90539752],
                          [36.10919371, 28.42887220, 33.68714437],
                          [31.81003830, 25.50702812, 30.04173816]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_buoyfls(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_buoyfls(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_buoyfls(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_buoyflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[38.41602527, 29.98454817, 35.45099283],
                          [33.93646240, 26.92597342, 31.63749229],
                          [38.41602527, 30.09636019, 35.57459211],
                          [33.93646240, 27.05292615, 31.77941140]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_buoyflx(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_buoyflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_buoyflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_frshflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[0.27428050, -17.10209951, -1.27575685],
                          [0.25283019, -17.11663761, -1.29367873],
                          [0.27428050, -17.10156642, -1.27517361],
                          [0.25283019, -17.11603581, -1.29301424]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_frshflx(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_frshflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_frshflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_heatflx(self):
        # cases: [jwarm, jcool]
        # total heat flux = latnflx + sensflx + rainflx + netlirr - netsirr
        latnf = np.array([[184.91334211, 133.43175366, 151.19456789],
                          [170.45205774, 123.62963458, 139.11084942],
                          [184.91334211, 133.78969897, 151.58612581],
                          [170.45205774, 124.03365974, 139.55690009]])

        sensf = np.array([[24.34435275, 19.84019748, 23.96742950],
                          [20.96527540, 17.52684365, 21.07171423],
                          [24.34435275, 19.92468342, 24.06116088],
                          [20.96527540, 17.62294816, 21.17960462]])

        rainf = np.array([[0.00000000, 110.67785202, 9.96157245],
                          [0.00000000, 110.67785202, 9.96157245],
                          [0.00000000, 110.90622850, 9.98146410],
                          [0.00000000, 110.96535859, 9.98688156]])

        netli = np.array([[41.43188924, 28.54298163, 42.15187511],
                          [39.19850664, 26.60211114, 39.95207621],
                          [41.43188924, 28.61328117, 42.22248731],
                          [39.19850664, 26.68338061, 40.03470156]])

        netsi_down = np.array([541.20150, 622.75500, 599.13000])
        netsi = np.tile(netsi_down, (4, 1))

        xpctd = latnf + sensf + rainf + netli - netsi

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_heatflx(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_heatflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_heatflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_latnflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[184.91334211, 133.43175366, 151.19456789],
                          [170.45205774, 123.62963458, 139.11084942],
                          [184.91334211, 133.78969897, 151.58612581],
                          [170.45205774, 124.03365974, 139.55690009]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_latnflx(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_latnflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_latnflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_mommflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[0.030974330742, 0.016090125876, 0.020303436138],
                          [0.030614459503, 0.015910513076, 0.020067006629],
                          [0.030974330742, 0.016096471089, 0.020310838724],
                          [0.030614459503, 0.015918213926, 0.020076101791]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_mommflx(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_mommflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_mommflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_netlirr(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[41.43188924, 28.54298163, 42.15187511],
                          [39.19850664, 26.60211114, 39.95207621],
                          [41.43188924, 28.61328117, 42.22248731],
                          [39.19850664, 26.68338061, 40.03470156]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_netlirr(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_netlirr(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_netlirr(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_rainflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[0.00000000, 110.67785202, 9.96157245],
                          [0.00000000, 110.67785202, 9.96157245],
                          [0.00000000, 110.90622850, 9.98146410],
                          [0.00000000, 110.96535859, 9.98688156]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_rainflx(*args_scalar)
        # for scalar inputs, rain_rate is by definition 0. set xpctd to 0.
        np.testing.assert_allclose(calc, xpctd[0:2, 0], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_rainflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_rainflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_sensflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[24.34435275, 19.84019748, 23.96742950],
                          [20.96527540, 17.52684365, 21.07171423],
                          [24.34435275, 19.92468342, 24.06116088],
                          [20.96527540, 17.62294816, 21.17960462]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_sensflx(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_sensflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_sensflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_sphum2m(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[19.42297996, 19.71138107, 19.47478118],
                          [19.41394979, 19.70380669, 19.46619271],
                          [19.42297996, 19.71166621, 19.47506527],
                          [19.41394979, 19.70411212, 19.46650610]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_sphum2m(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_sphum2m(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_sphum2m(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_stablty(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[-0.75134556, -1.53410291, -1.29155675],   # [00]
                          [-0.67713848, -1.40615623, -1.17706256],   # [01]
                          [-0.75134556, -1.53871067, -1.29521084],   # [10]
                          [-0.67713848, -1.41154284, -1.18138703]])  # [11]

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_stablty(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_stablty(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_stablty(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_tempa2m(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[28.36851533, 28.09024408, 27.81395793],
                          [28.35732615, 28.08200008, 27.80462176],
                          [28.36851533, 28.09054198, 27.81425524],
                          [28.35732615, 28.08234637, 27.80497552]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_tempa2m(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_tempa2m(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_tempa2m(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_tempskn(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[31.12500000, 31.03300000, 31.02200000],
                          [30.76398727, 30.71905805, 30.66606321],
                          [31.12500000, 31.04435295, 31.03340467],
                          [30.76398727, 30.73222317, 30.67945497]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_tempskn(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_tempskn(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_tempskn(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_wind10m(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[4.84933264, 3.42025464, 3.87210524],
                          [4.85069650, 3.42094844, 3.87301659],
                          [4.84933264, 3.42023139, 3.87207813],
                          [4.85069650, 3.42091723, 3.87297987]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_wind10m(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_vector_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_wind10m(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        # VECTOR CASES with time-vectorized sensor heights and algorithm switches
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, self.npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, self.npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_wind10m(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    """
    Tests for coolskin/warmlayer products on native METBK timebase.

    Test values were generated in matlab and implemented by Russell Desiderio.
    Vector cases with time-vectorized sensor heights and algorithm switches
    are not necessary because these data products are produced on the native
    time base of the METBK instrumentation (and the hourly versions have
    already been vetted).
    """
    def test_met_heatflx_minute(self):
        # cases: [jwarm, jcool]
        # total heat flux = latnflx + sensflx + rainflx + netlirr - netsirr
        latnf = np.array([[184.91334211, 133.43175366, 151.19456789],
                          [170.45205774, 123.62963458, 139.11084942],
                          [184.91334211, 133.50681019, 151.24665632],
                          [170.45205774, 123.71656274, 139.18576360]])

        sensf = np.array([[24.34435275, 19.84019748, 23.96742950],
                          [20.96527540, 17.52684365, 21.07171423],
                          [24.34435275, 19.85791503, 23.97990023],
                          [20.96527540, 17.54752327, 21.08983701]])

        rainf = np.array([[0.00000000, 110.67785202, 9.96157245],
                          [0.00000000, 110.67785202, 9.96157245],
                          [0.00000000, 110.72575756, 9.96421969],
                          [0.00000000, 110.73973748, 9.96582507]])

        netli = np.array([[41.43188924, 28.54298163, 42.15187511],
                          [39.19850664, 26.60211114, 39.95207621],
                          [41.43188924, 28.55772741, 42.16127194],
                          [39.19850664, 26.61960347, 39.96595888]])

        netsi_down = np.array([541.20150, 622.75500, 599.13000])
        netsi = np.tile(netsi_down, (4, 1))

        xpctd = latnf + sensf + rainf + netli - netsi

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_heatflx_minute(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_permin_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_heatflx_minute(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_latnflx_minute(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[184.91334211, 133.43175366, 151.19456789],
                          [170.45205774, 123.62963458, 139.11084942],
                          [184.91334211, 133.50681019, 151.24665632],
                          [170.45205774, 123.71656274, 139.18576360]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_latnflx_minute(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_permin_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_latnflx_minute(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_netlirr_minute(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[41.43188924, 28.54298163, 42.15187511],
                          [39.19850664, 26.60211114, 39.95207621],
                          [41.43188924, 28.55772741, 42.16127194],
                          [39.19850664, 26.61960347, 39.96595888]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_netlirr_minute(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_permin_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_netlirr_minute(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_sensflx_minute(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[24.34435275, 19.84019748, 23.96742950],
                          [20.96527540, 17.52684365, 21.07171423],
                          [24.34435275, 19.85791503, 23.97990023],
                          [20.96527540, 17.54752327, 21.08983701]])

        # SCALAR CASES [00] and [01]:
        calc = np.zeros(2)
        iwarm = 0
        for icool in range(2):
            args_scalar = self.args_scalar_inputs + (iwarm, icool)
            calc[icool] = mb.met_sensflx_minute(*args_scalar)
        np.testing.assert_allclose(calc, xpctd[0:2, self.kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASES
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            for icool in range(2):
                ctr = icool + iwarm * 2
                args_vector = self.args_permin_inputs + (iwarm, icool)
                calc[ctr, :] = mb.met_sensflx_minute(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    """
    Tests for products that require neither coolskin nor warmlayer corrections.

    Test values were generated in matlab and implemented by Russell Desiderio.
    """
    def test_met_barpres(self):
        pr_air_mbar = np.array([1005.1, 1006.2, 1007.3])
        xpctd = pr_air_mbar * 100.0

        # SCALAR CASE
        kk = 1
        calc = mb.met_barpres(pr_air_mbar[kk])
        np.testing.assert_allclose(calc, xpctd[kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASE
        calc = mb.met_barpres(pr_air_mbar)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_netsirr(self):
        xpctd = np.array([541.20150, 622.75500, 599.13000])
        Rshort_down = np.array([572.7, 659.0, 634.0])

        # SCALAR CASE
        kk = 1
        calc = mb.met_netsirr(Rshort_down[kk])
        np.testing.assert_allclose(calc, xpctd[kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASE
        calc = mb.met_netsirr(Rshort_down)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_netsirr_hourly(self):
        #
        #    This Vector Case differs from that for test_met_netsirr in that
        #    30 points of Rshort_down data spanning 3 hours (self.Rshort_down)
        #    are input instead of a vector of 3 values (Rshort_down).
        #
        xpctd = np.array([541.20150, 622.75500, 599.13000])
        Rshort_down = np.array([572.7, 659.0, 634.0])

        # SCALAR CASE
        kk = 1
        calc = mb.met_netsirr_hourly(Rshort_down[kk], 6540998.0)
        np.testing.assert_allclose(calc, xpctd[kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASE (input data span 3 hours)
        calc = mb.met_netsirr_hourly(self.Rshort_down, self.timestamp)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_rainrte(self):
        # SCALAR CASE
        xpctd_scalar = 0.0
        calc = mb.met_rainrte(42.0, 6540998.0)
        np.testing.assert_allclose(calc, xpctd_scalar, rtol=1.e-8, atol=0.0)

        # VECTOR CASE (input data span 3 hours)
        xpctd = np.array([0.0, 17.3, 1.5])
        calc = mb.met_rainrte(self.cumu_prcp, self.timestamp)

        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_salsurf(self):
        #
        #   Test values generated from ctd_functions.ctd_pracsal, using ztmpwat [m]
        #  as a proxy for pressure [db].
        #
        xpctd = np.array([20.80589579, 25.51852081, 29.42245514])
        cond_S_m = np.array([3.0, 4.0, 5.0])
        tC_sea = np.array([20.0, 25.0, 30.0])
        ztmpwat = np.array([1.0, 1.5, 2.0])

        # SCALAR CASE
        kk = 1
        calc = mb.met_salsurf(cond_S_m[kk], tC_sea[kk], ztmpwat[kk])
        np.testing.assert_allclose(calc, xpctd[kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASE
        calc = mb.met_salsurf(cond_S_m, tC_sea, ztmpwat)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_spechum(self):
        xpctd = np.array([19.12588067, 19.49340148, 19.23895340])
        tC_air = np.array([28.220, 27.961, 27.670])
        pr_air = np.array([1005.1, 1006.2, 1007.3])
        relhum = np.array([79.47, 82.30, 82.72])

        # SCALAR CASE
        kk = 1
        calc = mb.met_spechum(tC_air[kk], pr_air[kk], relhum[kk])
        np.testing.assert_allclose(calc, xpctd[kk], rtol=1.e-8, atol=0.0)

        # VECTOR CASE
        calc = mb.met_spechum(tC_air, pr_air, relhum)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    """
    Tests for current and wind products.

    Except as noted, test values were generated in matlab and implemented by
    Russell Desiderio.
    """
    def test_met_current_direction(self):
        #    Desiderio 13-Jul-2015.
        #    Revised unit tests to check action of use_velptmn_with_metbk switch

        # test data inputs - hit all quadrants.
        hf_rt3 = np.sqrt(3.0)/2.0
        vle = np.array([0.0, hf_rt3, 0.5, 0.0, -0.5, -hf_rt3,
                        -1.0, -hf_rt3, -0.5, 0.0, 0.5, hf_rt3])
        vln = np.array([0.0, 0.5, hf_rt3, 1.0, hf_rt3, 0.5,
                        0.0, -0.5, -hf_rt3, -1.0, -hf_rt3, -0.5])

        # expected output;  0 = North, 90 = East.
        xpctd = np.array([90.0, 60.0, 30.0, 0.0, 330.0, 300.0,
                          270.0, 240.0, 210.0, 180.0, 150.0, 120.0])

        # (0) no switch - calculated = nans
        calc = mb.met_current_direction(vle, vln)
        np.testing.assert_array_almost_equal(calc, xpctd*np.nan, decimal=4)

        # (1) scalar switch tests
        # (1a) all good data
        calc = mb.met_current_direction(vle, vln, 1)
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=4)
        # (1b) all bad data
        calc = mb.met_current_direction(vle, vln, 0)
        np.testing.assert_array_almost_equal(calc, xpctd*np.nan, decimal=4)

        # (2) vectorized switch tests
        npts = vle.shape[0]
        # (2a) all good data
        calc = mb.met_current_direction(vle, vln, np.ones(npts))
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=4)
        # (2b) all bad data
        calc = mb.met_current_direction(vle, vln, np.zeros(npts))
        np.testing.assert_array_almost_equal(calc, xpctd*np.nan, decimal=4)
        # (2b) mixed data
        use_metbk = np.ones(npts)
        for ii in [1, 6, 9, 10]:
            xpctd[ii] = np.nan
            use_metbk[ii] = 0
        calc = mb.met_current_direction(vle, vln, use_metbk)
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=4)

    def test_met_current_speed(self):
        #"""
        #    Test the surface current algorithm using test data generated in Matlab
        #    from the compass plot function example:
        #
        #    >> rng(0,'twister') % initialize random number generator
        #    >> M = randn(15,15);
        #    >> Z = eig(M);
        #    >> vle = real(Z);
        #    >> vln = imag(Z);
        #    >> crnt = sqrt(vle.^2 + vln.^2);
        #
        #    C. Wingard 2014-07-01
        #    R. Desiderio 13-Jul-2015.
        #       Revised unit tests to check action of use_velptmn_with_metbk switch
        #"""

        # test data inputs
        vle = np.array([-3.1330, -3.1330, -2.9908, 1.0666, 1.0666,
                        2.1770, 2.1770, -0.8572, -0.8572, -1.4486,
                        -1.4486, 0.2149, 0.2149, 1.4811, 0.7050])
        vln = np.array([1.1221, -1.1221, 0.0000, 2.8685, -2.8685,
                        1.6654, -1.6654, 2.2209, -2.2209, 0.9033,
                        -0.9033, 1.7152, -1.7152, 0.0000, 0.0000])

        # expected output
        crnt = np.array([3.3279, 3.3279, 2.9908, 3.0604, 3.0604,
                         2.7410, 2.7410, 2.3806, 2.3806, 1.7072,
                         1.7072, 1.7286, 1.7286, 1.4811, 0.7050])

        # no switch - calculated = nans
        out = mb.met_current_speed(vle, vln)
        np.testing.assert_array_almost_equal(out, crnt*np.nan, decimal=4)

        # compute the surface current - all good data
        out = mb.met_current_speed(vle, vln, 1)
        # and compare the expected to the calculated
        np.testing.assert_array_almost_equal(out, crnt, decimal=4)

        # compute the surface current - all bad data
        out = mb.met_current_speed(vle, vln, 0)
        # and compare the expected to the calculated
        np.testing.assert_array_almost_equal(out, crnt*np.nan, decimal=4)

        # vectorized switch tests
        npts = vle.shape[0]
        # all good data
        out = mb.met_current_speed(vle, vln, np.ones(npts))
        np.testing.assert_array_almost_equal(out, crnt, decimal=4)
        # all bad data
        out = mb.met_current_speed(vle, vln, np.zeros(npts))
        np.testing.assert_array_almost_equal(out, crnt*np.nan, decimal=4)
        # time vectorized, mixed good and bad data
        use_metbk = np.ones(vle.shape[0])
        for ii in [1, 6, 9, 10]:
            crnt[ii] = np.nan
            use_metbk[ii] = 0
        out = mb.met_current_speed(vle, vln, use_metbk)
        np.testing.assert_array_almost_equal(out, crnt, decimal=4)

    def test_met_relwind_direction(self):
        #
        # 2015-07-14. Russell Desiderio. Revised unit tests to check action of
        #                                use_velptmn_with_metbk switch
        #
        hf_rt3 = np.sqrt(3.0)/2.0
        wind_vle = np.array([0.0, hf_rt3, 0.5, 0.0, -0.5, -hf_rt3,
                            -1.0, -hf_rt3, -0.5, 0.0, 0.5, hf_rt3]) * 3.0
        wind_vln = np.array([0.0, 0.5, hf_rt3, 1.0, hf_rt3, 0.5,
                             0.0, -0.5, -hf_rt3, -1.0, -hf_rt3, -0.5]) * 3.0
        npts = wind_vle.shape[0]
        current_vle = np.tile(-0.5 * 2.0, npts)
        current_vln = np.tile(-hf_rt3 * 2.0, npts)

        # to calculate expected values, use met_current_direction, which has already
        # been checked (and not called by met_relwind_direction).
        xpctd = mb.met_current_direction(wind_vle-current_vle, wind_vln-current_vln, 1)

        # no current data and no switch - current data to be considered all bad, by default
        # expected result is xpctd*np.nan
        calc = mb.met_relwind_direction(wind_vle, wind_vln)
        np.testing.assert_array_almost_equal(calc, xpctd*np.nan, decimal=5)

        # no switch - current data to be considered all bad, by default
        # expected result is xpctd*np.nan
        calc = mb.met_relwind_direction(wind_vle, wind_vln, current_vle, current_vln)
        np.testing.assert_array_almost_equal(calc, xpctd*np.nan, decimal=5)

        # all good data: set the use_velptmn switch to 1 in the calling argument list.
        calc = mb.met_relwind_direction(wind_vle, wind_vln, current_vle, current_vln, 1)
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=5)

        # all bad data: set the use_velptmn switch to 0 in the calling argument list.
        # expected result is xpctd*np.nan
        calc = mb.met_relwind_direction(wind_vle, wind_vln, current_vle, current_vln, 0)
        np.testing.assert_array_almost_equal(calc, xpctd*np.nan, decimal=5)

        # time-vectorized switch cases
        # all good
        use_metbk = np.ones(npts)
        calc = mb.met_relwind_direction(wind_vle, wind_vln, current_vle, current_vln, use_metbk)
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=5)
        # all bad
        use_metbk = np.zeros(npts)
        calc = mb.met_relwind_direction(wind_vle, wind_vln, current_vle, current_vln, use_metbk)
        np.testing.assert_array_almost_equal(calc, xpctd*np.nan, decimal=5)
        # mixed good and bad data
        use_metbk = np.ones(npts)
        for ii in [1, 6, 9, 10]:
            xpctd[ii] = np.nan
            use_metbk[ii] = 0
        calc = mb.met_relwind_direction(wind_vle, wind_vln, current_vle, current_vln, use_metbk)
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=5)

    def test_met_relwind_speed(self):
        #
        #2015-07-14. Russell Desiderio. Revised unit tests to check action of
        #                               use_velptmn_with_metbk switch.
        #

        # These tests differ from those for relwind_direction, current_speed, and current_direction
        # in that bad current data should not result in nans for the relative windspeed; rather, the
        # windspeed itself should be calculated as the data product as if the current velocities are 0.

        hf_rt3 = np.sqrt(3.0)/2.0
        wind_vle = np.array([0.0, hf_rt3, 0.5, 0.0, -0.5, -hf_rt3,
                            -1.0, -hf_rt3, -0.5, 0.0, 0.5, hf_rt3]) * 3.0
        wind_vln = np.array([0.0, 0.5, hf_rt3, 1.0, hf_rt3, 0.5,
                             0.0, -0.5, -hf_rt3, -1.0, -hf_rt3, -0.5]) * 3.0
        npts = wind_vle.shape[0]
        current_vle = np.tile(-0.5 * 2.0, npts)
        current_vln = np.tile(-hf_rt3 * 2.0, npts)

        xpctd_goodcurrent = np.array([2.000000, 4.836559, 5.000000, 4.836559, 4.358899, 3.605551,
                                      2.645751, 1.614836, 1.000000, 1.614836, 2.645751, 3.605551])

        xpctd_badcurrent = np.hstack((0.0, np.tile(3.0, npts-1)))

        # no current data and no switch - current data to be considered all bad, by default
        # expected result is xpctd_badcurrent
        calc = mb.met_relwind_speed(wind_vle, wind_vln)
        np.testing.assert_array_almost_equal(calc, xpctd_badcurrent, decimal=5)

        # all bad current data: no switch
        calc = mb.met_relwind_speed(wind_vle, wind_vln, current_vle, current_vln)
        np.testing.assert_array_almost_equal(calc, xpctd_badcurrent, decimal=5)
        # all good current data: set the use_velptmn switch to 1 in the calling argument list.
        calc = mb.met_relwind_speed(wind_vle, wind_vln, current_vle, current_vln, 1)
        np.testing.assert_array_almost_equal(calc, xpctd_goodcurrent, decimal=5)
        # all bad current data: set the use_velptmn switch to 0 in the calling argument list.
        calc = mb.met_relwind_speed(wind_vle, wind_vln, current_vle, current_vln, 0)
        np.testing.assert_array_almost_equal(calc, xpctd_badcurrent, decimal=5)

        # time-vectorized switch cases
        # all good
        use_metbk = np.ones(npts)
        calc = mb.met_relwind_speed(wind_vle, wind_vln, current_vle, current_vln, use_metbk)
        np.testing.assert_array_almost_equal(calc, xpctd_goodcurrent, decimal=5)
        # all bad
        use_metbk = np.zeros(npts)
        calc = mb.met_relwind_speed(wind_vle, wind_vln, current_vle, current_vln, use_metbk)
        np.testing.assert_array_almost_equal(calc, xpctd_badcurrent, decimal=5)
        # mixed good and bad
        use_metbk = np.ones(npts)
        xpctd = np.copy(xpctd_goodcurrent)
        for ii in [1, 6, 9, 10]:
            xpctd[ii] = xpctd_badcurrent[ii]
            use_metbk[ii] = 0
        calc = mb.met_relwind_speed(wind_vle, wind_vln, current_vle, current_vln, use_metbk)
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=5)

    def test_met_windavg(self):
        #"""
        #    Date, lat, and lon provided approximately equal magnetic
        #    declinations of -17 and +17 degrees
        #
        #    have to convert these dates to ntp timestamp (seconds since
        #    1900-01-01)
        #
        #    C. Wingard 2014-07-01
        #"""

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

        ve_cor = mb.met_windavg_mag_corr_east(ve, vn, lat, lon, date_ts)
        vn_cor = mb.met_windavg_mag_corr_north(ve, vn, lat, lon, date_ts)

        # test data was only given to 2 decimals (despite that the
        # function can calculate to better accuracy based on comparison
        # of this function to the DPS Matlab function). So here the test
        # data is only tested to 2 decimal places
        np.testing.assert_array_almost_equal(ve_cor, ve_expected, decimal=2)
        np.testing.assert_array_almost_equal(vn_cor, vn_expected, decimal=2)

    """
    Tests involving time.
    """
    def test_make_hourly_data(self):
        #"""
        #    Note that except for the first two elements, the expected values are the
        #    same as those given above in the "scalar testset".
        #        The first element contains the hourly averages of the cumulative
        #            precipitation; when differenced, this will give the expected
        #            rain rate array [0.0, 17.3, 1,5].
        #        The second element contains the hourly timestamps which have been
        #            coded to be at the midpoint of the bins.
        #    Note also that the output of make_hourly_data is a list of np.arrays.
        #
        #    Initial code: 2014-09-20, Russell Desiderio.
        #"""
        xpctd = [np.array([0., 17.3, 18.8]),
                 np.array([866200., 869800., 873400.]),
                 np.array([40.5, 40.6, 40.7]),
                 np.array([31.125, 31.033, 31.022]),
                 np.array([4.805, 3.396, 3.843]),
                 np.array([28.22, 27.961, 27.67]),
                 np.array([79.47, 82.3, 82.72]),
                 np.array([1005.1, 1006.2, 1007.3]),
                 np.array([572.7, 659., 634.]),
                 np.array([443.3, 456., 441.9]),
                 np.array([0.1, 1.1, 2.1])]

        args = [self.cumu_prcp, self.timestamp, self.lon, self.tC_sea, self.wnd, self.tC_air,
                self.relhum, self.pr_air, self.Rshort_down, self.Rlong_down, self.lat]

        args = mb.condition_data(*args)

        calc = mb.make_hourly_data(*args)

        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=1.e-8)

    def test_met_timeflx(self):
        """
            Uses the timestamp data as in test_make_hourly_data.

            Initial code: 2014-10-22, Russell Desiderio.
        """
        xpctd = np.array([866200., 869800., 873400.])

        calc = mb.met_timeflx(self.timestamp)

        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=1.e-8)

    def test_multiple_days(self):
        """
            Written to test nan results when a day's data does not start before 6AM.
            Use latent heat flux calculation as the test product to monitor.

            Initial code: 2014-10-27, Russell Desiderio.
        """
        # this test originally generated runtime warnings whenever logical indexing or
        # np.where were used on arrays containing nan values; the calculations however
        # gave the expected and desired results. these warnings can be turned off by
        # executing: np.seterr(invalid='ignore'); however, instead i directly trapped
        # out these instances in the code.

        # jwarm=jcool=1 for latent heat flux product
        xpctd_1day = np.array([170.45205774, 124.03365974, 139.55690009])
        xpctd_3day = np.tile(xpctd_1day, 3)

        # first check that a run of 3 consecutive days of good data gives the expected result.
        xpctd = np.copy(xpctd_3day)
        calc = mb.met_latnflx(*self.args_multiple_days)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=1.e-8)
        #timeflx = mb.met_timeflx(self.args_multiple_days[4])
        #print timeflx

        # now check to see what happens if a day of data does not start before 6AM;
        # 3 tests, one each for adding 5 hours to one of the day's timestamps.
        #
        # note that the test values will change and also a different number of
        # binned points can result if a non-integral number of hours is used.
        for ii in range(3):
            # if the np.copy operations are not used, then the tests pass, but for the
            # 'wrong' reasons; for example, in the last iteration all xpctd values become
            # nans and all local times will be later than 6AM (instead of only the 3rd day's
            # times shifting to later than 6AM).
            xpctd = np.copy(xpctd_3day)
            # set ii_th day's values to nans
            xpctd[ii*3:ii*3+3] = np.nan
            alt_time = np.copy(self.args_multiple_days[4])
            # add 5 hours to ii_th day's times
            alt_time[ii*30:ii*30+30] = alt_time[ii*30:ii*30+30] + 5 * 3600.0
            args = self.args_multiple_days[0:4] + (alt_time,) + self.args_multiple_days[5:]
            calc = mb.met_latnflx(*args)
            #timeflx = mb.met_timeflx(alt_time)  # make sure python is creating copies
            #print timeflx, xpctd
            np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=1.e-8)

    def test_time_calcs_with_actual_metbk_data(self):
        """
            Test manipulation of date and time by running a flux algorithm with actual
            raw (each_minute) METBK data (from Irminger Sea deployment).

            This code also provides a template on how to calculate L2 METBK data products
            requiring the warmlayer and coolskin algorithms.

            NOTE: These PRECIPM_L0 data exhibit systematic spiking. The repeating cycle
            is a spike of 3 points followed by 12 points of unaliased data. This pattern
            is evident when plots of these precipm data v. time are zoomed in.

            2014-12-19: Russell Desiderio. Initial code.
            2014-12-29: Russell Desiderio. Extended tests.
        """
        # read in the input L0 test data ...
        file = os.path.join(os.getcwd(),
                            'ion_functions/data/matlab_scripts/metbk/metbk_test_dat.txt')
        with open(file, 'r') as f:
            lines = f.readlines()

        # ... parse it into a data array
        txtdata = []
        for ii in lines:
            txtdata.append((ii.strip().split(",")))
        data = np.array(txtdata, dtype=float)

        # to ensure no loss of precision, the OOI CI timestamps were broken up into a
        # more significant and a less significant field; a timestamp of 3619382459.096
        # is represented by the fields 3619382 and 459.096.
        #
        # reconstitute OOI CI timestamps [sec since 1900-01-01; UT]
        timstmp = data[:, 0] * 1000 + data[:, 1]

        # document input L0 data by assigning DPS variable names to the data columns
        barpres = data[:, 2]
        relhumi = data[:, 3]
        tempair = data[:, 4]
        longirr = data[:, 5]
        precipm = data[:, 6]
        tempsrf = data[:, 7]
        condsrf = data[:, 8]    # not used in any METBK DPA algorithms
        shrtirr = data[:, 9]
        wndrawE = data[:, 10]
        wndrawN = data[:, 11]

        # the test data came from the Irminger Sea deployment, so
        lon = np.copy(barpres)
        lon[:] = -39.
        lat = np.copy(barpres)
        lat[:] = 60.

        # sensor heights
        ztmpwat = 1.5
        zwindsp = 8.0
        ztmpair = 5.0
        zhumair = 4.0

        # correct wind components for magnetic variation; this is necessary in the
        # general case because the L1 VELPT current data product values, required for
        # the relative wind calculation, have themselves been magnetically corrected.
        ve_cor = mb.met_windavg_mag_corr_east(wndrawE, wndrawN, lat, lon, timstmp)
        vn_cor = mb.met_windavg_mag_corr_north(wndrawE, wndrawN, lat, lon, timstmp)
        # for this test, however, assume there is no current.
        zero_1D = np.array([0.0])
        rel_wnd_spd = mb.met_relwind_speed(ve_cor, vn_cor, zero_1D, zero_1D)

        # construct argument tuple for flux data products
        args = (tempsrf, rel_wnd_spd, tempair, relhumi, timstmp, lon, ztmpwat, zwindsp,
                ztmpair, zhumair, lat, barpres, shrtirr, longirr, precipm)

        # latent heat flux, an hourly data product
        latnflx = mb.met_latnflx(*args)
        # calculate time base for latnflx (and all hourly data products)
        hourly_time_base = mb.met_timeflx(timstmp)

        # TESTS:
        # (1) verify that the first and last ooici timestamps correlate with
        #     the datetimes in the raw METBK text files.
        #
        #     to convert ooi ci timestamps to readable datetime strings, first
        #     convert them to posix (unix) time by subtracting the difference in
        #     epochs, then use a function from the datetime module.
        epoch_offset = 2208988800  # number of seconds between 1900 and 1970 epochs
        utc_rawdata_first = dt.datetime.utcfromtimestamp(timstmp[0] - epoch_offset)
        xpctd_first = '2014-09-11 00:00:59.096000'
        np.testing.assert_equal(utc_rawdata_first.isoformat(' '), xpctd_first)

        utc_rawdata_last = dt.datetime.utcfromtimestamp(timstmp[-1] - epoch_offset)
        xpctd_last = '2014-09-17 23:59:38.992000'
        np.testing.assert_equal(utc_rawdata_last.isoformat(' '), xpctd_last)

        # (2) there should be 168 hourly timestamps, because the rawdata spans 7 complete days.
        np.testing.assert_equal(hourly_time_base.shape[0], 168)

        # (3a) the 1st hourly timestamp will be 30 minutes later than the first rawstamp.
        utc_hourly_first = dt.datetime.utcfromtimestamp(hourly_time_base[0] - epoch_offset)
        xpctd_first = '2014-09-11 00:30:59.096000'
        np.testing.assert_equal(utc_hourly_first.isoformat(' '), xpctd_first)

        # (3b) the last hourly timestamp will be (168-1) hours later than the first hourly stamp.
        utc_hourly_last = dt.datetime.utcfromtimestamp(hourly_time_base[-1] - epoch_offset)
        xpctd_last = '2014-09-17 23:30:59.096000'
        np.testing.assert_equal(utc_hourly_last.isoformat(' '), xpctd_last)

        # (4) Local time check: because the longitude is 39 West, Irminger local time
        #     is 39 * 240 seconds = 156 minutes earlier than UTC. Therefore, the
        #     first 4 hourly timestamps are:
        #
        #         UTC           local
        #     00:30:59.096   21:54:59.096
        #     01:30:59.096   22:54:59.096
        #     02:30:59.096   23:54:59.096
        #     03:30:59.096   00:54:59.096
        #
        #     The test is that the first non-NaN data value should be the 4th, because the
        #     previous day's data did not start before 6:00 AM local as required by the
        #     warmlayer algorithm.
        calc = np.isnan(latnflx[0:4])
        xpctd = np.array([True, True, True, False])
        np.testing.assert_equal(calc, xpctd)

    def test_warmlayer_time_keys(self):
        """
        idx_warm, newday, nanmask = warmlayer_time_keys(localdate)

            where

        idx_warm = indices of data records to be processed by the warmlayer routine;
                   these are data for days for which there are data before a threshold
                   time value early in the morning (usually set to equatorial sunrise).
        newday = boolean array: true for the first record of a day, false otherwise.
        nanmask = boolean array: true for indices of data records not to be processed
                  by the warmlayer routine.
        localdate = local (not UTC) date and time [sec since 01-01-1900]

        Initial code: 2014-10-27, Russell Desiderio.
        """
        # test data: set local date and times for each day
        # no times before 6AM (21600 seconds)
        day1 = np.array([30000.0, 34000.0, 40000.0]) + 3 * 86400.0
        # 1st time before 6AM
        day2 = np.array([20000.0, 25000.0, 30000.0, 36000.0]) + 4 * 86400.0
        day3 = np.array([50000.0, 60000.0]) + 7 * 86400.0
        day4 = np.array([10000.0, 15000.0, 20000.0, 26000.0]) + 9 * 86400.0
        # entire test time record
        localdate = np.hstack((day1, day2, day3, day4))

        # expected
        xpctd_idx = np.array([3, 4, 5, 6, 9, 10, 11, 12])
        xpctd_new = np.array([1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0], dtype=bool)
        xpctd_nan = np.array([1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0], dtype=bool)
        xpctd = (xpctd_idx, xpctd_new, xpctd_nan)

        # calculated
        calc = mb.warmlayer_time_keys(localdate)

        np.testing.assert_equal(calc, xpctd)

    """
    Miscellaneous tests.
    """
    def test_rain_heat_flux(self):
        """
            Tests new formulation of rain heat flux, independent of coare bulk algorithms.
        """
        Tsea = 25.0
        Tair = 23.5
        RH = 80.0       # relative humidity [%]
        rainrate = 1.0  # [mm/hour]
        Pr = 1015.0     # atmospheric pressure [mbar]

        # check value is from the matlab code test_rain_heat_flux.m
        xpctd = 4.663360812871
        calc = mb.rain_heat_flux(rainrate, Tsea, Tair, RH, Pr)
        np.testing.assert_allclose(calc, xpctd, rtol=0, atol=1.e-12)

    def test_time_vectorized_heights_and_switches(self):
        """
        Description:

            Time-vectorization of the sensor heights (ztmpwat, ztmpair, zhumair, zwindsp,
            and zinvpbl) and algorithm switches (jwarm, jcool) are tested using the DPAs
            met_latnflx and met_wind10m. The check values are identical to the corresponding
            vector case values.

        Implemented by:

            2014-07-14: Russell Desiderio. Initial Code.
            2015-10-26: Russell Desiderio. Added test using met_wind10m. Minor DPA mods were
                                           made to several DPAs involving binning the invariant
                                           sensor heights as a fix for redmine ticket #8592 (pull
                                           request 232). This documentation added 2017-02-14.
        """
        xpctd = np.array([[184.91334211, 133.43175366, 151.19456789],
                          [170.45205774, 123.62963458, 139.11084942],
                          [184.91334211, 133.78969897, 151.58612581],
                          [170.45205774, 124.03365974, 139.55690009]])

        # latnflx: Time-vectorized cases only
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        npts = self.args_vector_ztmvec[0].shape  # number of time points
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_latnflx(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

        xpctd = np.array([[4.84933264, 3.42025464, 3.87210524],
                          [4.85069650, 3.42094844, 3.87301659],
                          [4.84933264, 3.42023139, 3.87207813],
                          [4.85069650, 3.42091723, 3.87297987]])

        # wind10m: Time-vectorized cases only
        # sensor height arrays in self.args_vector_ztmvec have already been time-vectorized
        npts = self.args_vector_ztmvec[0].shape  # number of time points
        calc = np.zeros((4, 3))
        for iwarm in range(2):
            jwarm_tmvec = np.tile(iwarm, npts)  # time vectorize jwarm
            for icool in range(2):
                ctr = icool + iwarm * 2
                jcool_tmvec = np.tile(icool, npts)  # time vectorize jcool
                args_vector = self.args_vector_ztmvec + (jwarm_tmvec, jcool_tmvec)
                calc[ctr, :] = mb.met_wind10m(*args_vector)
        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

