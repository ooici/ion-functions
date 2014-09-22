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

import ion_functions.data.met_functions as mb


@attr('UNIT', group='func')
class TestMetFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        """
            "SCALAR" TESTSET:

            These test values are adapted from array indices 1356, 1357, and
            1358 from the test dataset Revelle10minutesLeg3_r3.mat, from
            ftp://ftp.etl.noaa.gov/users/cfairall/bulkalg/cor3_5/.

            They were chosen because 1357 had significant Solar, IR, and Rainrate
            sensor values. the first time record was changed so that a local time
            (which is a function of longitude) earlier than 6:00 AM would result
            so that the warmlayer code would run and give non-trivial results; time
            records are spaced 600 seconds apart. The vectorized warmlayer tests do
            NOT use these test data.

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
        timestamp = np.array([2200.0, 2800.0, 3400.0])
        lon = np.array([80.5, 80.6, 80.7]) - 40.0
        ztmpwat = 1.5
        zwindsp = 8.0
        ztmpair = 5.0
        zhumair = 4.0
        lat = np.array([0.1, 1.1, 2.1])
        pr_air = np.array([1005.1, 1006.2, 1007.3])
        Rshort_down = np.array([572.7, 659.0, 634.0])
        Rlong_down = np.array([443.3, 456.0, 441.9])
        rain_rate = np.array([0.0, 17.3, 1.5])
        zinvpbl = 600.0

        """
            Even though the OOI METBK data products will always be computed with
            jwarm = jcool = 1, the more basic test cases (jwarm=0 especially) need
            to be checked first.

            Because jwarm=1 requires at least 2 values to run, no scalar tests can
            be done for this case.

            Package the input arguments into tuples for convenience in calling
            the routines to be tested.
        """
        self.kk = 0  # so that any one of the 3 records can be easily selected
        # NOTE that the correct results for data products involving rain will only
        #     be obtained for self.kk=0, because by definition if there is only one
        #     value for cumulative precipitation, then the rain rate is 0.

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

        """
            VECTOR TESTSET:

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
            testset data and produces the scalar testset data.
        """
        # roughly 6-minute input, to be made into hourly testset.
        # note that the time points are not evenly spaced; and, there
        # are 12 points in the 1st hour, 8 and 10 in the next two.
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

        # the next 3 sets of values are not currently used in the unit tests.
        self.relwinddir = np.array([122.5, 118.5, 120.5, 114.5, 119.5, 116.5, 123.5, 124.5, 115.5, 125.5,
                                    121.5, 117.5, 206.5, 212.5, 208.5, 213.5, 207.5, 210.5, 209.5, 211.5,
                                    296.5, 300.5, 295.5, 303.5, 304.5, 299.5, 298.5, 302.5, 301.5, 297.5])

        self.vle_water = np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                                   0.2, 0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
                                   0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4])

        self.vln_water = self.vle_water + 0.3

    """
    Description

        Tests for data products requiring coolskin and warmlayer corrections.

        The matlab program WarmCoolLayer_OOI_DPA_calculation.m was written to generate
        unit test data for these data products. The program compares values calculated
        from reference coare version 3.5 matlab code and my (Desiderio) refactored matlab
        code. The refactored code calculates two versions of the data products: (1) 'old',
        which are calculated the same as in the reference code; and (2) 'new', where the
        matlab code incorporated the changes to be made in the OOI DPA implementation in
        python (for example, the celsius to kelvin conversion constant was corrected from
        273.16 to 273.15). For each of these data products, the reference values generated
        by the original unchanged code were checked against the 'old' values given by the
        refactored code to make sure they were identical. The 'new' values calculated by
        the refactored code were then incorporated as the target test values for the python
        code to calculate.

        The code is available in the data product specification artifact section of
        alfresco referenced below.

            Main (calling) code:
                WarmCoolLayer_OOI_DPA_calculation_hourlydata.m

            Code to generate sub-hourly dataset to match hourly test values in Main
                make_hourly_dataset_for_metbk.m

            Reference code:
                coare35vnWarm_edson.m
                coare35vn_edson.m

            Refactored code:
                coare35vnWarm_rad_v4.m
                coare35vn_rad_v4.m

    Implemented by:

        2014-09-13: Russell Desiderio. Initial Code.
        2014-09-20: Russell Desiderio. Added tests of DPAs on sub-hourly testdata.

    References

        OOI (2014). 1341-00370_BULKFLX Artifacts. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> REFERENCE >> Data Product Specification Artifacts
            >> 1341-00370_BULKFLX
    """

    def test_met_stablty(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[-0.75134556, -1.53410291, -1.29155675],   # [00]
                          [-0.67713848, -1.40615623, -1.17706256],   # [01]
                          [-0.75134556, -1.53871067, -1.29598877],   # [10]
                          [-0.67713848, -1.41154284, -1.18240994]])  # [11]

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

    def test_met_rainflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[0.00000000, 84.63739038, 7.58665072],
                          [0.00000000, 79.07514491, 7.04086184],
                          [0.00000000, 84.80046816, 7.60371799],
                          [0.00000000, 79.26057671, 7.06089135]])

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

    def met_mommflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[0.030974330742, 0.016090125876, 0.020303436138],
                          [0.030614459503, 0.015910513076, 0.020067006629],
                          [0.030974330742, 0.016096471089, 0.020312413574],
                          [0.030614459503, 0.015918213926, 0.020078251168]])

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

    def test_met_sensflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[24.34435275, 19.84019748, 23.96742950],
                          [20.96527540, 17.52684365, 21.07171423],
                          [24.34435275, 19.92468342, 24.08112584],
                          [20.96527540, 17.62294816, 21.20514191]])

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

    def test_met_latnflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[184.91334211, 133.43175366, 151.19456789],
                          [170.45205774, 123.62963458, 139.11084942],
                          [184.91334211, 133.78969897, 151.66954491],
                          [170.45205774, 124.03365974, 139.66250129]])

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

    def test_met_buoyflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[38.41602527, 29.98454817, 35.45099283],
                          [33.93646240, 26.92597342, 31.63749229],
                          [38.41602527, 30.09636019, 35.60092032],
                          [33.93646240, 27.05292615, 31.81300501]])

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

    def test_met_buoyfls(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[36.10919371, 28.32153986, 33.56844146],
                          [31.81003830, 25.38513247, 29.90539752],
                          [36.10919371, 28.42887220, 33.71242942],
                          [31.81003830, 25.50702812, 30.07401106]])

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

    def test_met_netlirr(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[41.43188924, 28.54298163, 42.15187511],
                          [39.19850664, 26.60211114, 39.95207621],
                          [41.43188924, 28.61328117, 42.23752326],
                          [39.19850664, 26.68338061, 40.05425084]])

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

    def test_met_tempa2m(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[28.36851533, 28.09024408, 27.81395793],
                          [28.35732615, 28.08200008, 27.80462176],
                          [28.36851533, 28.09054198, 27.81431854],
                          [28.35732615, 28.08234637, 27.80505917]])

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

    def test_met_sphum2m(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[19.42297996, 19.71138107, 19.47478118],
                          [19.41394979, 19.70380669, 19.46619271],
                          [19.42297996, 19.71166621, 19.47512583],
                          [19.41394979, 19.70411212, 19.46658035]])

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

    def test_met_wind10m(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[4.84933264, 3.42025464, 3.87210524],
                          [4.85069650, 3.42094844, 3.87301659],
                          [4.84933264, 3.42023139, 3.87207237],
                          [4.85069650, 3.42091723, 3.87297121]])

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

    def test_met_frshflx(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[0.27428050, -17.10209951, -1.27575685],
                          [0.25283019, -17.11663761, -1.29367873],
                          [0.27428050, -17.10156642, -1.27504935],
                          [0.25283019, -17.11603581, -1.29285692]])

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

    def test_met_tempskn(self):
        # cases: [jwarm, jcool]
        xpctd = np.array([[31.12500000, 31.03300000, 31.02200000],
                          [30.76398727, 30.71905805, 30.66606321],
                          [31.12500000, 31.04435295, 31.03583298],
                          [30.76398727, 30.73222317, 30.68262322]])

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

    def test_met_heatflx(self):
        # cases: [jwarm, jcool]
        # total heat flux = latnflx + sensflx + rainflx + netlirr - netsirr
        latnf = np.array([[184.91334211, 133.43175366, 151.19456789],
                          [170.45205774, 123.62963458, 139.11084942],
                          [184.91334211, 133.78969897, 151.66954491],
                          [170.45205774, 124.03365974, 139.66250129]])

        sensf = np.array([[24.34435275, 19.84019748, 23.96742950],
                          [20.96527540, 17.52684365, 21.07171423],
                          [24.34435275, 19.92468342, 24.08112584],
                          [20.96527540, 17.62294816, 21.20514191]])

        rainf = np.array([[0.00000000, 84.63739038, 7.58665072],
                          [0.00000000, 79.07514491, 7.04086184],
                          [0.00000000, 84.80046816, 7.60371799],
                          [0.00000000, 79.26057671, 7.06089135]])

        netli = np.array([[41.43188924, 28.54298163, 42.15187511],
                          [39.19850664, 26.60211114, 39.95207621],
                          [41.43188924, 28.61328117, 42.23752326],
                          [39.19850664, 26.68338061, 40.05425084]])

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

    """
        Description

        Tests for products that require neither coolskin nor warmlayer corrections.

        Except as noted, test values were generated in matlab and implemented by
        Russell Desiderio.
    """
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

    def test_met_rainrte(self):
        # SCALAR CASE
        xpctd_scalar = 0.0
        calc = mb.met_rainrte(42.0, 6540998.0)
        np.testing.assert_allclose(calc, xpctd_scalar, rtol=1.e-8, atol=0.0)

        # VECTOR CASE
        xpctd = np.array([0.0, 17.3, 1.5])
        calc = mb.met_rainrte(self.cumu_prcp, self.timestamp)

        np.testing.assert_allclose(calc, xpctd, rtol=1.e-8, atol=0.0)

    def test_met_salsurf(self):
        """
            Test values generated from ctd_functions.ctd_pracsal, using ztmpwat [m]
            as a proxy for pressure [db].
        """
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

    def test_met_current_direction(self):
        # test data inputs - hit all quadrants.
        hf_rt3 = np.sqrt(3.0)/2.0
        vle = np.array([0.0, hf_rt3, 0.5, 0.0, -0.5, -hf_rt3,
                        -1.0, -hf_rt3, -0.5, 0.0, 0.5, hf_rt3])
        vln = np.array([0.0, 0.5, hf_rt3, 1.0, hf_rt3, 0.5,
                        0.0, -0.5, -hf_rt3, -1.0, -hf_rt3, -0.5])

        # expected output
        xpctd = np.array([90.0, 60.0, 30.0, 0.0, 330.0, 300.0,
                          270.0, 240.0, 210.0, 180.0, 150.0, 120.0])

        # compute the surface current directions:
        # 0 = North, 90 = East.
        calc = mb.met_current_direction(vle, vln)

        # and compare the expected to the calculated
        np.testing.assert_array_almost_equal(calc, xpctd, decimal=4)

    def test_met_relative_wind(self):
        hf_rt3 = np.sqrt(3.0)/2.0
        wind_vle = np.array([0.0, hf_rt3, 0.5, 0.0, -0.5, -hf_rt3,
                            -1.0, -hf_rt3, -0.5, 0.0, 0.5, hf_rt3]) * 3.0
        wind_vln = np.array([0.0, 0.5, hf_rt3, 1.0, hf_rt3, 0.5,
                             0.0, -0.5, -hf_rt3, -1.0, -hf_rt3, -0.5]) * 3.0
        current_vle = -0.5 * 2.0
        current_vln = -hf_rt3 * 2.0

        xpctd_relspeed = np.array([2.000000, 4.836559, 5.000000, 4.836559, 4.358899, 3.605551,
                                   2.645751, 1.614836, 1.000000, 1.614836, 2.645751, 3.605551])

        calc = mb.met_relwind_speed(wind_vle, wind_vln, current_vle, current_vln)
        np.testing.assert_array_almost_equal(calc, xpctd_relspeed, decimal=5)
        #
        #
        xpctd_reldir = mb.met_current_direction(wind_vle-current_vle, wind_vln-current_vln)
        calc = mb.met_relwind_direction(wind_vle, wind_vln, current_vle, current_vln)
        np.testing.assert_array_almost_equal(calc, xpctd_reldir, decimal=5)

    def test_metbk_windavg(self):
        """
            Date, lat, and lon provided approximately equal magnetic
            declinations of -17 and +17 degrees

            have to convert these dates to ntp timestamp (seonds since
            1900-01-01)

            C. Wingard 2014-07-01
        """

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

    def test_met_current_speed(self):
        """
            Test the surface current algorithm using test data generated in Matlab
            from the compass plot function example:

            >> rng(0,'twister') % initialize random number generator
            >> M = randn(15,15);
            >> Z = eig(M);
            >> vle = real(Z);
            >> vln = imag(Z);
            >> crnt = sqrt(vle.^2 + vln.^2);

            C. Wingard 2014-07-01
        """

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

        # compute the surface current
        out = mb.met_current_speed(vle, vln)

        # and compare the expected to the calculated
        np.testing.assert_array_almost_equal(out, crnt, decimal=4)

    def test_make_hourly_data(self):
        """
            Note that except for the first two elements, the expected values are the
            same as those given above in the "scalar testset".
                The first element contains the hourly averages of the cumulative
                    precipitation; when differenced, this will give the expected
                    rain rate array [0.0, 17.3, 1,5].
                The second element contains the hourly timestamps which have been
                    coded to be at the midpoint of the bins.

            Initial code: 2014-09-20, Russell Desiderio.
        """
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
