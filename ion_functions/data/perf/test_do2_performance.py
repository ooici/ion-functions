#!/usr/bin/env python

from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data.do2_functions import do2_dofst_frequency, do2_dofst_volt

import numpy as np

# number of elements in the array to run
NeLEMENTS = 10000


class TestDo2Performance(PerformanceTestCase):
    """
    Performance tests for the DO2 family of functions to compare to the
    optimization benchmark of running the function for 10000 data
    records in under 10 seconds.

    Implemented by:
        2014-03-04: Stuart Pearce. Initial code.
    """
    def setUp(self):
        # create Nx1 arrays (standard is 10000)
        self.lat = np.ones(NeLEMENTS) * 45.0  # latitude
        self.lon = np.ones(NeLEMENTS) * -125.0  # longitude
        self.pres = np.ones(NeLEMENTS) * 112.1  # pressure
        self.temp = np.ones(NeLEMENTS) * 20.2  # temperature
        self.salt = np.ones(NeLEMENTS) * 33.4  # salinity

    def test_dofst_volt(self):
        """
        Performance test for the do2_dofst_volt function for the DOFST
        (SBE 43) instrument data processing algorithm.
        """
        # N voltage counts
        volt_counts = np.ones(NeLEMENTS, dtype=np.int) * 16384
        # calibration coefficients
        A = -3.1867e-3
        B = 1.7749e-4
        C = -3.5718e-6
        E = 0.036
        Voffset = -0.5186
        Soc = 0.4396

        # run performance test
        stats = []
        self.profile(
            stats, do2_dofst_volt,
            volt_counts, Voffset, Soc, A, B, C, E,
            self.pres, self.temp, self.salt, self.lat, self.lon)

    def test_dofst_frequency(self):
        """
        Performance test for the do2_dofst_volt function for the DOFST
        (SBE 43F) instrument data processing algorithm.
        """
        freq = np.ones(NeLEMENTS, dtype=np.int) * 4591
        A = -4.1168e-3
        B = 2.4818e-4
        C = -3.8820e-6
        E = 0.036
        Foffset = -839.55
        Soc = 2.9968e-4

        stats = []
        self.profile(
            stats, do2_dofst_frequency,
            freq, Foffset, Soc, A, B, C, E,
            self.pres, self.temp, self.salt, self.lat, self.lon)
