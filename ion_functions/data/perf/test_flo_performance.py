"""
@package ion_functions.data.perf.test_adcp_performance
@file ion_functions/data/perf/test_adcp_performance.py
@author Christopher Wingard
@brief Performance tests for adcp_functions module
"""

import numpy as np
from nose.plugins.attrib import attr

from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data import flo_functions as fl


@attr('PERF', group='func')
class TestADCPPerformance(PerformanceTestCase):

    def setUp(self):
        # set test inputs
        ### optical backscatter ###
        self.scat_counts = 55
        self.scat_dark = 47
        self.scat_scale = 3.058e-6
        self.beta = fl.flo_beta(self.scat_counts, self.scat_dark, self.scat_scale)
        self.degC = 20.0
        self.psu = 32.0
        self.theta = 140.0
        self.wavelength = 700.0
        self.chi_factor = 1.096

        ### chla ###
        self.chla_counts = 55
        self.chla_dark = 45
        self.chla_scale = 0.0121

        ### cdom ###
        self.cdom_counts = 55
        self.cdom_dark = 48
        self.cdom_scale = 0.0848

    def test_flo_bback_total(self):
        stats = []

        beta = np.repeat(self.beta, 1000000)
        degC = np.repeat(self.degC, 1000000)
        psu = np.repeat(self.psu, 1000000)
        theta = np.repeat(self.theta, 1000000)
        wavelength = np.repeat(self.wavelength, 1000000)
        chi_factor = np.repeat(self.chi_factor, 1000000)

        self.profile(stats, fl.flo_bback_total, beta, degC, psu, theta, wavelength, chi_factor)

    def test_flo_beam(self):
        stats = []

        counts = np.repeat(self.scat_counts, 1000000)
        dark = np.repeat(self.scat_dark, 1000000)
        scale = np.repeat(self.scat_scale, 1000000)

        self.profile(stats, fl.flo_beta, counts, dark, scale)

    def test_flo_cdom(self):
        stats = []

        counts = np.repeat(self.cdom_counts, 1000000)
        dark = np.repeat(self.cdom_dark, 1000000)
        scale = np.repeat(self.cdom_scale, 1000000)

        self.profile(stats, fl.flo_cdom, counts, dark, scale)

    def test_flo_chla(self):
        stats = []

        counts = np.repeat(self.chla_counts, 1000000)
        dark = np.repeat(self.chla_dark, 1000000)
        scale = np.repeat(self.chla_scale, 1000000)

        self.profile(stats, fl.flo_chla, counts, dark, scale)
