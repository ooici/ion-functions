#!/usr/bin/env python

"""
@package ion_functions.test.hyd_functions
@file ion_functions/test/test_hyd_functions.py
@author Christopher Wingard
@brief Unit tests for hyd_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import hyd_functions as hydfunc


@attr('UNIT', group='func')
class TestHYDFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        # setup the test values for HYDAPBB
        self.gain = 6.
        self.wav = np.array([-2.40000, -0.31200, 0.01110, 0.00442])
        self.tsv = np.atleast_2d(np.array([-3.60855, -0.46911, 0.01669, 0.00665]))

        # setup the test values for HYDAPLF
        self.raw = np.array([0, 512, 1024, 2048, 3072, 4096])
        self.hydaplf = np.array([0.00000E+00, 1.63840E-03, 3.27680E-03,
                                 6.55360E-03, 9.83040E-03, 1.31072E-02])

    def test_hyd_bb_acoustic_pwaves(self):
        """
        Test hyd_bb_acoustic_pwaves function.

        Values based on those described in DPS, with modification to account
        for errors and inconsistencies in the DPS.

        OOI (2013). Data Product Specification for Acoustic Pressure Waves.
            Document Control Number 1341-00820.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00820_Data_Product_SPEC_HYDAPBB_OOI.pdf)

        Implemented by Christopher Wingard, May 2014
        """
        # calculate the broadband acoustic pressure wave
        out = hydfunc.hyd_bb_acoustic_pwaves(self.wav, self.gain)

        # How'd we do?
        np.testing.assert_allclose(out, self.tsv, rtol=1e-5, atol=1e-5)

    def test_hyd_lf_acoustic_pwaves(self):
        """
        Test hyd_lf_acoustic_pwaves function.

        Values based on a test data set created using Excel by the test author
        (DPS does not provide any test data).

        Implemented by Christopher Wingard, July 2014
        """
        # calculate the low frequency acoustic pressure wave
        out = hydfunc.hyd_lf_acoustic_pwaves(self.raw)

        # How'd we do?
        np.testing.assert_allclose(out, self.hydaplf, rtol=1e-6, atol=1e-6)
