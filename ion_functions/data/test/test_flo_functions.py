#!/usr/bin/env python
"""
@package ion_functions.test.flo_functions
@file ion_functions/test/flo_functions.py
@author Craig Risien
@brief Unit tests for flo_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import flo_functions as flofunc


@attr('UNIT', group='func')
class TestSFLFunctionsUnit(BaseUnitTestCase):

    def test_flo_bback_total(self):
        """
        Test flo_bback_total function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Optical Backscatter (Red
            Wavelengths). Document Control Number 1341-00540.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00540_Data_Product_SPEC_FLUBSCT_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """
        counts_output = np.array([
            55, 57, 55, 56, 54, 54, 55, 54, 55, 56, 55
        ])
        counts_dark = np.ones(11) * 47
        scale_factor = np.ones(11) * 3.058e-6
        beta = flofunc.flo_beta(counts_output, counts_dark, scale_factor)

        # set known outputs (adjusted to account for incorrect DPS wherein test
        # values were computed using different formula and Chi (X) Factor).
        bback = np.array([0.000171, 0.000213, 0.000171, 0.000192,
                          0.000151, 0.000151, 0.000171, 0.000151,
                          0.000171, 0.000192, 0.000171])

        # set default constants
        degC = np.ones(11) * 20.0
        psu = np.ones(11) * 32.0
        #theta = np.ones(11) * 117.0
        #wlngth = np.ones(11) * 700.0
        #xfactor = np.ones(11) * 1.08

        # calculate the total optical backscatter
        bbout = flofunc.flo_bback_total(beta, degC, psu)

        # How'd we do?
        np.testing.assert_allclose(bbout, bback, rtol=1e-6, atol=1e-6)

    def test_flo_chla_cdom(self):
        """
        Test the flo_chla and flo_cdom functions.

        Values based on test data in DPSs and available on Alfresco.

        OOI (2012). Data Product Specification for Fluorometric Chlorophyll-a
            Concentration. Document Control Number 1341-00530.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00530_Data_Product_SPEC_CHLAFLO_OOI.pdf)

        OOI (2012). Data Product Specification for Fluorometric CDOM
            Concentration. Document Control Number 1341-00550.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00550_Data_Product_SPEC_CDOMFLO_OOI.pdf)

        Implemented by Craig Risien, February 2014
        """
        # test inputs
        chla_counts_dark = np.ones(16) * 45
        chla_scale_factor = np.ones(16) * 0.0121
        cdom_counts_dark = np.ones(16) * 48
        cdom_scale_factor = np.ones(16) * 0.0848

        chla_counts_output = np.array([54, 52, 51, 52, 52, 51, 51, 52,
                                       52, 54, 51, 51, 51, 50, 50, 52])
        cdom_counts_output = np.array([51, 50, 50, 52, 50, 51, 51, 51,
                                       51, 51, 49, 52, 52, 53, 52, 51])

        # expected outputs
        chla_expected = np.array([0.1089, 0.0847, 0.0726, 0.0847,
                                  0.0847, 0.0726, 0.0726, 0.0847,
                                  0.0847, 0.1089, 0.0726, 0.0726,
                                  0.0726, 0.0605, 0.0605, 0.0847])
        cdom_expected = np.array([0.2544, 0.1696, 0.1696, 0.3392,
                                  0.1696, 0.2544, 0.2544, 0.2544,
                                  0.2544, 0.2544, 0.0848, 0.3392,
                                  0.3392, 0.424, 0.3392, 0.2544])

        # compute chla and cdom values
        chla_calc = flofunc.flo_chla(chla_counts_output, chla_counts_dark, chla_scale_factor)
        cdom_calc = flofunc.flo_cdom(cdom_counts_output, cdom_counts_dark, cdom_scale_factor)

        # compare calculated results to expected
        np.testing.assert_allclose(chla_calc, chla_expected, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(cdom_calc, cdom_expected, rtol=1e-6, atol=1e-6)
