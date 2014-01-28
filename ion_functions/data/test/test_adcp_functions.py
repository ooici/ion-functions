"""
@package ion_functions.test.adcp_functions
@file ion_functions/test/test_adcp_functions.py
@author Christopher Wingard
@brief Unit tests for adcp_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import adcp_functions as adcpfunc
from ion_functions.data.wmm import WMM
from ion_functions.data.generic_functions import wmm_model


@attr('UNIT', group='func')
class TestCTDFunctionsUnit(BaseUnitTestCase):

    def test_adcp_beam2earth(self):
        """
        Test adcp_beam2ins and adcp_ins2earth functions.

        Values based on those defined in DPS:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """
        # set test inputs
        b1 = np.array([-0.0300, -0.2950, -0.5140, -0.2340, -0.1880,
                       0.2030, -0.3250,  0.3050])
        b2 = np.array([0.1800, -0.1320,  0.2130,  0.3090,  0.2910,
                       0.0490,  0.1880,  0.3730])
        b3 = np.array([-0.3980, -0.4360, -0.1310, -0.4730, -0.4430,
                       0.1880, -0.1680,  0.2910])
        b4 = np.array([-0.2160, -0.6050, -0.0920, -0.0580,  0.4840,
                       -0.0050,  0.3380,  0.1750])
        heading = 98.4100
        pitch = 0.6900
        roll = -2.5400
        orient = 1

        # compute beam to instrument transform
        u, v, w, e = adcpfunc.adcp_beam2ins(b1, b2, b3, b4)

        # compute instrument to Earth transform
        got_uu, got_vv, got_ww = adcpfunc.adcp_ins2earth(
            u, v, w, heading, pitch, roll, orient)

        # set expected results
        uu = np.array([0.2175, -0.2814, -0.1002,  0.4831,  1.2380,
                       -0.2455,  0.6218, -0.1807])
        vv = np.array([-0.3367, -0.1815, -1.0522, -0.8676, -0.8919,
                       0.2585, -0.8497, -0.0873])
        ww = np.array([0.1401,  0.3977,  0.1870,  0.1637,  0.0091,
                       -0.1290,  0.0334, -0.3017])

        # test the transform
        self.assertTrue(np.allclose(got_uu, uu, rtol=1e4, atol=0))
        self.assertTrue(np.allclose(got_vv, vv, rtol=1e4, atol=0))
        self.assertTrue(np.allclose(got_ww, ww, rtol=1e4, atol=0))

    #
    #**** NO LONGER USED BECAUSE OF WMM.VELOCITY_CORRECTIONS ****
    #
    #def test_adcp_magvar(self):
    #    """
    #    Test adcp_beam2ins and adcp_ins2earth functions.
    #
    #    Values based on those defined in DPS:
    #
    #    OOI (2012). Data Product Specification for Velocity Profile and Echo
    #        Intensity. Document Control Number 1341-00750.
    #        https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
    #        >> Controlled >> 1000 System Level >>
    #        1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)
    #
    #    Implemented by Christopher Wingard, April 2013
    #    """
    #    # apply the correction, assumption in adcpfunc.adcp_magvar is the
    #    # magnetic declination would be a scalar and uu and vv would be
    #    # velocity profiles (e.g. arrays)
    #    wmm = WMM(wmm_model)
    #    uu_cor, vv_cor = wmm(16.9604, np.array([0.4413]),
    #                                          np.array([0.1719]))
    #
    #    # test the transform
    #    self.assertTrue(np.allclose(uu_cor, 0.4722, rtol=1e4, atol=0))
    #    self.assertTrue(np.allclose(vv_cor, 0.0357, rtol=1e4, atol=0))
