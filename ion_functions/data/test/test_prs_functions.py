#!/usr/bin/env python

"""
@package ion_functions.test.prs_functions
@file ion_functions, est/prs_functions.py
@author Christopher Wingard
@brief Unit tests for prs_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import prs_functions as prsfunc
from ion_functions.utils import fill_value


@attr('UNIT', group='func')
class TestPRSFunctionsUnit(BaseUnitTestCase):

    # setup up a list with the test data used by the three test functions
    lily = [
        [-150.000, 110.000, 0.15, 'N9651', 183, 186.0, -36.3, 129],
        [-140.000, 100.000, 71.99, 'N9651', 170, 172.0, -35.5, 116],
        [-130.000, 90.000, 144.21, 'N9651', 92, 158.1, -34.7, 37],
        [-120.000, 80.000, 215.53, 'N9651', 317, 144.2, -33.7, 261],
        [-110.000, 70.000, 288.12, 'N9651', 204, 130.4, -32.5, 146],
        [-100.000, 60.000, 359.76, 'N9651', 183, 116.6, -31.0, 124],
        [-90.000, 50.000, 0.15, 'N9676', 200, 103.0, -29.1, 139],
        [-80.000, 40.000, 71.99, 'N9676', 132, 89.4, -26.6, 69],
        [-70.000, 30.000, 144.21, 'N9676', 53, 76.2, -23.2, 346],
        [0.000, 20.000, 215.53, 'N9676', 332, 20.0, 90.0, 332],
        [0.000, 10.000, 288.12, 'N9676', 263, 10.0, 90.0, 263],
        [0.000, 0.000, 359.76, 'N9676', 200, 0.0, 0.0, 290],
        [0.000, -10.000, 0.15, 'N9656', 166, 10.0, -90.0, 346],
        [0.000, -20.000, 71.99, 'N9656', 130, 20.0, -90.0, 310],
        [10.000, -30.000, 144.21, 'N9656', 87, 31.6, -71.6, 249],
        [20.000, -40.000, 215.53, 'N9656', 351, 44.7, -63.4, 144],
        [30.000, -50.000, 288.12, 'N9656', 240, 58.3, -59.0, 29],
        [40.000, -60.000, 359.76, 'N9656', 166, 72.1, -56.3, 312],
        [50.000, -70.000, 0.15, 'N9652', 173, 86.0, -54.5, 317],
        [60.000, -80.000, 71.99, 'N9652', 133, 100.0, -53.1, 276],
        [70.000, -90.000, 144.21, 'N9652', 82, 114.0, -52.1, 224],
        [80.000, -100.000, 215.53, 'N9652', 347, 128.1, -51.3, 128],
        [90.000, -110.000, 288.12, 'N9652', 243, 142.1, -50.7, 24],
        [100.000, -120.000, 359.76, 'N9652', 173, 156.2, -50.2, 313],
        [110.000, -130.000, 0.15, 'N9655', 173, 170.3, -49.8, 313],
        [120.000, 160.000, 71.99, 'N9655', 152, 200.0, 53.1, 189],
        [130.000, 150.000, 144.21, 'N9655', 98, 198.5, 49.1, 139],
        [140.000, 140.000, 215.53, 'N9655', 337, 198.0, 45.0, 22],
        [150.000, 130.000, 288.12, 'N9655', 222, 198.5, 40.9, 271],
        [160.000, 120.000, 359.76, 'N9655', 173, 200.0, 36.9, 226]
    ]

    def test_prs_bottilt_ccmp(self):
        """
        Test prs_bottilt_ccmp function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)

        Note, DPS does not specify a test for this function. Values from CCMP
        lookup table in DPS used to create this test by the function and
        test_function author.

        Implemented by Christopher Wingard, July 2013
        """
        # set known inputs
        scmp = np.atleast_1d([row[2] for row in self.lily])
        snum = np.atleast_1d([row[3] for row in self.lily])

        # set known output for the corrected compass direction
        ccmp = np.atleast_1d([row[4] for row in self.lily])

        # calculate the corrected compass direction
        out = prsfunc.prs_bottilt_ccmp(scmp, snum)

        # How'd we do?
        np.testing.assert_array_equal(out, ccmp)

    def test_prs_bottilt_tmag(self):
        """
        Test prs_bottilt_tmag function.

        Note, DPS specifies a test data set that does not vary the X-tilt and
        Y-tilt inputs to this function. They are 330 and -330 microradians,
        respectively, for the entirety of the dataset. Test values below were
        created by function and test_function author to ensure function
        operates as expected over a range of potential values.

        OOI (2012). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)

        Implemented by Christopher Wingard, July 2013
        """
        # set inputs
        xtilt = np.atleast_1d([row[0] for row in self.lily])
        ytilt = np.atleast_1d([row[1] for row in self.lily])

        # set known output for the tilt magnitude
        tmag = np.atleast_1d([row[5] for row in self.lily])

        # calculate the tilt magnitude
        out = prsfunc.prs_bottilt_tmag(xtilt, ytilt)

        # How'd we do?
        np.testing.assert_allclose(out, tmag, rtol=0.1, atol=0.1)

    def test_prs_bottilt_tdir(self):
        """
        Test prs_bottilt_tdir function.

        Note, DPS specifies a test data set that does not vary the X-tilt and
        Y-tilt inputs to this function. They are 330 and -330 microradians,
        respectively, for the entirety of the dataset. Test values below were
        created by function and test_function author to ensure function
        operates as expected over a range of potential values.

        OOI (2012). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)

        Implemented by Christopher Wingard, July 2013
        """
        # set inputs
        xtilt = np.atleast_1d([row[0] for row in self.lily])
        ytilt = np.atleast_1d([row[1] for row in self.lily])
        ccmp = np.atleast_1d([row[4] for row in self.lily])

        # set known output
        tdir = np.atleast_1d([row[7] for row in self.lily])

        # calculate the tilt direction
        out = prsfunc.prs_bottilt_tdir(xtilt, ytilt, ccmp)

        # How'd we do?
        np.testing.assert_array_equal(out, tdir)
