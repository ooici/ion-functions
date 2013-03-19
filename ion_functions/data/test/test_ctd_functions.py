#!/usr/bin/env python

"""
@package ion_functions.test.ctd_functions
@file ion_functions/test/ctd_functions.py
@author Christopher Wingard
@brief Unit tests for ctd_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import ctd_functions as ctdfunc

@attr('UNIT', group='func')
class TestCTDFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_ctd_pracsal(self):
        """
        Test ctd_pracsal function.

        Values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Salinty. Document Control
            Number 1341-00040. https://alfresco.oceanobservatories.org/ (See: 
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_PRACSAL_OOI.pdf)    
        """

        c = np.array([5.407471, 5.407880, 5.041008, 3.463402, 3.272557, 3.273035])
        t = np.array([28, 28, 20, 6, 3, 2])
        p = np.array([0, 10, 150, 800, 2500, 5000])

        output = ctdfunc.ctd_pracsal(c,t,p)
        
        """
        Note, DPS rounds off output values to the tenth decimal place. For test
        to work, these were recalculated using the GSW Toolbox, Version 3.02 in
        Matlab R2013a and output to the 6th decimal place.
        
        >> sprintf('%.6f\t',gsw_SP_from_C(c*10,t,p))
        ans =
        33.495229	33.495224	36.995774	34.898526	34.999244	34.999494
        """
        check_values = np.array((33.495229,
                                 33.495224,
                                 36.995774,
                                 34.898526,
                                 34.999244,
                                 34.999494))
        self.assertTrue(np.allclose(output, check_values))
    
    def test_ctd_density(self):
        """
        Test ctd_pracsal function.

        Values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Density. Document Control
            Number 1341-00050. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf)
        """
        
        SP = np.array([33.5, 33.5, 37, 34.9, 35, 35])
        t = np.array([28, 28, 20, 6, 3, 2])
        p = np.array([0, 10, 150, 800, 2500, 5000])
        lat = 15.00
        lon = -55.00
        
        output = ctdfunc.ctd_density(SP,t,p,lat,lon)

        check_values = np.array((1021.26851,
                                 1021.31148,
                                 1026.94422,
                                 1031.13498,
                                 1039.28768,
                                 1050.30616))
        self.assertTrue(np.allclose(output, check_values))
