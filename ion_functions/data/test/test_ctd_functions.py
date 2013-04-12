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

    def test_ctd_sbe16plus_tempwat(self):
        """
        Test ctd_sbe16plus_tempwat function.

        Values based on those described in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)
            
        Implemented by Christopher Wingard, April 2013
        """
        t0 = 248471
        a0 = 1.281651e-3
        a1 = 2.706002e-4
        a2 = -1.027561e-6
        a3 = 1.749446e-7
        
        t = ctdfunc.ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3)
        self.assertTrue(np.allclose(t, 22.544681, rtol=1e-6, atol=0))
    

    def test_ctd_sbe16plus_preswat(self):
        """
        Test ctd_pracsal function.

        Values based on those defined in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
            
        Implemented by Christopher Wingard, April 2013
        """
        p0 = 528418
        therm0 = 24303
        ptempa0 = -6.87701000e+001
        ptempa1 = 5.05406200e+001
        ptempa2 = -2.15672900e-001
        ptca0 = 5.24965500e+005
        ptca1 = 7.23620100e+000
        ptca2 = -9.94485900e-002
        ptcb0 = 2.51220000e+001
        ptcb1 = -2.00000000e-004
        ptcb2 = 0.00000000e+000
        pa0 = 1.73472300e+000
        pa1 = 1.57475000e-002
        pa2 = -6.51927800e-010
        
        p = ctdfunc.ctd_sbe16plus_preswat(p0, therm0, ptempa0, ptempa1, ptempa2,
                          ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                          pa0, pa1, pa2)
        self.assertTrue(np.allclose(p, 27.282116, rtol=1e-6, atol=0))


    def test_ctd_sbe16plus_condwat(self):
        """
        Test ctd_pracsal function.

        Values based on those defined in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_PRACSAL_OOI.pdf)
            
        Implemented by Christopher Wingard, March 2013
        """
        c0 = 1673175
        t1 = 22.544681
        p1 = 27.282116
        g = -9.72193700e-001
        h = 1.38675900e-001
        i = -1.08398500e-004
        j = 2.63219300e-005
        cpcor = -9.57000000e-008
        ctcor = 3.2500e-006

        c = ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor)
        self.assertTrue(np.allclose(c, 4.969069, rtol=1e-6, atol=0))


    def test_ctd_pracsal(self):
        """
        Test ctd_pracsal function.

        Values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Salinty. Document Control
            Number 1341-00040. https://alfresco.oceanobservatories.org/ (See: 
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_PRACSAL_OOI.pdf)
            
        Implemented by Christopher Wingard, March 2013
        """

        c = np.array([5.407471, 5.407880, 5.041008, 3.463402, 3.272557, 3.273035])
        t = np.array([28, 28, 20, 6, 3, 2])
        p = np.array([0, 10, 150, 800, 2500, 5000])

        output = ctdfunc.ctd_pracsal(c,t,p)
        
        """
        Note, DPS rounds off output values to %.1f. For test to work, these were
        recalculated using the GSW Toolbox, Version 3.02 in Matlab R2013a and
        output using %.6f (see Matlab code snippet below). The DPS will be
        editted to correctly specify the higher precision.
        
        >> sprintf('%.6f\t',gsw_SP_from_C(c*10,t,p))
        ans =
        33.495229	33.495224	36.995774	34.898526	34.999244	34.999494
        """
        check_values = np.array([33.495229,
                                 33.495224,
                                 36.995774,
                                 34.898526,
                                 34.999244,
                                 34.999494])
        self.assertTrue(np.allclose(output, check_values))

    
    def test_ctd_density(self):
        """
        Test ctd_pracsal function.

        Values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Density. Document Control
            Number 1341-00050. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf)
            
        Implemented by Christopher Wingard, March 2013
        """
        
        SP = np.array([33.5, 33.5, 37, 34.9, 35, 35])
        t = np.array([28, 28, 20, 6, 3, 2])
        p = np.array([0, 10, 150, 800, 2500, 5000])
        lat = 15.00
        lon = -55.00
        
        output = ctdfunc.ctd_density(SP,t,p,lat,lon)

        check_values = np.array([1021.26851,
                                 1021.31148,
                                 1026.94422,
                                 1031.13498,
                                 1039.28768,
                                 1050.30616])
        self.assertTrue(np.allclose(output, check_values))
